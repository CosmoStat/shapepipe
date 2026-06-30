#!/usr/bin/env python3

"""
Initialise the directory structure for a v2.0 ShapePipe run.
Creates tile and exposure staging directories, config symlinks,
and a tile number list.

Config file is loaded first (if provided), then overridden by command-line arguments.
"""

import argparse
import sys
from pathlib import Path
import os

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML not installed. Install with: pip install pyyaml", file=sys.stderr)
    sys.exit(1)


class InitRun:
    def __init__(self):
        self.sp_root = Path.home() / "shapepipe"
        self.quiet = False
        self.defaults = {
            "type": "data",
            "subdir": "",
            "base_dir": ".",
            "config_dir": "",
            "params_src": str(Path.home() / "sp_validation" / "notebooks" / "params.py"),
            "tiles_src": "",
            "tile_ids_list": None,
            "sample": "",
            "input_sims_base_dir": "",
            "input_tiles": "",
            "input_exp": "",
            "sims_type": "",
            "num": "",
        }
        self.config = {}
        self.args = {}

    def load_config(self, config_file):
        """Load YAML config file."""
        config_path = Path(config_file)
        if not config_path.exists():
            print(f"ERROR: Config file not found: {config_file}", file=sys.stderr)
            sys.exit(2)

        with open(config_path) as f:
            self.config = yaml.safe_load(f) or {}

        # Map config keys to script parameters
        if "base" in self.config:
            self.config["base_dir"] = self.config.pop("base")
        if "type" in self.config:
            # Keep type from config (data or image_sims)
            self.config["type"] = self.config["type"]
        if "tile_IDs" in self.config:
            tile_ids = self.config.pop("tile_IDs")
            # tile_IDs can be a file path (string) or a list of IDs
            if isinstance(tile_ids, str):
                self.config["tiles_src"] = tile_ids
            elif isinstance(tile_ids, list):
                # Will be written to file in run directory
                self.config["tile_ids_list"] = tile_ids
        if "params_py" in self.config:
            self.config["params_src"] = self.config.pop("params_py")
        # Keep input_tiles and input_exp if present (they may contain placeholders)
        # They will be substituted before creating symlinks

    def merge_config(self):
        """Merge: defaults <- config <- command-line args."""
        params = dict(self.defaults)
        params.update(self.config)
        params.update(self.args)
        return params

    def substitute_placeholders(self, value, params):
        """Substitute {key} placeholders in path strings using params dict."""
        if not isinstance(value, str):
            return value
        result = value
        for key, val in params.items():
            if val is not None and val != "":
                # Convert to string and strip trailing slashes
                clean_val = str(val).rstrip('/')
                result = result.replace(f"{{{key}}}", clean_val)
        # Clean up any double slashes that may result
        while '//' in result:
            result = result.replace('//', '/')
        return result

    def construct_default_subdir(self, params):
        """Construct default subdir from sims_type and num if not explicitly set."""
        if not params["subdir"] and params.get("sims_type") and params.get("num"):
            params["subdir"] = f"1p2z_{params['sims_type']}_{params['num']}"

    def set_type_defaults(self, params):
        """Set type-specific defaults for config_dir and tiles_src."""
        if params["type"] == "data":
            config_dir = self.sp_root / "example" / "cfis"
            if not params["tiles_src"]:
                params["tiles_src"] = str(
                    self.sp_root / "auxdir" / "CFIS" / "tiles_202604" / "tiles_r.txt"
                )
        elif params["type"] == "image_sims":
            config_dir = self.sp_root / "example" / "cfis_image_sims"
            if not params["tiles_src"]:
                # Try to use num-specific tile file first
                if params.get("num"):
                    num_tile_file = self.sp_root / "auxdir" / "CFIS" / "im_sims_202606" / f"tile_numbers_{params['num']}.txt"
                    if num_tile_file.exists():
                        params["tiles_src"] = str(num_tile_file)
                    else:
                        params["tiles_src"] = str(
                            self.sp_root / "auxdir" / "CFIS" / "im_sims_202606" / "numbers.txt"
                        )
                else:
                    params["tiles_src"] = str(
                        self.sp_root / "auxdir" / "CFIS" / "im_sims_202606" / "numbers.txt"
                    )
        else:
            print(f"ERROR: Invalid type '{params['type']}'. Must be 'data' or 'image_sims'", file=sys.stderr)
            sys.exit(3)

        # Use config_dir from config, or fall back to type-specific default
        if not params["config_dir"]:
            params["config_dir"] = str(config_dir)
        return params

    def create_symlink(self, link_name, target, force=False):
        """Create symlink, return True if created."""
        link_path = Path(link_name)

        if link_path.exists() or link_path.is_symlink():
            if force:
                link_path.unlink()
            else:
                if link_path.is_symlink() and not self.quiet:
                    print(f"{link_name} symlink already exists, skipping")
                elif not self.quiet:
                    print(f"WARNING: {link_name} exists as a regular file, not creating symlink")
                return False

        link_path.symlink_to(target)
        if not self.quiet:
            print(f"Created symlink: {link_name} -> {target}")
        return True

    def write_tile_ids(self, tile_ids_list):
        """Write tile IDs list to tile_numbers.txt."""
        tile_file = Path("tile_numbers.txt")
        with open(tile_file, "w") as f:
            for tile_id in tile_ids_list:
                f.write(f"{tile_id}\n")
        if not self.quiet:
            print(f"Created tile_numbers.txt with {len(tile_ids_list)} tiles")

    def copy_params(self, params, force=False):
        """Copy and adapt params.py."""
        params_dest = Path("params.py")

        if params_dest.exists() and not force:
            if not self.quiet:
                print("params.py already exists, skipping")
            return

        params_src = Path(params["params_src"])
        if not params_src.exists():
            print(f"ERROR: Source params.py not found: {params_src}", file=sys.stderr)
            sys.exit(4)

        content = params_src.read_text()
        # Remove IMAFLAGS_ISO
        content = content.replace("'IMAFLAGS_ISO', ", "").replace(", 'IMAFLAGS_ISO'", "")
        content = content.replace('"IMAFLAGS_ISO", ', "").replace(', "IMAFLAGS_ISO"', "")

        if params["type"] == "image_sims":
            # Adapt for image_sims
            lines = content.split("\n")
            adapted = []
            for line in lines:
                if line.startswith("name = "):
                    line = f"name = '{params['subdir']}'"
                elif line.startswith("star_cat_path = "):
                    line = "star_cat_path = None"
                elif line.startswith("output_format = "):
                    line = "output_format = '.hdf5'"
                adapted.append(line)
            content = "\n".join(adapted)
            if not self.quiet:
                print("Copied and adapted: params.py <- " + params["params_src"])
                print(f"  name = '{params['subdir']}', star_cat_path = None, output_format = '.hdf5'")
        else:
            if not self.quiet:
                print(f"Copied: params.py <- {params['params_src']}")

        params_dest.write_text(content)

    def run(self, config_file=None, quiet=False, **kwargs):
        """Run initialization."""
        self.quiet = quiet

        # Load config first
        if config_file:
            self.load_config(config_file)

        # Apply command-line overrides (exclude quiet)
        self.args = {k: v for k, v in kwargs.items() if v is not None and k != "quiet"}

        # Merge all sources
        params = self.merge_config()

        # Construct default subdir from sims_type and num if needed
        self.construct_default_subdir(params)

        # Validate that we have a subdir
        if not params["subdir"]:
            print("ERROR: subdir must be provided (-s) or derivable from sims_type and num in config", file=sys.stderr)
            sys.exit(4)

        # Set type-specific defaults (including num-based tile_IDs for image_sims)
        params = self.set_type_defaults(params)

        # Expand placeholders in tiles_src (e.g. {num})
        if params.get("tiles_src"):
            params["tiles_src"] = self.substitute_placeholders(params["tiles_src"], params)

        # Use pwd if base_dir still empty
        if not params["base_dir"]:
            params["base_dir"] = str(Path.cwd())

        base_dir = Path(params["base_dir"]) / params["subdir"]

        if not self.quiet:
            print(f"Initialising ShapePipe run directory: {params['base_dir']}")
            print()

        # Create base directory
        base_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(base_dir)

        # Image sims specific: input symlinks
        if params["type"] == "image_sims":
            # input_tiles and input_exp symlinks (with placeholder substitution)
            if params.get("input_tiles"):
                input_tiles_path = self.substitute_placeholders(params["input_tiles"], params)
                self.create_symlink("input_tiles", input_tiles_path, params["force"])
            if params.get("input_exp"):
                input_exp_path = self.substitute_placeholders(params["input_exp"], params)
                self.create_symlink("input_exp", input_exp_path, params["force"])

        # Image sims specific: mask config symlink
        if params["type"] == "image_sims" and params.get("sample"):
            mask_src = (
                Path.home()
                / "sp_validation"
                / "config"
                / "calibration"
                / f"mask_v1.X.{params['sample']}_im_sim.yaml"
            )
            self.create_symlink("config_mask.yaml", str(mask_src), params["force"])

        # Create directories
        if not self.quiet:
            print("Creating tiles/ directory...")
        Path("tiles").mkdir(exist_ok=True)
        if not self.quiet:
            print("  tiles/ created (subdirs added at download time)")

        if not self.quiet:
            print("Creating exp/ directory...")
        Path("exp").mkdir(exist_ok=True)
        if not self.quiet:
            print("  exp/ created (subdirs added at download time)")

        if not self.quiet:
            print("Creating log and debug directories...")
        Path("logs").mkdir(exist_ok=True)
        Path("debug").mkdir(exist_ok=True)

        # Config symlink
        self.create_symlink("cfis", params["config_dir"], params["force"])

        # Tile numbers: either symlink to file or write list
        if params.get("tile_ids_list"):
            if not self.quiet:
                print("Creating tile_numbers.txt from tile_IDs list...")
            self.write_tile_ids(params["tile_ids_list"])
        else:
            if not self.quiet:
                print("Creating tile_numbers.txt symlink...")
            self.create_symlink("tile_numbers.txt", params["tiles_src"], params["force"])
            n_tiles = sum(1 for _ in open("tile_numbers.txt"))
            if not self.quiet:
                print(f"  {n_tiles} tiles")

        # params.py
        if not self.quiet:
            print("Creating params.py...")
        self.copy_params(params, params["force"])

        # Print summary
        print(f"Done: {params['subdir']}")

        if not self.quiet:
            print()
            print("Directory structure:")
            print(f"  {base_dir}")
            print("  ├── tiles/")
            print("  ├── exp/")
            print("  ├── logs/")
            print("  ├── debug/")
            print(f"  ├── cfis  ->  {params['config_dir']}")

            # Show tile_numbers.txt info
            if params.get("tile_ids_list"):
                print(f"  ├── tile_numbers.txt  ({len(params['tile_ids_list'])} tiles)")
            else:
                print(f"  ├── tile_numbers.txt  ->  {params['tiles_src']}")

            if params["type"] == "image_sims":
                if params.get("sample"):
                    mask_src = (
                        Path.home()
                        / "sp_validation"
                        / "config"
                        / "calibration"
                        / f"mask_v1.X.{params['sample']}_im_sim.yaml"
                    )
                    print(f"  ├── config_mask.yaml  ->  {mask_src}")
                if params.get("input_tiles"):
                    input_tiles_path = self.substitute_placeholders(params["input_tiles"], params)
                    print(f"  ├── input_tiles  ->  {input_tiles_path}")
                if params.get("input_exp"):
                    input_exp_path = self.substitute_placeholders(params["input_exp"], params)
                    print(f"  ├── input_exp  ->  {input_exp_path}")

            print(f"  └── params.py  <-  {params['params_src']}")


def main():
    parser = argparse.ArgumentParser(
        description="Initialise ShapePipe v2.0 run directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Read all from config file (subdir auto-generated from sims_type + num)
  init_run_v2.0.py -c config.yaml -t image_sims

  # Override config values with command-line options
  init_run_v2.0.py -c config.yaml -t image_sims -d /custom/base -S 6

  # Explicit subdir (overrides auto-generation)
  init_run_v2.0.py -c config.yaml -t image_sims -s my_subdir

  # No config file, use defaults and flags
  init_run_v2.0.py -t image_sims -d /tmp/test -s grid_1
        """,
    )

    parser.add_argument(
        "-c", "--config",
        help="Config file (YAML); entries can be overridden by following flags"
    )
    parser.add_argument(
        "-t", "--type",
        choices=["data", "image_sims"],
        help="Input type"
    )
    parser.add_argument(
        "-d", "--dir",
        dest="base_dir",
        help="Base run directory"
    )
    parser.add_argument(
        "-C", "--config-dir",
        dest="config_dir",
        help="Config directory (overrides type-specific default)"
    )
    parser.add_argument(
        "-s", "--subdir",
        help="Subdirectory for simulations (auto-generated from sims_type+num if not provided)"
    )
    parser.add_argument(
        "-T", "--tiles",
        dest="tiles_src",
        help="Tile numbers file"
    )
    parser.add_argument(
        "-i", "--input-sims",
        dest="input_sims_base_dir",
        help="Input simulations base directory (image_sims only)"
    )
    parser.add_argument(
        "--input-tiles",
        dest="input_tiles",
        help="Input tiles directory path (supports {subdir} and {input_sims_base_dir} placeholders)"
    )
    parser.add_argument(
        "--input-exp",
        dest="input_exp",
        help="Input exposures directory path (supports {subdir} and {input_sims_base_dir} placeholders)"
    )
    parser.add_argument(
        "-S", "--sample",
        help="Sample version for mask config (e.g., 6)"
    )
    parser.add_argument(
        "-P", "--params",
        dest="params_src",
        help="Source params.py to copy"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recreate existing symlinks and parameter files"
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress all output except final 'Done' message"
    )

    args = parser.parse_args()

    init = InitRun()
    init.run(
        config_file=args.config,
        type=args.type,
        base_dir=args.base_dir,
        config_dir=args.config_dir,
        subdir=args.subdir,
        tiles_src=args.tiles_src,
        input_sims_base_dir=args.input_sims_base_dir,
        input_tiles=args.input_tiles,
        input_exp=args.input_exp,
        sample=args.sample,
        params_src=args.params_src,
        force=args.force,
        quiet=args.quiet,
    )


if __name__ == "__main__":
    main()
