"""Build isolated campaigns and inspect Snakemake's resolved jobs.

No executor runs: the API's ``printdag`` operation resolves input functions
and job wildcards, including the compute parse's real SQLite index build.
Only the final access to ``WorkflowApi._workflow`` is private; Snakemake's
public API exposes graph printing but not the resolved Job objects.
"""

import io
import os
import sys
from contextlib import contextmanager, redirect_stdout
from dataclasses import dataclass
from pathlib import Path

import yaml
from snakemake.api import SnakemakeApi
from snakemake.settings.enums import RerunTrigger
from snakemake.settings.types import (
    DAGSettings,
    DeploymentSettings,
    ResourceSettings,
    WorkflowSettings,
)
from snakemake_interface_executor_plugins.settings import DeploymentMethod
from workflow.scripts import build_index

REPO = Path(__file__).resolve().parents[2]
MODES = [("data", "psfex"), ("image_sims", "fake")]
# psf_model=mccd is refused at parse time (refuse_unpersistable_psf), so it
# has no DAG to resolve; test_mccd_is_refused_during_parse pins the refusal.


@dataclass
class Campaign:
    """A campaign with two ready tiles and campaign-wide consumer edges."""

    root: Path
    input_type: str
    psf_model: str
    name: str = "dag-campaign"

    def __post_init__(self):
        """Write exposure lists, prepared manifests, and campaign history."""
        self.root.mkdir(parents=True, exist_ok=True)
        self.ready = {
            "123.456": ("2243881", "2243882"),
            "124.456": ("2243882",),
        }
        self.exposures = ("2243881", "2243882")
        self.unready = "125.456"
        self.outside = "126.456"
        self.ignored = "127.456"
        self.run_dir = self.root / "scratch" / self.name
        # The basename deliberately differs from run: to catch name sniffing.
        self.products_dir = self.root / "persistent" / self.name / "products"
        self.index_db = self.products_dir / "index" / "run_index.sqlite"
        self.tile_list = self.root / "tiles.txt"
        self.config_path = self.root / "run.yaml"
        self.state_dir = self.root / "state"
        self.image = self.root / "planning-only.sif"
        # Image resolution checks existence; DAG inspection never opens it.
        self.image.touch()
        self.state_dir.mkdir()
        self.tile_list.write_text("\n".join([
            *self.ready, self.unready, next(iter(self.ready)), "",
        ]))
        for tile, exposures in {
            **self.ready, self.outside: ("2243881",),
            self.ignored: ("2243882",),
        }.items():
            path = build_index.exp_list_path(self.run_dir, tile)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("".join(f"{exp}p\n" for exp in exposures))
        for tile in self.ready:
            for stage in ("tile_get_images", "tile_uncompress",
                          "tile_find_exposures"):
                self._manifest(tile, stage)
        self._manifest(self.outside, "tile_vignets")
        # Seed out-of-scope history; the compute parse indexes ready tiles.
        build_index.build(
            [self.outside, self.ignored], self.run_dir, self.index_db,
        )
        machine_defaults = {
            "tile_list": "$base_dir/tiles.txt",
            "retrieve": "symlink",
            "container": str(self.image),
            "inputs": {
                "tiles": "$base_dir/inputs/$run/tiles",
                "exposures": "$base_dir/inputs/$run/exposures",
            },
            "outputs": {
                "run_dir": "$base_dir/scratch/$run",
                "products_dir": "$base_dir/persistent/$run/products",
                "index_db": (
                    "$base_dir/persistent/$run/products/index/run_index.sqlite"
                ),
            },
        }
        self.config = {
            "run": self.name,
            "machine": "candide",
            "input_type": self.input_type,
            "psf_model": self.psf_model,
            "psf_dict": str(self.root / "psf_dict.pickle")
            if self.psf_model == "fake" else "",
            "clean": True,
            "clean_tiles": True,
            "clean_ignore_tiles": [self.ignored],
            "machines": {"candide": {
                "base_dir": str(self.root),
                "data": machine_defaults,
                "image_sims": machine_defaults,
            }},
        }
        self.write_config()

    def _manifest(self, tile, stage):
        path = self.tile_manifest(tile, stage)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}\n")

    def write_config(self):
        """Write the fixture's run configuration."""
        self.config_path.write_text(yaml.safe_dump(self.config))

    def omit_run(self):
        """Omit run: while keeping every path explicit and resolvable."""
        self.config.pop("run")
        for mode in ("data", "image_sims"):
            defaults = self.config["machines"]["candide"][mode]
            # Resolve only fixture templates, not the production config reader.
            text = yaml.safe_dump(defaults)
            text = text.replace("$base_dir", str(self.root))
            text = text.replace("$run", self.name)
            self.config["machines"]["candide"][mode] = yaml.safe_load(text)
        self.write_config()

    def tile_manifest(self, tile, stage):
        """Return an independently specified scratch manifest path."""
        return (self.run_dir / "tiles" / tile[:2] / tile
                / "manifests" / f"{stage}.json")

    def final_cat(self, tile):
        """Return the expected persistent catalogue for one tile."""
        return (self.products_dir / "tiles" / tile[:2] / tile
                / f"final_cat-{tile}.fits")

    def persist_manifest(self, exp):
        """Return the expected persistent manifest for one exposure."""
        return (self.products_dir / "exp" / exp[:2] / exp
                / "manifests" / "exp_persist.json")


@dataclass
class ResolvedDAG:
    """A live API context's workflow, campaign, and graph rendering."""

    workflow: object
    campaign: Campaign
    dot: str

    @property
    def namespace(self):
        """Return the parsed Snakefile's namespace."""
        return self.workflow.globals

    @property
    def jobs(self):
        """Return resolved jobs, including already-satisfied dependencies."""
        return tuple(self.workflow.dag.jobs)

    @property
    def rule_names(self):
        """Return only rules with jobs in this campaign's compute DAG."""
        return {job.rule.name for job in self.jobs}

    @property
    def declared_rule_names(self):
        """Return every parsed rule, including rules without requested jobs."""
        return {rule.name for rule in self.workflow.rules}

    def jobs_for(self, rule):
        """Return jobs of a named rule in stable wildcard order."""
        return sorted(
            (job for job in self.jobs if job.rule.name == rule),
            key=lambda job: sorted(job.wildcards_dict.items()),
        )


def load_profile(name):
    """Read the committed cluster profile without invoking its executor."""
    path = REPO / "profiles" / name / "config.yaml"
    return yaml.safe_load(path.read_text())


@contextmanager
def resolve(campaign, monkeypatch):
    """Resolve ``all`` in an isolated state directory, without running jobs."""
    from snakemake import workflow as sm_workflow

    scripts = REPO / "workflow" / "scripts"
    module_names = {path.stem for path in scripts.glob("*.py")}
    # Snakefile declarations share Snakemake's module-global namespace.
    namespace = sm_workflow.__dict__
    original_namespace = dict(namespace)
    with monkeypatch.context() as patch:
        patch.syspath_prepend(str(scripts))
        for name in module_names:
            patch.delitem(sys.modules, name, raising=False)
        for name in list(os.environ):
            if name.startswith("SP_") or name == "SNAKEMAKE_PROFILE":
                patch.delenv(name)
        for name, value in {
            "SP_PHASE": "compute",
            "SP_PROFILE": "candide",
            "SP_RUN_CONFIG": campaign.config_path,
            "SP_STATE_DIR": campaign.state_dir,
            "SP_CONTAINER": campaign.image,
            "SP_CACHE_DIR": campaign.root / "cache",
            "SP_SANDBOX": campaign.root / "no-sandbox",
            "SP_MISSING_THRESHOLD": "0.34",
            "XDG_CACHE_HOME": campaign.root / "cache",
        }.items():
            patch.setenv(name, str(value))
        profile = load_profile("candide")
        try:
            with SnakemakeApi() as api:
                workflow_api = api.workflow(
                    snakefile=REPO / "workflow" / "Snakefile",
                    workdir=campaign.state_dir,
                    resource_settings=ResourceSettings(
                        # Cluster jobs retain each rule's full thread count.
                        nodes=profile["jobs"],
                        default_resources=profile["default-resources"],
                        overwrite_resources=profile.get("set-resources", {}),
                    ),
                    deployment_settings=DeploymentSettings(
                        deployment_method={DeploymentMethod.APPTAINER},
                        apptainer_args=profile["apptainer-args"],
                    ),
                    workflow_settings=WorkflowSettings(
                        runtime_source_cache_path=(
                            campaign.root / "source-cache"
                        ),
                    ),
                )
                dag_api = workflow_api.dag(DAGSettings(
                    targets={"all"},
                    # Snakemake 9's graph renderer compares the CLI string.
                    print_dag_as="dot",
                    rerun_triggers=RerunTrigger.parse_choices_set(
                        profile["rerun-triggers"],
                    ),
                ))
                stream = io.StringIO()
                with redirect_stdout(stream):
                    dag_api.printdag()
                yield ResolvedDAG(
                    workflow_api._workflow, campaign, stream.getvalue(),
                )
        finally:
            # Bare script imports carry module-level environment constants.
            # Remove this parse's copies before restoring any caller's modules.
            for name in module_names:
                sys.modules.pop(name, None)
            for name in namespace.keys() - original_namespace.keys():
                del namespace[name]
            namespace.update(original_namespace)
