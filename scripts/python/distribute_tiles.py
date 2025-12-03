#!/arc/home/kilbinger/.conda/envs/shapepipe/bin/python3.10
"""
distribute_tiles.py

Uses canfar.helpers.distributed.chunk to automatically distribute tiles
across replicas, then processes each tile assigned to this replica.
"""

import os
import sys
import subprocess
from multiprocessing import Pool
import fcntl
from canfar.helpers.distributed import chunk


def parse_arguments(args):
    """
    Extract the file IDs from command line arguments.

    Args:
        args (list of str): Original command line arguments (e.g., sys.argv[1:])

    Returns:
        str or None: The value passed with -f, or None if not present
    """
    # Default values
    parsed = {"dry_run": 0, "parallel_jobs": 1}

    i = 0
    while i < len(args):
        if args[i] == '-f' and i + 1 < len(args):
            parsed['file_ids'] = args[i + 1]
            i += 2
        elif args[i] == '--batch_num' and i + 1 < len(args):
            parsed['batch_num'] = int(args[i + 1])
            i += 2
        elif args[i] == '--batch_tot' and i + 1 < len(args):
            parsed['batch_tot'] = int(args[i + 1])
            i += 2
        elif args[i] == '--batch_size' and i + 1 < len(args):
            parsed['batch_size'] = int(args[i + 1])
            i += 2
        elif args[i] in ('-n', '--dry_run') and i + 1 < len(args):
            parsed['dry_run'] = int(args[i + 1])
            i += 2
        elif args[i] == '--parallel_jobs' and i + 1 < len(args):
            parsed['parallel_jobs'] = int(args[i + 1])
            i += 2
        elif args[i] == "--debug_out" and i + 1 < len(args):
            parsed["debug_out"] = args[i + 1]
            i += 2
        else:
            i += 1

    return parsed


def get_my_tiles(all_tiles, batch_num, batch_tot, batch_size):
      """
      Calculate which tiles this replica should process based on global position.
      Args:
          all_tiles: List of all tile IDs
          batch_num: Current batch number (1-indexed)
          batch_tot: Total number of batches
          batch_size: Number of replicas per batch
      Returns:
          list: Tiles assigned to this replica

      """
      # Get local replica info from environment
      local_replica_id = int(os.environ.get('REPLICA_ID', 1))

      start_idx = (batch_num - 1) * batch_size
      end_idx = min(batch_num * batch_size, len(all_tiles))
      
      batch_tiles = all_tiles[start_idx:end_idx]
      print(f"Batch {batch_num}/{batch_tot}, local replica {local_replica_id}/{batch_size}")
      print(f"Batch processes tiles {start_idx + 1} to {end_idx} ({len(batch_tiles)} tiles)")
      
      # Use chunk() to distribute this batch's tiles among local replicas
      # chunk() will use REPLICA_ID and REPLICA_COUNT automatically
      my_tiles = list(chunk(batch_tiles))
      
      return my_tiles

def get_tile_list(file_ids):
    """Read all tile IDs from file.
    
    Args:
        file_ids: Path to file containing tile IDs
        
    Returns:
        list: List of tile ID strings
    """
    with open(file_ids, 'r') as f:
        tiles = [line.strip() for line in f if line.strip()]
    return tiles


def build_process_command(tile_id, original_args):
    """
    Build command for a single tile by replacing -f <file_ids> with -e <tile_id>.

    Args:
        tile_id (str): Tile ID to process
        original_args (list of str): Original CLI arguments

    Returns:
        list of str: Command ready for subprocess
    """
    cmd = [f"{os.environ['HOME']}/shapepipe/scripts/sh/init_run_exclusive_canfar.sh"]

    # Transform arguments: replace -f <file_ids> with -e <tile_id>,
    # skip batch arguments, parallel_jobs, and dry_run if not 0, 1
    i = 0
    while i < len(original_args):
        if original_args[i] == '-f' and i + 1 < len(original_args):
            cmd.extend(['-e', tile_id])
            i += 2
        elif original_args[i] in ('--batch_num', '--batch_tot', '--batch_size', '--parallel_jobs'):
            i += 2
        elif original_args[i] in ('-n', '--dry_run'):
            if original_args[i] in (0, 1):
                cmd.extend(["-n", original_args[i]])
            i += 2
        else:
            cmd.append(original_args[i])
            i += 1

    return cmd


def process_single_tile(args_tuple):
    """
    Process a single tile (designed for multiprocessing).

    Args:
        args_tuple: Tuple of (tile_id, tile_num, total_tiles, original_args, dry_run, msg_batch)

    Returns:
        tuple: (tile_id, success, error_message)
    """
    tile_id, tile_num, total_tiles, original_args, dry_run, msg_batch = args_tuple

    print(f"{'='*5} Processing {msg_batch}tile num/total={tile_num}/{total_tiles}: ID={tile_id} {'='*5}")

    # Build command
    cmd = build_process_command(tile_id, original_args)
    print(f"Command: {' '.join(cmd)}")

    # Execute command
    try:
        if dry_run != 2:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"✓ Successfully processed tile {tile_id}")
            return (tile_id, True, None)
        else:
            print(f"dry_run=2")
            return (tile_id, True, None)
    except subprocess.CalledProcessError as e:
        error_msg = f"Failed to process tile {tile_id}: {e}"
        print(f"✗ {error_msg}", file=sys.stderr)
        return (tile_id, False, error_msg)
    except Exception as e:
        error_msg = f"Unexpected error processing tile {tile_id}: {e}"
        print(f"✗ {error_msg}", file=sys.stderr)
        return (tile_id, False, error_msg)


def print_debug(pat, tile_list, out_path, verbose=False):

    local_replica_id = int(os.environ.get('REPLICA_ID', 1))

    with open(out_path, "a") as f:
        # Exclusive lock for save parallel use
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            # Build output string
            output = f"{pat} distribute_tiles REPLICA_ID={local_replica_id}, tiles="
            output += " ".join(tile_list) + "\n"

            # Write to file
            f.write(output)

            # Optionally write to stdout
            if verbose:
                print(output, end="")
        finally:
            # Unlock
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def main():

    # Debug file line pattern
    pat = "- "

    # Parse arguments
    args = parse_arguments(sys.argv[1:])

    if not "file_ids" in args:
        print("Error: -f <file_ids> must be provided", file=sys.stderr)
        sys.exit(1)
    else:
        file_ids = args["file_ids"]
    
    # Read all tiles
    print(f"Reading tile list from {file_ids}")
    all_tiles = get_tile_list(file_ids)
    print(f"Total tiles in file: {len(all_tiles)}")
    
    print(
        "REPLICA_ID, REPLICA_COUNT=",
        os.environ.get('REPLICA_ID'),
        os.environ.get('REPLICA_COUNT')
    )


    # Check if we're in multi-batch mode
    if "batch_num" in args and "batch_tot" in args and "batch_size" in args:
        # Multi-batch mode: calculate global distribution
        my_tile_list = get_my_tiles(
            all_tiles,
            args['batch_num'],
            args['batch_tot'],
            args['batch_size']
        )
        msg_batch = f"batch {args['batch_num']}/{args['batch_tot']} "

    else:
        # Use chunk() to get tiles for this replica
        print("Using chunk() to determine tiles for this replica...")
        my_tile_list = list(chunk(all_tiles))
        msg_batch = ""
    
    if "debug_out" in args:
        print_debug(f"{pat}{msg_batch}", my_tile_list, args["debug_out"], verbose=True)

    print(f"This replica assigned {len(my_tile_list)} tiles")
    print(f"Parallel jobs: {args['parallel_jobs']}")

    if len(my_tile_list) > 0:
        print(f"First tile: {my_tile_list[0]}")
        if len(my_tile_list) > 1:
            print(f"Last tile: {my_tile_list[-1]}")
        else:
            print("Only one tile")

    # Process each tile assigned to this replica
    success_count = 0
    failure_count = 0

    # Prepare arguments for parallel processing
    process_args = [
        (tile_id, i, len(my_tile_list), sys.argv[1:], args["dry_run"], msg_batch)
        for i, tile_id in enumerate(my_tile_list, 1)
    ]

    if args['parallel_jobs'] > 1:
        # Parallel processing
        print(f"Processing tiles in parallel with {args['parallel_jobs']} workers")
        with Pool(processes=args['parallel_jobs']) as pool:
            results = pool.map(process_single_tile, process_args)
    else:
        # Sequential processing (original behavior)
        print("Processing tiles sequentially")
        results = [process_single_tile(arg) for arg in process_args]

    # Count successes and failures
    for tile_id, success, error_msg in results:
        if success:
            success_count += 1
        else:
            failure_count += 1

    # Summary
    print(
        f"{'='*5} Processing completed tot/suc/fail="
        + f"{len(my_tile_list)}/{success_count}/{failure_count}"
    )

    # Exit with error if any failures
    if failure_count > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
