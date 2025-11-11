#!/arc/home/kilbinger/.conda/envs/shapepipe/bin/python3.10
"""
distribute_tiles.py

Uses canfar.helpers.distributed.chunk to automatically distribute tiles
across replicas, then processes each tile assigned to this replica.
"""

import os
import sys
import subprocess
from canfar.helpers.distributed import chunk


def parse_arguments(args):
    """
    Extract the file IDs from command line arguments.

    Args:
        args (list of str): Original command line arguments (e.g., sys.argv[1:])

    Returns:
        str or None: The value passed with -f, or None if not present
    """
    parsed = {"dry_run": 0}

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
      print(f"Batch {batch_num}/{batch_tot}, Local replica {local_replica_id}/{batch_size}")
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
    # skip batch arguments and dry_run if not 0, 1
    i = 0
    while i < len(original_args):
        if original_args[i] == '-f' and i + 1 < len(original_args):
            cmd.extend(['-e', tile_id])
            i += 2
        elif original_args[i] in ('--batch_num', '--batch_tot', '--batch_size'):
            i += 2
        elif original_args[i] in ('-n', '--dry_run'):
            if original_args[i] in (0, 1):
                cmd.extend(["-n", original_args[i]])
            i += 2
        else:
            cmd.append(original_args[i])
            i += 1 

    return cmd


def main():
    # Parse arguments
    args = parse_arguments(sys.argv[1:])

    print("MKDEBUG", args["dry_run"])

    if not "file_ids" in args:
        print("Error: -f <file_ids> must be provided", file=sys.stderr)
        sys.exit(1)
    else:
        file_ids = args["file_ids"]
    
    # Read all tiles
    print(f"Reading tile list from {file_ids}")
    all_tiles = get_tile_list(file_ids)
    print(f"Total tiles in file: {len(all_tiles)}")
    
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
    
    print(f"This replica assigned {len(my_tile_list)} tiles")
    
    if len(my_tile_list) > 0:
        print(f"First tile: {my_tile_list[0]}")
        if len(my_tile_list) > 1:
            print(f"Last tile: {my_tile_list[-1]}")
    
    # Process each tile assigned to this replica
    success_count = 0
    failure_count = 0
    
    for i, tile_id in enumerate(my_tile_list, 1):
        print(f"{'='*5} Processing {msg_batch}tile {i}/{len(my_tile_list)}: {tile_id} {'='*5}")
        
        # Build command
        cmd = build_process_command(tile_id, sys.argv[1:])
        print(f"Command: {' '.join(cmd)}")
        
        # Execute command
        try:
            if args["dry_run"] != 2:
                result = subprocess.run(cmd, check=True)
                print(f"✓ Successfully processed tile {tile_id}")
            else:
                print(f"Not running {cmd} (dry_run=2)")
            success_count += 1
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to process tile {tile_id}: {e}", file=sys.stderr)
            failure_count += 1
            # Continue processing other tiles even if one fails
        except Exception as e:
            print(f"✗ Unexpected error processing tile {tile_id}: {e}", file=sys.stderr)
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