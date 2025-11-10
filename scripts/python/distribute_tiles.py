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
    if '-f' in args:
        index = args.index('-f')
        if index + 1 < len(args):
            return args[index + 1]
    return None


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
    
    # Transform arguments: replace -f <file_ids> with -e <tile_id>
    i = 0
    while i < len(original_args):
        if original_args[i] == '-f' and i + 1 < len(original_args):
            cmd.extend(['-e', tile_id])
            i += 2  # skip -f and its value
        else:
            cmd.append(original_args[i])
            i += 1
    
    return cmd


def main():
    # Parse arguments
    #args = parse_arguments(sys.argv[1:])
    file_ids = args = parse_arguments(sys.argv[1:])

    #if not args['file_ids']:
    if not file_ids:
        print("Error: -f <file_ids> must be provided", file=sys.stderr)
        sys.exit(1)
    
    # Read all tiles
    print(f"Reading tile list from {file_ids}")
    all_tiles = get_tile_list(file_ids)
    print(f"Total tiles in file: {len(all_tiles)}")
    
    # Use chunk() to get tiles for this replica
    print("Using chunk() to determine tiles for this replica...")
    my_tile_list = list(chunk(all_tiles))
    print(f"This replica assigned {len(my_tile_list)} tiles")
    
    if len(my_tile_list) > 0:
        print(f"First tile: {my_tile_list[0]}")
        if len(my_tile_list) > 1:
            print(f"Last tile: {my_tile_list[-1]}")
    
    # Process each tile assigned to this replica
    success_count = 0
    failure_count = 0
    
    for i, tile_id in enumerate(my_tile_list, 1):
        print(f"\n{'='*60}")
        print(f"Processing tile {i}/{len(my_tile_list)}: {tile_id}")
        print(f"{'='*60}")
        
        # Build command
        cmd = build_process_command(tile_id, sys.argv[1:])
        print(f"Command: {' '.join(cmd)}")
        
        # Execute command
        try:
            result = subprocess.run(cmd, check=True)
            print(f"✓ Successfully processed tile {tile_id}")
            success_count += 1
        except subprocess.CalledProcessError as e:
            print(f"✗ Failed to process tile {tile_id}: {e}", file=sys.stderr)
            failure_count += 1
            # Continue processing other tiles even if one fails
        except Exception as e:
            print(f"✗ Unexpected error processing tile {tile_id}: {e}", file=sys.stderr)
            failure_count += 1
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Processing complete")
    print(f"{'='*60}")
    print(f"Total tiles assigned to this replica: {len(my_tile_list)}")
    print(f"Successfully processed: {success_count}")
    print(f"Failed: {failure_count}")
    
    # Exit with error if any failures
    if failure_count > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()