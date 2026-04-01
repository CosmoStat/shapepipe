#!/usr/bin/env bash

# Script to remove duplicate tile runs, keeping only the newest
# Run from directory containing P?/tile_runs subdirectories

count=0
total_removed=0

echo "Starting duplicate tile run removal..."
echo "Looking for patches in P?/tile_runs structure..."
echo ""

# Process each patch directory
for patch_dir in P*/tile_runs; do
  if [ ! -d "$patch_dir" ]; then
    continue
  fi

  patch_name=$(dirname "$patch_dir")
  echo "Processing patch: $patch_name"

  # Process each tile directory within the patch
  for dir in "$patch_dir"/*; do
    if [ ! -d "$dir" ]; then
      continue
    fi

    # Check each of the three run types
    for run_type in run_sp_tile_PsViSmVi run_sp_tile_Mc run_sp_tile_Ms; do
      if compgen -G "$dir/output/${run_type}*" > /dev/null; then
        n=$(ls -dt "$dir/output/${run_type}"* 2>/dev/null | wc -l)
        if [ "$n" -gt 1 ]; then
          ((n_remove=n-1))
          echo "  Found $n copies of $run_type in $(basename "$dir"), removing $n_remove oldest"

          # Remove all but the newest
          ls -dt "$dir/output/${run_type}"* | tail -n "$n_remove" | while read old_run; do
            echo "    Removing: $(basename "$old_run")"
            rm -rf "$old_run"
            ((total_removed++))
          done

          ((count++))
        fi
      fi
    done
  done
  echo ""
done

echo "========================================"
echo "Summary:"
echo "  Processed tile directories: $count"
echo "  Total run directories removed: $total_removed"
echo "========================================"
