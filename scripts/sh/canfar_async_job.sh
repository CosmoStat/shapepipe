#!/bin/bash

source $HOME/shapepipe/scripts/sh/functions.sh


# 1. Extract the line from -f file
# Save original arguments
ARGS=("$@")

# Find -f argument and get its value; remove -f and -S <anything>
FILE_IDS=""
NEW_ARGS=()
i=0
while [[ $i -lt $# ]]; do
    arg="${ARGS[i]}"
    case "$arg" in
        -f)
            ((i++))
            FILE_IDS="${ARGS[i]}"
            ;;
        -S)
            ((i++))  # skip next argument (value of -S)
            ;;
        *)
            NEW_ARGS+=("$arg")
            ;;
    esac
    ((i++))
done

if [[ -z "$FILE_IDS" ]]; then
    echo "Error: -f <file_IDs> must be provided"
    exit 1
fi

# Get the line for this replica
LINE=$(get_line_from_file -f "$FILE_IDS")

# 2. Call init_run.sh with remaining args + -e line
$HOME/shapepipe/scripts/sh/init_run_exclusive_canfar.sh "${NEW_ARGS[@]}" -e "$LINE"
