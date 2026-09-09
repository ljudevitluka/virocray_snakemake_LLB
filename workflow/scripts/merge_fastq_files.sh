#!/bin/bash

# Set the base directory (current directory by default)
base_dir="${1:-.}"

# Traverse subdirectories
for subdir in "$base_dir"/*/ ; do
    echo "Checking $subdir"
    
    # Find matching _1.fq.gz and _2.fq.gz files
    fq1_files=("$subdir"*"_1.fq.gz")
    fq2_files=("$subdir"*"_2.fq.gz")

    # Count matching files
    num_fq1=${#fq1_files[@]}
    num_fq2=${#fq2_files[@]}

    # Merge only if there are multiple matching files
    if [[ "$num_fq1" -gt 1 && "$num_fq2" -gt 1 ]]; then
        echo "Merging files in $subdir"

        # Define output files
        merged1="${subdir%/}/merged_1.fq.gz"
        merged2="${subdir%/}/merged_2.fq.gz"

        # Merge using gzip stream (to avoid decompressing to disk)
        cat "${fq1_files[@]}" > "$merged1"
        cat "${fq2_files[@]}" > "$merged2"

        echo " -> Created $merged1 and $merged2"
    else
        echo " -> Skipping: only one or no paired files"
    fi
done
