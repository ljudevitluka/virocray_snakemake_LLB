#!/bin/bash

# Get current working directory
parent_dir=$(pwd)

echo "Running in: $parent_dir"

# Loop over subdirectories (assumed to be sample IDs)
for sample_dir in "$parent_dir"/*; do
    if [[ -d "$sample_dir" ]]; then
        sample_id=$(basename "$sample_dir")
        echo "Processing sample: $sample_id"

        # Find matching _1 and _2 files (first match only)
        fq1_file=("$sample_dir"/*_1.fq.gz)
        fq2_file=("$sample_dir"/*_2.fq.gz)

        # New filenames
        new_fq1="$sample_dir/${sample_id}_1.fq.gz"
        new_fq2="$sample_dir/${sample_id}_2.fq.gz"

        # Check existence and avoid overwriting
        if [[ -f "${fq1_file[0]}" && -f "${fq2_file[0]}" ]]; then
            echo "  Found:"
            echo "    ${fq1_file[0]}"
            echo "    ${fq2_file[0]}"

            # Rename if target names don't exist already
            if [[ ! -f "$new_fq1" && ! -f "$new_fq2" ]]; then
                mv "${fq1_file[0]}" "$new_fq1"
                mv "${fq2_file[0]}" "$new_fq2"
                echo "  Renamed to:"
                echo "    $new_fq1"
                echo "    $new_fq2"
            else
                echo "  Skipping: Target filenames already exist"
            fi
        else
            echo "  Skipping: _1 or _2 file not found in $sample_dir"
        fi
    fi
done
