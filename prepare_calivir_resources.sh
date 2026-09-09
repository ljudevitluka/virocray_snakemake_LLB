#!/usr/bin/env bash
set -euo pipefail

SOURCE_DIR="/scratch/lukabostjancic/ViroCray/WP1/CaliVir/01.RawData"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST_DIR="$PROJECT_DIR/resources"

dry_run="${DRY_RUN:-0}"
declare -A seen

while IFS= read -r -d '' source_file; do
    filename="$(basename "$source_file")"
    sample="$(basename "$(dirname "$source_file")")"

    case "$filename" in
        *_1.fq.gz) read_number=1 ;;
        *_2.fq.gz) read_number=2 ;;
        *) continue ;;
    esac

    destination="$DEST_DIR/${sample}_${read_number}.fq.gz"

    if [[ -n "${seen[$destination]+x}" ]]; then
        printf 'ERROR: multiple source files map to %s\n' "$destination" >&2
        exit 1
    fi
    seen["$destination"]="$source_file"

    if [[ "$dry_run" == 1 ]]; then
        printf 'would copy: %s -> %s\n' "$source_file" "$destination"
        continue
    fi

    mkdir -p "$DEST_DIR"
    if [[ -e "$destination" ]]; then
        printf 'ERROR: destination already exists: %s\n' "$destination" >&2
        exit 1
    fi
    cp -- "$source_file" "$destination"
    printf 'copied: %s -> %s\n' "$source_file" "$destination"
done < <(find "$SOURCE_DIR" -type f -name '*.fq.gz' -print0 | sort -z)

for sample in FI_GILL FI_GUT FI_HEP CS SS WS; do
    for read_number in 1 2; do
        expected="$DEST_DIR/${sample}_${read_number}.fq.gz"
        if [[ "$dry_run" != 1 && ! -f "$expected" ]]; then
            printf 'ERROR: missing expected file: %s\n' "$expected" >&2
            exit 1
        fi
    done
done

printf 'Resource preparation complete.\n'