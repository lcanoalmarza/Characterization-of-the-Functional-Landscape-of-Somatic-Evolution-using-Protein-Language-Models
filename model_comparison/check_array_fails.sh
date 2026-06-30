#!/bin/bash

# Usage: ./check_array.sh <input_txt> <log_prefix>
# Example: ./check_array.sh my_seqs.txt logs/ESM2_250_32838420

INPUT_FILE=$1
LOG_PREFIX=$2

if [[ -z "$INPUT_FILE" || -z "$LOG_PREFIX" ]]; then
    echo "Usage: $0 <input_txt> <log_prefix>"
    echo "Example: $0 data.txt logs/my_job_12345"
    exit 1
fi

i=1
fail_count=0

echo "Checking logs for array elements 1 to $(wc -l < "$INPUT_FILE")..."
echo "----------------------------------------------------"

while read -r fasta; do
    # This matches the pattern: logs/name_of_array_X.err
    log_file="${LOG_PREFIX}_${i}.err"

    if [[ -f "$log_file" ]]; then
        # -s checks if the file is NOT empty (contains error text)
        if [[ -s "$log_file" ]]; then
            echo "FAILED | Task ID: $i | FASTA: $fasta"
            echo "       | Log: $log_file"
            ((fail_count++))
        fi
    else
        echo "MISSING| Task ID: $i | FASTA: $fasta (File not found)"
    fi

    ((i++))
done < "$INPUT_FILE"

echo "----------------------------------------------------"
echo "Scan complete. Total failed/non-empty logs: $fail_count"
