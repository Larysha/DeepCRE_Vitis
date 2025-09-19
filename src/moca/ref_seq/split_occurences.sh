#!/bin/bash

# Split large BLAMM occurrences.txt files into manageable chunks
# Preserves original file and creates numbered chunks for parallel processing
# Usage: ./split_occurrences.sh [project_id] [lines_per_chunk] [input_file]

set -euo pipefail

# Default parameters
DEFAULT_PROJECT="vitis_stress_02"
DEFAULT_CHUNK_SIZE="1000000"
DEFAULT_INPUT="../../../out/moca_results/mo_proj/genomic_vitis_vinifera_PN40024_1000bp_pt0.0001/occurrences.txt"

# Command line arguments with defaults
PROJECT_ID="${1:-$DEFAULT_PROJECT}"
CHUNK_SIZE="${2:-$DEFAULT_CHUNK_SIZE}"
INPUT_FILE="${3:-$DEFAULT_INPUT}"

# Validate chunk size is numeric
if ! [[ "$CHUNK_SIZE" =~ ^[0-9]+$ ]]; then
    echo "Error: Chunk size must be a positive integer, got: $CHUNK_SIZE" >&2
    exit 1
fi

# Start timing and resource monitoring
start_time=$(date +%s)
start_resources=$(ps -o pid,%cpu,%mem,vsz,rss,tty,stat,start_time --no-headers $$)

echo "=== BLAMM Occurrences File Splitter ==="
echo "Project ID: $PROJECT_ID"
echo "Input file: $INPUT_FILE"
echo "Lines per chunk: $CHUNK_SIZE"
echo ""

# Validate input file exists and is readable
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found" >&2
    echo "Expected to find BLAMM output file in current directory" >&2
    exit 1
fi

if [[ ! -r "$INPUT_FILE" ]]; then
    echo "Error: Cannot read input file '$INPUT_FILE'" >&2
    exit 1
fi

# Get file stats before processing
file_size=$(stat -c%s "$INPUT_FILE" 2>/dev/null || stat -f%z "$INPUT_FILE" 2>/dev/null || echo "unknown")
total_lines=$(wc -l < "$INPUT_FILE")

echo "Input file stats:"
echo "  Size: $file_size bytes"
echo "  Total lines: $total_lines"
echo "  Estimated chunks: $((total_lines / CHUNK_SIZE + 1))"
echo ""

# Create output directory structure in same location as input file
INPUT_DIR=$(dirname "$INPUT_FILE")
OUTPUT_DIR="${INPUT_DIR}/occ${PROJECT_ID}"
mkdir -p "$OUTPUT_DIR"

echo "Creating chunks in directory: $OUTPUT_DIR"

# Split the file with descriptive naming
# Using numeric suffixes starting from 001 for better sorting
split --lines="$CHUNK_SIZE" \
     --numeric-suffixes=1 \
     --suffix-length=3 \
     --verbose \
     "$INPUT_FILE" \
     "${OUTPUT_DIR}/chunk_"

# Verify split operation
chunk_count=$(find "$OUTPUT_DIR" -name "chunk_*" -type f | wc -l)

if [[ $chunk_count -eq 0 ]]; then
    echo "Error: No chunks were created" >&2
    exit 1
fi

echo ""
echo "=== SPLIT SUMMARY ==="
echo "Created $chunk_count chunk files:"

# Show chunk file details
for chunk_file in "${OUTPUT_DIR}"/chunk_*; do
    if [[ -f "$chunk_file" ]]; then
        chunk_lines=$(wc -l < "$chunk_file")
        chunk_name=$(basename "$chunk_file")
        printf "  %-12s: %8d lines\n" "$chunk_name" "$chunk_lines"
    fi
done

# Verify total lines are preserved
total_chunk_lines=0
for chunk_file in "${OUTPUT_DIR}"/chunk_*; do
    if [[ -f "$chunk_file" ]]; then
        chunk_lines=$(wc -l < "$chunk_file")
        total_chunk_lines=$((total_chunk_lines + chunk_lines))
    fi
done

echo ""
echo "Verification:"
echo "  Original lines: $total_lines"
echo "  Chunk lines:    $total_chunk_lines"

if [[ $total_lines -eq $total_chunk_lines ]]; then
    echo "  Line count verified"
else
    echo "  Line count mismatch!" >&2
    exit 1
fi

# Calculate runtime and resource usage
end_time=$(date +%s)
end_resources=$(ps -o pid,%cpu,%mem,vsz,rss,tty,stat,start_time --no-headers $$)
runtime=$((end_time - start_time))

echo ""
echo "=== PERFORMANCE SUMMARY ==="
echo "Runtime: ${runtime} seconds"
echo "Throughput: $((total_lines / runtime)) lines/second"
echo ""
echo "Resource usage:"
echo "$start_resources" | awk '{print "  Start: " $0}'
echo "$end_resources" | awk '{print "  End:   " $0}'

echo ""
echo "=== NEXT STEPS ==="
echo "1. Process chunks independently for parallel analysis"
echo "2. Use chunk files for distributed computing if needed"
echo "3. Original file preserved at: $INPUT_FILE"
echo "4. Chunk directory: $OUTPUT_DIR"
echo ""
echo "Example parallel processing:"
echo "  for chunk in ${OUTPUT_DIR}/chunk_*; do"
echo "    process_motifs.sh \"\$chunk\" &"
echo "  done"
echo "  wait"