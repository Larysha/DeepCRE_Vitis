#!/bin/bash

# Extract gene sequences with flanking regions and rich annotation headers
# Handles both GFF3 and GTF formats
# Creates FASTA headers with gene IDs, GENOMIC coordinates, and strand information
# Adapted for moca pipeline
# Usage: ./extract_range_to_fasta_genomic.sh [species_model_name] [flank_size]

set -euo pipefail 

# Default parameters
SPECIES_MODEL="${1:-vitis_vinifera_PN40024}"
FLANK_SIZE="${2:-1000}"

# File paths - search in standard locations
GENOME_DIRS=(
    "../../../vitis_data/genome"
    "../../vitis_data/genome"
    "/home/rish/phd_2025/deepcre_vitis/vitis_cre/vitis_data/genome"
)
ANNOTATION_DIRS=(
    "../../../vitis_data/gene_models"
    "../../vitis_data/gene_models"
    "/home/rish/phd_2025/deepcre_vitis/vitis_cre/vitis_data/gene_models"
)

# Function to find file in multiple directories
find_file() {
    local filename="$1"
    shift
    local dirs=("$@")
    
    for dir in "${dirs[@]}"; do
        if [[ -f "$dir/$filename" ]]; then
            echo "$dir/$filename"
            return 0
        fi
    done
    return 1
}

# Find genome and annotation files
echo "Searching for input files..."

GENOME_FILE="${SPECIES_MODEL}_dna.fa"
if ! GENOME_PATH=$(find_file "$GENOME_FILE" "${GENOME_DIRS[@]}"); then
    echo "Error: Could not find genome file $GENOME_FILE in directories: ${GENOME_DIRS[*]}" >&2
    exit 1
fi

# Try both GFF3 and GTF extensions  
ANNOTATION_FILE=""
for ext in "gff3" "gtf"; do
    test_file="${SPECIES_MODEL}.${ext}"
    if ANNOTATION_PATH=$(find_file "$test_file" "${ANNOTATION_DIRS[@]}"); then
        ANNOTATION_FILE="$test_file"
        break
    fi
done

if [[ -z "$ANNOTATION_FILE" ]]; then
    echo "Error: Could not find annotation file ${SPECIES_MODEL}.gff3 or ${SPECIES_MODEL}.gtf in directories: ${ANNOTATION_DIRS[*]}" >&2
    exit 1
fi

echo "Found genome: $GENOME_PATH"
echo "Found annotation: $ANNOTATION_PATH"
echo "Using flank size: ${FLANK_SIZE}bp"


OUTPUT_DIR="../../../out/moca_results/ref_seq"
mkdir -p "$OUTPUT_DIR"

RANGES_FILE="genomic_coords_gene_ranges_${FLANK_SIZE}bp.txt"
OUTPUT_FASTA="${OUTPUT_DIR}/genomic_coords_${SPECIES_MODEL}_${FLANK_SIZE}bp_flanks.fa" 

# Main task starts here
# Start timer
start_time=$(date +%s)
echo "Extracting gene information with annotations..."

# Extract gene information with rich metadata
# Creates both ranges file and metadata for header construction
awk -v flank="$FLANK_SIZE" '
BEGIN { OFS="\t" }
# Skip comments and headers
/^#/ { next }
# Process gene features (column 3 contains "gene")
$3 ~ /gene/ {
    chr = $1
    gene_start = $4
    gene_end = $5
    strand = $7
    attributes = $9
    
    # Calculate flanked coordinates
    flank_start = gene_start - flank
    flank_end = gene_end + flank
    
    # Bounds checking - dont go below 1
    if (flank_start < 1) flank_start = 1
    
    # Get gene ID from attributes (works for both GFF3 and GTF)
    gene_id = ""
    if (match(attributes, /ID=([^;]+)/, arr)) {
        gene_id = arr[1]
    } else if (match(attributes, /gene_id[= ]"?([^;"]+)"?/, arr)) {
        gene_id = arr[1]
    } else {
        # Fallback: create ID from coordinates
        gene_id = chr "_" gene_start "_" gene_end
    }
    
    # Create range for samtools
    range = chr ":" flank_start "-" flank_end

    # Header information
    # Format: >GeneID Chr Strand GenomicStart-GenomicEnd
    header = gene_id "\t" chr "\t" strand "\t" flank_start "-" flank_end
    
    # Print range for samtools and header info for post-processing
    print range "\t" header
}' "$ANNOTATION_PATH" > "$RANGES_FILE"

# Split ranges and headers for processing
cut -f1 "$RANGES_FILE" > "${RANGES_FILE}.ranges"
cut -f2- "$RANGES_FILE" > "${RANGES_FILE}.headers"

# Count extracted ranges
gene_count=$(wc -l < "${RANGES_FILE}.ranges")
echo "Extracted $gene_count gene ranges with annotations"

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is required but not found in PATH" >&2
    exit 1
fi

echo "Extracting sequences with samtools..."

# Extract sequences using samtools
samtools faidx "$GENOME_PATH" \
    --region-file "${RANGES_FILE}.ranges" \
    --output "${OUTPUT_FASTA}.tmp"

# Post-process FASTA to add rich headers (switching over to Python for flexibility)
echo "Adding gene annotation to FASTA headers..."
python3 - << EOF
import sys

# Read headers file
with open("${RANGES_FILE}.headers", 'r') as f:
    headers = [line.strip().split('\t') for line in f]

# Read and process FASTA
with open("${OUTPUT_FASTA}.tmp", 'r') as infile, open("${OUTPUT_FASTA}", 'w') as outfile:
    header_idx = 0
    for line in infile:
        if line.startswith('>'):
            # Replace header with rich annotation
            gene_id, chr_name, strand, coords = headers[header_idx]
            new_header = f">{gene_id}\t{chr_name}\t{strand}\t{coords}\n"
            outfile.write(new_header)
            header_idx += 1
        else:
            outfile.write(line)
EOF

# Verify output
if [[ ! -f "$OUTPUT_FASTA" ]]; then
    echo "Error: Failed to create output FASTA file" >&2
    exit 1
fi

seq_count=$(grep -c "^>" "$OUTPUT_FASTA" || true)
echo "Successfully extracted $seq_count sequences with rich headers"

echo "Example header format:"
head -n 1 "$OUTPUT_FASTA"

# Runtime summary
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Completed in ${runtime} seconds"

# Cleanup intermediate files
rm "$RANGES_FILE" "${RANGES_FILE}.ranges" "${RANGES_FILE}.headers" "${OUTPUT_FASTA}.tmp"

echo "Output written to: $OUTPUT_FASTA"
echo ""
echo "FASTA header format: >GeneID Chr Strand GenomicCoordinates"
echo "This format provides full traceability to genome positions"