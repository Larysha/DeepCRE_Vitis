#!/bin/bash

# BLAMM motif mapping pipeline for moca_blue
# Maps learned motifs (PFMs) back to extracted gene sequences
# 
# REQUIREMENTS: BLAMM must be installed and available in PATH (system-wide installation)
# Test with: blamm --help
# 
# WHY PFMs: BLAMM requires Position Frequency Matrices (PFMs) - raw count matrices with 
# positive integers. BLAMM internally converts these to PWMs for scoring and generates
# proper p-value thresholds. CWMs/PWMs with negative log-odds values are NOT suitable.
#
# BLAMM WORKFLOW: Uses manifest files (.mf) that point to FASTA files
# Adapted for moca_blue pipeline
# Usage: ./blamm_motif_mapping.sh [species_model] [p_value_threshold] [flank_size]

set -euo pipefail  

# Parameters
SPECIES_MODEL="${1:-vitis_vinifera_PN40024}"
P_THRESHOLD="${2:-0.0001}"
FLANK_SIZE="${3:-1000}"

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

INPUT_SEQ_DIRS=("results/moca_blue/ref_seq" "vitis_cre/src/results/moca_blue/ref_seq")
MOTIFS_DIRS=("results/moca_blue/mo_nom" "vitis_cre/src/results/moca_blue/mo_nom")
OUTPUT_DIR="vitis_cre/src/results/moca_blue/mo_proj"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find input files
SEQUENCE_FILE=""
if ! SEQUENCE_FILE=$(find_file "${SPECIES_MODEL}_${FLANK_SIZE}bp_flanks.fa" "${INPUT_SEQ_DIRS[@]}"); then
    echo "Error: Could not find sequence file ${SPECIES_MODEL}_${FLANK_SIZE}bp_flanks.fa in directories: ${INPUT_SEQ_DIRS[*]}" >&2
    exit 1
fi

MOTIFS_FILE=""
if ! MOTIFS_FILE=$(find_file "rdf5_PFM_patternvitisssr_pfm-motifs.jaspar" "${MOTIFS_DIRS[@]}"); then
    echo "Error: Could not find motifs file rdf5_PFM_patternvitisssr_pfm-motifs.jaspar in directories: ${MOTIFS_DIRS[*]}" >&2
    exit 1
fi

PROJECT_ID="${SPECIES_MODEL}_${FLANK_SIZE}bp_pt${P_THRESHOLD}"

# Validate inputs
echo "Validating input files..."
echo "Found sequence file: $SEQUENCE_FILE"
echo "Found motifs file (PFM format): $MOTIFS_FILE"

# Check BLAMM is available in PATH
if ! command -v blamm &> /dev/null; then
    echo "Error: BLAMM not found in PATH" >&2
    echo "Please ensure BLAMM is installed system-wide and accessible via 'blamm' command" >&2
    echo "Test with: blamm --help" >&2
    exit 1
fi

echo "Starting BLAMM motif mapping pipeline..."
echo "Sequence file: $SEQUENCE_FILE"
echo "Motifs file (PFM format): $MOTIFS_FILE"
echo "P-value threshold: $P_THRESHOLD"
echo "Project ID: $PROJECT_ID"
echo ""

start_time=$(date +%s)

echo "Step 0: Creating BLAMM manifest file..."
MANIFEST_FILE="sequences.mf"

# Create manifest file pointing to FASTA
# Format: group_name<TAB>path_to_fasta
# Convert relative path to absolute path to avoid issues
ABSOLUTE_SEQUENCE_FILE=$(realpath "$SEQUENCE_FILE")
echo -e "gene_sequences\t$ABSOLUTE_SEQUENCE_FILE" > "$MANIFEST_FILE"

echo "Created manifest file: $MANIFEST_FILE"
echo "Manifest contents:"
cat "$MANIFEST_FILE"
echo ""

echo "Step 1: Creating sequence dictionary..."
blamm dict "$MANIFEST_FILE"

echo "Step 2: Generating empirical score distributions from PFMs..."
echo "BLAMM will convert PFMs to PWMs internally and generate proper thresholds"
blamm hist -e "$MOTIFS_FILE" "$MANIFEST_FILE"

echo "Step 3: Scanning sequences for motif occurrences..."
echo "Using reverse complement scanning (-rc) with p-value threshold $P_THRESHOLD"
blamm scan -rc -pt "$P_THRESHOLD" "$MOTIFS_FILE" "$MANIFEST_FILE"

echo "Step 4: Organizing outputs..."
PROJECT_OUTPUT_DIR="${OUTPUT_DIR}/${PROJECT_ID}"
mkdir -p "$PROJECT_OUTPUT_DIR"

# Move BLAMM outputs to project directory
if [[ -f "occurrences.txt" ]]; then
    mv occurrences.txt "$PROJECT_OUTPUT_DIR/"
    echo "Motif occurrences saved to: ${PROJECT_OUTPUT_DIR}/occurrences.txt"
else
    echo "Warning: No occurrences.txt file generated" >&2
fi

if [[ -f "PWMthresholds.txt" ]]; then
    mv PWMthresholds.txt "$PROJECT_OUTPUT_DIR/"
fi

# Move histogram files
for hist_file in hist_*; do
    if [[ -f "$hist_file" ]]; then
        mv "$hist_file" "$PROJECT_OUTPUT_DIR/"
    fi
done

# Move manifest and dictionary files for record keeping
if [[ -f "$MANIFEST_FILE" ]]; then
    mv "$MANIFEST_FILE" "$PROJECT_OUTPUT_DIR/"
fi

if [[ -f "${MANIFEST_FILE}.dict" ]]; then
    mv "${MANIFEST_FILE}.dict" "$PROJECT_OUTPUT_DIR/"
fi

# Generate summary statistics
if [[ -f "${PROJECT_OUTPUT_DIR}/occurrences.txt" ]]; then
    total_occurrences=$(wc -l < "${PROJECT_OUTPUT_DIR}/occurrences.txt")
    unique_motifs=$(cut -f3 "${PROJECT_OUTPUT_DIR}/occurrences.txt" | sort -u | wc -l)
    unique_sequences=$(cut -f1 "${PROJECT_OUTPUT_DIR}/occurrences.txt" | sort -u | wc -l)
    
    echo ""
    echo "=== MAPPING SUMMARY ==="
    echo "Total motif occurrences: $total_occurrences"
    echo "Unique motifs found: $unique_motifs"
    echo "Sequences with motifs: $unique_sequences"
    echo "Results directory: $PROJECT_OUTPUT_DIR"
    echo ""
    echo "Output format (GTF-like):"
    echo "sequence_id tool motif_id start end score strand . ."
    echo ""
    echo "First few results:"
    head -5 "${PROJECT_OUTPUT_DIR}/occurrences.txt" 2>/dev/null || echo "No results to display"
fi

# Calculate runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Pipeline completed in ${runtime} seconds"

echo ""
echo "=== NEXT STEPS ==="
echo "1. Coordinates in occurrences.txt are relative to extracted sequences (not genome)"
echo "2. To get genomic coordinates, add back the flanking offset (${FLANK_SIZE}bp)"
echo "3. Use header metadata in FASTA file to map back to original gene coordinates"
echo "4. Consider filtering results by score/significance before downstream analysis"
echo "5. Results are in GTF format - compatible with many genomics tools"
echo "6. Run R scripts from mo_proj folder to apply biological filters and analysis"
echo ""
echo "=== FILES CREATED ==="
echo "- ${PROJECT_OUTPUT_DIR}/occurrences.txt (main results - ready for R filtering)"
echo "- ${PROJECT_OUTPUT_DIR}/PWMthresholds.txt (score thresholds used)"
echo "- ${PROJECT_OUTPUT_DIR}/sequences.mf (manifest file)"
echo "- ${PROJECT_OUTPUT_DIR}/sequences.mf.dict (sequence dictionary)"
echo "- ${PROJECT_OUTPUT_DIR}/hist_* (score distribution histograms)"