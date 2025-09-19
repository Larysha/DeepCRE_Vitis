#!/bin/bash

# BLAMM motif mapping pipeline for genomic coordinate FASTA files
# Maps learned motifs (PFMs) back to gene sequences with rich annotation
# Creates BLAMM manifest file and handles coordinate mapping
# 
# REQUIREMENTS: BLAMM must be installed and available in PATH (system-wide installation)
# Test with: blamm --help
# 
# WHY PFMs: BLAMM requires Position Frequency Matrices (PFMs) - raw count matrices with 
# positive integers. BLAMM internally converts these to PWMs for scoring and generates
# proper p-value thresholds. CWMs/PWMs with negative log-odds values are NOT suitable.
#
# BLAMM WORKFLOW: Uses manifest files (.mf) that point to FASTA files
# a P value of 0.0001 = 1 in 10,000 chance that a motif hit occurred by random chance
# this is quite stringent - only keeping high confidence motif hits
#
# The p-value is calculated against the background nucleotide 
# composition BLAMM detected: 32.2% A/T, 17.8% C/G—typical for plant genomes.
#
# Adapted for moca pipeline
# Usage: ./blamm_motif_mapping_genomic.sh [species_model] [p_value_threshold] [flank_size]

set -euo pipefail  # Exit on error, undefined vars, pipe failures

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

# Input/Output directories
INPUT_SEQ_DIRS=("../../../out/moca_results/ref_seq")
MOTIFS_DIRS=("../../../out/moca_results/mo_nom")
OUTPUT_DIR="../../../out/moca_results/mo_proj"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find input files - use genomic coordinate version
SEQUENCE_FILE=""
if ! SEQUENCE_FILE=$(find_file "genomic_coords_${SPECIES_MODEL}_${FLANK_SIZE}bp_flanks.fa" "${INPUT_SEQ_DIRS[@]}"); then
    echo "Error: Could not find genomic coordinate sequence file genomic_coords_${SPECIES_MODEL}_${FLANK_SIZE}bp_flanks.fa" >&2
    echo "Run extract_range_to_fasta_genomic.sh first to generate annotated sequences" >&2
    exit 1
fi

# Find PFM motifs file
MOTIFS_FILE=""
if ! MOTIFS_FILE=$(find_file "rdf5_PFM_patternvitisssr_pfm-motifs.jaspar" "${MOTIFS_DIRS[@]}"); then
    echo "Error: Could not find PFM motifs file in directories: ${MOTIFS_DIRS[*]}" >&2
    echo "Expected file: rdf5_PFM_patternvitisssr_pfm-motifs.jaspar" >&2
    echo "Run mo_nom/get_matrices.R to generate PFM motifs" >&2
    exit 1
fi

PROJECT_ID="genomic_${SPECIES_MODEL}_${FLANK_SIZE}bp_pt${P_THRESHOLD}"

# Validate inputs
echo "Validating input files..."
echo "Found sequence file: $SEQUENCE_FILE"
echo "Found PFM motifs file: $MOTIFS_FILE"

# Check BLAMM is available in PATH
if ! command -v blamm &> /dev/null; then
    echo "Error: BLAMM not found in PATH" >&2
    echo "Please ensure BLAMM is installed system-wide and accessible via 'blamm' command" >&2
    echo "Test with: blamm --help" >&2
    exit 1
fi

echo "Starting BLAMM motif mapping pipeline with genomic coordinates..."
echo "Sequence file: $SEQUENCE_FILE"
echo "PFM motifs file: $MOTIFS_FILE"
echo "P-value threshold: $P_THRESHOLD"
echo "Project ID: $PROJECT_ID"
echo ""

# Start timing
start_time=$(date +%s)

# Step 0: Create BLAMM manifest file (sequences.mf)
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

# Step 1: Create sequence dictionary
echo "Step 1: Creating sequence dictionary..."
blamm dict "$MANIFEST_FILE"

# Step 2: Generate empirical PWM score distributions
echo "Step 2: Generating empirical score distributions from PFMs..."
echo "BLAMM will convert PFMs to PWMs internally and generate proper thresholds"
blamm hist -e "$MOTIFS_FILE" "$MANIFEST_FILE"

# Step 3: Scan sequences for motif occurrences
echo "Step 3: Scanning sequences for motif occurrences..."
echo "Using reverse complement scanning (-rc) with p-value threshold $P_THRESHOLD"
blamm scan -rc -pt "$P_THRESHOLD" "$MOTIFS_FILE" "$MANIFEST_FILE"

# Step 4: Organize outputs
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
    
    # Show example occurrences with genomic coordinates
    echo ""
    echo "Output format (GTF-like):"
    echo "sequence_id tool motif_id start end score strand . ."
    echo ""
    echo "Example occurrences (first 5 lines):"
    head -n 5 "${PROJECT_OUTPUT_DIR}/occurrences.txt" 2>/dev/null || echo "No results to display"
fi

# Calculate runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo ""
echo "Pipeline completed in ${runtime} seconds"

echo ""
echo "=== GENOMIC COORDINATE ADVANTAGES ==="
echo "✓ Occurrences include gene IDs, chromosomes, and strand information"
echo "✓ Coordinates are relative to extracted sequences (add to sequence start for genome position)"
echo "✓ Rich FASTA headers enable direct traceability to genome features"
echo "✓ Results ready for gene annotation and functional analysis"
echo "✓ Distance-based filtering can use gene features instead of arbitrary thresholds"
echo ""
echo "=== NEXT STEPS ==="
echo "1. Analyze ${PROJECT_OUTPUT_DIR}/occurrences.txt for motif-gene associations"
echo "2. Convert relative coordinates to absolute genomic positions if needed"
echo "3. Apply distance-based filtering using gene features (promoter/terminator regions)"
echo "4. Perform GO enrichment analysis on genes containing specific motifs"
echo "5. Compare motif occurrence patterns across different gene sets"
echo ""
echo "=== FILES CREATED ==="
echo "- ${PROJECT_OUTPUT_DIR}/occurrences.txt (main results - ready for R filtering)"
echo "- ${PROJECT_OUTPUT_DIR}/PWMthresholds.txt (score thresholds used)"
echo "- ${PROJECT_OUTPUT_DIR}/sequences.mf (manifest file)"
echo "- ${PROJECT_OUTPUT_DIR}/sequences.mf.dict (sequence dictionary)"
echo "- ${PROJECT_OUTPUT_DIR}/hist_* (score distribution histograms)"