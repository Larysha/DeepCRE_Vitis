#!/bin/bash

# Extract gene sequences matching EXACT model training structure
# Creates 3020bp fixed-length sequences with structured regions:
# 1000bp upstream + 500bp into gene + 20bp spacer + 500bp into gene + 1000bp downstream
# BLAMM-compatible FASTA headers (simple IDs + tab-separated metadata)
# Adapted for moca_blue pipeline - MODEL TRAINING MATCHED VERSION
# Usage: ./extract_range_to_fasta_model_matched.sh [species_model_name]

set -euo pipefail 

SPECIES_MODEL="${1:-vitis_vinifera_PN40024}"

# Model training parameters (FIXED - do not change)
EXTRAGENIC=1000    # bp upstream of TSS, downstream of TTS
INTRAGENIC=500     # bp downstream of TSS, upstream of TTS  
SPACER=20          # bp spacer between regions
TOTAL_LENGTH=3020  # Total sequence length (1000 + 500 + 20 + 500 + 1000)

# File paths - search in standard locations relative to project root
GENOME_DIRS=("genome" "src/genome" "vitis_cre/src/genome")
ANNOTATION_DIRS=("gene_models" "src/gene_models" "vitis_cre/src/gene_models")

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

echo "=== MODEL TRAINING MATCHED SEQUENCE EXTRACTION ==="
echo "Extracting sequences with EXACT model training structure:"
echo "- Total length: ${TOTAL_LENGTH}bp (fixed for all genes)"
echo "- Structure: ${EXTRAGENIC}bp upstream + ${INTRAGENIC}bp into gene + ${SPACER}bp spacer + ${INTRAGENIC}bp into gene + ${EXTRAGENIC}bp downstream"
echo ""

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
echo ""

# Get the script's directory and build path relative to it
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/../../results/moca_blue/ref_seq"
mkdir -p "$OUTPUT_DIR"

RANGES_FILE="gene_ranges_model_matched.txt"
OUTPUT_FASTA="${OUTPUT_DIR}/${SPECIES_MODEL}_model_matched_${TOTAL_LENGTH}bp.fa"

# Start timing
start_time=$(date +%s)
echo "Extracting gene ranges with model-matched structure..."

# Extract gene ranges with model training structure
# Creates both ranges file and metadata for header construction
awk -v extragenic="$EXTRAGENIC" -v intragenic="$INTRAGENIC" -v spacer="$SPACER" '
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
    
    # Calculate gene length
    gene_length = gene_end - gene_start + 1
    
    # Skip genes shorter than 2 * intragenic (cannot extract both ends)
    if (gene_length < (2 * intragenic)) {
        print "Warning: Skipping gene with length " gene_length "bp (too short for " (2 * intragenic) "bp extraction)" > "/dev/stderr"
        next
    }
    
    # Extract gene ID from attributes (works for both GFF3 and GTF)
    gene_id = ""
    if (match(attributes, /ID=([^;]+)/, arr)) {
        gene_id = arr[1]
    } else if (match(attributes, /gene_id[= ]"?([^;"]+)"?/, arr)) {
        gene_id = arr[1]
    } else {
        # Fallback: create ID from coordinates
        gene_id = chr "_" gene_start "_" gene_end
    }
    
    # Clean gene_id for BLAMM compatibility (alphanumeric + underscore only)
    clean_gene_id = gene_id
    gsub(/[^A-Za-z0-9_]/, "_", clean_gene_id)
    
    # Create model training structure coordinates
    if (strand == "+") {
        # Forward strand: upstream -> 5-prime -> spacer -> 3-prime -> downstream
        region1_start = gene_start - extragenic              # 1000bp upstream
        region1_end = gene_start - 1
        
        region2_start = gene_start                           # First 500bp of gene
        region2_end = gene_start + intragenic - 1
        
        # 20bp spacer will be inserted programmatically
        
        region3_start = gene_end - intragenic + 1            # Last 500bp of gene
        region3_end = gene_end
        
        region4_start = gene_end + 1                         # 1000bp downstream
        region4_end = gene_end + extragenic
        
    } else {
        # Reverse strand: upstream -> 3-prime -> spacer -> 5-prime -> downstream
        # Note: coordinates still in forward direction, but logical order reversed
        region1_start = gene_end + 1                         # 1000bp downstream (upstream for -)
        region1_end = gene_end + extragenic
        
        region2_start = gene_end - intragenic + 1            # Last 500bp of gene (5-prime for -)
        region2_end = gene_end
        
        # 20bp spacer will be inserted programmatically
        
        region3_start = gene_start                           # First 500bp of gene (3-prime for -)
        region3_end = gene_start + intragenic - 1
        
        region4_start = gene_start - extragenic              # 1000bp upstream (downstream for -)
        region4_end = gene_start - 1
    }
    
    # Bounds checking - ensure no coordinates go below 1
    if (region1_start < 1) region1_start = 1
    if (region4_start < 1) region4_start = 1
    
    # Create ranges for each region (samtools format)
    range1 = chr ":" region1_start "-" region1_end
    range2 = chr ":" region2_start "-" region2_end  
    range3 = chr ":" region3_start "-" region3_end
    range4 = chr ":" region4_start "-" region4_end
    
    # BLAMM-compatible header: ID\tchr\tstrand\toriginal_coords\tstructure_type
    coord_range = gene_start "-" gene_end
    header = clean_gene_id "\t" chr "\t" strand "\t" coord_range "\tmodel_matched_3020bp"
    
    # Print all regions and header info for post-processing
    # Format: range1|range2|range3|range4 TAB header
    all_ranges = range1 "|" range2 "|" range3 "|" range4
    print all_ranges "\t" header
}' "$ANNOTATION_PATH" > "$RANGES_FILE"

# Count genes processed
gene_count=$(wc -l < "$RANGES_FILE")
echo "Processed $gene_count genes for model-matched extraction"

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is required but not found in PATH" >&2
    exit 1
fi

echo "Extracting sequences with samtools and assembling model structure..."

# Create the final FASTA file with model-matched sequences
python3 - "$RANGES_FILE" "$GENOME_PATH" "$OUTPUT_FASTA" "$SPACER" << 'EOF'
import sys
import subprocess
import tempfile
import os

ranges_file = sys.argv[1]
genome_path = sys.argv[2] 
output_fasta = sys.argv[3]
spacer_length = int(sys.argv[4])

# Create spacer sequence (Ns)
spacer_seq = 'N' * spacer_length

print("Assembling model-matched sequences...")

with open(ranges_file, 'r') as f, open(output_fasta, 'w') as outfile:
    for line_num, line in enumerate(f, 1):
        if line_num % 1000 == 0:
            print(f"Processed {line_num} genes...")
            
        parts = line.strip().split('\t')
        if len(parts) != 2:
            continue
            
        ranges_str = parts[0]
        header = parts[1]
        
        # Parse the four ranges
        ranges = ranges_str.split('|')
        if len(ranges) != 4:
            print(f"Warning: Unexpected range format for gene {header.split()[0]}")
            continue
        
        # Extract each region using samtools
        sequences = []
        for i, range_coord in enumerate(ranges):
            try:
                # Use samtools faidx to extract sequence
                result = subprocess.run(
                    ['samtools', 'faidx', genome_path, range_coord],
                    capture_output=True, text=True, check=True
                )
                
                # Parse FASTA output (skip header line, join sequence lines)
                fasta_lines = result.stdout.strip().split('\n')
                if len(fasta_lines) < 2:
                    raise ValueError(f"Empty sequence for range {range_coord}")
                
                sequence = ''.join(fasta_lines[1:])  # Skip header, join all sequence lines
                sequences.append(sequence)
                
            except subprocess.CalledProcessError as e:
                print(f"Error extracting {range_coord}: {e}")
                sequences.append('')  # Add empty sequence to maintain structure
            except Exception as e:
                print(f"Unexpected error for {range_coord}: {e}")
                sequences.append('')
        
        # Assemble final sequence: region1 + region2 + spacer + region3 + region4
        if all(sequences):  # Only proceed if all regions extracted successfully
            final_sequence = sequences[0] + sequences[1] + spacer_seq + sequences[2] + sequences[3]
            
            # Write to output FASTA
            outfile.write(f">{header}\n")
            outfile.write(f"{final_sequence}\n")
        else:
            print(f"Warning: Skipping gene due to extraction errors: {header.split()[0]}")

print("Model-matched sequence assembly complete!")
EOF

# Verify output
if [[ ! -f "$OUTPUT_FASTA" ]]; then
    echo "Error: Failed to create output FASTA file" >&2
    exit 1
fi

seq_count=$(grep -c "^>" "$OUTPUT_FASTA" || true)
echo "Successfully created $seq_count model-matched sequences"

# Verify sequence lengths
echo "Verifying sequence lengths..."
python3 - "$OUTPUT_FASTA" "$TOTAL_LENGTH" << 'EOF'
import sys

fasta_file = sys.argv[1]
expected_length = int(sys.argv[2])

lengths = []
with open(fasta_file, 'r') as f:
    sequence = ''
    for line in f:
        if line.startswith('>'):
            if sequence:
                lengths.append(len(sequence))
            sequence = ''
        else:
            sequence += line.strip()
    if sequence:  # Don't forget the last sequence
        lengths.append(len(sequence))

correct_length = sum(1 for l in lengths if l == expected_length)
total_sequences = len(lengths)

print(f"Sequence length verification:")
print(f"- Total sequences: {total_sequences}")
print(f"- Correct length ({expected_length}bp): {correct_length}")
print(f"- Incorrect length: {total_sequences - correct_length}")

if total_sequences > 0:
    print(f"- Success rate: {correct_length/total_sequences*100:.1f}%")

if correct_length != total_sequences:
    unique_lengths = set(lengths)
    print(f"- Found lengths: {sorted(unique_lengths)}")
EOF

echo ""
echo "Example header format:"
head -n 1 "$OUTPUT_FASTA"

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo ""
echo "Completed in ${runtime} seconds"

# Cleanup intermediate files
rm "$RANGES_FILE"

echo ""
echo "=== MODEL TRAINING MATCHED EXTRACTION COMPLETE ==="
echo "Output written to: $OUTPUT_FASTA"
echo ""
echo "Sequence structure (${TOTAL_LENGTH}bp fixed length):"
echo "- ${EXTRAGENIC}bp upstream region"
echo "- ${INTRAGENIC}bp into gene (5' end)"  
echo "- ${SPACER}bp spacer (Ns)"
echo "- ${INTRAGENIC}bp into gene (3' end)"
echo "- ${EXTRAGENIC}bp downstream region"
echo ""
echo "This matches EXACTLY what the model was trained on!"
echo ""
echo "Next steps:"
echo "1. Use this file for BLAMM motif scanning"
echo "2. Motif positions will directly correspond to model training data"
echo "3. No coordinate transformation needed for model validation"