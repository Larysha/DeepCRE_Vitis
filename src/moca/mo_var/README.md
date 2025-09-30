# Genotype/Variety Variance Analysis (mo_var)

## Overview

This module analyzes **motif conservation and mutation patterns across multiple genotypes or varieties** to understand the relationship between regulatory element variation and differential gene expression.

**Refactored from:** `mo_genotype_variance.v1.4.R` (original moca_blue project)
**Original Author:** Dr. Simon M. Zumkeller (2023-08-30)
**Reference:** https://www.nature.com/articles/s41467-024-47744-0

## Biological Question

**Do regulatory element variations explain expression divergence between genotypes?**

### Key Hypotheses:
1. **Differentially expressed (DE) genes** should show enrichment for **MUTATED motifs** (variable across genotypes)
2. **Uniformly expressed (UE) genes** should show enrichment for **CONSERVED motifs** (present in all genotypes)
3. Motif variance patterns correlate with expression phenotypes

## Analysis Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│ STEP 0: Multi-Genotype Data Preparation                        │
│ - Obtain expression data for multiple varieties                │
│ - Identify DE and UE gene sets                                 │
│ - In silico perturbation predictions (future)                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 1: Preprocessing                                          │
│ Script: 00_prepare_genotype_gene_lists.R                       │
│ - Extract DE and UE gene lists                                 │
│ - Create GFF annotation subsets                                │
│ - Prepare coordinate files for FASTA extraction                │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 2: Sequence Extraction (MANUAL)                           │
│ bash ../ref_seq/extract_range_to_fasta_genomic.sh              │
│ - Extract regulatory regions for DE genes (all genotypes)      │
│ - Extract regulatory regions for UE genes (all genotypes)      │
│ - Merge multi-genotype FASTAs                                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 3: Motif Mapping (MANUAL)                                 │
│ blamm <pwms> <merged_fasta> > occurrences.txt                  │
│ - Map motifs to all genotype sequences                         │
│ - Filter with: ../mo_proj/filter_blamm_occ_00.R                │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 4: Variance Analysis                                      │
│ Script: 01_genotype_variance_analysis.R                        │
│ - Count motif occurrences per gene per genotype                │
│ - Classify motifs as conserved or mutated                      │
│ - Statistical comparison: DE vs UE genes                       │
│ - Bootstrap resampling for confidence intervals                │
│ - Generate visualizations                                      │
└─────────────────────────────────────────────────────────────────┘
```

## Scripts

### 00_prepare_genotype_gene_lists.R
**Purpose:** Preprocessing - prepare gene lists and annotations

**Inputs:**
- Multi-genotype expression data (FUTURE - placeholder structure defined)
- Gene annotation GFF file
- Expression classification criteria

**Outputs:**
- `*_differential_expressed_locations.txt` - DE gene coordinates
- `*_uniformly_expressed_locations.txt` - UE gene coordinates
- Gene-location mapping files

**Status:** ⚠️ TEMPLATE - awaiting multi-genotype data

---

### 01_genotype_variance_analysis.R
**Purpose:** Main analysis - motif conservation patterns

**Inputs:**
- Filtered motif occurrences (from BLAMM + filter_blamm_occ_00.R)
- DE gene annotations
- UE gene annotations

**Analysis Steps:**
1. **Section 1-2:** Load gene annotations (DE and UE)
2. **Section 3:** Load and filter motif occurrences
3. **Section 4:** Merge motifs with DE genes, classify conservation
4. **Section 5:** Merge motifs with UE genes (with bootstrap resampling)
5. **Section 6:** Statistical comparison of conservation rates
6. **Section 7:** Bootstrap analysis for confidence intervals
7. **Section 8:** Create classified gene lists
8. **Section 9:** Generate visualizations
9. **Section 10:** Write output files

**Outputs:**
- `*_combined_gene_lists.csv` - Genes classified by expression and motif status
- `*_variance_comparison.csv` - Conservation rate statistics
- `*_conservation_rates.pdf` - Bar plot visualization
- `*_gene_conservation_boxplot.pdf` - Bootstrap distribution plot

**Status:** ⚠️ TEMPLATE - awaiting input data

## Key Concepts

### Conservation Classification

**CONSERVED motif:**
- Present in **ALL** genotypes at a gene locus
- Mathematical criterion: `n_occurrences / n_genotypes = integer`
- Example: 6 occurrences across 3 genotypes → 6/3 = 2 (conserved)

**MUTATED motif:**
- Present in **SOME** but not all genotypes
- Mathematical criterion: `n_occurrences / n_genotypes ≠ integer`
- Example: 5 occurrences across 3 genotypes → 5/3 = 1.67 (mutated)

### Bootstrap Resampling

Because UE gene sets are typically much larger than DE sets:
1. Randomly sample UE genes to match DE set size
2. Repeat 100-1000 times
3. Calculate average statistics
4. Generate confidence intervals

This ensures fair statistical comparison between gene sets of different sizes.

## Input Data Requirements

### When Data Available:

**1. Multi-Genotype Expression Data**
```
gene_id       | genotype    | expression_tpm | expression_class | differential_status
Vitvi01g00010 | PN40024     | 125.3          | high            | uniform
Vitvi01g00010 | Cabernet    | 118.7          | high            | uniform
Vitvi01g00020 | PN40024     | 5.2            | low             | differential
Vitvi01g00020 | Cabernet    | 95.3           | high            | differential
```

**2. Motif Occurrences** (from BLAMM)
- Must include genotype identifier in location string
- Format: `genotype_chr:start-end`
- Filtered to gene proximity regions (±1500bp)

**3. Gene Annotations**
- GFF3 format
- With flanking regions (±1000bp)

## Future Data Sources

### In Silico Perturbation (Planned)
- Train CNN models on variety-specific data
- Generate in silico perturbed sequences
- Predict expression for sequences with/without specific motifs
- Compare predictions across varieties

### Multi-Variety RNA-seq
- Differential expression across Vitis varieties
- Drought-tolerant vs sensitive cultivars
- PN40024, Cabernet Sauvignon, Chardonnay, etc.

## Interpretation Guide

### Expected Results:

**High motif mutation rate in DE genes + Low mutation rate in UE genes**
→ Regulatory variation drives expression divergence

**Similar mutation rates in both gene sets**
→ Expression differences likely due to other factors (trans-acting, epigenetic)

**Specific motif types (p0m/p1m) with differential conservation**
→ Motif class-specific evolutionary constraints

## Usage Examples

### When Data Available:

```bash
# Step 1: Preprocessing
Rscript 00_prepare_genotype_gene_lists.R \
  path/to/multi_genotype_expression.csv \
  path/to/gene_annotation.gff3 \
  path/to/output_dir

# Step 2: Extract sequences (manual)
bash ../ref_seq/extract_range_to_fasta_genomic.sh \
  vitis_vinifera_PN40024 \
  output_dir/20250930_vitis_genotypes_differential_expressed_locations.txt

# Repeat for UE genes and all genotypes

# Step 3: Merge FASTAs
cat genotype1.fasta genotype2.fasta genotype3.fasta > merged_genotypes.fasta

# Step 4: BLAMM mapping
blamm path/to/motif_pwms.txt merged_genotypes.fasta > occurrences.txt

# Step 5: Filter occurrences
Rscript ../mo_proj/filter_blamm_occ_00.R occurrences.txt

# Step 6: Run variance analysis
Rscript 01_genotype_variance_analysis.R \
  path/to/filtered_occurrences.txt \
  path/to/de_gene_locations.txt \
  path/to/ue_gene_locations.txt \
  path/to/output_dir
```

## Current Status

⚠️ **FRAMEWORK COMPLETE - AWAITING DATA**

Both scripts are fully implemented with:
- ✅ Detailed section-by-section annotations
- ✅ Complete analysis logic (pseudocode where data unavailable)
- ✅ Statistical methods defined
- ✅ Visualization functions ready
- ✅ Pattern-based file finding
- ✅ Compatible with existing pipeline structure

**Next Steps:**
1. Obtain/generate multi-genotype expression data
2. Run preprocessing script
3. Execute pipeline with actual data

## Notes

- All scripts use **automatic file finding** (no hardcoded paths)
- Compatible with existing `mo_proj` pipeline outputs
- Follows vitis_cre naming conventions (p0m/p1m format)
- Extensively annotated for future reference
- No modifications made to original `deepcre_original` repository

## References

1. **Original Implementation:** moca_blue `mo_genotype_variance.v1.4.R`
2. **Methods Paper:** https://www.nature.com/articles/s41467-024-47744-0
3. **Related Pipelines:**
   - `mo_proj/` - Single-genome motif analysis
   - `mo_go/` - Functional enrichment analysis

## Contact

For questions about this module, refer to:
- Original author: Dr. Simon M. Zumkeller
- Refactored implementation: LR (2025)
- Pipeline documentation: `../README.md`