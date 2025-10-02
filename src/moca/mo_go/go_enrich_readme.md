# MO_GO Analysis Pipeline

This pipeline analyzes motif-gene-GO relationships to identify functional enrichment of regulatory motifs in Vitis vinifera drought response.

## Pipeline Overview

```
01_gene_mapper → 02_motif_predictability → 03_go_enrichment → 04_visualization (future)
```

## Scripts

### Script 01: Gene-Motif-GO Mapper
**File**: `01_mo_gene_mapper_vitis.R`
**Purpose**: Integrate motif mappings with GO annotations and expression/prediction data

**Inputs**:
- Filtered motif mappings (`*_filtered.csv`)
- CNN predictions
- Expression targets (TPM-based)
- GO annotations (GMT format)

**Outputs**:
- `gene_mapping.csv` - Gene-level summary with expression/prediction classes
- `detailed_motif_go_associations.csv` - All gene-motif-GO triplets

### Script 02: Motif Predictability Analysis
**File**: `02_motif_predictability_vitis.R`
**Purpose**: Analyze overall motif performance (NOT broken down by GO category)

**Follows**: Original moca_blue methodology (motif_predictabilityV1.5.R)

**Inputs**:
- Gene mapping from script 01
- Filtered motif mappings

**Outputs**:
- `motif_predictability_performance.csv`

**Key Metrics** (per motif):
- `expr_class_low/high` - Expression class distribution
- `prob_class_low/high` - Prediction class distribution
- `perf_FALSE/TRUE` - Correct vs incorrect predictions
- `mo_metacluster` - Motif type (0=p0m activator, 1=p1m repressor)
- `prob_rat01` - Directional ratio test (does motif predict expected direction?)
- `expr_rat01` - Expression ratio test
- `TPR_p0_expr`, `TPR_p1_expr` - True positive rates for expression
- `TPR_p0_prob`, `TPR_p1_prob` - True positive rates for predictions
- `TPR_TF` - Overall prediction accuracy
- `chi_expr_class`, `chi_prob_class`, `chi_TF` - Chi-squared tests

**Note**: Original moca_blue loaded GO data but never performed GO-specific analysis (line 139 commented out). This script matches that behavior.

### Script 03: GO Enrichment Analysis
**File**: `03_go_enrichment_analysis.R`
**Purpose**: Test which GO terms are enriched in correctly predicted genes per motif

**NEW FUNCTIONALITY**: Original scripts never performed this analysis

**Inputs**:
- Motif performance from script 02
- Detailed associations from script 01
- Gene mapping from script 01

**Analysis**:
For each motif:
1. Identify genes with correct predictions
2. For each GO term, test enrichment using **hypergeometric test**:
   - Null hypothesis: Genes randomly distributed across GO terms
   - Alternative: GO term over-represented in motif's correctly predicted genes
3. Apply Benjamini-Hochberg FDR correction
4. Calculate fold enrichment (observed/expected)

**Filters**:
- Minimum 3 genes with motif in GO term
- FDR < 0.05 for significance
- Fold enrichment ≥ 1.5 for biological relevance

**Outputs**:
- `all_results.csv` - All motif-GO enrichment tests
- `significant_enrichments.csv` - FDR < 0.05
- `relevant_enrichments.csv` - FDR < 0.05 AND fold ≥ 1.5
- `motif_summary.csv` - Per-motif: how many enriched GO terms
- `go_term_summary.csv` - Per-GO: which motifs enrich it
- `top_per_motif.csv` - Top 10 enriched terms per motif

**Key Columns**:
- `genes_in_term` - Observed count
- `expected` - Expected under null
- `fold_enrichment` - Observed/expected
- `p_value` - Raw hypergeometric p-value
- `FDR` - Benjamini-Hochberg corrected p-value

## Running the Pipeline

```bash
cd vitis_cre/src/moca/mo_go

# Step 1: Map genes to motifs and GO categories
Rscript 01_mo_gene_mapper_vitis.R

# Step 2: Analyze motif performance (overall)
Rscript 02_motif_predictability_vitis.R

# Step 3: GO enrichment analysis
Rscript 03_go_enrichment_analysis.R
```

## Biological Interpretation

### What this pipeline tells you:

1. **Script 02**: Which motifs predict expression well overall?
   - High `TPR_TF` = motif is predictive
   - `prob_rat01 = 1` = motif shows expected directional bias
   - Significant `chi_TF` = performance differs from random

2. **Script 03**: For well-performing motifs, what biological functions are they regulating?
   - Significant FDR = GO term enriched beyond chance
   - High fold enrichment = strong biological signal
   - Example: "Motif epm_vitis_ssr_p0m15_M01_F predicts drought response genes enriched in GO:0009414 (response to water deprivation), FDR=0.0001, 3.5-fold enriched"

### Example biological question answered:

> "Are p0m activator motifs binding to genes involved in water stress response?"

**Answer from pipeline**:
- Script 02: p0m motifs have high TPR_TF (predictive)
- Script 03: p0m motifs enriched in GO terms like:
  - GO:0009414 (response to water deprivation)
  - GO:0006979 (response to oxidative stress)
  - GO:0009651 (response to salt stress)
- **Conclusion**: Yes, these motifs regulate drought-related pathways

## Comparison to Original moca_blue

### What original scripts did:
✅ Motif performance analysis (contingency tables, TPR, chi-squared)
✅ Loaded GO annotations
❌ Did NOT perform GO enrichment tests (contingency_table_GO commented out)

### What this pipeline adds:
✅ All original functionality preserved in script 02
✅ NEW: Hypergeometric GO enrichment tests (script 03)
✅ NEW: FDR correction
✅ NEW: Effect size (fold enrichment)
✅ NEW: Multiple output summaries

## Notes

- **Statistical test**: Hypergeometric (standard for GO enrichment, mathematically equivalent to Fisher's exact for 2x2 tables)
- **Multiple testing**: Benjamini-Hochberg FDR correction across all motif-GO tests
- **Effect size**: Report both p-value AND fold enrichment (avoid "statistically significant but biologically boring" results)
- **Minimum threshold**: 3 genes required per motif-GO combination
