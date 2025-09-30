# Quick Reference: Genotype Variance Analysis

## Purpose
Analyze motif conservation/mutation patterns across genotypes to explain expression divergence.

## Key Question
**Do regulatory differences explain why the same gene is expressed differently in different varieties?**

## Critical Insight from Paper
**This is NOT a traditional differential expression analysis.**

The paper uses **CNN model predictions** to classify genes, not direct RNA-seq data:
1. Train DeepCRE models on multiple varieties
2. Generate CNN probability predictions (0-1 scale) for each gene
3. Calculate **variance in predictions** across varieties
4. Variance **>0.005** = differential expression | **≤0.005** = uniform expression

**Why?** This validates that CNN predictions correlate with actual regulatory variation. Genes with predicted differential expression should show enrichment for mutated motifs.

---

## Data You Need (Future)

### 1. CNN Prediction Data for Multiple Genotypes (Paper Methodology)
```
gene_id       | genotype   | cnn_prob | class | variance | status
Vitvi01g00010 | PN40024    | 0.62     | high  | 0.003    | uniform
Vitvi01g00010 | Cabernet   | 0.58     | high  | 0.003    | uniform
Vitvi01g00020 | PN40024    | 0.35     | low   | 0.012    | differential
Vitvi01g00020 | Cabernet   | 0.78     | high  | 0.012    | differential
```

**Where this comes from:**
- Train DeepCRE MSR models on multiple Vitis varieties
- Generate CNN probability predictions (0-1 scale) for each gene
- Calculate variance in predictions across varieties
- **Key threshold:** variance >0.005 = differential, ≤0.005 = uniform

**Note:** Paper uses CNN predictions, not RNA-seq data directly. This validates that predicted differential expression correlates with regulatory variation.

### 2. Gene Classifications
- **DE genes:** Variance in CNN predictions >0.005 across genotypes
- **UE genes:** Variance in CNN predictions ≤0.005 across genotypes

---

## Analysis Logic (Simplified)

### Step 1: Map Motifs to All Genotypes
```
Gene X in Genotype A: [p0m02, p1m05, p0m07]
Gene X in Genotype B: [p0m02, p0m07]          ← p1m05 MISSING
Gene X in Genotype C: [p0m02, p1m05, p0m07]
```

### Step 2: Classify Each Motif
```
p0m02: Present 3/3 genotypes → CONSERVED
p1m05: Present 2/3 genotypes → MUTATED (variable)
p0m07: Present 3/3 genotypes → CONSERVED
```

### Step 3: Count Gene-Level Patterns
```
Gene X: 2 conserved motifs, 1 mutated motif
Gene Y: 3 conserved motifs, 0 mutated motifs
Gene Z: 0 conserved motifs, 2 mutated motifs
```

### Step 4: Compare Gene Sets
```
DE genes: 40% have mutated motifs  ← HIGH
UE genes: 15% have mutated motifs  ← LOW

→ Conclusion: Regulatory variation drives expression differences
```

---

## Scripts Overview

### 00_prepare_genotype_gene_lists.R
**What it does:** Identifies DE and UE genes from expression data

**Inputs:**
- Multi-genotype expression data
- Gene annotations (GFF)

**Outputs:**
- DE gene coordinate file
- UE gene coordinate file

**When to run:** Before FASTA extraction

---

### 01_genotype_variance_analysis.R
**What it does:** Main analysis comparing motif patterns

**Inputs:**
- Filtered motif occurrences (from BLAMM)
- DE gene coordinates
- UE gene coordinates

**Outputs:**
- Gene lists by conservation status
- Statistical comparison tables
- Visualization plots

**When to run:** After motif mapping complete

---

## Manual Steps Required

Between the scripts, you need to:

1. **Extract FASTA sequences**
   ```bash
   bash extract_range_to_fasta_genomic.sh genome_name gene_coordinates.txt
   ```

2. **Merge genotype FASTAs**
   ```bash
   cat genotype1.fasta genotype2.fasta > merged.fasta
   ```

3. **Run BLAMM** (with paper parameters)
   ```bash
   blamm -e 0.0001 -w 14 motif_pwms.txt merged.fasta > occurrences.txt
   ```
   - e-value <0.0001 (paper specification)
   - word_size 14bp

4. **Filter occurrences**
   ```bash
   Rscript ../mo_proj/filter_blamm_occ_00.R occurrences.txt
   ```
   - Then filter to EPM preferred positional ranges (from mo_range outputs)
   - Critical: motifs outside preferred ranges should be excluded

---

## Key Statistics Explained

### Conservation Rate
```
conserved_rate = conserved_motifs / total_motifs
```
Higher in UE genes (expected)

### Mutation Rate (Gene-Level)
```
mutation_rate = genes_with_any_mutated / total_genes
```
Higher in DE genes (expected)

### Bootstrap Confidence
```
Repeat analysis 1000 times with random sampling
→ Get distribution of results
→ Calculate mean and confidence intervals
```
Paper uses: 100 genes per sample, 1000 iterations

### Statistical Tests (Paper Methodology)
- **Fisher's exact test** for gene set intersections
- **Chi-squared test** for conservation rate differences
- **Significance threshold:** p <0.0001
- **Example from paper:** 15 genotypes, 314 DE genes, 27,993 UE genes

---

## Expected Results

### If Hypothesis is TRUE:

```
┌─────────────────┬──────────┬──────────┐
│ Gene Set        │ Conserved│ Mutated  │
├─────────────────┼──────────┼──────────┤
│ DE genes        │ 60%      │ 40%   ← HIGH
│ UE genes        │ 85%      │ 15%   ← LOW
└─────────────────┴──────────┴──────────┘

→ Motif variation explains expression differences
```

### If Hypothesis is FALSE:

```
┌─────────────────┬──────────┬──────────┐
│ Gene Set        │ Conserved│ Mutated  │
├─────────────────┼──────────┼──────────┤
│ DE genes        │ 70%      │ 30%
│ UE genes        │ 75%      │ 25%   ← SIMILAR
└─────────────────┴──────────┴──────────┘

→ Expression differences due to other factors
```

---

## Script Annotations Guide

Both scripts have **detailed section-by-section annotations** with:

### Section Headers
```R
######################################################################################
# SECTION X: DESCRIPTIVE NAME
######################################################################################
```

### Description Blocks
```R
#
# DESCRIPTION:
# What this section does
#
# BIOLOGICAL INTERPRETATION:
# Why we do this
#
# ORIGINAL LOGIC:
# Reference to original implementation
#
```

### Code Comments
```R
# STEP 1: Detailed explanation of step
# STEP 2: What happens next
```

### Pseudocode
```R
# PSEUDOCODE - run when data available
# actual_results <- function_call(data)

cat("  PLACEHOLDER: Function defined\n")
cat("  TODO: Run when data available\n")
```

---

## Current Status

**FRAMEWORK READY - AWAITING DATA**

### What's Complete:
- Full analysis logic implemented
- All functions defined
- Statistical methods ready
- Visualization code prepared
- Extensive documentation

### What's Needed:
- Multi-genotype expression data
- In silico perturbation predictions
- Or cross-variety RNA-seq results

---

## Tips for Future Work

### When You Get Data:

1. **Check data format matches expected structure**
   - Column names correct?
   - Genotype identifiers consistent?
   - Expression classifications present?

2. **Start with small test**
   - Use 2 genotypes first
   - Test with 100 genes
   - Verify logic before full run

3. **Monitor bootstrap progress**
   - Scripts print iteration counts
   - Check for convergence
   - Adjust iterations if needed

4. **Validate results**
   - Do conservation rates make sense?
   - Are p-values reasonable?
   - Check a few genes manually

### Common Issues:

**Location string format:**
- Must be: `genotype_chr:start-end`
- Script extracts genotype prefix
- Make sure underscores are correct

**Sample size imbalance:**
- UE sets typically much larger
- Bootstrap handles this automatically (100 genes, 1000 iterations)
- But verify sample sizes are appropriate

**Motif length filtering:**
- Default minimum = 14bp (paper specification)
- Adjust WORD_SIZE parameter if needed
- Check motif length distribution first

**Positional filtering:**
- **Critical:** Must filter to EPM preferred ranges (from mo_range outputs)
- Example: epmVitis-M006-p0m02 might prefer -500 to -100 from TSS
- Paper excluded motifs outside preferred ranges
- This is essential for accurate results

---

## Related Documentation

- **Full Details:** `README.md` in this directory
- **Original Script:** `../../../deepcre_original/moca_blue/mo_proj/mo_genotype_variance.v1.4.R`
- **Pipeline Overview:** `../../README.md`
- **Methods Paper:** https://www.nature.com/articles/s41467-024-47744-0

---

## Quick Troubleshooting

| Problem | Solution |
|---------|----------|
| No file found errors | Update file path patterns in script header |
| Empty results | Check location string format matches between files |
| Bootstrap too slow | Reduce iterations (lines with `n_iterations` parameter) |
| Plots look weird | Check data has both conserved and mutated categories |
| Gene counts don't match | Verify GFF gene IDs match expression data IDs |

---

## Contact

Questions? Check:
1. Section-by-section annotations in the scripts
2. README.md in this directory
3. Original script comments
4. Methods section of reference paper