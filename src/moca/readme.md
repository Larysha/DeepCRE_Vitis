# MoCa: Motif Characterisation & Annotation from Deep Learning

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/)

A comprehensive toolbox for analysing DNA motifs derived from deep learning model feature extraction with [TF-MoDISco](https://github.com/kundajelab/tfmodisco). Transforms discovered regulatory elements into biologically interpretable modules through statistical validation, clustering, and genomic mapping.

**Original pipeline**: Dr Simon M. Zumkeller ([DeepCRE](https://github.com/NAMlab/DeepCRE))  
Adapted here for a V. vinifera drought stress study

## Quick Start

**Input**: TF-MoDISco HDF5 files from model interpretation  
**Output**: Ranked regulatory motifs with genomic locations and functional annotations

**Workflow**: Extract motifs → Analyse positions and importance → Cluster similarities → Map to genome → Validate motif and model performance → Functional enrichment → Cross-genotype analysis

## Project Structure

```
moca/
├── mo_nom/          # Extract motifs and assign names
├── mo_imp/          # Visualise saliency maps
├── mo_range/        # Analyse positional preferences
├── mo_clu/          # Cluster similar motifs
├── ref_seq/         # Map motifs to genome sequences
├── mo_proj/         # Filter and validate performance
├── mo_go/           # Functional enrichment and GO analysis
├── mo_var/          # Cross-genotype motif variance analysis (in progress)
└── utils.R          # Shared utilities
```

## Module Overview

### 1. `mo_nom/` - Motif Extraction
**Purpose**: Extract motifs from TF-MoDISco and create systematic nomenclature

**Key script**: `get_matrices.R`

**Outputs**:
- JASPAR format motifs for genome scanning
- Summary statistics with importance scores
- Systematic naming: `epm_vitis_ssr_p0m05` (Expression Predictive Motif)

**Motif types**:
- **PFM**: Raw count matrices for BLAMM scanning
- **CWM**: Neural network contributions preserving regulatory direction

---

### 2. `mo_imp/` - Saliency Analysis
**Purpose**: Visualise position-wise importance and rank genes by regulatory impact

**Key script**: `plot_saliency_maps.R`

**Outputs**:
- Gene importance rankings with percentile scores
- Position-wise saliency profiles across regulatory regions
- Sequence logos for top regulatory genes

**Application**: Identify exemplar genes with strongest regulatory signals for detailed study

---

### 3. `mo_range/` - Positional Preferences  
**Purpose**: Characterise where motifs occur relative to gene structure (TSS vs TTS)

**Key scripts**: `get_seql_per_pattern.R` → `motif_ranges_meta.R`

**Outputs**:
- Q10/Q90 positional quartiles for filtering
- TSS vs TTS enrichment patterns
- Coordinate-standardised motif distributions

---

### 4. `mo_clu/` - Motif Clustering
**Purpose**: Group similar motifs to reduce redundancy and reveal regulatory families

**Key script**: `cluster_motifs.R`

**Method**: Sandelin-Wassermann algorithm with 95th percentile similarity threshold

**Outputs**:
- Unique vs redundant motif classifications
- Phylogenetic trees showing regulatory relationships
- Pairwise similarity matrices for custom analysis

**Application**: Focus downstream analysis on distinct regulatory mechanisms

---

### 5. `ref_seq/` - Genome Mapping
**Purpose**: Map discovered motifs to actual genomic sequences with statistical significance

**Key scripts**: 
- `extract_range_to_fasta_genomic.sh` - Extract gene sequences with flanking regions
- `blamm_motif_mapping.sh` - Statistical motif scanning with empirical p-values

**Method**: BLAMM (BLAS Accelerated Motif Matching) with p < 0.0001 threshold

**Outputs**: Statistically significant motif occurrences with genomic coordinates

---

### 6. `mo_proj/` - Performance Analysis
**Purpose**: Filter, validate, and assess predictive performance of discovered motifs

#### Core Scripts:

**`filter_blamm_occ.R`**: Convert BLAMM results to genomic coordinates and filter to training regions

**`filter_genomic_features.R`**: Multi-level filtering pipeline:
- Q10/Q90 positional filtering (core distribution zones)
- Weight region filtering (≥20% preference threshold)  
- Word size filtering (≥14bp minimum)

**`motif_enrichment.R`**: Statistical validation core
- Tests motif enrichment in expected expression classes
- Calculates importance scores: log2(fold-change) with metacluster-specific directions
- Chi-square significance testing with multiple comparison awareness

**`motif_performance.R`**: Visualisation and reporting
- Performance plots by metacluster
- Top performer rankings
- Model validation metrics (accuracy, F1, precision, recall)

---

### 7. `mo_go/` - Functional Enrichment Analysis
**Purpose**: Integrate GO annotations with motif analysis to understand biological context and identify functionally coherent regulatory modules

#### Core Scripts:

**`convert_gene_ids.R`**: One-time gene ID standardisation
- Converts GO annotation IDs to match motif analysis naming conventions
- Builds hash table lookup from genome GFF3 for efficient conversion
- Input: GO GMT format + genome annotation
- Output: Converted GMT with standardised gene IDs

**`mo_gene_mapper_vitis.R`**: Comprehensive gene-motif-function integration
- Cross-references genes across expression data, model predictions, motif occurrences, and GO annotations
- Calculates motif-GO enrichment statistics (observed/expected ratios with significance thresholds)
- Identifies enriched (≥2-fold, ≥3 genes) and depleted (≤0.5-fold, ≥3 genes) associations
- Generates gene sets organised by motif combinations and GO categories (≥3-5 genes per set)
- Creates high-confidence gene sets with multiple lines of evidence

**`motif_predictability_vitis.R`**: Context-dependent motif performance analysis
- Evaluates motif predictive accuracy within specific GO categories (10-1000 genes per combination)
- Compares GO-specific performance to baseline via fold change
- Identifies enhanced (≥1.2-fold) and reduced (≤0.8-fold) performance contexts
- Calculates True Positive Rates for expected expression classes (p0m→high, p1m→low)
- Retains top 10,000 combinations for tractable analysis

**`functional_enrichment.R`**: Statistical association testing
- Creates contingency tables between motifs and GO categories
- Integrates existing enrichment and performance results
- Applies confidence filtering: |prob - 0.5| ≥ 0.1 for reliable predictions
- Motif frequency threshold: ≥10 genes for statistical power
- Calculates True Positive and True Negative prediction counts per motif-GO combination

**Outputs**:
- Gene-level integration tables with expression, predictions, motifs, and GO annotations
- Motif-GO enrichment statistics with fold changes and significance flags
- Context-dependent performance metrics (enhanced/reduced categories)
- Gene sets for pathway analysis and experimental validation

**Application**: Reveals whether motifs have context-dependent regulatory activity. A motif performing well globally but poorly in specific functions may indicate missing co-factors or chromatin state differences. Enhanced performance in specific functions strengthens biological relevance.

---

### 8. `mo_var/` - Cross-Genotype Variance Analysis
**Purpose**: Analyse motif conservation and mutation patterns across multiple genotypes/varieties to understand regulatory variation underlying expression divergence

**Status**: IN PROGRESS - Framework complete, awaiting multi-variety data

#### Core Scripts:

**`prepare_genotype_gene_lists.R`**: Data preparation and gene classification
- Classifies genes based on CNN prediction variance across varieties (threshold: >0.005 = differential, ≤0.005 = uniform)
- Extracts coordinates for differentially expressed (DE) and uniformly expressed (UE) gene sets
- Prepares files for FASTA extraction and BLAMM mapping across genotypes

**`genotype_variance_analysis.R`**: Conservation pattern analysis
- Counts motif occurrences per gene per genotype
- Classifies motifs as conserved (present in all genotypes) or mutated (variable across genotypes)
- Bootstrap resampling for robust comparison (100 genes, 1000 iterations)
- Statistical tests: Fisher's exact and Chi-squared (p < 0.0001)

**Key Methodology** (from paper):
- Uses CNN model predictions, not direct RNA-seq data
- Validates that predicted differential expression correlates with regulatory variation
- BLAMM parameters: e-value < 0.0001, word size 14bp
- Critical filtering: EPMs must be within preferred positional ranges (from mo_range outputs)

**Expected Results**:
- DE genes should show enrichment for mutated motifs (regulatory divergence)
- UE genes should show enrichment for conserved motifs (regulatory stability)

**Application**: Identifies genotype-specific regulatory elements that may explain expression differences between varieties. Useful for understanding drought tolerance mechanisms across cultivars.

---

## Key Concepts


**Metaclusters** (Based on TF-MoDISco conventional output):
- **p0m (metacluster 0)**: Associated with high expression (activators)
- **p1m (metacluster 1)**: Associated with low expression (repressors)

**Validation Strategy** (based on DeepCRE model):
- Uses held-out validation genes not seen during CNN training
- Tests whether motifs show expected expression associations in independent data
- Distinguishes genuine regulatory patterns from spurious correlations

**Quality Indicators**:
- **High importance score** (|score| > 1.0): ≥2-fold enrichment, biologically meaningful
- **Statistical significance** (p < 0.05): Pattern unlikely due to chance
- **High seqlet count** (>100): Robust pattern with multiple supporting instances

## Requirements

**R packages**: `dplyr`, `TFBSTools`, `universalmotif`, `rhdf5`, `reticulate`  
**External tools**: BLAMM, samtools  
**Input format**: TF-MoDISco HDF5 files with standard metacluster structure

## Usage

1. Set working directory to project root
2. Run modules sequentially: `mo_nom → mo_range → mo_clu → ref_seq → mo_proj → mo_go`
3. Optional: `mo_var` for cross-genotype analysis (requires multi-variety data)
4. Scripts auto-detect files based on naming conventions, but can also be specified as command line arguments (refer to script comments to see other customisable parameters)
5. Outputs saved to `../../out/moca_results/[module_name]/`

## Output Interpretation

**Top-performing motifs** (high importance + significance) represent:
- Key drought response elements for experimental validation
- Targets for transcription factor binding studies (DAP-seq, EMSA, ChIP-seq)
- Candidates for transgenic promoter analysis in *Vitis vinifera*

**Performance visualisations** reveal:
- Whether CNN learnt distinct regulatory strategies (metacluster separation)
- Model calibration quality (dual importance correlation)
- Overall predictive accuracy on independent validation data

The pipeline transforms abstract neural network features into ranked, validated regulatory candidates ready for experimental follow-up.

## Citation

If you use moca in your research, please cite:

- **Original DeepCRE framework & moca pipeline:** Peleke et al. (2024). "Deep learning the cis-regulatory code for gene expression in selected model plants." Nature Communications. https://doi.org/10.1038/s41467-024-47744-0
- **Original moca scripts:** Dr Simon M. Zumkeller, NAMlab (https://github.com/NAMlab/DeepCRE/tree/main/moca_blue)
- **TF-MoDISco:** Shrikumar et al. (2020). Nature Machine Intelligence

---

**Note:** This pipeline is under active development for *Vitis vinifera* drought stress analysis