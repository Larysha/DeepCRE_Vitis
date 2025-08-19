# moca_blue_reimplemented

**MOtif Characterization & Annotation from DEEP LEARNING feature enrichment**

## Attribution & Adaptation

**Original pipeline:** All scripts originally developed by [Dr Simon M. Zumkeller](https://github.com/NAMlab/DeepCRE/tree/main/moca_blue) for the DeepCRE model published in [Peleke et al. 2024, Nature Communications](https://www.nature.com/articles/s41467-024-47744-0).

**This repository:** Contains refactored and organised versions of the original moca_blue scripts, adapted for water stress studies in *Vitis vinifera* gene regulatory networks.

**Scope of changes:** improved error handling and reproducibility and adaptation for *Vitis* genomic resources. Core algorithms and methodology remain unchanged from the original implementation for now.

---

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/)

## Overview

moca_blue integrates with the [DeepCRE](https://github.com/NAMlab/DeepCRE) framework to process TF-MoDISco outputs from deep learning models trained on gene expression data. The pipeline enables systematic characterisation of regulatory motifs discovered through neural network interpretation, bridging the gap between model features and biological understanding.

**Input**: HDF5 files from TF-MoDISco feature extraction  
**Output**: Biologically interpretable regulatory modules with genomic annotations

### Key Features

- **Motif extraction** from neural network attribution scores
- **Positional preference analysis** relative to gene features
- **Phylogenetic clustering** using established similarity metrics
- **Genomic mapping** with statistical significance testing
- **Cross-species compatibility** and modular architecture

## Installation

### Prerequisites

**System Requirements:**
- R ≥ 4.0
- Python ≥ 3.7
- BLAMM motif scanner
- samtools ≥ 1.9

**R Dependencies:**
```r
install.packages(c("dplyr", "purrr", "ape", "ggtree"))
BiocManager::install(c("TFBSTools", "universalmotif", "motifStack", "rhdf5"))
```

**External Tools:**
```bash
# BLAMM installation (required for motif mapping)
git clone https://github.com/biointec/blamm.git
cd blamm && make
# Ensure blamm is in PATH

# samtools (for sequence extraction)
conda install samtools
```

### Setup

1. Clone the repository:
```bash
git clone https://github.com/Larysha/DeepCRE_Vitis.git
cd DeepCRE_Vitis/src
```

## Quick Start

### Basic Pipeline

```bash
# Set working directory to project root
cd /path/to/deepcre_project/src

# 1. Extract motifs from TF-MoDISco output
    # ssr stands for single species reference trained model (in this case vitis)
Rscript moca_blue/mo_nom/get_matrices.R vitis ssr

# 2. Analyse positional preferences
Rscript moca_blue/mo_range/get_seql_per_pattern.R
Rscript moca_blue/mo_range/motif_ranges_meta.R

# 3. Cluster similar motifs
Rscript moca_blue/mo_clu/cluster_motifs.R vitis ssr

# 4. Extract reference sequences
bash moca_blue/ref_seq/extract_range_to_fasta.sh vitis_vinifera_PN40024

# 5. Map motifs to genome
bash moca_blue/ref_seq/blamm_motif_mapping.sh vitis_vinifera_PN40024

# 6. Filter and analyse mappings
Rscript moca_blue/mo_proj/occ_filter_v1.1.R
```

### Project Structure

```
project_root/
├── moca_blue/                 # Pipeline modules
│   ├── mo_nom/               # Motif extraction & nomenclature
│   ├── mo_imp/               # Saliency analysis
│   ├── mo_range/             # Positional preferences
│   ├── mo_clu/               # Motif clustering
│   ├── ref_seq/              # Reference mapping
│   ├── mo_proj/              # Genomic projection
│   └── utils.R               # Shared utilities
├── results
|   ├── moca_blue/            # Pipeline outputs
|   ├── modisco/              # TF-MoDISco inputs
|   ├── shap/                 # SHAP interpretation files
└── genome/                   # Reference sequences
```

## Pipeline Modules

### 1. mo_nom: Motif Extraction

Extracts position frequency matrices (PFMs), position weight matrices (PWMs), and contribution weight matrices (CWMs) from TF-MoDISco HDF5 outputs.

**Key Script:** `get_matrices.R`

**Matrix Types:**
- **PFM**: Raw count matrices for BLAMM input
- **PWM**: Frequency matrices (legacy compatibility)
- **CWM**: Neural network attribution scores for clustering

**Output Format:**
```
epm_<species>_<model>_<pattern_id>_<metacluster>_<strand>_<seqlet_count>
```

### 2. mo_range: Positional Analysis

Characterises spatial distribution of motifs relative to transcription start sites (TSS) and transcription termination sites (TTS).

**Scripts:**
- `get_seql_per_pattern.R`: Extract seqlet coordinates
- `motif_ranges_meta.R`: Analyse positional preferences

**Windows Analysed:**
- TSS-proximal: 1–1500 bp
- TTS-proximal: 1520–3000 bp

### 3. mo_clu: Motif Clustering

Groups motifs using the Sandelin-Wassermann similarity algorithm to identify redundant patterns and build phylogenetic relationships.

**Key Script:** `cluster_motifs.R`

**Features:**
- Sandelin-Wassermann 2004 similarity metric
- 95th percentile threshold for redundancy detection
- Phylogenetic tree construction
- Information content analysis

### 4. ref_seq: Reference Mapping

Extracts gene sequences and maps discovered motifs using BLAMM with statistical significance testing.

**Scripts:**
- `extract_range_to_fasta.sh`: Sequence extraction with flanking regions
- `blamm_motif_mapping.sh`: Statistical motif mapping

**Requirements:**
- Genome FASTA with samtools index
- Gene annotation (GFF3/GTF format)
- PFM matrices in JASPAR format

### 5. mo_proj: Genomic Projection

Filters mapping results and provides performance analysis of discovered regulatory modules.

**Key Features:**
- Statistical significance filtering
- Gene-motif association networks
- Expression prediction evaluation
- Functional enrichment analysis

## Input Requirements

### Essential Files

1. **TF-MoDISco output** (HDF5 format)
   - Generated by DeepCRE interpretation pipeline
   - Contains contribution scores and seqlet coordinates

2. **Reference genome** (FASTA format)
   - Indexed with samtools faidx
   - Chromosome naming consistent with annotation

3. **Gene annotation** (GFF3/GTF format)
   - Standard gene model features
   - Coordinate system matching genome assembly

4. **Expression metadata** (CSV format)
   - Gene identifiers matching annotation
   - Expression class labels (high/low)

### Naming Conventions

**HDF5 files:** `<species>_<model>_modisco.hdf5`  
**Genome files:** `<species>_<assembly>_dna.fa`  
**Annotation:** `<species>_<assembly>.gff3`

## Output Interpretation

### Matrix Types

**CWM (Contribution Weight Matrix):**
- Represents neural network importance scores
- Used for clustering analysis
- Scaled by seqlet count for comparability

**PFM (Position Frequency Matrix):**
- Raw nucleotide counts
- Required for statistical motif mapping
- Input format for BLAMM scanner

### Motif Nomenclature

```
epm_vitis_ssr_0_0_F_3432
│   │     │   │ | │  │
│   │     │   │ | │  └─ Supporting evidence (no. "seqlets")
│   │     │   │ | └─ Strand orientation  
│   │     │   │ |─ Expression association
│   │     │   └─ Pattern identifier
│   │     └─ Model name
│   └─ Species identifier
└─ Expression Pattern Module prefix
```

**Metaclusters:**
- `metacluster_0`: High expression associated
- `metacluster_1`: Low expression associated

### Performance Metrics

**Motif Quality Indicators:**
- `contrib_score_sum`: Total neural network contribution
- `seqlet_count`: Number of supporting instances
- Information content (bits)
- Consensus sequence confidence

## Advanced Usage

### Custom Species

```bash
# Modify parameters for new species
SPECIES="arabidopsis"
MODEL="msr" # note I haven't tested the msr model yet

Rscript moca_blue/mo_nom/get_matrices.R $SPECIES $MODEL
# Continue with remaining pipeline steps...
```

### Parameter Tuning

**BLAMM p-value threshold:**
```bash
bash moca_blue/ref_seq/blamm_motif_mapping.sh species_name 0.001 1500
```

**Clustering similarity threshold:**
```r
# Edit cluster_motifs.R
similarity_threshold <- 0.90  # Default: 0.95
```

### Integration with DeepCRE

moca_blue seamlessly integrates with the DeepCRE pipeline:

1. Train models using `train_models.py`
2. Generate interpretations with `deepcre_interpret.py`
3. Run TF-MoDISco feature extraction
4. Process outputs with moca_blue pipeline

## Troubleshooting

### Common Issues

**BLAMM parsing errors:**
- Verify JASPAR format compliance
- Check for non-integer values in PFM matrices
- Ensure tab-separated headers

**Memory limitations:**
- Use file chunking for large BLAMM outputs
- Monitor disk usage during motif mapping
- Consider conservative p-value thresholds

**Coordinate mismatches:**
- Verify consistent genome assemblies
- Check chromosome naming conventions
- Validate GFF3 coordinate systems

### Performance Optimisation

**Large datasets:**
```bash
# Process BLAMM output in chunks
bash moca_blue/ref_seq/split_file.sh occurrences.txt 1000000
```

## Citation

If you use moca_blue in your research, please cite:

- **Original DeepCRE framework & moca_blue pipeline:** Peleke et al. (2024). "Predictive models of gene expression reveal mechanisms of transcriptional regulation in rice and tomato." Nature Communications. https://doi.org/10.1038/s41467-024-47744-0
- **Original moca_blue scripts:** Dr Simon M. Zumkeller, NAMlab (https://github.com/NAMlab/DeepCRE/tree/main/moca_blue)
- **TF-MoDISco:** Shrikumar et al. (2020). Nature Machine Intelligence




## Acknowledgements

- **Original moca_blue scripts:** Dr Simon M. Zumkeller, NAMlab
- **DeepCRE framework:** Peleke et al., Nature Communications (2024)
- **BLAMM motif scanner:** Biointec Lab, Ghent University
- **TF-MoDISco:** Kundaje Lab, Stanford University

---

**Note:** This pipeline is under development. Please check for updates and report any issues through the GitHub repository.
**Contact**: larysha@sun.ac.za
