# moca_blue: Motif Characterisation & Annotation from Deep Learning Feature Enrichment

*A comprehensive toolkit for extracting, processing, and analysing DNA motifs derived from deep learning model feature extraction, with a focus on plant gene regulatory networks.*

**Developed by:** Dr. Simon M. Zumkeller  
**Enhanced by:** Larysha Rothmann   
**Context:** Deep learning drought adaptation in *Vitis vinifera's* gene regulatory network

---

## Overview

The moca_blue suite provides a complete pipeline for converting deep learning model outputs (specifically MoDisco HDF5 files) into interpretable, standardised motif representations. This toolkit bridges the gap between complex neural network feature extraction and traditional bioinformatics motif analysis workflows.

The pipeline extracts three complementary matrix representations from MoDisco results:
- **Position Frequency Matrices (PFMs)**: Raw nucleotide frequencies across seqlets
- **Position Weight Matrices (PWMs)**: Log-odds scores relative to background frequencies  
- **Contribution Weight Matrices (CWMs)**: Model-derived importance scores scaled by seqlet abundance

These matrices capture different aspects of regulatory element function and can be used together to provide a comprehensive view of *cis*-regulatory logic in plant genomes.

## Scientific Context

This toolkit was developed as part of the DeepCRE project investigating drought adaptation mechanisms in *Vitis vinifera* gene regulatory networks. The work builds upon recent advances in applying deep learning to plant genomics, particularly the DeepCRE framework by Peleke et al. (2024), which predicts gene expression classes from promoter sequences.

### Why Multiple Matrix Types Matter

Understanding gene regulation requires multiple complementary perspectives on DNA sequence-function relationships:

**Position Frequency Matrices (PFMs)** represent the observed nucleotide frequencies at each position within aligned seqlets. These provide the most direct representation of what sequences were actually found by the model, weighted by how many examples (seqlets) support each pattern. PFMs are essential for understanding the raw binding preferences and sequence diversity captured by the model.

**Position Weight Matrices (PWMs)** convert these frequencies into information-theoretic measures by calculating log-odds scores relative to background nucleotide frequencies. This transformation reveals which positions and nucleotides are most informative for distinguishing functional binding sites from random sequence. PWMs are particularly valuable for motif comparison, clustering, and database searching because they normalise for genome composition differences.

**Contribution Weight Matrices (CWMs)** preserve the actual importance scores learned by the deep learning model, reflecting not just binding preferences but the functional impact of each position on gene expression predictions. Unlike traditional motifs, CWMs capture the model's learned understanding of how sequence variation translates to regulatory outcome. The scaling approach used here weights these scores by seqlet abundance, ensuring that motifs are ranked by their aggregate importance rather than being dominated by rare but extreme examples.

### Integration with Gene Regulatory Networks

The extracted motifs serve as building blocks for reconstructing gene regulatory networks (GRNs) under environmental stress conditions. When combined with expression data, chromatin accessibility measurements (ATAC-seq), and transcription factor binding profiles (DAP-seq), these motifs enable:

- Identification of water stress-responsive regulatory modules
- Prediction of transcription factor binding sites across cultivars
- Assessment of regulatory impact from genomic variants (SNPs, indels)
- Integration with co-expression networks to validate functional relationships

This multi-layered approach is essential for understanding how *Vitis vinifera* cultivars adapt to water stress through regulatory network modifications.

---

## Installation & Dependencies

### R Dependencies
```r
# Core packages
install.packages(c("dplyr", "purrr"))

# Bioconductor packages  
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")
```

### System Requirements
- R >= 4.0.0
- Sufficient memory for loading large HDF5 files (typically 4-16 GB depending on model size)
- For large-scale analyses: multi-core system recommended

---

## Pipeline Architecture

The moca_blue pipeline follows a modular architecture designed for reproducibility and extensibility:

```
INPUT DIRECTORY/
└── 0MOTIFS/
    ├── modisco.hdf5         # MoDisco feature extraction results
    └── ref_seq/             # Reference genome data
        ├── *.fasta          # Genome sequences
        ├── *.gff            # Gene annotations  
        └── meta-data/       # Experimental metadata

PROCESSING MODULES/
├── mo_nom/                  # Motif extraction & nomenclature
├── mo_range/                # Positional preference analysis
├── mo_clu/                  # Motif clustering & comparison
└── mo_proj/                 # Genome mapping & interpretation

OUTPUT/
├── JASPAR format files      # Standardised motif representations
├── Summary statistics       # Motif importance metrics
└── Visualisation plots      # Logo plots & diagnostic figures
```

### Data Flow

The pipeline processes MoDisco HDF5 files through several key stages:

1. **Extraction**: Raw matrices are extracted from the nested HDF5 structure, handling both forward and reverse strand patterns across multiple metaclusters
2. **Processing**: Different mathematical transformations are applied based on matrix type (frequency scaling for PFMs, log-odds calculation for PWMs, importance weighting for CWMs)
3. **Nomenclature**: Systematic naming ensures traceability back to original patterns while encoding key metadata (species, model, strand, seqlet count)
4. **Quality Control**: Summary statistics help identify high-confidence motifs and potential processing artifacts
5. **Standardisation**: Output in JASPAR format ensures compatibility with downstream motif analysis tools

---

## Usage

### Basic Usage

The enhanced motif extraction script processes all three matrix types in a single run:

```bash
# Using defaults (recommended for most analyses)
Rscript get_motifs.R

# With custom parameters
Rscript get_motifs.R vitis ssr 20241201 ./data/modisco/vitis_drought_modisco.hdf5 ./results/

# Command line arguments:
# 1. Species identifier (default: "vitis")  
# 2. Model identifier (default: "ssr")
# 3. Date string (default: current date)
# 4. Input HDF5 file path (auto-detected if not specified)
# 5. Output directory (default: "./results/modisco/stats")
```

### Advanced Configuration

#### Custom Background Frequencies for PWM Calculation

The default PWM calculation assumes equal nucleotide frequencies (0.25 each), but plant genomes often show AT/GC bias that should be accounted for:

```r
# Example: AT-rich genome (modify in script)
background_freq <- c(0.29, 0.21, 0.21, 0.29)  # A, C, G, T

# Genome-derived backgrounds (recommended approach)
# Calculate from your reference FASTA:
# A_freq <- count_A / total_bases, etc.
```

For accurate PWM calculations across different cultivars or species, we recommend deriving background frequencies from the actual reference genome being used. This ensures that motif significance scores reflect genuine binding preferences rather than genome composition effects.

#### Handling Variable-Length Motifs

Unlike the original scripts, the enhanced version automatically handles variable-length motifs without hard-coded dimension assumptions. However, for downstream compatibility, you may want to pad or trim motifs to standard lengths:

```r
# This functionality can be added to the processing functions if needed
# standardise_motif_length <- function(matrix, target_length = 15) { ... }
```

---

## Output Files

### JASPAR Format Motifs

Three files are generated, each containing motifs in standard JASPAR format:

- `rdf5_PFM_pattern[SPEC][MODEL]_pfm-motifs.jaspar`: Position frequency matrices
- `rdf5_epm[SPEC][MODEL]_pwm-motifs.jaspar`: Position weight matrices  
- `rdf5_epm[SPEC][MODEL]_cwm-motifs.jaspar`: Contribution weight matrices

### Summary Statistics

The comprehensive summary table (`[DATE]_[SPEC][MODEL]_motif_summary.csv`) includes:

- **motif**: Systematic motif identifier
- **matrix_type**: PFM, PWM, or CWM
- **score_sum**: Total importance/frequency score
- **score_max/min**: Range of values within the motif
- **score_range**: Dynamic range of the motif

This summary enables rapid identification of high-confidence motifs and comparison of motif strength across conditions or cultivars.

### Motif Nomenclature System

The systematic naming convention encodes essential metadata:

```
epm_<species>_<model>_<pattern_id>_<metacluster>_<strand>_<seqlet_count>

Example: epm_vitis_ssr_0_metacluster_0_F_1247
- Species: vitis (Vitis vinifera)
- Model: ssr (single-species reference)  
- Pattern: 0 (first pattern discovered)
- Metacluster: 0 (associated with high expression)
- Strand: F (forward)
- Seqlets: 1247 (number of supporting examples)
```

This naming system ensures full traceability and enables systematic comparison across different model runs, conditions, or species.

---

## Integration with Downstream Analysis

### Motif Database Searching

The JASPAR format output is compatible with major motif analysis tools:

```bash
# MEME Suite
tomtom motifs.jaspar JASPAR2022_CORE_plants_non-redundant_pfms_meme.txt

# FIMO for genome-wide scanning
fimo --oc ./fimo_results motifs.jaspar genome.fasta

# Comparison with PlantTFDB
# Import JASPAR files into R for programmatic analysis
```

### Gene Regulatory Network Construction

The extracted motifs serve as foundational components for GRN reconstruction:

1. **Motif-to-TF Assignment**: Compare extracted motifs to known plant transcription factor binding profiles using tools like TOMTOM or STAMP
2. **Genome Scanning**: Use FIMO or similar tools to identify putative binding sites across gene promoters
3. **Expression Integration**: Combine binding predictions with RNA-seq data to infer active regulatory relationships
4. **Network Validation**: Use chromatin accessibility (ATAC-seq) and binding data (ChIP-seq/DAP-seq) to validate predicted interactions

### Cultivar-Specific Analysis

For multi-cultivar studies, the pipeline supports comparative analysis:

```r
# Compare motif conservation across cultivars
# This analysis framework can be extended for SNP impact assessment:
# 1. Extract motifs from reference cultivar
# 2. Apply variants from other cultivars to regulatory sequences  
# 3. Re-score motifs with modified sequences
# 4. Quantify predicted regulatory impact
```

---

## Troubleshooting

### Common Issues

**Memory errors with large HDF5 files**: Increase available RAM or process metaclusters separately. The script can be modified to handle metaclusters iteratively rather than loading all data simultaneously.

**Missing patterns in output**: Check that the HDF5 file structure matches expected MoDisco format. Some patterns may lack contribution scores or sequence data, triggering warnings but not errors.

**Unexpected motif dimensions**: Variable-length motifs are handled automatically, but downstream tools may require consistent motif lengths. Consider standardisation if compatibility issues arise.

**Background frequency optimization for PWMs**: If PWM values seem extreme or unrealistic, verify that background frequencies match your genome composition. Consider calculating empirical backgrounds from your reference FASTA.

### Validation Approaches

To ensure motif extraction quality:

1. **Visual inspection**: Generate sequence logos for high-scoring motifs and compare to known plant TF binding profiles
2. **Cross-validation**: Extract motifs from subsets of your data and assess consistency  
3. **Literature comparison**: Compare extracted motifs to published binding sites for known stress-responsive transcription factors
4. **Functional validation**: Test predicted binding sites using reporter assays or ChIP-qPCR

---

## Citation

If you use moca_blue in your research, please cite:

- **Original moca_blue framework**: Zumkeller, S.M. et al. (in preparation)
- **DeepCRE methodology**: Peleke et al. (2024). Deep learning of *cis*-regulatory elements from sequence explains widespread variation in gene expression across plant species and tissues. *Plant Cell* 36: 482-505.
- **MoDisco algorithm**: Shrikumar, A., et al. (2020). Technical note on transcription factor motif discovery from importance scores (TF-MoDISco) version 0.5.6.5. *arXiv preprint* arXiv:1811.00416.

---

## Development Roadmap

### Planned Enhancements

- **Automated background frequency calculation** from reference FASTA files
- **Motif clustering and family assignment** using sequence similarity metrics
- **Integration with PlantTFDB** for automated transcription factor annotation  
- **Batch processing capabilities** for large-scale comparative studies
- **Web interface development** for interactive motif exploration and visualisation

### Contributing

This toolkit is part of an active research project. We welcome contributions, bug reports, and feature requests. Please contact the development team or submit issues through the project repository.

For questions related to drought adaptation research in grapevine or deep learning applications in plant genomics, please see our related publications or contact the research team directly.

---

## Acknowledgments

This work is supported by the Deep learning drought adaptation in *Vitis vinifera's* gene regulatory network project at Stellenbosch University. We thank the broader plant genomics community for developing the foundational tools and databases that make this research possible.

Special acknowledgment to the developers of MoDisco, the MEME Suite, and the R/Bioconductor ecosystem, whose tools form the backbone of this pipeline.