# Exome-Seq and Metabolomics Integration Analysis for Early-Onset Schizophrenia

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![Bash](https://img.shields.io/badge/Bash-%3E%3D4.0-green.svg)](https://www.gnu.org/software/bash/)

> **Reproducible analysis pipelines for integrative exome sequencing (WES) and LC-MS/MS metabolomics in early-onset schizophrenia (EOS).**

This repository contains the complete computational workflow used in our study, encompassing **variant calling & CNV detection**, **differential metabolite profiling**, and **multi-omics association analyses**.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Workflow](#workflow)
- [Requirements](#requirements)
- [Usage](#usage)
- [Data Availability](#data-availability)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

Schizophrenia is a complex psychiatric disorder with strong genetic and metabolic components. This repository implements a three-module integrative framework:

1. **Genomic Module**: SNV/CNV detection and functional annotation from whole-exome sequencing (WES) data.
2. **Metabolomic Module**: Peak picking, quality control, and differential abundance analysis from LC-MS/MS data.
3. **Integrative Module**: Variant-symptom, metabolite-symptom, and variant-metabolite association analyses, culminating in a protein-metabolite network.

---

## Repository Structure

| File | Description | Module |
|------|-------------|--------|
| `1.1.Variant_detection.sh` | WES alignment & variant calling (BWA + GATK) | Genomic |
| `1.2.VariantQC.sh` | Variant quality control & Hail filtering | Genomic |
| `1.3.Variant_function_prediction_and_OR_filtering.sh` | VEP annotation & OR filtering | Genomic |
| `1.4.CNV_detection_and_filtering.sh` | XHMM-based CNV detection | Genomic |
| `2.LC-MSMS.data_process_and_statistics.R` | XCMS peak picking -> Masscleaner QC -> Statistics | Metabolomic |
| `3.1.Variant_symptom_assoc.R` | Gene variant x PANSS association | Integrative |
| `3.2.Metabolite_symptom_corr.R` | Metabolite x PANSS correlation | Integrative |
| `3.3.Variant_metabolite_assoc.R` | Variant x metabolite association | Integrative |
| `3.4.Metabolite_metabolite_correlation.R` | Metabolite correlation network | Integrative |
| `3.5.Protein_metabolite_network.R` | Protein-metabolite network construction | Integrative |

---

## Workflow

![Workflow Diagram](workflow_diagram.svg)

**Three modules:**

| Module | Input | Key Tools | Output |
|--------|-------|-----------|--------|
| **Genomic** | Raw FASTQ (WES) | BWA, Picard, GATK, XHMM, VEP | Filtered variants & CNVs |
| **Metabolomic** | Raw .mzML (LC-MS/MS) | XCMS, Masscleaner, massstat | Differential metabolites (DEMs) |
| **Integrative** | Variants + DEMs + Clinical | R (tidyverse, limma, igraph) | Association networks |

---

## Requirements

### System Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| [BWA](http://bio-bwa.sourceforge.net/) | >= 0.7.17 | Sequence alignment |
| [Picard](https://broadinstitute.github.io/picard/) | >= 2.26 | BAM processing |
| [GATK](https://gatk.broadinstitute.org/) | >= 4.2 | Variant calling & VQSR |
| [Bedtools](https://bedtools.readthedocs.io/) | >= 2.30 | Coverage analysis |
| [XHMM](http://atgu.mgh.harvard.edu/xhmm/) | >= 1.0 | CNV detection |
| [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) | >= 104 | Variant annotation |
| [bcftools](https://samtools.github.io/bcftools/) | >= 1.14 | VCF manipulation |
| [Hail](https://hail.is/) | >= 0.2 | Large-scale QC |
| [Python](https://www.python.org/) | >= 3.8 | Hail scripting |

### R Packages

```r
# Core data manipulation
install.packages(c("tidyverse", "magrittr", "broom", "scales"))

# Omics analysis
BiocManager::install(c("xcms", "MSnbase", "limma", "SummarizedExperiment", "S4Vectors"))

# Tidymass ecosystem
remotes::install_github("tidymass/masscleaner")
remotes::install_github("tidymass/massdataset")
remotes::install_github("tidymass/massstat")
remotes::install_github("tidymass/massqc")

# Visualization & statistics
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "corrplot", 
                   "ggrepel", "ggsci", "VennDiagram", "igraph", "ggsignif", 
                   "cowplot", "car", "pls", "ropls", "Hmisc", "corrr"))
```

> **Note:** Please modify all hardcoded paths (`~/EOSCZ/...`) in each script to match your local environment before execution.

---

## Usage

### Genomic Pipeline

```bash
# Step 1: Align and call variants (run separately for SCZ and NOR groups)
bash 1.1.Variant_detection.sh

# Step 2: Merge groups and perform QC
bash 1.2.VariantQC.sh

# Step 3: Annotate and filter by OR > 1
bash 1.3.Variant_function_prediction_and_OR_filtering.sh

# Step 4: Detect CNVs
bash 1.4.CNV_detection_and_filtering.sh
```

### Metabolomic Pipeline

```bash
# Phase 1-3: Peak picking, QC, and statistics (positive mode example)
Rscript 2.LC-MSMS.data_process_and_statistics.R positive
```

### Association Analyses

```r
# Run each association script in R or via command line
Rscript 3.1.Variant_symptom_assoc.R
Rscript 3.2.Metabolite_symptom_corr.R
Rscript 3.3.Variant_metabolite_assoc.R
Rscript 3.4.Metabolite_metabolite_correlation.R
Rscript 3.5.Protein_metabolite_network.R
```

> **Important:** Scripts `3.1`-`3.5` depend on outputs from `1.*` and `2.*`. Ensure upstream steps are completed before running integrative analyses.

---

## Data Availability

Raw sequencing and metabolomics data are available under controlled access. Processed variant lists and metabolite intensity tables are provided as supplementary material to the publication (see [Citation](#citation)).

---

## Citation

If you use this code in your research, please cite:

> **Yang, X.**, et al. (2026). *Integrative exome sequencing and metabolite profiling reveals novel genetic-metabolic interactions in early-onset schizophrenia.* [Journal Name], [Volume](Issue), pp-pp. DOI: [10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx)

*(Citation will be updated upon publication.)*

---

## Contact

For questions, bug reports, or collaboration inquiries, please:

- **Open an issue** on this GitHub repository
- **Email:** [yxx1123@smu.edu.cn](mailto:yxx1123@smu.edu.cn)

---

## License

This project is licensed under the MIT License.

---

<p align="center">
  <i>Built for reproducible psychiatric genomics.</i>
</p>
