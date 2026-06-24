🧬 Exome-Seq & Metabolomics Integration Analysis for Early-Onset Schizophrenia
https://opensource.org/licenses/MIT
https://www.r-project.org/
https://www.gnu.org/software/bash/
Reproducible analysis pipelines for integrative exome sequencing (WES) and LC-MS/MS metabolomics in early-onset schizophrenia (EOS).
This repository contains the complete computational workflow used in our study, encompassing variant calling & CNV detection, differential metabolite profiling, and multi-omics association analyses.
📋 Table of Contents
Overview
Repository Structure
Workflow
Requirements
Usage
Data Availability
Citation
Contact
🔬 Overview
Schizophrenia is a complex psychiatric disorder with strong genetic and metabolic components. This repository implements a three-module integrative framework:
Genomic Module: SNV/CNV detection and functional annotation from whole-exome sequencing (WES) data.
Metabolomic Module: Peak picking, quality control, and differential abundance analysis from LC-MS/MS data.
Integrative Module: Variant-symptom, metabolite-symptom, and variant-metabolite association analyses, culminating in a protein-metabolite network.
📁 Repository Structure
plain
复制
.
├── 1.1.Variant_detection.sh              # WES alignment & variant calling (BWA-GATK)
├── 1.2.VariantQC.sh                      # Variant quality control & Hail filtering
├── 1.3.Variant_function_prediction_and_OR_filtering.sh  # VEP annotation & OR filtering
├── 1.4.CNV_detection_and_filtering.sh  # XHMM-based CNV detection
├── 2.LC-MSMS.data_process_and_statistics.R              # XCMS → Masscleaner → Statistics
├── 3.1.Variant_symptom_assoc.R           # Gene variant × PANSS association
├── 3.2.Metabolite_symptom_corr.R       # Metabolite × PANSS correlation
├── 3.3.Variant_metabolite_assoc.R      # Variant × metabolite association
├── 3.4.Metabolite_metabolite_correlation.R              # Metabolite correlation network
├── 3.5.Protein_metabolite_network.R    # Protein-metabolite network construction
├── README.md
└── (optional) config.R / config.sh     # Shared variables (recommended)
Naming convention: Scripts are numbered by module (1.* = Genomic, 2.* = Metabolomic, 3.* = Integrative).
🔄 Workflow
Mermaid
全屏 
下载 
复制
代码
预览
🔗 Integrative Module

⚗️ Metabolomic Module

🧬 Genomic Module

BWA MEM

Picard

GATK

CombineGVCFs

VQSR

VEP

R/OR filter

XHMM

XCMS

Masscleaner

massstat

3.1

3.2

3.3

3.4

3.5

Raw FASTQ
Aligned BAM
Mark Duplicates
G VCF
Genotyped VCF
Filtered VCF
Annotated Variants
Pathogenic Variants
CNV Calls
Raw .mzML
Peak Picking
QC & Normalization
Differential Metabolites
Variant × Symptom
Metabolite × Symptom
Variant × Metabolite
Metabolite Correlation
Network Assembly
Protein-MetaboliteNetwork
🛠️ Requirements
System Dependencies
表格
Tool	Version	Purpose
BWA	≥0.7.17	Sequence alignment
Picard	≥2.26	BAM processing
GATK	≥4.2	Variant calling & VQSR
Bedtools	≥2.30	Coverage analysis
XHMM	≥1.0	CNV detection
VEP	≥104	Variant annotation
bcftools	≥1.14	VCF manipulation
Hail	≥0.2	Large-scale QC
Python	≥3.8	Hail scripting
R Packages
r
复制
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
Note: Please modify all hardcoded paths (~/EOSCZ/...) in each script to match your local environment before execution.
🚀 Usage
Genomic Pipeline
bash
复制
# Step 1: Align and call variants (run separately for SCZ and NOR groups)
bash 1.1.Variant_detection.sh

# Step 2: Merge groups and perform QC
bash 1.2.VariantQC.sh

# Step 3: Annotate and filter by OR > 1
bash 1.3.Variant_function_prediction_and_OR_filtering.sh

# Step 4: Detect CNVs
bash 1.4.CNV_detection_and_filtering.sh
Metabolomic Pipeline
bash
复制
# Phase 1–3: Peak picking, QC, and statistics (positive mode example)
Rscript 2.LC-MSMS.data_process_and_statistics.R positive
Association Analyses
r
复制
# Run each association script in R or via command line
Rscript 3.1.Variant_symptom_assoc.R
Rscript 3.2.Metabolite_symptom_corr.R
Rscript 3.3.Variant_metabolite_assoc.R
Rscript 3.4.Metabolite_metabolite_correlation.R
Rscript 3.5.Protein_metabolite_network.R
⚠️ Important: Scripts 3.1–3.5 depend on outputs from 1.* and 2.*. Ensure upstream steps are completed before running integrative analyses.
📊 Data Availability
Raw sequencing and metabolomics data are available under controlled access. Processed variant lists and metabolite intensity tables are provided as supplementary material to the publication (see Citation).
📖 Citation
If you use this code in your research, please cite:
Yang, X., et al. (2025). Integrative exome sequencing and metabolite profiling reveals novel genetic-metabolic interactions in early-onset schizophrenia. [Journal Name], Volume, pp–pp. DOI: 10.xxxx/xxxxx
(Citation will be updated upon publication.)
📬 Contact
For questions, bug reports, or collaboration inquiries, please:
Open an issue on this GitHub repository
Email: yxx1123@smu.edu.cn
⚖️ License
This project is licensed under the MIT License — see the LICENSE file for details.
<p align="center">
  <i>Built with 🧬 + 🧪 + 💻 for reproducible psychiatric genomics.</i>
</p>
