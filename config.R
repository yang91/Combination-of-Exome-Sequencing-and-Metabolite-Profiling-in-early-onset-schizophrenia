# =============================================================================
# config.R — Shared Configuration for EOSCZ Analysis Pipeline
# =============================================================================
# Source this file at the top of every analysis script:
#   source("config.R")
#
# To customize paths without modifying this file, set environment variables
# before sourcing (e.g., in .Renviron or at the R console):
#   Sys.setenv(EOSCZ_HOME = "/your/custom/path")
# =============================================================================

#-------------------------------------------------------------------------------
# 0. Project Root (customizable via environment variable)
#-------------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("EOSCZ_HOME", unset = "~/EOSCZ")

#-------------------------------------------------------------------------------
# 1. Genomic Analysis Paths
#-------------------------------------------------------------------------------
EXOME_DIR    <- file.path(PROJECT_ROOT, "ExomeSeq-Result")
VAR_DIR      <- file.path(EXOME_DIR, "variant")
BWA_DIR      <- file.path(EXOME_DIR, "bwa")
PICARD_DIR   <- file.path(EXOME_DIR, "picard")
GATK_DIR     <- file.path(EXOME_DIR, "GATK")
REF_INDEX_DIR <- file.path(PROJECT_ROOT, "Ref_and_Index")

# Reference files (adjust filenames to match your local setup)
REF_FA       <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "Homo_sapiens_assembly38.fa")
DBSNP        <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "dbsnp_144.hg38.withchr.vcf")
KNOWN_INDEL  <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "Homo_sapiens_assembly38.known_indels.vcf")
GOLDEN_INDEL <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "Mills_and_1000G_gold_standard.indels.hg38.vcf")
HAPMAP       <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "hapmap_3.3.hg38.vcf.gz")
OMNI         <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "1000G_omni2.5.hg38.vcf.gz")
OKG          <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "1000G_phase1.snps.high_confidence.hg38.vcf.gz")
AXIOM        <- file.path(REF_INDEX_DIR, "GATK-hg38bundle", "Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz")
BWA_INDEX    <- file.path(REF_INDEX_DIR, "BWA_INDEX")
BED_FILE     <- file.path(REF_INDEX_DIR, "T086V4_MT.merged.success.liftover.to.hg38.bed")

# VEP
VEP_CACHE    <- file.path(REF_INDEX_DIR, ".vep")

# Gene conversion table
GENE_CONV_FILE <- file.path(REF_INDEX_DIR, "combine.ensemb_104.orgHsegdb.ncbi_name_searching.manually_checked.txt")

# Gene score file
GENE_SCORE_FILE <- file.path(EXOME_DIR, "gene_score.high-risk.txt")

#-------------------------------------------------------------------------------
# 2. Metabolomic Analysis Paths
#-------------------------------------------------------------------------------
META_ROOT    <- file.path(PROJECT_ROOT, "MetaSeq-Result")
SEQDATA_DIR  <- file.path(PROJECT_ROOT, "MetaSeq-Data", "convert_data")
XCMS_OUT     <- file.path(META_ROOT, "XCMS")
MASS_DIR     <- file.path(META_ROOT, "tidymass")
STAT_DIR     <- file.path(META_ROOT, "statistics")
META_STAT_DIR <- STAT_DIR  # alias for backward compatibility

# Annotation files (adjust to match your local annotation file names)
POS_ANNT_FILE <- file.path(SEQDATA_DIR, "positive-all-identification.csv")
NEG_ANNT_FILE <- file.path(SEQDATA_DIR, "negative-all-identification.csv")

#-------------------------------------------------------------------------------
# 3. Clinical & Network Paths
#-------------------------------------------------------------------------------
CLIN_DIR     <- file.path(PROJECT_ROOT, "Clinical_Info")
NET_DIR      <- file.path(PROJECT_ROOT, "Network_file")

# Clinical data files
SPINFO_FILE  <- file.path(CLIN_DIR, "sample_information.csv")
SAMPLE_LIST  <- file.path(PROJECT_ROOT, "sample_list.txt")

#-------------------------------------------------------------------------------
# 4. Variant Filtering Files (outputs from upstream scripts)
#-------------------------------------------------------------------------------
VEP_ANNOT_FILE         <- file.path(VAR_DIR, "Maf_0.05.potential_harmful.ORgt1_variants.txt")
DISEASE_VARIANT_FILE   <- file.path(VAR_DIR, "disease_causing_ORgt1_variants_in_GSgt4_gene.txt")
VEP_SERVER_FILE        <- file.path(VAR_DIR, "VEP_server.txt")

# Metabolite result files
POS_STAT_FILE <- file.path(STAT_DIR, "positive.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")
NEG_STAT_FILE <- file.path(STAT_DIR, "negative.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")

# Network intermediate files
VAR_MET_ASSOC_FILE <- file.path(NET_DIR, "Variant-metabolite.association.RData")
DEM_COR_FILE       <- file.path(NET_DIR, "DEM_NonDEM.correlation.RData")

# Output file templates
PANSS_ASSOC_OUTPUT     <- file.path(CLIN_DIR, "High_risk_genes_associated_with_PANSS.txt")
METABOLITE_COR_OUTPUT  <- file.path(CLIN_DIR, "MultiVar_LM.DEMs_correlated_with_PANSS.txt")
NETWORK_EDGE_OUTPUT    <- file.path(NET_DIR, "HighRisk_gene.Annotated_metabolites.network_edge.txt")
NETWORK_NODE_OUTPUT    <- file.path(NET_DIR, "HighRisk_gene.Annotated_metabolites.network_node.txt")

#-------------------------------------------------------------------------------
# 5. Analysis Parameters
#-------------------------------------------------------------------------------
RSD_THRESHOLD <- 0.3        # Relative standard deviation threshold for QC
N_CORES       <- 4          # Number of cores for parallel processing
RUN_MODE      <- "both"     # Default: "positive", "negative", or "both"

# XCMS parameters
CENT_WAVE_PARAMS  <- list(
  snthresh = 6,
  ppm = 45,
  peakwidth = c(5, 25),
  mzdiff = 0.01
)
MERGE_PARAMS <- list(ppm = 45)

# Statistical thresholds
P_THRESHOLD    <- 0.05      # Standard p-value threshold
PADJ_THRESHOLD <- 0.05      # FDR-adjusted p-value threshold
MAF_THRESHOLD  <- 0.05      # Minor allele frequency threshold
COR_THRESHOLD  <- 0.3       # Minimum correlation coefficient
Q_THRESHOLD    <- 0.05      # q-value (BH-adjusted) threshold
COR_ABS_THRESHOLD <- 0.9    # Minimum absolute correlation for network edges

# Sample filtering
MIN_CALL_RATE  <- 0.8       # Minimum sample call rate for Hail QC
MIN_DP_MEAN    <- 10        # Minimum mean depth
MAX_DP_MEAN    <- 1000      # Maximum mean depth
MIN_GQ_MEAN    <- 25        # Minimum mean genotype quality

# CADD threshold for SAV scoring
CADD_PHRED_THRESHOLD <- 32

# VEP skip lines (header lines to skip when reading VEP output)
VEP_SKIP_LINES <- 44

#-------------------------------------------------------------------------------
# 6. Color Schemes & Visualization Settings
#-------------------------------------------------------------------------------
GROUP_COLORS <- c(
  SCZ = "#FBABA5",
  NOR = "#33CCD0",
  QC  = "#63AA83",
  case     = "#FBABA5",
  control  = "#33CCD0",
  Patient  = "#FBABA5",
  Normal   = "#33CCD0"
)

MV_PLOT_COLORS <- c("#63AA83", "#33CCD0", "#FBABA5")

# Plot dimensions (inches)
PLOT_WIDTH  <- 8
PLOT_HEIGHT <- 6

# PDF device settings
PDF_WIDTH  <- 10
PDF_HEIGHT <- 8

#-------------------------------------------------------------------------------
# 7. Impact Level Mapping (for VEP annotation)
#-------------------------------------------------------------------------------
EFFECT_LEVELS <- c(
  HIGH    = 4,
  MODERATE = 3,
  LOW     = 2,
  MODIFIER = 1
)

#-------------------------------------------------------------------------------
# 8. Helper: Create output directories if they don't exist
#-------------------------------------------------------------------------------
init_directories <- function() {
  dirs <- c(EXOME_DIR, VAR_DIR, BWA_DIR, PICARD_DIR, GATK_DIR,
            META_ROOT, XCMS_OUT, MASS_DIR, STAT_DIR,
            CLIN_DIR, NET_DIR)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      message("Created directory: ", d)
    }
  }
  invisible(dirs)
}

#-------------------------------------------------------------------------------
# 9. Helper: Validate that critical input files exist
#-------------------------------------------------------------------------------
check_input_files <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing required input file(s):\n  ", paste(missing, collapse = "\n  "))
  }
  message("All required input files found.")
  invisible(TRUE)
}

#-------------------------------------------------------------------------------
# 10. Session info logging for reproducibility
#-------------------------------------------------------------------------------
log_session <- function(log_file = NULL) {
  si <- capture.output(sessionInfo())
  if (!is.null(log_file)) {
    writeLines(si, log_file)
    message("Session info written to: ", log_file)
  } else {
    message(paste(si, collapse = "\n"))
  }
  invisible(si)
}

#-------------------------------------------------------------------------------
# 11. Print configuration summary (useful for debugging)
#-------------------------------------------------------------------------------
print_config <- function() {
  cat("\n========== EOSCZ Pipeline Configuration ==========\n")
  cat("Project root:", PROJECT_ROOT, "\n")
  cat("Exome dir:  ", EXOME_DIR, "\n")
  cat("Meta dir:   ", META_ROOT, "\n")
  cat("Clinical:   ", CLIN_DIR, "\n")
  cat("Network:    ", NET_DIR, "\n")
  cat("Reference:  ", REF_FA, "\n")
  cat("Sample info:", SPINFO_FILE, "\n")
  cat("Cores:      ", N_CORES, "\n")
  cat("Run mode:   ", RUN_MODE, "\n")
  cat("================================================\n\n")
}

# Optionally print config on load (can be commented out)
# print_config()
