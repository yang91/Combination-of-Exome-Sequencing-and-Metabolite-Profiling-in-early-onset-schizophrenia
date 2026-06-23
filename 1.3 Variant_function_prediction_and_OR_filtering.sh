#!/bin/bash
#===============================================================================
# ExomeSeq Variant Annotation Pipeline
# Before running, ensure bcftools and vep are in your PATH.
# All seqdata, reference annotations and output results should be in the same 
# project directory (use symlinks if needed).
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
RESULT_DIR="${HOME}/ExomeSeq-Result/variant"
GENOME="${HOME}/Ref_and_Index/GATK-hg38bundle/Homo_sapiens_assembly38.fa"
REF_DIR="${HOME}/Ref_and_Index"
VEP_CACHE="${REF_DIR}/.vep"

# Create output directories
mkdir -p "${RESULT_DIR}"
mkdir -p "${RESULT_DIR}/maf_005"

#-------------------------------------------------------------------------------
# Step 0: Split multi-allelic variants (optional, uncomment if needed)
#-------------------------------------------------------------------------------
bcftools norm -m -both -f "${GENOME}" -o "${RESULT_DIR}/final.QC.vcf" \
     "${RESULT_DIR}/final.QC.vcf.bgz"

#-------------------------------------------------------------------------------
# Step 1: Run VEP command-line annotation
#-------------------------------------------------------------------------------
echo "[1/3] Running VEP command-line annotation..."

vep -i "${RESULT_DIR}/final.QC.vcf" \
    --plugin "ExACpLI,${VEP_CACHE}/Plugins/ExACpLI_values.txt" \
    --plugin "MPC,${VEP_CACHE}/Plugins/fordist_constraint_official_mpc_values_v2.txt.gz" \
    --offline \
    --individual all \
    --assembly GRCh38 \
    --cache \
    --dir_cache "${VEP_CACHE}" \
    -o "${RESULT_DIR}/VEP_command.txt"

echo "VEP command-line annotation completed."

#-------------------------------------------------------------------------------
# Step 2: Run VEP webserver annotation (manual step - assumed already done)
# The annotated file VEP_server.txt should be placed in ${RESULT_DIR}
#-------------------------------------------------------------------------------

if [[ ! -f "${RESULT_DIR}/VEP_server.txt" ]]; then
    echo "WARNING: ${RESULT_DIR}/VEP_server.txt not found. Please upload final.QC.vcf to VEP webserver and save the result as VEP_server.txt in ${RESULT_DIR}"
    exit 1
fi

#-------------------------------------------------------------------------------
# Step 3: R analysis - variant filtering and OR calculation
#-------------------------------------------------------------------------------
echo "[2/3] Running R variant analysis..."

Rscript --vanilla - <<'RSCRIPT' \
    --args "${RESULT_DIR}" "${REF_DIR}"

# R code starts here
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript script.R <RESULT_DIR> <REF_DIR>")
}
RESULT_DIR <- args[1]
REF_DIR <- args[2]

library(magrittr)

#-------------------------------------------------------------------------------
# Step 1: Import VEP command annotation
#-------------------------------------------------------------------------------
cat("[R] Importing VEP command results...\n")

vep_cmd_file <- file.path(RESULT_DIR, "VEP_command.txt")
if (!file.exists(vep_cmd_file)) {
    stop("VEP command output not found: ", vep_cmd_file)
}

var_cmd <- read.table(vep_cmd_file, header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, quote = "", 
                       comment.char = "", check.names = FALSE)

# Add group labels
var_cmd$Group <- ifelse(grepl("NOR", var_cmd$Sample), "Normal",
                 ifelse(grepl("SCZ", var_cmd$Sample), "Patient", NA))

#-------------------------------------------------------------------------------
# Step 2: Calculate odds ratios using Fisher's exact test
#-------------------------------------------------------------------------------
cat("[R] Calculating odds ratios...\n")

# Count samples per group
n_scz_total <- length(unique(var_cmd$Sample[var_cmd$Group == "Patient"]))
n_nor_total <- length(unique(var_cmd$Sample[var_cmd$Group == "Normal"]))

cat("  SCZ samples:", n_scz_total, "\n")
cat("  NOR samples:", n_nor_total, "\n")

variants <- unique(var_cmd$Uploaded_variation)
or_results <- lapply(variants, function(v) {
    tmp <- var_cmd[var_cmd$Uploaded_variation == v, 
                   c("Uploaded_variation", "Sample", "Group")] %>% unique()

    n1 <- sum(tmp$Group == "Patient")
    n2 <- sum(tmp$Group == "Normal")

    # Construct 2x2 contingency table
    #           Variant+  Variant-
    # SCZ         n1      n_scz_total - n1
    # NOR         n2      n_nor_total - n2

    mat <- matrix(c(n1, n_scz_total - n1,
                    n2, n_nor_total - n2),
                  nrow = 2, byrow = TRUE)

    # Handle edge cases where Fisher test may fail
    if (any(mat < 0)) {
        return(data.frame(
            Uploaded_variation = v,
            SCZ_num = n1,
            NOR_num = n2,
            OR = NA,
            p_value = NA,
            stringsAsFactors = FALSE
        ))
    }

    # Only run fisher.test if table is valid
    if (sum(mat[1, ]) == 0 || sum(mat[2, ]) == 0 || 
        mat[1,1] == sum(mat[1,]) || mat[2,1] == sum(mat[2,])) {
        # All or none have variant in a group
        return(data.frame(
            Uploaded_variation = v,
            SCZ_num = n1,
            NOR_num = n2,
            OR = ifelse(n1 > 0 && n2 == 0, Inf, 
                 ifelse(n1 == 0 && n2 > 0, 0, 1)),
            p_value = NA,
            stringsAsFactors = FALSE
        ))
    }

    sta <- fisher.test(mat)
    data.frame(
        Uploaded_variation = v,
        SCZ_num = n1,
        NOR_num = n2,
        OR = as.numeric(sta$estimate),
        p_value = sta$p.value,
        stringsAsFactors = FALSE
    )
})

or_df <- do.call(rbind, or_results)
var_cmd <- merge(var_cmd, or_df, by = "Uploaded_variation", all.x = TRUE)

#-------------------------------------------------------------------------------
# Step 3: Import VEP webserver annotation
#-------------------------------------------------------------------------------
cat("[R] Importing VEP webserver results...\n")

vep_web_file <- file.path(RESULT_DIR, "VEP_server.txt")
var_web <- read.table(vep_web_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "",
                       comment.char = "", check.names = FALSE)

# Select relevant annotation columns
keep_cols <- c(
    "X.Uploaded_variation", "Location", "Consequence", "IMPACT", 
    "Gene", "Feature_type", "Feature", "BIOTYPE", "cDNA_position",
    "CDS_position", "Protein_position", "Amino_acids", "Codons",
    "DISTANCE", "STRAND", "FLAGS", "SIFT", "PolyPhen",
    "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS", "DisGeNET", "CADD_PHRED", "CADD_RAW",
    "Aloft_Confidence", "Aloft_Fraction_transcripts_affected", "Aloft_pred",
    "Aloft_prob_Dominant", "Aloft_prob_Recessive", "Aloft_prob_Tolerant",
    "MutationTaster_AAE", "MutationTaster_converted_rankscore",
    "MutationTaster_model", "MutationTaster_pred", "MutationTaster_score",
    "clinvar_MedGen_id", "clinvar_OMIM_id", "clinvar_Orphanet_id",
    "clinvar_clnsig", "clinvar_hgvs", "clinvar_id", "clinvar_review",
    "clinvar_trait", "clinvar_var_source"
)

# Only keep columns that exist in the data
keep_cols <- intersect(keep_cols, colnames(var_web))
var_web_sel <- var_web[, keep_cols, drop = FALSE]

#-------------------------------------------------------------------------------
# Step 4: Merge annotations and filter by OR > 1
#-------------------------------------------------------------------------------
cat("[R] Merging annotations and filtering...\n")

# Merge: keep only variants with OR info from command results
var_merged <- merge(
    var_cmd[, c("Uploaded_variation", "OR", "p_value")],
    var_web_sel,
    by.x = "Uploaded_variation",
    by.y = "X.Uploaded_variation",
    all.x = TRUE
)

# Filter OR > 1
var_or1 <- subset(var_merged, OR > 1)

# Handle CADD_PHRED
var_or1$CADD_PHRED <- suppressWarnings(as.numeric(var_or1$CADD_PHRED))
var_or1$CADD_PHRED[is.na(var_or1$CADD_PHRED)] <- 0

#-------------------------------------------------------------------------------
# Step 5: Separate PTV and SAV
#-------------------------------------------------------------------------------
cat("[R] Classifying variants...\n")

# Protein Truncating Variants (HIGH impact)
var_ptv <- subset(var_or1, IMPACT == "HIGH")
cat("  PTV variants:", length(unique(var_ptv$Uploaded_variation)), "\n")

# Single Amino acid Variants (MODERATE impact, excluding TFBS variants)
var_sav <- subset(var_or1, IMPACT == "MODERATE")
var_sav <- var_sav[!grepl("TFBS_ablation", var_sav$Consequence, fixed = TRUE), ]
cat("  SAV variants:", length(unique(var_sav$Uploaded_variation)), "\n")

#-------------------------------------------------------------------------------
# Step 6: Rank scoring for SAVs
#-------------------------------------------------------------------------------
cat("[R] Scoring SAVs...\n")

# Initialize scoring columns
var_sav$sift_jud <- 0
var_sav$polyphen_jud <- 0
var_sav$cadd_jud <- 0
var_sav$mutationtaster_jud <- 0

# SIFT: deleterious(0) = damaging
var_sav$sift_jud[grepl("deleterious\\(0\\)", var_sav$SIFT, perl = TRUE)] <- 1

# PolyPhen: probably_damaging(1) = damaging
var_sav$polyphen_jud[grepl("probably_damaging\\(1\\)", var_sav$PolyPhen, perl = TRUE)] <- 1

# CADD: PHRED >= 32
var_sav$cadd_jud[var_sav$CADD_PHRED >= 32] <- 1

# MutationTaster: A (disease_causing_automatic) or D (disease_causing) = damaging
# N (polymorphism) = benign
mt_pred <- var_sav$MutationTaster_pred
var_sav$mutationtaster_jud[grepl("[AD]", mt_pred, perl = TRUE)] <- 1
var_sav$mutationtaster_jud[grepl("N", mt_pred, perl = TRUE)] <- 0

# Total score
var_sav$jud_total <- var_sav$sift_jud + var_sav$polyphen_jud + 
                     var_sav$cadd_jud + var_sav$mutationtaster_jud

# Keep SAVs with score >= 3
var_sav_rank <- subset(var_sav, jud_total >= 3)
cat("  SAVs with score >= 3:", nrow(var_sav_rank), "\n")

#-------------------------------------------------------------------------------
# Step 7: Convert ENSG to Entrez ID
#-------------------------------------------------------------------------------
cat("[R] Converting gene IDs...\n")

cvt_file <- file.path(REF_DIR, "combine.ensemb_104.orgHsegdb.ncbi_name_searching.manually_checked.txt")
if (!file.exists(cvt_file)) {
    stop("Gene conversion file not found: ", cvt_file)
}

cvt <- read.table(cvt_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(cvt)[1] <- "Gene"
colnames(cvt)[6] <- "ENTREZID"

# Merge with conversion table
var_sav_cvt <- merge(var_sav_rank, cvt[, c("Gene", "ENTREZID")], 
                      by.x = "Gene", by.y = "Gene", all.x = TRUE)
var_ptv_cvt <- merge(var_ptv, cvt[, c("Gene", "ENTREZID")], 
                      by.x = "Gene", by.y = "Gene", all.x = TRUE)

#-------------------------------------------------------------------------------
# Step 8: Combine and output
#-------------------------------------------------------------------------------
cat("[R] Writing output...\n")

# Ensure same columns before rbind
common_cols <- intersect(colnames(var_sav_cvt), colnames(var_ptv_cvt))
var_sav_cvt <- var_sav_cvt[, common_cols, drop = FALSE]
var_ptv_cvt <- var_ptv_cvt[, common_cols, drop = FALSE]

var_selected <- rbind(var_sav_cvt, var_ptv_cvt) %>% unique()

output_file <- file.path(RESULT_DIR, "Maf_0.05.potential_harmful.ORgt1_variants.txt")

write.table(var_selected, file = output_file,
            row.names = FALSE, col.names = TRUE, 
            sep = "\t", quote = FALSE)

cat("[R] Done. Output:", output_file, "\n")
cat("[R] Total selected variants:", length(unique(var_selected$Uploaded_variation)), "\n")

# R code ends here
RSCRIPT

echo "[3/3] Pipeline completed successfully."
