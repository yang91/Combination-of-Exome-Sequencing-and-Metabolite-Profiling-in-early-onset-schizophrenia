# =============================================================================
# 3.1.Variant_symptom_assoc.R
# Variant-Symptom Association Analysis
# =============================================================================
# Associates high-risk gene variants with PANSS symptom scores using t-tests.
#
# Input:  VEP annotations, disease-causing variants, gene scores, sample info
# Output: PANSS_ASSOC_OUTPUT (significant gene-PANSS associations)
# =============================================================================

source("config.R")
source("utils.R")

library(ggsignif)
library(cowplot)

# =============================================================================
# 1. Load Data
# =============================================================================

# Gene scores (used for filtering high-risk genes)
gene_scores <- read.csv(GENE_SCORE_FILE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Disease-causing variants with OR > 1 in genes with GS > 4
disease_variants <- read.csv(DISEASE_VARIANT_FILE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# VEP annotation results (skip header lines)
vep_results <- read.csv(VEP_ANNOT_FILE, sep = '\t', stringsAsFactors = FALSE, skip = VEP_SKIP_LINES, row.names = 1)

# Sample information
sample_info <- load_sample_info()
sample_names <- rownames(sample_info)
n_samples <- length(sample_names)

# =============================================================================
# 2. Variant Data Processing (Variant-Level)
# =============================================================================

# Filter VEP results to retain only target genes and variants
filtered_vep <- vep_results %>%
  dplyr::filter(
    Gene %in% disease_variants$Gene,
    X.Uploaded_variation %in% disease_variants$X.Uploaded_variation
  )

# Build per-variant presence (0/1) and effect (impact level 1-4) matrices
variant_mats <- build_variant_matrices(filtered_vep, sample_names)
var_presence <- variant_mats$presence
var_effect   <- variant_mats$effect

# =============================================================================
# 3. Variant Data Aggregation (Gene-Level)
# =============================================================================

gene_list <- unique(disease_variants$Gene)
n_genes <- length(gene_list)

gene_presence <- matrix(0, nrow = n_genes, ncol = n_samples)
gene_effect   <- matrix(0, nrow = n_genes, ncol = n_samples)
colnames(gene_presence) <- colnames(gene_effect) <- sample_names
rownames(gene_presence) <- rownames(gene_effect) <- gene_list

for (i in seq_len(n_genes)) {
  gene_name <- gene_list[i]
  gene_variants <- disease_variants %>%
    dplyr::filter(Gene == gene_name) %>%
    dplyr::pull(X.Uploaded_variation)
  
  for (j in seq_len(n_samples)) {
    sample_id <- sample_names[j]
    sample_variants <- gene_variants[gene_variants %in% rownames(var_presence)]
    
    if (length(sample_variants) > 0) {
      presence_values <- as.numeric(var_presence[sample_variants, sample_id])
      effect_values <- as.numeric(var_effect[sample_variants, sample_id])
      
      gene_presence[i, j] <- min(sum(presence_values), 1)
      gene_effect[i, j]   <- max(effect_values)
    }
  }
}

# Add sample information column for each variant
disease_variants$sample <- apply(disease_variants, 1, function(row) {
  variant_id <- row[1]
  vep_records <- vep_results[vep_results$X.Uploaded_variation == variant_id, ]
  sample_ids <- extract_vep_samples(vep_records)$sample
  valid_samples <- sample_ids[sample_ids %in% sample_names]
  paste(unique(valid_samples), collapse = ',')
})

# =============================================================================
# 4. Load PANSS Score Data
# =============================================================================

panss <- sample_info[sample_info$with_panss == 'Y',
                     c('Sample', 'Name', 'Onset_age', 'PANSS.total',
                       'BPRS', 'PANSS.positive', 'PANSS.negative', 'Source')]

# =============================================================================
# 5. Compare PANSS Scores Between Mutated and Non-Mutated Samples
# =============================================================================

# Map each gene to its mutated samples
gene2sample <- tapply(disease_variants$sample, disease_variants$Gene, function(s) {
  unique(unlist(strsplit(s, ',')))
})

# Convert to uniform data frame for downstream analysis
uniq_genes <- data.frame(
  Gene   = names(gene2sample),
  Sample = sapply(gene2sample, paste, collapse = ','),
  row.names = NULL
)

# Test each gene individually
mut_stat <- apply(uniq_genes, 1, function(x) {
  gene_id  <- x["Gene"]
  mut_sp   <- unlist(strsplit(x["Sample"], ','))
  gene_sym <- gene_scores[gene_scores$ENSEMBL == gene_id, "SYMBOL"]
  if (length(gene_sym) == 0) gene_sym <- NA
  
  # ---- Prepare PANSS data frame ----
  panss_dt <- data.frame(
    sample         = panss$Sample,
    panss_total    = as.numeric(panss[["PANSS.total"]]),
    panss_positive = as.numeric(panss[["PANSS.positive"]]),
    panss_negative = as.numeric(panss[["PANSS.negative"]]),
    mutype = factor(rep("UnMut", nrow(panss)), levels = c("UnMut", "Mut"))
  )
  panss_dt$mutype[panss_dt$sample %in% mut_sp] <- "Mut"
  
  # ---- Statistical tests ----
  if (sum(panss_dt$mutype == "Mut") > 1) {
    tt <- t.test(panss_total    ~ mutype, data = panss_dt)
    tp <- t.test(panss_positive ~ mutype, data = panss_dt)
    tn <- t.test(panss_negative ~ mutype, data = panss_dt)
    
    # Calculate group means
    mean_tbl <- panss_dt %>%
      dplyr::group_by(mutype) %>%
      dplyr::summarise(dplyr::across(
        c(panss_total, panss_positive, panss_negative),
        \(x) mean(x, na.rm = TRUE)
      ))
    
    c(
      gene_sym, gene_id, paste(mut_sp, collapse = ","),
      tt$p.value, tp$p.value, tn$p.value,
      mean_tbl$panss_total[mean_tbl$mutype == "Mut"],
      mean_tbl$panss_total[mean_tbl$mutype == "UnMut"],
      mean_tbl$panss_positive[mean_tbl$mutype == "Mut"],
      mean_tbl$panss_positive[mean_tbl$mutype == "UnMut"],
      mean_tbl$panss_negative[mean_tbl$mutype == "Mut"],
      mean_tbl$panss_negative[mean_tbl$mutype == "UnMut"]
    )
  } else {
    c(
      gene_sym, gene_id, paste(mut_sp, collapse = ","),
      NA_real_, NA_real_, NA_real_,
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_
    )
  }
}) %>%
  t() %>%
  as.data.frame()

colnames(mut_stat) <- c(
  "Gene_symbol", "Ensembl", "Mut_samples",
  "P_total", "P_positive", "P_negative",
  "Mean_Total_Mut", "Mean_Total_UnMut",
  "Mean_Positive_Mut", "Mean_Positive_UnMut",
  "Mean_Negative_Mut", "Mean_Negative_UnMut"
)

# Convert numeric columns
num_cols <- c("P_total", "P_positive", "P_negative",
              "Mean_Total_Mut", "Mean_Total_UnMut",
              "Mean_Positive_Mut", "Mean_Positive_UnMut",
              "Mean_Negative_Mut", "Mean_Negative_UnMut")
mut_stat[num_cols] <- lapply(mut_stat[num_cols], as.numeric)

# Extract significant results (p < 0.05 in any domain)
sig_mut <- mut_stat %>%
  dplyr::filter(dplyr::if_any(c(P_total, P_positive, P_negative), ~ .x < 0.05))

write.table(sig_mut, file = PANSS_ASSOC_OUTPUT,
            row.names = FALSE, col.names = TRUE, quote = TRUE, sep = '\t')

message("Variant-symptom association analysis complete.")
message("Results saved to: ", PANSS_ASSOC_OUTPUT)
message("Total genes tested: ", n_genes)
message("Significant associations: ", nrow(sig_mut))
