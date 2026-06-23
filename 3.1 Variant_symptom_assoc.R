library(ggsignif)
library(cowplot)

# =============================================================================
# 0. Configuration and Path Setup
# =============================================================================

EXOME_DIR <- "~/EOSCZ/ExomeSeq-Result/"
VAR_DIR <- "~/EOSCZ/ExomeSeq-Result/variant"
CLIN_DIR <- "~/EOSCZ/Clinical_Info/"


# ============================================================================
# 1. Load and Process Mutation Data
# ============================================================================

# Gene score evaluation model from an in-house script of unpublished data
gs_file <- file.path(EXOME_DIR, "gene_score.high-risk.txt")
gene_scores <- read.csv(gs_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Load disease-causing variants with OR > 1 in genes with GS > 4
var_file <- file.path(VAR_DIR, "disease_causing_ORgt1_variants_in_GSgt4_gene.txt")
disease_variants <- read.csv(var_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Load VEP annotation results (skip first 44 header lines)
vep_file <- file.path(VAR_DIR, "Maf_0.05.potential_harmful.ORgt1_variants.txt")
vep_results <- read.csv(vep_file, sep = '\t', stringsAsFactors = FALSE, skip = 44, row.names = 1)

# Load sample information
sp_file <- file.path(CLIN_DIR, "sample_information.csv")
sample_info <- read.csv(sp_file, header = TRUE, stringsAsFactors = FALSE)

# Exclude QC samples and sort by group
sample_info <- sample_info %>%
  dplyr::filter(group != 'QC') %>%
  dplyr::arrange(group)

rownames(sample_info) <- sample_info$sample.name
sample_names <- sample_info$sample.name
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

variant_list <- unique(filtered_vep$X.Uploaded_variation)
n_variants <- length(variant_list)

# Initialize matrices:
# var_presence: 1 = sample carries the variant, 0 = does not carry
# var_effect: variant impact level (4 = HIGH, 3 = MODERATE, 2 = LOW, 1 = MODIFIER)
var_presence <- matrix(0, nrow = n_variants, ncol = n_samples)
var_effect   <- matrix(0, nrow = n_variants, ncol = n_samples)
colnames(var_presence) <- colnames(var_effect) <- sample_names
rownames(var_presence) <- rownames(var_effect) <- variant_list

# Impact level mapping
EFFECT_LEVELS <- c('HIGH' = 4, 'MODERATE' = 3, 'LOW' = 2, 'MODIFIER' = 1)

# Process each variant
for (i in seq_len(n_variants)) {
  variant_id <- variant_list[i]
  variant_records <- dplyr::filter(filtered_vep, X.Uploaded_variation == variant_id)
  
  # Extract sample IDs and impact levels
  variant_info <- apply(variant_records, 1, function(row) {
    fields <- unlist(strsplit(row, ';'))
    
    # Extract sample ID
    sample_field <- fields[grepl('IND=', fields)]
    sample_id <- gsub('IND=', '', sample_field)
    
    # Extract impact level
    impact_field <- fields[grepl('IMPACT=', fields)]
    impact_level <- gsub('IMPACT=', '', impact_field)
    
    data.frame(sample = sample_id, impact = impact_level, stringsAsFactors = FALSE)
  })
  
  # Combine all records
  variant_df <- do.call(rbind, variant_info)
  variant_df$impact_num <- EFFECT_LEVELS[variant_df$impact]
  
  # Get samples carrying this variant (deduplicated)
  carrier_samples <- unique(variant_df$sample)
  carrier_samples <- carrier_samples[carrier_samples %in% sample_names]
  
  # Fill matrices: record maximum impact level per sample
  for (sample_id in carrier_samples) {
    sample_effects <- variant_df$impact_num[variant_df$sample == sample_id]
    var_presence[i, sample_id] <- 1
    var_effect[i, sample_id] <- max(sample_effects, na.rm = TRUE)
  }
}

rownames(var_effect) <- rownames(var_presence) <- variant_list


# =============================================================================
# 3. Variant Data Aggregation (Gene-Level)
# =============================================================================

gene_list <- unique(disease_variants$Gene)
n_genes <- length(gene_list)

# Initialize gene-level matrices
gene_presence <- matrix(0, nrow = n_genes, ncol = n_samples)  # whether sample carries any variant in gene
gene_effect   <- matrix(0, nrow = n_genes, ncol = n_samples)  # maximum impact level in gene
colnames(gene_presence) <- colnames(gene_effect) <- sample_names
rownames(gene_presence) <- rownames(gene_effect) <- gene_list

for (i in seq_len(n_genes)) {
  gene_name <- gene_list[i]
  gene_variants <- disease_variants %>%
    dplyr::filter(Gene == gene_name) %>%
    dplyr::pull(X.Uploaded_variation)
  
  # Calculate variant count and maximum impact per sample
  for (j in seq_len(n_samples)) {
    sample_id <- sample_names[j]
    sample_variants <- gene_variants[gene_variants %in% variant_list]
    
    if (length(sample_variants) > 0) {
      presence_values <- as.numeric(var_presence[sample_variants, sample_id])
      effect_values <- as.numeric(var_effect[sample_variants, sample_id])
      
      gene_presence[i, j] <- min(sum(presence_values), 1)  # cap at 1 (presence/absence)
      gene_effect[i, j] <- max(effect_values)
    }
  }
}

# Add sample information column for each variant
disease_variants$sample <- apply(disease_variants, 1, function(row) {
  variant_id <- row[1]
  vep_records <- vep_results[vep_results$X.Uploaded_variation == variant_id, ]
  
  # Extract sample IDs
  sample_ids <- apply(vep_records, 1, function(vep_row) {
    fields <- unlist(strsplit(vep_row[14], split = ';'))
    gsub('IND=', '', fields[1])
  })
  
  # Filter valid samples and deduplicate
  valid_samples <- sample_ids[sample_ids %in% sample_names]
  paste(unique(valid_samples), collapse = ',')
})


# ============================================================================
# 4. Load PANSS Score Data
# ============================================================================

sp_file <- file.path(CLIN_DIR, "sample_information.csv")
sample_info <- read.csv(sp_file, header = TRUE, stringsAsFactors = FALSE)
panss <- sample_info[,c('Sample', 'Name', 'Onset_age', 'PANSS.total','BPRS', 'PANSS.positive', 'PANSS.negative', 'Source')]

# Keep only samples present in metabolite data
panss <- panss[panss$Sample %in% colnames(dem_intensity), ]


# ============================================================================
# 5. Compare PANSS Scores Between Mutated and Non-Mutated Samples
# ============================================================================

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
  mut_sp   <- unlist(strsplit(x["Sample"], ','))  # already deduplicated
  gene_sym <- gs_file[gs_file$ENSEMBL == gene_id, "SYMBOL"]
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
      group_by(mutype) %>%
      summarise(across(
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

# Extract significant results (p < 0.05 in any domain)
sig_mut <- mut_stat %>%
  filter(if_any(c(P_total, P_positive, P_negative), ~ .x < 0.05))

write.table(sig_mut, file = paste0(CLIN_DIR, '/High_risk_genes_associated_with_PANSS.txt',
  row.names = FALSE, col.names = TRUE, quote = TRUE, sep = '\t')
