# =============================================================================
# PANSS-Metabolite Correlation Analysis
# =============================================================================

library(magrittr)
library(ggplot2)
library(dplyr)
library(ggrepel)

# =============================================================================
# 0. Configuration and Path Setup
# =============================================================================

META_STAT_DIR  <- "~/EOSCZ/MetaSeq-Result/statistics"
CLIN_DIR <- "~/EOSCZ/Clinical_Info/"

# =============================================================================
# 1. Load Metabolite Intensity Data
# =============================================================================

# --- 1.1 Load positive and negative mode metabolite data ---
META_STAT_DIR
pos_file <- file.path(META_STAT_DIR, "positive.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")
neg_file <- file.path(META_STAT_DIR, "negative.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")
pos_metabolites <- read.csv(pos_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
neg_metabolites <- read.csv(neg_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Combine and select key columns
all_metabolites <- rbind(pos_metabolites, neg_metabolites)
selected_columns <- c(1, 4, 13, 34, 39:106)  # XCMS.id, model, annotation, p-value, sample intensities
all_metabolites <- all_metabolites[, selected_columns]

# --- 1.2 Select differentially expressed metabolites (DEM) ---
dem_metabolites <- all_metabolites %>%
  dplyr::filter(All_sample.LM.padj < 0.05)

# --- 1.3 Merge metabolites with duplicate annotations (average intensities) ---
merge_duplicate_metabolites <- function(metabolite_df) {
  # Find metabolites with duplicate annotations
  duplicate_names <- metabolite_df$MS2Metabolite[duplicated(metabolite_df$MS2Metabolite)] %>%
    unique() %>%
    .[-1]  # remove NA or first element
  
  if (length(duplicate_names) == 0) {
    return(metabolite_df)
  }
  
  # Separate duplicate and unique metabolites
  dup_metabolites <- dplyr::filter(metabolite_df, MS2Metabolite %in% duplicate_names)
  unique_metabolites <- dplyr::filter(metabolite_df, !MS2Metabolite %in% duplicate_names)
  
  # Merge duplicate metabolites
  merged_metabolites <- data.frame()
  for (i in seq_along(duplicate_names)) {
    met_name <- duplicate_names[i]
    met_group <- dplyr::filter(dup_metabolites, MS2Metabolite == met_name)
    
    # Columns 1-4: concatenate or take first value
    new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(metabolite_df)))
    colnames(new_row) <- colnames(metabolite_df)
    
    for (j in 1:4) {
      new_row[1, j] <- paste(met_group[, j], collapse = ';')
    }
    
    # Columns 5-23: take first value
    new_row[1, 5:23] <- met_group[1, 5:23]
    
    # Column 24: concatenate
    new_row[1, 24] <- paste(met_group[, 24], collapse = ';')
    
    # Columns 25-72 (sample intensities): calculate mean
    intensity_cols <- 25:72
    new_row[1, intensity_cols] <- apply(met_group[, intensity_cols], 2, function(x) {
      mean(as.numeric(x), na.rm = TRUE)
    })
    
    merged_metabolites <- rbind(merged_metabolites, new_row)
  }
  
  rbind(unique_metabolites, merged_metabolites)
}

dem_metabolites <- merge_duplicate_metabolites(dem_metabolites)

# --- 1.4 Build analysis data frames ---
prepare_intensity_data <- function(metabolite_df, sample_names) {
  intensity_data <- metabolite_df[, sample_names]
  intensity_data$XCMS.id <- paste0(metabolite_df$XCMS.id, '_', metabolite_df$model)
  intensity_data$model <- metabolite_df$model
  
  # Reorder columns: XCMS.id, samples, model
  n_sample_cols <- length(sample_names)
  intensity_data <- intensity_data[, c(n_sample_cols + 1, 1:n_sample_cols, n_sample_cols + 2)]
  
  return(intensity_data)
}

dem_intensity <- prepare_intensity_data(dem_metabolites, sample_names)

# =============================================================================
# 2. Load PANSS Clinical Scores
# =============================================================================

sp_file <- file.path(CLIN_DIR, "sample_information.csv")
sample_info <- read.csv(sp_file, header = TRUE, stringsAsFactors = FALSE)
panss <- sample_info[sample_info$with_panss=='Y',
                     c('Sample', 'Name', 'Onset_age', 'PANSS.total','BPRS', 'PANSS.positive', 'PANSS.negative', 'Source')]

# =============================================================================
# 3. Prepare PANSS Score Vectors
# =============================================================================

ps_total <- as.numeric(panss$PANSS.total)      # Total PANSS score
ps_positive <- as.numeric(panss$PANSS.positive)  # Positive symptom score
ps_negative <- as.numeric(panss$PANSS.negative)  # Negative symptom score

# Assign sample names and sort consistently
names(ps_total) <- names(ps_positive) <- names(ps_negative) <- panss$Sample
ps_total <- ps_total[order(names(ps_total))]
ps_positive <- ps_positive[order(names(ps_positive))]
ps_negative <- ps_negative[order(names(ps_negative))]

# =============================================================================
# 4. Helper Function: Flatten Correlation Matrix
# =============================================================================

flatten_corr_matrix <- function(cormat, pmat) {
  upper_triangle <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[upper_triangle]],
    column = rownames(cormat)[col(cormat)[upper_triangle]],
    cor = cormat[upper_triangle],
    p = pmat[upper_triangle]
  )
}

# =============================================================================
# 5. Calculate Correlation Between All DEMs and PANSS Scores
# =============================================================================

dem_panss_sp <- dem_intensity
rownames(dem_panss_sp) <- dem_intensity$XCMS.id

dem_panss_sp_cor <- data.frame()

for (i in seq_len(nrow(dem_panss_sp))) {
  # Extract metabolite intensities for patient samples only
  x <- dem_panss_sp[i, names(ps_total)] %>% as.numeric()
  names(x) <- names(ps_total)
  
  # Calculate Pearson correlations with three PANSS subscores
  r_total <- cor.test(x, ps_total)
  r_positive <- cor.test(x, ps_positive)
  r_negative <- cor.test(x, ps_negative)
  
  # Store results
  dem_panss_sp_cor[i, 1:3] <- c(
    dem_panss_sp[i, 'XCMS.id'],
    dem_panss_sp[i, 'MS2Metabolite'],
    dem_panss_sp[i, 'MS2superclass']
  )
  dem_panss_sp_cor[i, 4:9] <- c(
    r_total$estimate, r_total$p.value,
    r_positive$estimate, r_positive$p.value,
    r_negative$estimate, r_negative$p.value
  )
}

# Set column names
colnames(dem_panss_sp_cor) <- c(
  'XCMS.id', 'MS2Metabolite', 'MS2superclass',
  'r_total.r', 'r_total.p', 'r_positive.r', 'r_positive.p',
  'r_negative.r', 'r_negative.p'
)

# Convert correlation and p-value columns to numeric
dem_panss_sp_cor$r_total.r <- as.numeric(dem_panss_sp_cor$r_total.r)
dem_panss_sp_cor$r_positive.r <- as.numeric(dem_panss_sp_cor$r_positive.r)
dem_panss_sp_cor$r_negative.r <- as.numeric(dem_panss_sp_cor$r_negative.r)
dem_panss_sp_cor$r_total.p <- as.numeric(dem_panss_sp_cor$r_total.p)
dem_panss_sp_cor$r_positive.p <- as.numeric(dem_panss_sp_cor$r_positive.p)
dem_panss_sp_cor$r_negative.p <- as.numeric(dem_panss_sp_cor$r_negative.p)

# =============================================================================
# 6. Filter Significant Correlations (|r| >= 0.3 and p < 0.1)
# =============================================================================

# Total PANSS score
select_cor_total <- rbind(
  subset(dem_panss_sp_cor, as.numeric(r_total.r) >= 0.3),
  subset(dem_panss_sp_cor, as.numeric(r_total.r) <= -0.3)
) %>%
  dplyr::filter(r_total.p < 0.1) 

# Positive symptom score
select_cor_positive <- rbind(
  subset(dem_panss_sp_cor, as.numeric(r_positive.r) >= 0.3),
  subset(dem_panss_sp_cor, as.numeric(r_positive.r) <= -0.3)
) %>%
  dplyr::filter(r_positive.p < 0.1) 

# Negative symptom score
select_cor_negative <- rbind(
  subset(dem_panss_sp_cor, as.numeric(r_negative.r) >= 0.3),
  subset(dem_panss_sp_cor, as.numeric(r_negative.r) <= -0.3)
) %>%
  dplyr::filter(r_negative.p < 0.1)  # 16 metabolites

# Find metabolites correlated with all three PANSS subscores
select_cor_total %>%
  dplyr::filter(XCMS.id %in% select_cor_positive$XCMS.id) %>%
  dplyr::filter(XCMS.id %in% select_cor_negative$XCMS.id)

# =============================================================================
# 7. Save Significant Correlation Results
# =============================================================================

select_cor_total$type <- 'PANSS_total_cor'
select_cor_positive$type <- 'PANSS_positive_cor'
select_cor_negative$type <- 'PANSS_negative_cor'

sig_cor <- unique(rbind(select_cor_total, select_cor_positive, select_cor_negative))

write.table(sig_cor, file = paste0(CLIN_DIR, '/MultiVar_LM.DEMs_correlated_with_PANSS.txt'),
  row.names = FALSE, col.names = TRUE, quote = TRUE, sep = '\t')
