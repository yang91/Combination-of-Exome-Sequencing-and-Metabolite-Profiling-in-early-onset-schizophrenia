# =============================================================================
# 3.2.Metabolite_symptom_corr.R
# Metabolite-Symptom Correlation Analysis
# =============================================================================
# Calculates Pearson correlations between differentially expressed metabolites
# (DEMs) and PANSS symptom scores. Filters for |r| >= 0.3 and p < 0.1.
#
# Input:  Positive/negative mode metabolite statistics, sample info with PANSS
# Output: METABOLITE_COR_OUTPUT (DEM-PANSS correlation results)
# =============================================================================

source("config.R")
source("utils.R")

library(magrittr)
library(ggplot2)
library(dplyr)
library(ggrepel)

# =============================================================================
# 1. Load Metabolite Data
# =============================================================================

all_metabolites <- load_metabolite_data()

# Select key columns used in original analysis (kept for compatibility)
select_cols <- c("XCMS.id", "model", "MS2Metabolite", "MS2superclass",
                 "MS2class", "All_sample.LM.padj")
available_cols <- intersect(select_cols, colnames(all_metabolites))
all_metabolites <- all_metabolites[, available_cols, drop = FALSE]

# Filter DEMs and merge duplicates
dem_metabolites <- filter_dem(all_metabolites)
dem_metabolites <- merge_duplicate_metabolites(dem_metabolites)

# Prepare intensity data
dem_intensity <- prepare_intensity_data(dem_metabolites)

# =============================================================================
# 2. Load PANSS Clinical Scores
# =============================================================================

spinfo <- load_sample_info(with_panss_only = TRUE)
panss <- spinfo[, c('Sample', 'Name', 'Onset_age', 'PANSS.total',
                    'BPRS', 'PANSS.positive', 'PANSS.negative', 'Source')]

ps_total    <- setNames(as.numeric(panss$PANSS.total),    panss$Sample)
ps_positive <- setNames(as.numeric(panss$PANSS.positive), panss$Sample)
ps_negative <- setNames(as.numeric(panss$PANSS.negative), panss$Sample)

# Sort consistently
ps_total    <- ps_total[order(names(ps_total))]
ps_positive <- ps_positive[order(names(ps_positive))]
ps_negative <- ps_negative[order(names(ps_negative))]

# =============================================================================
# 3. Calculate Correlation Between All DEMs and PANSS Scores
# =============================================================================

dem_panss_sp <- dem_intensity
rownames(dem_panss_sp) <- dem_intensity$XCMS.id

# Match PANSS samples with metabolite sample columns
sample_names <- names(ps_total)
common_samples <- intersect(sample_names, colnames(dem_panss_sp))
if (length(common_samples) == 0) {
  stop("No common samples between PANSS and metabolite data.")
}

# Pre-allocate results data frame
dem_panss_sp_cor <- data.frame(
  XCMS.id = character(nrow(dem_panss_sp)),
  MS2Metabolite = character(nrow(dem_panss_sp)),
  MS2superclass = character(nrow(dem_panss_sp)),
  r_total.r = numeric(nrow(dem_panss_sp)),
  r_total.p = numeric(nrow(dem_panss_sp)),
  r_positive.r = numeric(nrow(dem_panss_sp)),
  r_positive.p = numeric(nrow(dem_panss_sp)),
  r_negative.r = numeric(nrow(dem_panss_sp)),
  r_negative.p = numeric(nrow(dem_panss_sp)),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(dem_panss_sp))) {
  x <- as.numeric(dem_panss_sp[i, common_samples])
  names(x) <- common_samples
  
  ct <- cor.test(x, ps_total[common_samples])
  cp <- cor.test(x, ps_positive[common_samples])
  cn <- cor.test(x, ps_negative[common_samples])
  
  dem_panss_sp_cor[i, ] <- c(
    dem_panss_sp[i, 'XCMS.id'],
    dem_panss_sp[i, 'MS2Metabolite'],
    dem_panss_sp[i, 'MS2superclass'],
    ct$estimate, ct$p.value,
    cp$estimate, cp$p.value,
    cn$estimate, cn$p.value
  )
}

# Convert correlation and p-value columns to numeric
num_cols <- c("r_total.r", "r_total.p", "r_positive.r", "r_positive.p",
              "r_negative.r", "r_negative.p")
dem_panss_sp_cor[num_cols] <- lapply(dem_panss_sp_cor[num_cols], as.numeric)

# =============================================================================
# 4. Filter Significant Correlations (|r| >= 0.3 and p < 0.1)
# =============================================================================

filter_sig_cor <- function(df, r_col, p_col, r_thresh = 0.3, p_thresh = 0.1) {
  rbind(
    subset(df, df[[r_col]] >= r_thresh),
    subset(df, df[[r_col]] <= -r_thresh)
  ) %>% dplyr::filter(!!rlang::sym(p_col) < p_thresh)
}

select_cor_total    <- filter_sig_cor(dem_panss_sp_cor, "r_total.r",    "r_total.p")
select_cor_positive <- filter_sig_cor(dem_panss_sp_cor, "r_positive.r", "r_positive.p")
select_cor_negative <- filter_sig_cor(dem_panss_sp_cor, "r_negative.r", "r_negative.p")

# Find metabolites correlated with all three PANSS subscores
select_cor_total %>%
  dplyr::filter(XCMS.id %in% select_cor_positive$XCMS.id) %>%
  dplyr::filter(XCMS.id %in% select_cor_negative$XCMS.id)

# =============================================================================
# 5. Save Significant Correlation Results
# =============================================================================

select_cor_total$type    <- 'PANSS_total_cor'
select_cor_positive$type <- 'PANSS_positive_cor'
select_cor_negative$type <- 'PANSS_negative_cor'

sig_cor <- unique(rbind(select_cor_total, select_cor_positive, select_cor_negative))

write.table(sig_cor, file = METABOLITE_COR_OUTPUT,
            row.names = FALSE, col.names = TRUE, quote = TRUE, sep = '\t')

message("Metabolite-symptom correlation analysis complete.")
message("Results saved to: ", METABOLITE_COR_OUTPUT)
message("Significant correlations: ", nrow(sig_cor))
