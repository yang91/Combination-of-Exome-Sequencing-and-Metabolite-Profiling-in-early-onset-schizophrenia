# =============================================================================
# Metabolite-metabolite correlation Analysis
# =============================================================================

library(magrittr)
library(Hmisc)
library(dplyr)

# =============================================================================
# 0. Configuration and Path Setup
# =============================================================================

META_STAT_DIR  <- "~/EOSCZ/MetaSeq-Result/statistics/"
NET_DIR <- "~/EOSCZ/Network_file/"

# =============================================================================
# 1. Metabolite Data Loading and Processing
# =============================================================================

# --- 1.1 Load positive and negative mode metabolite data ---
pos_file <- file.path(META_STAT_DIR, "positive.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")
neg_file <- file.path(META_STAT_DIR, "negative.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")
pos_metabolites <- read.csv(pos_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
neg_metabolites <- read.csv(neg_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Combine and select key columns
all_metabolites <- rbind(pos_metabolites, neg_metabolites)
selected_columns <- c(1, 4, 13, 34, 39:106)  # XCMS.id, model, annotation, p-value, sample intensities
all_metabolites <- all_metabolites[, selected_columns]

# --- 1.2 Separate differentially expressed metabolites (DEM) and non-DEM ---
dem_metabolites <- all_metabolites %>%
  dplyr::filter(All_sample.LM.padj < 0.05)

non_dem_metabolites <- all_metabolites %>%
  dplyr::filter(!XCMS.id %in% dem_metabolites$XCMS.id)

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
non_dem_intensity <- prepare_intensity_data(non_dem_metabolites, sample_names)

# =============================================================================
# 2. Calculate metabolite-metabolite correlation
# =============================================================================

# Function to convert correlation matrix into correlation paris
flatternCorrMatrix <- function(cormat, pmat){
  ut <- upper.tri(cormat)
  data.frame(row = rownames(cormat)[row(cormat)[ut]],
             column = rownames(cormat)[col(cormat)[ut]],
             cor = cormat[ut],
             p = pmat[ut])
}

# --- Calculate the correlation ---
all_intensity <- rbind(dem_intensity, non_dem_intensity)
rownames(all_intensity) <- all_intensity$XCMS.id
dem_cor <- rcorr(as.matrix(t(all_intensity[,2:49])), type = 'spearman')
dem_cor_flat <- flatternCorrMatrix(dem_cor$r, dem_cor$P)
dem_cor_flat$q <- p.adjust(dem_cor_flat$p, method = 'BH')
save(dem_cor_flat,'./MultiVar.LM.dem_cor_flat.rds')

# =============================================================================
# 3. Save Results
# =============================================================================

output_file <- file.path(NET_DIR, "DEM_NonDEM.correlation.RData")
save(dem_cor_flat, file = output_file)
