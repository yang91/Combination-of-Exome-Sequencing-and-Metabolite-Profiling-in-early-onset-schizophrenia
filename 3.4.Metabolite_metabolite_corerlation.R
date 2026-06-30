# =============================================================================
# 3.4.Metabolite_metabolite_corerlation.R
# Metabolite-Metabolite Correlation Analysis
# =============================================================================
# Calculates Spearman correlations among all annotated metabolites.
#
# Input:  Positive/negative mode metabolite statistics
# Output: DEM_COR_FILE (flattened correlation matrix)
# =============================================================================

source("config.R")
source("utils.R")

library(magrittr)
library(Hmisc)
library(dplyr)

# =============================================================================
# 1. Load Metabolite Data
# =============================================================================

all_metabolites <- load_metabolite_data()

# Separate DEMs and non-DEMs
dem_metabolites <- filter_dem(all_metabolites)
non_dem_metabolites <- all_metabolites %>%
  dplyr::filter(!XCMS.id %in% dem_metabolites$XCMS.id)

# Merge duplicate annotations
dem_metabolites <- merge_duplicate_metabolites(dem_metabolites)
non_dem_metabolites <- merge_duplicate_metabolites(non_dem_metabolites)

# Prepare intensity data
dem_intensity <- prepare_intensity_data(dem_metabolites)
non_dem_intensity <- prepare_intensity_data(non_dem_metabolites)

# =============================================================================
# 2. Calculate Metabolite-Metabolite Correlation
# =============================================================================

all_intensity <- rbind(dem_intensity, non_dem_intensity)
rownames(all_intensity) <- all_intensity$XCMS.id

# Extract sample columns dynamically
sample_cols <- setdiff(colnames(all_intensity), c("XCMS.id", "model"))

dem_cor <- rcorr(as.matrix(t(all_intensity[, sample_cols])), type = 'spearman')
dem_cor_flat <- flatten_corr_matrix(dem_cor$r, dem_cor$P)
dem_cor_flat$q <- p.adjust(dem_cor_flat$p, method = 'BH')

# =============================================================================
# 3. Save Results
# =============================================================================

save(dem_cor_flat, file = DEM_COR_FILE)

message("Metabolite-metabolite correlation analysis complete.")
message("Results saved to: ", DEM_COR_FILE)
