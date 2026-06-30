# =============================================================================
# 3.3.Variant_metabolite_assoc.R
# Variant-Metabolite Association Analysis
# =============================================================================
# Tests metabolite differences between variant carriers and non-carriers.
#
# Input:  VEP annotations, disease variants, metabolite statistics
# Output: VAR_MET_ASSOC_FILE (dem_association, non_dem_association)
# =============================================================================

source("config.R")
source("utils.R")

library(magrittr)
library(Hmisc)
library(dplyr)

# =============================================================================
# 1. Load Variant Data
# =============================================================================

gene_scores <- read.csv(GENE_SCORE_FILE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
disease_variants <- read.csv(DISEASE_VARIANT_FILE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
vep_results <- read.csv(VEP_ANNOT_FILE, sep = '\t', stringsAsFactors = FALSE, skip = VEP_SKIP_LINES, row.names = 1)

sample_info <- load_sample_info()
sample_names <- rownames(sample_info)
n_samples <- length(sample_names)

# Filter VEP to target genes and variants
filtered_vep <- vep_results %>%
  dplyr::filter(
    Gene %in% disease_variants$Gene,
    X.Uploaded_variation %in% disease_variants$X.Uploaded_variation
  )

# Build variant presence and effect matrices
variant_mats <- build_variant_matrices(filtered_vep, sample_names)
var_presence <- variant_mats$presence
var_effect <- variant_mats$effect

# Gene-level aggregation
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
      gene_presence[i, j] <- min(sum(var_presence[sample_variants, sample_id]), 1)
      gene_effect[i, j]   <- max(var_effect[sample_variants, sample_id])
    }
  }
}

# Add sample information to disease variants
disease_variants$sample <- apply(disease_variants, 1, function(row) {
  variant_id <- row[1]
  vep_records <- vep_results[vep_results$X.Uploaded_variation == variant_id, ]
  sample_ids <- extract_vep_samples(vep_records)$sample
  valid_samples <- sample_ids[sample_ids %in% sample_names]
  paste(unique(valid_samples), collapse = ',')
})

# =============================================================================
# 2. Metabolite Data Loading and Processing
# =============================================================================

all_metabolites <- load_metabolite_data()

# Separate DEMs and non-DEMs
dem_metabolites <- filter_dem(all_metabolites)
non_dem_metabolites <- all_metabolites %>%
  dplyr::filter(!XCMS.id %in% dem_metabolites$XCMS.id)

# Merge duplicate annotations
dem_metabolites <- merge_duplicate_metabolites(dem_metabolites)
non_dem_metabolites <- merge_duplicate_metabolites(non_dem_metabolites)

# Select reliable annotations
dem_metabolites_annt <- select_annt(dem_metabolites)
non_dem_metabolites_annt <- select_annt(non_dem_metabolites)

# Prepare intensity data
dem_intensity <- prepare_intensity_data(dem_metabolites, sample_names = sample_names)
non_dem_intensity <- prepare_intensity_data(non_dem_metabolites, sample_names = sample_names)

# =============================================================================
# 3. Core Analysis: Compare Carrier vs Non-Carrier
# =============================================================================

compare_metabolite_difference <- function(variant_info, metabolite_data) {
  carrier_samples <- unlist(strsplit(unique(variant_info$sample), split = ','))
  scz_samples <- sample_names[grep('SCZ', sample_names)]
  nor_samples <- sample_names[grep('NOR', sample_names)]
  
  scz_carriers <- carrier_samples[grep('SCZ', carrier_samples)]
  scz_non_carriers <- scz_samples[!(scz_samples %in% scz_carriers)]
  
  if (any(grepl('NOR', variant_info$sample))) {
    nor_carriers <- carrier_samples[grep('NOR', carrier_samples)]
    nor_non_carriers <- nor_samples[!(nor_samples %in% nor_carriers)]
  } else {
    nor_non_carriers <- nor_samples
  }
  
  # Dynamic sample columns from metabolite_data
  met_sample_names <- setdiff(colnames(metabolite_data), c("XCMS.id", "model"))
  
  results <- data.frame(
    XCMS.id = character(nrow(metabolite_data)),
    scz_var_vs_nor_no_var.pvalue = numeric(nrow(metabolite_data)),
    scz_var_vs_scz_no_var.pvalue = numeric(nrow(metabolite_data)),
    stringsAsFactors = FALSE
  )
  
  for (j in seq_len(nrow(metabolite_data))) {
    test1 <- t.test(
      as.numeric(metabolite_data[j, scz_carriers]),
      as.numeric(metabolite_data[j, nor_non_carriers])
    )
    
    test2 <- t.test(
      as.numeric(metabolite_data[j, scz_carriers]),
      as.numeric(metabolite_data[j, scz_non_carriers])
    )
    
    results[j, ] <- c(
      metabolite_data[j, 'XCMS.id'],
      test1$p.value,
      test2$p.value
    )
  }
  
  results$scz_var_vs_nor_no_var.pvalue <- as.numeric(results$scz_var_vs_nor_no_var.pvalue)
  results$scz_var_vs_scz_no_var.pvalue <- as.numeric(results$scz_var_vs_scz_no_var.pvalue)
  
  return(results)
}

# =============================================================================
# 4. Filter Genes with Variants in >=2 SCZ Samples
# =============================================================================

get_multisample_genes <- function(disease_variants) {
  duplicated_genes <- disease_variants$Gene[duplicated(disease_variants$Gene)]
  multi_sample_genes <- c()
  
  for (gene in unique(duplicated_genes)) {
    gene_samples <- disease_variants[disease_variants$Gene == gene, 'sample']
    has_no_nor <- !any(grepl('NOR', gene_samples))
    has_unique_samples <- !any(duplicated(gene_samples))
    if (has_no_nor && has_unique_samples) {
      multi_sample_genes <- c(multi_sample_genes, gene)
    }
  }
  
  multi_variant_genes <- c()
  multi_variant_records <- disease_variants[grepl(',', disease_variants$sample), ]
  for (i in seq_len(nrow(multi_variant_records))) {
    if (!grepl('NOR', multi_variant_records$sample[i])) {
      multi_variant_genes <- c(multi_variant_genes, multi_variant_records$Gene[i])
    }
  }
  
  unique(c(multi_sample_genes, multi_variant_genes))
}

target_genes <- get_multisample_genes(disease_variants)

# =============================================================================
# 5. Calculate Ranks
# =============================================================================

calculate_ranks <- function(intensity_data) {
  met_sample_names <- setdiff(colnames(intensity_data), c("XCMS.id", "model"))
  rank_data <- as.data.frame(t(apply(intensity_data[, met_sample_names], 1, rank)))
  rank_data$XCMS.id <- intensity_data$XCMS.id
  rank_data <- rank_data[, c(ncol(rank_data), 1:(ncol(rank_data) - 1))]
  return(rank_data)
}

dem_ranks <- calculate_ranks(dem_intensity)
non_dem_ranks <- calculate_ranks(non_dem_intensity)

# =============================================================================
# 6. Run Association Analysis
# =============================================================================

run_association_analysis <- function(genes, variant_data, intensity_data, rank_data) {
  results_list <- list()
  
  for (i in seq_along(genes)) {
    gene_name <- genes[i]
    gene_variants <- variant_data %>%
      dplyr::filter(Gene == gene_name)
    
    combined_samples <- paste(gene_variants$sample, collapse = ',')
    gene_variants$sample[1] <- combined_samples
    variant_record <- gene_variants[1, ]
    
    # Rank-based tests
    rank_results <- compare_metabolite_difference(variant_record, rank_data)
    rank_results <- rank_results[, 1:3]
    colnames(rank_results)[2:3] <- c('rank.scz_var_vs_nor_no_var.pvalue', 
                                      'rank.scz_var_vs_scz_no_var.pvalue')
    rank_results$rank.scz_var_vs_nor_no_var.padj <- p.adjust(
      rank_results$rank.scz_var_vs_nor_no_var.pvalue, method = 'BH'
    )
    rank_results$rank.scz_var_vs_scz_no_var.padj <- p.adjust(
      rank_results$rank.scz_var_vs_scz_no_var.pvalue, method = 'BH'
    )
    
    # Intensity-based tests
    intensity_results <- compare_metabolite_difference(variant_record, intensity_data)
    intensity_results <- intensity_results[, 1:3]
    colnames(intensity_results)[2:3] <- c('intensity.scz_var_vs_nor_no_var.pvalue',
                                          'intensity.scz_var_vs_scz_no_var.pvalue')
    intensity_results$intensity.scz_var_vs_nor_no_var.padj <- p.adjust(
      intensity_results$intensity.scz_var_vs_nor_no_var.pvalue, method = 'BH'
    )
    intensity_results$intensity.scz_var_vs_scz_no_var.padj <- p.adjust(
      intensity_results$intensity.scz_var_vs_scz_no_var.pvalue, method = 'BH'
    )
    
    combined_results <- merge(intensity_results, rank_results, by = 'XCMS.id')
    results_list[[i]] <- combined_results
    names(results_list)[i] <- gene_name
  }
  
  return(results_list)
}

# Run analysis separately for DEM and non-DEM
dem_association <- run_association_analysis(
  target_genes, disease_variants, dem_intensity, dem_ranks
)

non_dem_association <- run_association_analysis(
  target_genes, disease_variants, non_dem_intensity, non_dem_ranks
)

# =============================================================================
# 7. Result Summary and Statistics
# =============================================================================

summarize_significant_results <- function(dem_results, non_dem_results) {
  for (i in seq_along(dem_results)) {
    gene_name <- names(dem_results)[i]
    dem_res <- dem_results[[i]]
    non_dem_res <- non_dem_results[[i]]
    
    stats <- list(
      intensity_scz_vs_scz = list(
        dem = sum(dem_res$intensity.scz_var_vs_scz_no_var.padj < 0.05, na.rm = TRUE),
        non_dem = sum(non_dem_res$intensity.scz_var_vs_scz_no_var.padj < 0.05, na.rm = TRUE)
      ),
      intensity_scz_vs_nor = list(
        dem = sum(dem_res$intensity.scz_var_vs_nor_no_var.padj < 0.05, na.rm = TRUE),
        non_dem = sum(non_dem_res$intensity.scz_var_vs_nor_no_var.padj < 0.05, na.rm = TRUE)
      ),
      rank_scz_vs_scz = list(
        dem = sum(dem_res$rank.scz_var_vs_scz_no_var.padj < 0.05, na.rm = TRUE),
        non_dem = sum(non_dem_res$rank.scz_var_vs_scz_no_var.padj < 0.05, na.rm = TRUE)
      ),
      rank_scz_vs_nor = list(
        dem = sum(dem_res$rank.scz_var_vs_nor_no_var.padj < 0.05, na.rm = TRUE),
        non_dem = sum(non_dem_res$rank.scz_var_vs_nor_no_var.padj < 0.05, na.rm = TRUE)
      )
    )
    
    dem_shared <- dem_res %>%
      dplyr::filter(
        intensity.scz_var_vs_scz_no_var.padj < 0.05,
        intensity.scz_var_vs_nor_no_var.padj < 0.05,
        rank.scz_var_vs_scz_no_var.padj < 0.05,
        rank.scz_var_vs_nor_no_var.padj < 0.05
      )
    
    non_dem_shared <- non_dem_res %>%
      dplyr::filter(
        intensity.scz_var_vs_scz_no_var.padj < 0.05,
        intensity.scz_var_vs_nor_no_var.padj < 0.05,
        rank.scz_var_vs_scz_no_var.padj < 0.05,
        rank.scz_var_vs_nor_no_var.padj < 0.05
      )
    
    message(sprintf(
      "%s | Intensity(SCZ-mut vs SCZ-non): DEM=%d, NonDEM=%d | Intensity(SCZ-mut vs NOR): DEM=%d, NonDEM=%d | Rank(SCZ-mut vs SCZ-non): DEM=%d, NonDEM=%d | Rank(SCZ-mut vs NOR): DEM=%d, NonDEM=%d | Shared: DEM=%d, NonDEM=%d",
      gene_name,
      stats$intensity_scz_vs_scz$dem, stats$intensity_scz_vs_scz$non_dem,
      stats$intensity_scz_vs_nor$dem, stats$intensity_scz_vs_nor$non_dem,
      stats$rank_scz_vs_scz$dem, stats$rank_scz_vs_scz$non_dem,
      stats$rank_scz_vs_nor$dem, stats$rank_scz_vs_nor$non_dem,
      nrow(dem_shared), nrow(non_dem_shared)
    ))
  }
}

summarize_significant_results(dem_association, non_dem_association)

# =============================================================================
# 8. Save Results
# =============================================================================

save(dem_association, non_dem_association, file = VAR_MET_ASSOC_FILE)

write.table(dem_metabolites,           file = file.path(META_STAT_DIR, 'DEMs_with_annt.txt'),           quote = FALSE, row.names = FALSE)
write.table(non_dem_metabolites,       file = file.path(META_STAT_DIR, 'Non-DEMs_with_annt.txt'),       quote = FALSE, row.names = FALSE)
write.table(dem_metabolites_annt,      file = file.path(META_STAT_DIR, 'DEMs_with_annt.filter.txt'),  quote = FALSE, row.names = FALSE)
write.table(non_dem_metabolites_annt,  file = file.path(META_STAT_DIR, 'Non-DEMs_with_annt.filter.txt'), quote = FALSE, row.names = FALSE)

message("Variant-metabolite association analysis complete.")
message("Results saved to: ", VAR_MET_ASSOC_FILE)
