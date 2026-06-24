# =============================================================================
# Variant-Metabolite Association Analysis
# =============================================================================

library(magrittr)
library(Hmisc)
library(dplyr)

# =============================================================================
# 0. Configuration and Path Setup
# =============================================================================

EXOME_DIR <- "~/EOSCZ/ExomeSeq-Result/"
VAR_DIR <- "~/EOSCZ/ExomeSeq-Result/variant/"
META_STAT_DIR  <- "~/EOSCZ/MetaSeq-Result/statistics/"
CLIN_DIR <- "~/EOSCZ/Clinical_Info/"
NET_DIR <- "~/EOSCZ/Network_file/"

# =============================================================================
# 1. Data Loading and Preprocessing
# =============================================================================

# --- 1.1 Load gene score data ---
# Gene score evaluation model from an in-house script of unpublished data
gs_file <- file.path(EXOME_DIR, "gene_score.high-risk.txt")
gene_scores <- read.csv(gs_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Load disease-causing variants with OR > 1 in genes with GS > 4
var_file <- file.path(VAR_DIR, "disease_causing_ORgt1_variants_in_GSgt4_gene.txt")
disease_variants <- read.csv(var_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Load VEP annotation results (skip first 44 header lines)
vep_file <- file.path(VAR_DIR, "Maf_0.05.potential_harmful.ORgt1_variants.txt")
vep_results <- read.csv(vep_file, sep = '\t', stringsAsFactors = FALSE, skip = 44, row.names = 1)

# --- 1.2 Load sample information ---
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

# =============================================================================
# 4. Metabolite Data Loading and Processing
# =============================================================================

# --- 4.0 Function to select metabolites with accurate annotation ---
select_annt <- function(metabolite){
  ms1_ant <- rbind(metabolite %>% dplyr::filter(NumberMS1hmdb>0),
                 metabolite %>% dplyr::filter(NumberMS1kegg>0)) %>% unique()
              
  ms2_ant <- metabolite %>% dplyr::filter(MS2Metabolite!='-')
  ms2_ant_only <- ms2_ant %>% dplyr::filter(NumberMS1hmdb==0) %>% dplyr::filter(NumberMS1kegg==0) 

  ms1_2_ant <- rbind(ms1_ant %>% dplyr::filter(MS2Metabolite!='-'),
                   ms2_ant %>% dplyr::filter(NumberMS1hmdb>0),
                   ms2_ant %>% dplyr::filter(NumberMS1kegg>0)) %>% unique()
  ms1_2_ant_match <- vector()
  n <- vector()
  for(i in 1:dim(ms1_2_ant)[1]){
    tmp <- ms1_2_ant[i,]
    j <- 0
    if(grepl("tmp$MS2Metabolite", tmp$MS1hmdbName, perl = T)){
      j <- j+1
    }
    if(grepl("tmp$MS2Metabolite", tmp$MS1keggName, perl = T)){
      j <- j+1
    }
    if(tmp$MS2hmd!='-' && tmp$MS2hmd!='' && !is.na(tmp$MS2hmd) && grepl(tmp$MS2hmd, tmp$MS1hmdbID)){
      j <- j+1
    }
    if(tmp$MS2kegg!='-' && tmp$MS2kegg!='' && !is.na(tmp$MS2kegg) && grepl(tmp$MS2kegg, tmp$MS1hmdbTokegg)){
      j <- j+1
    }
    if(tmp$MS2kegg!='-' &&  tmp$MS2kegg!='' && !is.na(tmp$MS2kegg) &&grepl(tmp$MS2kegg, tmp$MS1keggID)){
      j <- j+1
    }
  
    if(j>0){  n[i] <- 'TRUE'  }
    else{ n[i] <- 'FALSE'  }
  }
  ms1_2_ant_match <- ms1_2_ant[as.logical(n),]
  res <- rbind(ms2_ant_only, ms1_2_ant_match)
  return(res)
}


# --- 4.1 Load positive and negative mode metabolite data ---
META_STAT_DIR
pos_file <- file.path(META_STAT_DIR, "positive.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")
neg_file <- file.path(META_STAT_DIR, "negative.Two_group.Univariate-t.Multivariat_lm.with_annt.txt")
pos_metabolites <- read.csv(pos_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
neg_metabolites <- read.csv(neg_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Combine and select key columns
all_metabolites <- rbind(pos_metabolites, neg_metabolites)
selected_columns <- c(1, 4, 13, 34, 39:106)  # XCMS.id, model, annotation, p-value, sample intensities
all_metabolites <- all_metabolites[, selected_columns]

# --- 4.2 Separate differentially expressed metabolites (DEM) and non-DEM ---
dem_metabolites <- all_metabolites %>%
  dplyr::filter(All_sample.LM.padj < 0.05)

non_dem_metabolites <- all_metabolites %>%
  dplyr::filter(!XCMS.id %in% dem_metabolites$XCMS.id)

# --- 4.3 Merge metabolites with duplicate annotations (average intensities) ---
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
dem_metabolites_annt <- select_annt(dem_metabolites)
non_dem_metabolites_annt <- select_annt(non_dem_metabolites)

# --- 4.4 Build analysis data frames ---
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
# 5. Core Analysis Function: Compare Metabolite Differences
#    Between Variant Carriers and Non-Carriers
# =============================================================================

#' Compare metabolite differences between variant carriers and non-carriers
#' 
#' @param variant_info Data frame containing variant sample information
#' @param metabolite_data Metabolite data (intensity or rank values)
#' @return Data frame with statistical test results
compare_metabolite_difference <- function(variant_info, metabolite_data) {
  
  # Extract sample groups
  carrier_samples <- unlist(strsplit(unique(variant_info$sample), split = ','))
  scz_samples <- sample_names[grep('SCZ', sample_names)]
  nor_samples <- sample_names[grep('NOR', sample_names)]
  
  # Define sample groups
  scz_carriers <- carrier_samples[grep('SCZ', carrier_samples)]
  scz_non_carriers <- scz_samples[!(scz_samples %in% scz_carriers)]
  
  if (any(grepl('NOR', variant_info$sample))) {
    nor_carriers <- carrier_samples[grep('NOR', carrier_samples)]
    nor_non_carriers <- nor_samples[!(nor_samples %in% nor_carriers)]
  } else {
    nor_non_carriers <- nor_samples
  }
  
  # Perform statistical tests
  results <- data.frame(
    XCMS.id = character(nrow(metabolite_data)),
    scz_var_vs_nor_no_var.pvalue = numeric(nrow(metabolite_data)),
    scz_var_vs_scz_no_var.pvalue = numeric(nrow(metabolite_data)),
    stringsAsFactors = FALSE
  )
  
  for (j in seq_len(nrow(metabolite_data))) {
    # Comparison 1: SCZ variant carriers vs healthy controls (non-carriers)
    test1 <- t.test(
      as.numeric(metabolite_data[j, scz_carriers]),
      as.numeric(metabolite_data[j, nor_non_carriers])
    )
    
    # Comparison 2: SCZ variant carriers vs SCZ non-carriers
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
  
  # Convert numeric columns
  results$scz_var_vs_nor_no_var.pvalue <- as.numeric(results$scz_var_vs_nor_no_var.pvalue)
  results$scz_var_vs_scz_no_var.pvalue <- as.numeric(results$scz_var_vs_scz_no_var.pvalue)
  
  return(results)
}

# =============================================================================
# 6. Filter Target Genes (Genes with variants in >= 2 SCZ samples)
# =============================================================================

# Identify genes carrying variants in multiple SCZ samples
get_multisample_genes <- function(disease_variants) {
  
  # Case 1: Same gene appears in multiple distinct SCZ samples (non-duplicated samples)
  duplicated_genes <- disease_variants$Gene[duplicated(disease_variants$Gene)]
  multi_sample_genes <- c()
  
  for (gene in unique(duplicated_genes)) {
    gene_samples <- disease_variants[disease_variants$Gene == gene, 'sample']
    
    # Check if only in SCZ samples and samples are non-duplicated
    has_no_nor <- !any(grepl('NOR', gene_samples))
    has_unique_samples <- !any(duplicated(gene_samples))
    
    if (has_no_nor && has_unique_samples) {
      multi_sample_genes <- c(multi_sample_genes, gene)
    }
  }
  
  # Case 2: Single sample carries multiple variants in same gene (comma-separated sample names)
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
# 7. Association Analysis: Gene Variants and Metabolites
# =============================================================================

# --- 7.1 Calculate ranks (for non-parametric comparison) ---
calculate_ranks <- function(intensity_data) {
  rank_data <- as.data.frame(t(apply(intensity_data[, sample_names], 1, rank)))
  rank_data$XCMS.id <- intensity_data$XCMS.id
  rank_data <- rank_data[, c(ncol(rank_data), 1:(ncol(rank_data) - 1))]
  return(rank_data)
}

dem_ranks <- calculate_ranks(dem_intensity)
non_dem_ranks <- calculate_ranks(non_dem_intensity)

# --- 7.2 Run association analysis ---
run_association_analysis <- function(genes, variant_data, intensity_data, rank_data) {
  
  results_list <- list()
  
  for (i in seq_along(genes)) {
    gene_name <- genes[i]
    
    # Get all variant samples for this gene
    gene_variants <- variant_data %>%
      dplyr::filter(Gene == gene_name)
    
    # Combine all samples (comma-separated)
    combined_samples <- paste(gene_variants$sample, collapse = ',')
    gene_variants$sample[1] <- combined_samples
    
    # Use first row as representative for testing
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
    intensity_results <- compare_metabolite_difference(variant_record, intensity_data[, 1:49])
    intensity_results <- intensity_results[, 1:3]
    colnames(intensity_results)[2:3] <- c('intensity.scz_var_vs_nor_no_var.pvalue',
                                            'intensity.scz_var_vs_scz_no_var.pvalue')
    intensity_results$intensity.scz_var_vs_nor_no_var.padj <- p.adjust(
      intensity_results$intensity.scz_var_vs_nor_no_var.pvalue, method = 'BH'
    )
    intensity_results$intensity.scz_var_vs_scz_no_var.padj <- p.adjust(
      intensity_results$intensity.scz_var_vs_scz_no_var.pvalue, method = 'BH'
    )
    
    # Merge results
    combined_results <- merge(
      intensity_results, rank_results,
      by = 'XCMS.id'
    )
    
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
# 8. Result Summary and Statistics
# =============================================================================

summarize_significant_results <- function(dem_results, non_dem_results) {
  
  for (i in seq_along(dem_results)) {
    gene_name <- names(dem_results)[i]
    dem_res <- dem_results[[i]]
    non_dem_res <- non_dem_results[[i]]
    
    # Count significant results per condition (padj < 0.05)
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
    
    # Count metabolites significant in all four conditions
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
    
    # Output summary
    message(sprintf(
      paste0(
        "%s | ",
        "Intensity(SCZ-mut vs SCZ-non): DEM=%d, NonDEM=%d | ",
        "Intensity(SCZ-mut vs NOR): DEM=%d, NonDEM=%d | ",
        "Rank(SCZ-mut vs SCZ-non): DEM=%d, NonDEM=%d | ",
        "Rank(SCZ-mut vs NOR): DEM=%d, NonDEM=%d | ",
        "Shared: DEM=%d, NonDEM=%d"
      ),
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
# 9. Save Results
# =============================================================================
output_file <- file.path(NET_DIR, "Variant-metabolite.association.RData")
save(dem_association, non_dem_association, file = output_file)

dem_file <- file.path(META_STAT_DIR,'DEMs_with_annt.txt')
write.table(dem_metabolites, file = dem_file, quote = F, row.names = F)

non_dem_file <- file.path(META_STAT_DIR,'Non-DEMs_with_annt.txt')
write.table(non_dem_metabolites, file = non_dem_file, quote = F, row.names = F)

dem_annt_file <- file.path(META_STAT_DIR,'DEMs_with_annt.filter.txt')
write.table(dem_metabolites_annt, file = dem_annt_file, quote = F, row.names = F)

non_dem_annt_file <- file.path(META_STAT_DIR,'Non-DEMs_with_annt.filter.txt')
write.table(non_dem_metabolites_annt, file = non_dem_annt_file, quote = F, row.names = F)
