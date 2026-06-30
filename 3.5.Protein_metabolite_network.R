# =============================================================================
# 3.5.Protein_metabolite_network.R
# Protein-Metabolite Network Construction
# =============================================================================
# Combines variant-metabolite associations and metabolite-metabolite correlations.
#
# Input:  VAR_MET_ASSOC_FILE, DEM_COR_FILE, gene scores, metabolite annotations
# Output: NETWORK_EDGE_OUTPUT, NETWORK_NODE_OUTPUT
# =============================================================================

source("config.R")
source("utils.R")

library(magrittr)
library(dplyr)

# =============================================================================
# 1. Load Data
# =============================================================================

# Variant-metabolite associations (dem_association, non_dem_association)
load(VAR_MET_ASSOC_FILE)

# Metabolite correlation (dem_cor_flat)
load(DEM_COR_FILE)

# Select significant metabolite correlation pairs
sig_cor_pair <- subset(dem_cor_flat, q < 0.05 & abs(cor) >= COR_ABS_THRESHOLD)

# Select pairs containing at least one DEM
# dem_association is a list; extract all XCMS.id from all gene results
all_dem_ids <- unique(unlist(lapply(dem_association, function(x) x$XCMS.id)))
dem_cor_pair <- rbind(
  subset(sig_cor_pair, row %in% all_dem_ids),
  subset(sig_cor_pair, column %in% all_dem_ids)
) %>% unique()

# Gene scores
gene_scores <- read.csv(GENE_SCORE_FILE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Metabolite annotations (output from 3.3)
dem_metabolites_annt <- read.csv(file.path(META_STAT_DIR, 'DEMs_with_annt.filter.txt'),
                                  sep = '\t', stringsAsFactors = FALSE, row.names = 1)
non_dem_metabolites_annt <- read.csv(file.path(META_STAT_DIR, 'Non-DEMs_with_annt.filter.txt'),
                                      sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# NOTE: The following objects are expected to exist in the R environment
# or be loaded from prior analysis steps. They are not defined in this script
# because they were created interactively during the original pipeline.
#   special_id_3   : mapping of concatenated XCMS.id to MS2Metabolite
#   dem_with_annt  : DEM annotation table (used for name resolution)
#   non_dem_with_annt : non-DEM annotation table (used for name resolution)
# If these are not in the environment, load them from saved files before running.

# =============================================================================
# 2. Extract Significant Variant-Metabolite Associations
# =============================================================================

var_met_pair <- data.frame()
gene_list <- names(dem_association)

for (i in seq_along(gene_list)) {
  symbol <- subset(gene_scores, ENSEMBL == gene_list[i])$SYMBOL
  
  # Extract variant-DEM association pairs
  dem_edge <- dem_association[[gene_list[i]]]
  dem_edge$reg <- apply(dem_edge, 1, function(x) {
    y <- ifelse(as.numeric(x[4]) > as.numeric(x[5]), 'Up', 'Down')
    return(y)
  })
  
  # NOTE: Column names here match the output of 3.3.Variant_metabolite_assoc.R
  # Original 3.5 had a bug using "scz_varvsnor_novar" (no underscores).
  # Fixed to use "scz_var_vs_nor_no_var" (with underscores) to match 3.3 output.
  dem_edge_kpt <- dem_edge %>% 
    dplyr::filter(intensity.scz_var_vs_nor_no_var.padj < 0.05) %>%
    dplyr::filter(intensity.scz_var_vs_scz_no_var.padj < 0.05) %>%
    dplyr::filter(rank.scz_var_vs_nor_no_var.padj < 0.05) %>%
    dplyr::filter(rank.scz_var_vs_scz_no_var.padj < 0.05)
  
  dem_node <- dem_edge_kpt$'XCMS.id'
  if (length(dem_node) > 0) {
    dem_reg <- dem_edge_kpt$'reg'
    dem_net <- data.frame(
      node1 = rep(symbol, length(dem_node)), 
      node1_type = rep('MutGene', length(dem_node)),
      node2 = dem_node, 
      node2_type = rep('DEM', length(dem_node)),
      edge_type = 'PMI-var', 
      regulation = dem_reg
    )
  } else {
    dem_net <- data.frame()
  }
  
  # Extract variant-nonDEM association pairs
  nondem_edge <- non_dem_association[[i]]
  nondem_edge$reg <- apply(nondem_edge, 1, function(x) {
    y <- ifelse(as.numeric(x[4]) > as.numeric(x[5]), 'Up', 'Down')
    return(y)
  })
  
  nondem_edge_kpt <- nondem_edge %>% 
    dplyr::filter(intensity.scz_var_vs_nor_no_var.padj < 0.05) %>%
    dplyr::filter(intensity.scz_var_vs_scz_no_var.padj < 0.05) %>%
    dplyr::filter(rank.scz_var_vs_nor_no_var.padj < 0.05) %>%
    dplyr::filter(rank.scz_var_vs_scz_no_var.padj < 0.05) %>%
    dplyr::filter(XCMS.id %in% unique(c(dem_cor_pair$row, dem_cor_pair$column)))
  
  nondem_node <- nondem_edge_kpt$'XCMS.id'
  if (length(nondem_node) > 0) {
    nondem_reg <- nondem_edge_kpt$reg
    nondem_net <- data.frame(
      node1 = rep(symbol, length(nondem_node)), 
      node1_type = rep('MutGene', length(nondem_node)),
      node2 = nondem_node,
      node2_type = rep('NonDEM', length(nondem_node)),
      edge_type = 'PMI-var',
      regulation = nondem_reg
    )
  } else {
    nondem_net <- data.frame()
  }
  
  var_met_pair <- rbind(var_met_pair, dem_net, nondem_net) %>% unique()
}

# =============================================================================
# 3. Build Metabolite-Metabolite Pairs
# =============================================================================

met_met_pair <- data.frame(
  node1 = dem_cor_pair$row,
  node1_type = rep('-', nrow(dem_cor_pair)),
  node2 = dem_cor_pair$column,
  node2_type = rep('-', nrow(dem_cor_pair)),
  edge_type = 'MMR-cor',
  stringsAsFactors = FALSE
)

# Classify node types using utils.R helper
dem_ids <- paste0(dem_metabolites_annt$XCMS.id, '_', dem_metabolites_annt$model)
nondem_ids <- paste0(non_dem_metabolites_annt$XCMS.id, '_', non_dem_metabolites_annt$model)

met_met_pair$node1_type <- classify_metabolite_node(met_met_pair$node1, dem_ids, nondem_ids)
met_met_pair$node2_type <- classify_metabolite_node(met_met_pair$node2, dem_ids, nondem_ids)

met_met_pair$regulation <- apply(dem_cor_pair, 1, function(x) {
  ifelse(as.numeric(x[3]) > 0, 'Pos', 'Neg')
})

# =============================================================================
# 4. Combine Network and Resolve Names
# =============================================================================

net <- rbind(var_met_pair, met_met_pair)

for (i in seq_len(nrow(net))) {
  tmp <- net[i, ]
  
  # Node 1
  if (grepl('Gene', tmp$node1_type)) {
    node1_id <- subset(gene_scores, SYMBOL == tmp$node1)$ENTREZID
  } else {
    node1_id <- tmp$node1
    if (grepl(';', node1_id)) {
      net[i, 'node1'] <- subset(special_id_3, XCMS.id == node1_id)$MS2Metabolite
    } else {
      x <- substr(tmp$node1, 1, nchar(tmp$node1) - 9)
      y <- substr(tmp$node1, nchar(tmp$node1) - 8 + 1, nchar(tmp$node1))
      net[i, 'node1'] <- subset(dem_with_annt, XCMS.id == x & model == y)$MS2Metabolite
    }
  }
  
  # Node 2
  node2_id <- tmp$node2
  if (grepl(';', node2_id)) {
    net[i, 'node2'] <- subset(special_id_3, XCMS.id == node2_id)$MS2Metabolite
  } else {
    m <- substr(tmp$node2, 1, nchar(tmp$node2) - 9)
    n <- substr(tmp$node2, nchar(tmp$node2) - 8 + 1, nchar(tmp$node2))
    if (tmp$node2_type == 'DEM') {
      net[i, 'node2'] <- subset(dem_with_annt, XCMS.id == m & model == n)$MS2Metabolite
    } else {
      net[i, 'node2'] <- subset(non_dem_with_annt, XCMS.id == m & model == n)$MS2Metabolite
    }
  }
  net[i, 7:8] <- c(node1_id, node2_id)
}

colnames(net)[c(1, 3, 7:8)] <- c('node1_name', 'node2_name', 'node1_id', 'node2_id')
net <- net[, c(7:8, 1, 3, 2, 4:6)]

# =============================================================================
# 5. Build Node Metadata
# =============================================================================

net_node <- data.frame(
  node_id = c(net$node1_id, net$node2_id),
  node_name = c(net$node1_name, net$node2_name),
  node_type = c(net$node1_type, net$node2_type)
) %>% unique()

net_node$metabolite_class <- net_node$super_class <- net_node$value <- '-'

add <- apply(net_node, 1, function(x) {
  if (x[3] == 'MutGene') {
    y <- subset(gene_scores, SYMBOL == x[2])$GS
    m <- '-'
    n <- '-'
  } else {
    id <- substr(x[1], 1, nchar(x[1]) - 9)
    mod <- substr(x[1], nchar(x[1]) - 8 + 1, nchar(x[1]))
    y <- ifelse(x[3] == 'DEM',
                subset(dem_metabolites_annt, XCMS.id == id & model == mod)$'SCZ_allvsNOR.FC',
                subset(non_dem_metabolites_annt, XCMS.id == id & model == mod)$'SCZ_allvsNOR.FC')
    m <- ifelse(x[3] == 'DEM',
                subset(dem_metabolites_annt, XCMS.id == id & model == mod)$'MS2superclass',
                subset(non_dem_metabolites_annt, XCMS.id == id & model == mod)$'MS2superclass')
    n <- ifelse(x[3] == 'DEM',
                subset(dem_metabolites_annt, XCMS.id == id & model == mod)$'MS2class',
                subset(non_dem_metabolites_annt, XCMS.id == id & model == mod)$'MS2class')
  }
  return(c(y, m, n))
})

net_node[, 4:6] <- t(add)
net_node <- as.data.frame(net_node)

# =============================================================================
# 6. Save Outputs
# =============================================================================

write.table(net,     file = NETWORK_EDGE_OUTPUT, row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(net_node, file = NETWORK_NODE_OUTPUT, row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

message("Network construction complete.")
message("Edges: ", nrow(net))
message("Nodes: ", nrow(net_node))
