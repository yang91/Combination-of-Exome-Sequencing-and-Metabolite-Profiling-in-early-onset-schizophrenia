
# =============================================================================
# 0. Configuration and Path Setup
# =============================================================================

META_STAT_DIR  <- "~/EOSCZ/MetaSeq-Result/statistics"
CLIN_DIR <- "~/EOSCZ/Clinical_Info/"
NET_DIR <- "~/EOSCZ/Network_file/"

# =============================================================================
# 1. Data Loading and Preprocessing
# =============================================================================
# Variant-metabolite association (dem_association, non_dem_association)
load(file.path(NET_DIR, "Variant-metabolite.association.RData"))

# Metabolite correlation (dem_cor_flat)
load(file.path(NET_DIR, "DEM_NonDEM.correlation.RData"))
sig_cor_pair <- subset(dem_cor_flat, q<0.05&abs(cor)>=0.9) # select significant metabolite correlation pairs
# Select the correlation pairs containing at least one DEM
dem_cor_pair <- rbind(subset(sig_cor_pair, row %in% dem_association), 
                      subset(sig_cor_pair, column %in% dem_association)) %>% unique()

# Gene score evaluation model from an in-house script of unpublished data
gs_file <- file.path(EXOME_DIR, "gene_score.high-risk.txt")
gene_scores <- read.csv(gs_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# DEMs and Non-DEMannotations with accurate annotation
dem_annt_file <- file.path(META_STAT_DIR,'DEMs_with_annt.filter.txt')
dem_metabolites_annt <- read.csv(dem_annt_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

non_dem_annt_file <- file.path(META_STAT_DIR,'Non-DEMs_with_annt.filter.txt')
non_dem_metabolites_annt <- read.csv(non_dem_annt_file, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# =============================================================================
# 2. Extract siginificant variant-metabolite association
# =============================================================================
var_met_pair <- data.frame()
gene_list <- names(dem_association)

for(i in 1:length(gene_list)){
  symbol <- subset(gene_scores,ENSEMBL==gene_list[i])$SYMBOL
  
  # Extract variant-DEM association pairs
  dem_edge <- dem_association[[gene_list[i]]]
  dem_edge$reg <- apply(dem_edge,1,function(x){
    y <- ifelse(x[4]>x[5],'Up','Down')
    return(y)
  })
  dem_edge_kpt <- dem_edge %>% 
    dplyr::filter(intensity.scz_varvsnor_novar.padj<0.05) %>%
    dplyr::filter(intensity.scz_varvsscz_novar.padj<0.05) %>%
    dplyr::filter(rank.scz_varvsnor_novar.padj<0.05) %>%
    dplyr::filter(rank.scz_varvsscz_novar.padj<0.05)
  dem_node <- dem_edge_kpt$'XCMS.id'
  if(length(dem_node)>0){
    dem_reg <- dem_edge_kpt$'reg'
    dem_net <- data.frame(node1 = rep(symbol, length(dem_node)), 
                          node1_type = rep('MutGene', length(dem_node)),
                          node2 = dem_node, 
                          node2_type = rep('DEM',length(dem_node)),
                          edge_type = 'PMI-var', 
                          regulation = dem_reg)
  }

  # Extract variant-nonDEM association pairs
  nondem_edge <- non_dem_association[[i]]
  nondem_edge$reg <- apply(nondem_edge,1,function(x){
    y <- ifelse(x[4]>x[5],'Up','Down')
    return(y)
    })
  nondem_edge_kpt <- nondem_edge %>% 
    dplyr::filter(intensity.scz_varvsnor_novar.padj<0.05) %>%
    dplyr::filter(intensity.scz_varvsscz_novar.padj<0.05) %>%
    dplyr::filter(rank.scz_varvsnor_novar.padj<0.05) %>%
    dplyr::filter(rank.scz_varvsscz_novar.padj<0.05) %>%
    dplyr::filter(XCMS.id %in% unique(c(dem_cor_pair$row, dem_cor_pair$column)))
  nondem_node <- nondem_edge_kpt$'XCMS.id'
  if(length(nondem_node)>0){
    nondem_reg <- nondem_edge_kpt$reg
    nondem_net <- data.frame(node1 = rep(symbol, length(nondem_node)), 
                           node1_type = rep('MutGene', length(nondem_node)),
                           node2 = nondem_node,
                           node2_type = rep('NonDEM', length(nondem_node)),
                           edge_type = 'PMI-var',
                           regulation = nondem_reg)
  }

  var_met_pair <- rbind(var_met_pair,dem_net,nondem_net) %>% unique()
}

met_met_pair <- data.frame(node1 = dem_cor_pair$row, node1_type = rep('-',dim(dem_cor_pair)[1]), 
                           node2 = dem_cor_pair$column, node2_type = rep('-',dim(dem_cor_pair)[1]),
                           edge_type = 'MMR-cor')
met_met_pair[met_met_pair$node1 %in% paste0(dem_metabolites_annt$XCMS.id,'_',dem_metabolites_annt$model),'node1_type'] <- 'DEM'
met_met_pair[met_met_pair$node1 %in% paste0(non_dem_metabolites_annt$XCMS.id,'_',non_dem_metabolites_annt$model),'node1_type'] <- 'NonDEM'
met_met_pair[met_met_pair$node2 %in% paste0(dem_metabolites_annt$XCMS.id,'_',dem_metabolites_annt$model),'node2_type'] <- 'DEM'
met_met_pair[met_met_pair$node2 %in% paste0(non_dem_metabolites_annt$XCMS.id,'_',non_dem_metabolites_annt$model),'node2_type'] <- 'NonDEM'
met_met_pair$regulation <- apply(dem_cor_pair,1,function(x){
  y <- ifelse(as.numeric(x[3])>0,'Pos','Neg')
  return(y)
})

net <- rbind(var_met_pair, met_met_pair)
for(i in 1:nrow(net)){
  tmp <- net[i,]
  node1_id <- vector()
  node2_id <- vector()
  
  # Nde1_type only have MutGene and DEM
  if(grepl('Gene',tmp$node1_type)){
    node1_id <- subset(gene_scores,SYMBOL==tmp$node1)$ENTREZID
  }  
  else{
    node1_id <- tmp$node1
    if(grepl(';',node1_id)){
      net[i,'node1'] <- subset(special_id_3,XCMS.id==node1_id)$MS2Metabolite
    }
    else{
      x <- substr(tmp$node1,1,nchar(tmp$node1)-9)
      y <- substr(tmp$node1,nchar(tmp$node1)-8+1,nchar(tmp$node1))
      net[i,'node1'] <- subset(dem_with_annt,XCMS.id==x&model==y)$MS2Metabolite
    }
  }
  
  # Node2_type only have DEM and NonDEM
  node2_id <- tmp$node2
  if(grepl(';',node2_id)){
    net[i,'node2'] <- subset(special_id_3,XCMS.id==node2_id)$MS2Metabolite
  }
  else{
    m <- substr(tmp$node2,1,nchar(tmp$node2)-9)
    n <- substr(tmp$node2,nchar(tmp$node2)-8+1,nchar(tmp$node2))
    if(tmp$node2_type=='DEM'){
      net[i,'node2'] <- subset(dem_with_annt,XCMS.id==m&model==n)$MS2Metabolite
    }
    else{
      net[i,'node2'] <- subset(non_dem_with_annt,XCMS.id==m&model==n)$MS2Metabolite
    }
  }
  net[i,7:8] <- c(node1_id, node2_id)
}
colnames(net)[c(1,3,7:8)] <- c('node1_name','node2_name','node1_id','node2_id')
net <- net[,c(7:8,1,3,2,4:6)]


net_node <- data.frame(node_id = c(net$node1_id,net$node2_id),
                        node_name = c(net$node1_name, net$node2_name),
                        node_type = c(net$node1_type, net$node2_type)) %>% unique()
net_node$metabolite_class <- net_node$super_class <- net_node$value <- '-'
add <- apply(net_node,1,function(x){
  if(x[3]=='MutGene'){
    y <- subset(gene_scores,SYMBOL==x[2])$GS
    m <- '-'
    n <- '-'
  }
  else{
    id <- substr(x[1],1,nchar(x[1])-9)
    mod <- substr(x[1],nchar(x[1])-8+1,nchar(x[1]))
    y <- ifelse(x[3]=='DEM',
                subset(dem_metabolites_annt,XCMS.id==id&model==mod)$'SCZ_allvsNOR.FC',
                subset(non_dem_metabolites_annt,XCMS.id==id&model==mod)$'SCZ_allvsNOR.FC')
    m <- ifelse(x[3]=='DEM',
                   subset(dem_metabolites_annt,XCMS.id==id&model==mod)$'MS2superclass',
                   subset(non_dem_metabolites_annt,XCMS.id==id&model==mod)$'MS2superclass')
    n <- ifelse(x[3]=='DEM',
                   subset(dem_metabolites_annt,XCMS.id==id&model==mod)$'MS2class',
                   subset(non_dem_metabolites_annt,XCMS.id==id&model==mod)$'MS2class')
  }
  return(c(y,m,n))
  })
net_node[,4:6] <- t(add)

table(net_node$node_type)
net_node <- as.data.frame(net_node)

net_file <- file.path(NET_DIR,'HighRisk_gene.Annotated_metabolites.network_edge.txt')
node_file <- file.path(NET_DIR,'HighRisk_gene.Annotated_metabolites.network_node.txt')
write.table(net,net_file, row.names = F, col.names = T, sep = '\t', quote = F)
write.table(net_node, node_file, row.names = F, col.names = T, sep = '\t', quote = F)
