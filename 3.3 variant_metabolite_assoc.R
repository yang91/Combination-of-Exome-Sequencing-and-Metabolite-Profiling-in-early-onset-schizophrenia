library(magrittr)
library(Hmisc)
setwd("~/EOSCZ/")

#### 0. Data processing ####
# read the gene score file. The gene score evaluation model was our in-house script from one unpublished data.
gs <- read.csv('./ExomeSeq-Result/gene_score.high-risk.txt', header = T, sep = '\t', stringsAsFactors = F)
# read the OR>1 disease causing variants in GS>4 genes
gs_var <- read.csv('/ExomeSeq-Result/disease_causing_ORgt1_variants_in_GSgt4_gene.txt', header = T, sep = '\t', stringsAsFactors = F)
# read the raw vep result to extract variants corresponding samples
vep <- read.csv('./ExomeSeq-Result/variant/Maf_0.05.potential_harmful.ORgt1_variants.txt', header = T, sep = '\t', skip = 44, stringsAsFactors = F)

# read early-onset schizophrenia samples 
sp_info <- read.csv('./clinical_info/sample_information.csv', header = T, stringsAsFactors = F)
sp_info <- sp_info[sp_info$group!='QC',]
sp_info <- sp_info[order(sp_info$group,decreasing = F),]
rownames(sp_info) <- sp_info$sample.name

# Process data based on variant-based level
vep_gs_var <- subset(vep, (Gene %in% gs_var$Gene)&(X.Uploaded_variation %in% gs_var$X.Uploaded_variation))
vep_list <- unique(vep_gs_var$X.Uploaded_variation)
var_ipt <- data.frame() # The data frame to store the samples containing variants
var_eft <- data.frame() # The data frame to store the prediction effect of variants. If one variant have several effects in different transcripts, record the most severe one.
var_ipt[1,1:48] <- var_eft[1,1:48] <- 0
colnames(var_ipt) <- colnames(var_eft) <- sp_info$sample.name
for(i in 1:length(vep_list)){
  var_ipt[i,1:48] <- '0'
  var_eft[i,1:48] <- '0'
  tmp <- subset(vep_gs_var,X.Uploaded_variation==vep_list[i])
  vt <- apply(tmp,1,function(x){
    y <- unlist(strsplit(x,';'))
    z1 <- y[grepl('IND',y)]
    z1 <- gsub('IND=','',z1)
    z2 <- y[grepl('IMPACT',y)]
    z2 <- gsub('IMPACT=','',z2)
    vt <- data.frame(sample = z1, conseq = z2)
    return(unlist(vt))
  })
  sp <- unique(vt[1,])
  vt[2,] <- gsub('HIGH',4,vt[2,])
  vt[2,] <- gsub('MODERATE',3,vt[2,])
  vt[2,] <- gsub('LOW',2,vt[2,])
  vt[2,] <- gsub('MODIFIER',1,vt[2,])
  if(length(sp)==1){
    if(sp %in% sp_info$sample.name){
      var_ipt[i,sp] <- 1
      var_eft[i,sp] <- max(vt[2,])
    }
  }
  else{
    for(j in 1:length(sp)){
      if(sp[j] %in% sp_info$sample.name){
        var_ipt[i,sp[j]] <- 1
        var_eft[i,sp[j]] <- max(vt[2,grepl(sp[j],vt[1,])])
      }
    }
  }
}
rownames(var_eft) <- rownames(var_ipt) <- vep_list

# Collapse the variant information to gene-based level
gene_ipt <- data.frame() # The data frame to store the samples containing either variants in corresponding gene
gene_eft <- data.frame() # The data frame to store the prediction effect of variants inn one gene (record the most severe one).
gene_list <- unique(gs_var$Gene)
for(i in 1:length(gene_list)){
  tmp <- subset(gs_var,Gene==gene_list[i])$X.Uploaded_variation
  gene_ipt[i,1:48] <- apply(var_ipt[tmp,],2,function(x){
    y <- sum(as.numeric(x));return(y)
  })
  gene_eft[i,1:48] <- apply(var_eft[tmp,],2,function(x){
    y <- sum(as.numeric(x));return(y)
  })
}
colnames(gene_eft) <- colnames(gene_ipt) <- colnames(var_eft)
rownames(gene_eft) <- rownames(gene_ipt) <- gene_list
for(i in 1:dim(gene_ipt)[2]){
  if(length(which(gene_ipt[,i]>1))>1){
    k <- which(gene_ipt[,i]>1)
    sp <- colnames(gene_ipt)[i]
    for(j in 1:length(k)){
      gn <- rownames(gene_ipt[k[j],i])
      tmp <- subset(gs_var,Gene==gn)$X.Uploaded_variation
      gene_ipt[k[j],i] <- 1
      gene_eft[k[j],i] <- max(as.numeric(var_eft[tmp,sp]))
    }
  }
  else{
    next;
  }
}
gs_var$sample <- apply(gs_var, 1, function(x){
  tmp <- vep[vep$X.Uploaded_variation==x[1],]
  sp <- apply(tmp,1,function(y){
    z <- strsplit(y[14], split = ';')
    z <- unlist(z)
    return(z[1])
  })
  sp <- gsub('IND=','',sp)
  sp <- sp[sp %in% sp_info$sample.name]
  sp <- paste0(unique(sp),collapse = ',')
  return(sp)
})

# Remove the redundant metbolite intensities. For the metabolite with several peaks, calculate the average intensity for further analysis.
# In this step, all metabolites were consider to compute the association between gene variants and metabolite intensities, no matter if this metabolite has accurate annotation.
pos <- read.csv('./MetaSeq-Result/statistics/positive.Two_group.Univariate-t.Multivariat_lm.with_annt.txt', header = T, sep = '\t', stringsAsFactors = F)
neg <- read.csv('./MetaSeq-Result/statistics/negative.Two_group.Univariate-t.Multivariat_lm.with_annt.txt', header = T, sep = '\t', stringsAsFactors = F)
all_met <- rbind(pos,neg)
all_met <- all_met[,c(1,4,13,34,39:106)]
dem <- subset(all_met,All_sample.LM.padj<0.05)
non_dem <- subset(all_met, !(XCMS.id %in% dem$XCMS.id))

dup <- dem$MS2Metabolite[duplicated(dem$MS2Metabolite)] %>% unique()
dup <- dup[-c(1)]
a <- subset(dem,MS2Metabolite %in% dup)
b <- subset(dem,!(MS2Metabolite %in% dup))
c <- data.frame()
for(i in 1:length(dup)){
  tmp <- subset(a,MS2Metabolite==dup[i])
  for(j in c(1:4)){
    c[i,j] <- paste0(tmp[,j], collapse = ';')
  }

  c[i,5:23] <- tmp[1,5:23]
  c[i,24] <- paste0(tmp[,24], collapse = ';')
  c[i,25:72] <- apply(tmp[,25:72],2,function(x){x <- as.numeric(x);return(mean(x))})
}
colnames(c) <- colnames(a)
dem <- rbind(b,c)

dem_int <- dem[,sp_info$sample.name]
dem_int$XCMS.id <- paste0(dem$XCMS.id,'_',dem$model)
dem_int$model <- dem$model
dem_int <- dem_int[,c(49,1:48,50)] # XCMS.id, sample name, model

non_dem_int <- non_dem[,sp_info$sample.name]
non_dem_int$XCMS.id <- paste0(non_dem$XCMS.id,'_',non_dem$model)
non_dem_int$model <- non_dem$model
non_dem_int <- non_dem_int[,c(49,1:48,50)]

#### 1. Compare the intensity difference between mutation carriers and non-mutation carriers ####
# Both intensity and rank based on intensity were considered in the comparison #
estimate_met_diff <- function(mut,met){
  res <- data.frame()
  sp <- unlist(strsplit(unique(mut$sample), split = ','))
  nor_sp <- colnames(met)[grep('NOR',colnames(met))]
  scz_sp <- colnames(met)[grep('SCZ',colnames(met))]
  if(grepl('NOR',mut$sample)){
    nor_var_sp <- sp[grep('NOR',sp)]
    nor_novar_sp <- nor_sp[!(nor_sp %in% nor_var_sp)]
  }else{
    nor_novar_sp <- nor_sp
  }
  scz_var_sp <- sp[grep('SCZ',sp)]
  scz_novar_sp <- scz_sp[!(scz_sp %in% scz_var_sp)]
  for(j in 1:dim(met)[1]){
    sta1 <- t.test(met[j,scz_var_sp], met[j,nor_novar_sp])
    sta2 <- t.test(met[j,scz_var_sp], met[j,scz_novar_sp])
    res[j,1:9] <- c(met[j,'XCMS.id'],sta1$p.value,sta2$p.value, 
                    sta1$estimate,sta2$estimate,sta1$stderr,sta2$stderr)
  }
  colnames(res) <- c('XCMS.id','scz_varvsnor_novar.p.value','scz_varvsscz_novar.p.value',
                     'mean of scz_var_sp met rank','mean of nor_novar_sp met rank','mean of scz_var_sp met rank','mean of scz_novar_sp met rank',
                     'scz_varvsnor_novar.stderr','scz_varvsscz_novar.stderr')
  return(res)
}

# Only consider the genes harboring variants in >= 2 SCZ samples.
# Two comparison were applied: Mutation carriers (SCZ) vs non-mutation carriers (SCZ); Mutation carriers (SCZ) vs healthy controls
##### 1.1 DEMs analysis #####
x <- vector()
for(i in 1:length(gs_var$Gene[duplicated(gs_var$Gene)])){
  gene <- gs_var$Gene[duplicated(gs_var$Gene)][i]
  sample <- gs_var[gs_var$Gene==gene,'sample']
  
  if((!grepl('NOR',paste0(sample,collapse = ';')))&(!duplicated(sample))){
    x <- c(x,gene) %>% unique()
  }
}
length(x)==length(unique(gs_var$Gene[duplicated(gs_var$Gene)]))）
y <- vector()
for(i in 1:length(gs_var[grepl(',',gs_var$sample),'Gene'])){
  tmp <- gs_var[grepl(',',gs_var$sample),][i,]
  if(!(grepl('NOR',tmp$sample))){
    y <- c(tmp$Gene,y) %>% unique
  }
}
multivar_onegene_id <- c(x,y) %>% unique()

# Rank every DEM based on intensity
dem_rank <- as.data.frame(t(apply(dem_int[,2:49],1,rank)))
dem_rank$XCMS.id <- dem_int$XCMS.id
dem_rank <- dem_rank[,c(49,1:48)]
met_var_relation <- list()
for(m in 1:length(multivar_onegene_id)){
  tmp <- subset(gs_var,Gene==multivar_onegene_id[m])
  tmp$sample[1] <- paste0(tmp$sample,collapse = ',')
  # Rank comparision
  lp_rank <- estimate_met_diff(tmp[1,],dem_rank)
  lp_rank <- lp_rank[,1:3]
  colnames(lp_rank)[2:3] <- c('rank.scz_varvsnor_novar.p.value','rank.scz_varvsscz_novar.p.value')
  lp_rank$'rank.scz_varvsnor_novar.padj' <- p.adjust(lp_rank$'rank.scz_varvsnor_novar.p.value', method = 'BH')
  lp_rank$'rank.scz_varvsscz_novar.padj' <- p.adjust(lp_rank$'rank.scz_varvsscz_novar.p.value', method = 'BH')
  # Normalized intensity comparison
  lp_real <- estimate_met_diff(tmp[1,],dem_int[,1:49])
  lp_real <- lp_real[,1:3]
  colnames(lp_real)[2:3] <- c('intensity.scz_varvsnor_novar.p.value',
                              'intensity.scz_varvsscz_novar.p.value')
  lp_real$'intensity.scz_varvsnor_novar.padj' <- p.adjust(lp_real$'intensity.scz_varvsnor_novar.p.value', method = 'BH')
  lp_real$'intensity.scz_varvsscz_novar.padj' <- p.adjust(lp_real$'intensity.scz_varvsscz_novar.p.value', method = 'BH')
  
  lp <- merge(lp_real,lp_rank, by.x = 'XCMS.id', by.y = 'XCMS.id')
  
  met_var_relation[[m]] <- lp
  names(met_var_relation)[m] <- multivar_onegene_id[m]

}

#####  1.2 Non-dem analysis #####
non_dem_rank <- as.data.frame(t(apply(non_dem_int[,2:49],1,rank)))
non_dem_rank$XCMS.id <- non_dem_int$XCMS.id
non_dem_rank <- non_dem_rank[,c(49,1:48)]
nondem_met_var_relation <- list()
for(n in 1:length(multivar_onegene_id)){
  tmp <- subset(gs_var,Gene==multivar_onegene_id[n])
  tmp$sample[1] <- paste0(tmp$sample,collapse = ',')
  lp_rank <- estimate_met_diff(tmp[1,],non_dem_rank)
  lp_rank <- lp_rank[,1:3]
  colnames(lp_rank)[2:3] <- c('rank.scz_varvsnor_novar.p.value','rank.scz_varvsscz_novar.p.value')
  lp_rank$'rank.scz_varvsnor_novar.padj' <- p.adjust(lp_rank$'rank.scz_varvsnor_novar.p.value', method = 'BH')
  lp_rank$'rank.scz_varvsscz_novar.padj' <- p.adjust(lp_rank$'rank.scz_varvsscz_novar.p.value', method = 'BH')
  lp_real <- estimate_met_diff(tmp[1,],non_dem_int[,1:49])
  lp_real <- lp_real[,1:3]
  colnames(lp_real)[2:3] <- c('intensity.scz_varvsnor_novar.p.value',
                              'intensity.scz_varvsscz_novar.p.value')
  lp_real$'intensity.scz_varvsnor_novar.padj' <- p.adjust(lp_real$'intensity.scz_varvsnor_novar.p.value', method = 'BH')
  lp_real$'intensity.scz_varvsscz_novar.padj' <- p.adjust(lp_real$'intensity.scz_varvsscz_novar.p.value', method = 'BH')
  
  lp <- merge(lp_real,lp_rank, by.x = 'XCMS.id', by.y = 'XCMS.id')
  
  nondem_met_var_relation[[n]] <- lp
  names(nondem_met_var_relation)[n] <- multivar_onegene_id[n]
}

##### 1.3 Count the number of variant-associated DEMs and non-DEMs #####
for(i in 1:length(met_var_relation)){
  dem_a <- subset(met_var_relation[[i]],(intensity.scz_varvsscz_novar.padj<0.05)) %>% dim()
  nondem_a <- subset(nondem_met_var_relation[[i]],(intensity.scz_varvsscz_novar.padj<0.05)) %>% dim()
  
  dem_b <- subset(met_var_relation[[i]],(intensity.scz_varvsnor_novar.padj<0.05)) %>% dim()
  nondem_b <- subset(nondem_met_var_relation[[i]],(intensity.scz_varvsnor_novar.padj<0.05)) %>% dim()
  
  dem_c <- subset(met_var_relation[[i]],(rank.scz_varvsscz_novar.padj<0.05)) %>% dim()
  nondem_c <- subset(nondem_met_var_relation[[i]],(rank.scz_varvsscz_novar.padj<0.05)) %>% dim()
  
  dem_d <- subset(met_var_relation[[i]],(rank.scz_varvsnor_novar.padj<0.05)) %>% dim()
  nondem_d <- subset(nondem_met_var_relation[[i]],(rank.scz_varvsnor_novar.padj<0.05)) %>% dim()
  
  dem_fin <- met_var_relation[[i]] %>% 
    dplyr::filter(intensity.scz_varvsscz_novar.padj<0.05) %>%
    dplyr::filter(intensity.scz_varvsnor_novar.padj<0.05) %>% 
    dplyr::filter(rank.scz_varvsscz_novar.padj<0.05) %>% 
    dplyr::filter(rank.scz_varvsnor_novar.padj<0.05)
  nondem_fin <- nondem_met_var_relation[[i]] %>% 
    dplyr::filter(intensity.scz_varvsscz_novar.padj<0.05) %>%
    dplyr::filter(intensity.scz_varvsnor_novar.padj<0.05) %>% 
    dplyr::filter(rank.scz_varvsscz_novar.padj<0.05) %>% 
    dplyr::filter(rank.scz_varvsnor_novar.padj<0.05)
  
  print(paste0(multivar_onegene_id[[i]],
               "  intensity.mut-SCZ vs nonmut-SCZ  ",
               'Sig DEMs: ',dem_a[1],' ; Sig NonDEMs: ',nondem_a[1],
               "  intensity.mut-SCZ vs NOR  ",
               'Sig DEMs: ',dem_b[1],' ; Sig NonDEMs: ',nondem_b[1],
               "  rank.mut-SCZ vs nonmut-SCZ  ",
               'Sig DEMs: ',dem_c[1],' ; Sig NonDEMs: ',nondem_c[1],
               "  rank.mut-SCZ vs NOR  ",
               'Sig DEMs: ',dem_d[1],' ; Sig NonDEMs: ',nondem_d[1],
              "  shared  ",
              'Sig DEMs: ', nrow(dem_fin),' ; Sig NonDEMs: ',nrow(nondem_fin))
  
}

#### 2. Save and output ####
save(met_var_relation, nondem_met_var_relation, file = './Association-analysis/variant-metabolite.association.RData')
