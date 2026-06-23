===============================================================================
# Complete Metabolomics Pipeline
# Phase 1: XCMS Peak Picking -> Phase 2: Masscleaner QC -> Phase 3: Statistics
# Supports: "positive" or "negative" modes via command line argument. 
# Usage: Rscript pipeline.R [pos|neg|both]
#===============================================================================

args <- commandArgs(trailingOnly = TRUE)
RUN_MODE <- ifelse(length(args) > 0, args[1], "both")
if (!RUN_MODE %in% c("positive","negative")) {
  stop("Usage: Rscript pipeline.R [positive|negative]")
}

#===============================================================================
library(S4Vectors)
library(xcms)
library(MSnbase)
library(magrittr)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)
library(masscleaner)
library(massdataset)
library(massstat)
library(massqc)
library(tidyverse)
library(scales)
library(corrplot)
library(ropls)
library(ggsci)
library(car)
library(pls)
library(VennDiagram)

#-------------------------------------------------------------------------------
# 0. Configuration
#-------------------------------------------------------------------------------
SEQDATA_DIR     <- file.path("~/EOSCZ/MetaSeq-Data", "convert_data") # The input files is the .mzML converted from .raw file format using ProteoWizard
RESULT_DIR      <- "~/EOSCZ/MetaSeq-Result/"
XCMS_OUT        <- file.path(RESULT_DIR, "XCMS")
MASS_DIR        <- file.path(RESULT_DIR, "tidymass")
STAT_DIR        <- file.path(RESULT_DIR, "statistics")

#-------------------------------------------------------------------------------
# 1. Ion-Mode-Specific Functions (positive or Negative)
#-------------------------------------------------------------------------------

#' Phase 1: XCMS workflow for ONE ion mode
run_xcms_for_mode <- function(ion_mode, file_list, cent_params, merge_params) {
  cat("\n========== XCMS:", toupper(ion_mode), "Mode ==========\n")

  raw <- readMSData(
    files = file_list$file_name,
    pdata = new("NAnnotatedDataFrame", file_list[, c("sample_name", "group")]),
    mode = "onDisk"
  )
  cat("  Samples:", length(file_list$file_name), 
      "| MS1:", sum(msLevel(raw) == 1), 
      "| MS2:", sum(msLevel(raw) == 2), "\n")

  # TIC
  tic <- chromatogram(raw, aggregationFun = "sum")
  pdf(file.path(XCMS_OUT, paste0(ion_mode, ".TIC.pdf")), width = 10, height = 6)
  plot(tic, col = GROUP_COLORS[raw$group], main = paste(toupper(ion_mode), "TIC"))
  dev.off()

  # BPC correlation heatmap
  bpc <- chromatogram(raw, aggregationFun = "max")
  bpis_bin <- bin(bpc, binSize = 1.18)
  cormat <- do.call(cbind, lapply(bpis_bin, intensity)) %>% log2() %>% cor()
  colnames(cormat) <- rownames(cormat) <- raw$sample_name
  ann <- data.frame(group = raw$group)
  rownames(ann) <- raw$sample_name

  pdf(file.path(XCMS_OUT, paste0(ion_mode, ".correlation_heatmap.pdf")), width = 8, height = 7)
  pheatmap(cormat, annotation = ann, annotation_color = list(group = GROUP_COLORS),
           main = paste(toupper(ion_mode), "Sample Correlation"))
  dev.off()

  # Peak detection
  cat("  Peak detection...\n")
  peaks <- findChromPeaks(raw, param = cent_params)
  cat("  Peaks detected:", nrow(chromPeaks(peaks)), "\n")

  # Refinement
  cat("  Refining peaks...\n")
  refined <- refineChromPeaks(peaks, param = merge_params)
  cat("  After refinement:", nrow(chromPeaks(refined)), "\n")

  # Peak summary
  summary_fun <- function(z) c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
  peak_summary <- do.call(rbind, lapply(
    split.data.frame(chromPeaks(refined), f = chromPeaks(refined)[, "sample"]),
    FUN = summary_fun
  ))
  rownames(peak_summary) <- basename(fileNames(refined))
  write.csv(peak_summary, file.path(XCMS_OUT, paste0(ion_mode, ".peak_summary.csv")))

  # RT alignment
  cat("  RT alignment...\n")
  aligned <- adjustRtime(refined, param = ObiwarpParam(binSize = 0.1))

  adj_bpc <- chromatogram(aligned, aggregationFun = "max", include = "none")
  pdf(file.path(XCMS_OUT, paste0(ion_mode, ".RT_alignment.pdf")), width = 10, height = 8)
  par(mfrow = c(2, 1), mar = c(4.5, 4.2, 1, 0.5))
  plot(adj_bpc, col = GROUP_COLORS[aligned$group], main = paste(toupper(ion_mode), "Aligned BPC"))
  plotAdjustedRtime(aligned, col = GROUP_COLORS[aligned$group])
  dev.off()

  # Correspondence
  cat("  Feature grouping...\n")
  pdp <- PeakDensityParam(sampleGroups = aligned$group, minFraction = 0.5, binSize = 0.015, bw = 5)
  grouped <- groupChromPeaks(aligned, param = pdp, msLevel = 1L)

  cat("  Gap filling...\n")
  filled <- fillChromPeaks(grouped, param = ChromPeakAreaParam())

  # Export
  res <- quantify(filled, value = "into", method = "sum")
  peak_table <- featureValues(filled) %>% as.data.frame()
  sample_map <- setNames(colData(res)$sample_name, colnames(peak_table))
  colnames(peak_table) <- sample_map[colnames(peak_table)]

  feat_info <- as.data.frame(rowData(res))
  peak_table <- cbind(feat_info, peak_table)
  keep_cols <- c("mz", "rt", setdiff(colnames(peak_table), colnames(feat_info)))
  peak_simp <- peak_table[, keep_cols]

  write.csv(peak_table, file.path(XCMS_OUT, paste0(ion_mode, ".Peak_table.csv")), quote = FALSE)
  write.csv(peak_simp, file.path(XCMS_OUT, paste0(ion_mode, ".Peak_table.simplify.csv")), quote = FALSE)

  cat("  Output:", nrow(peak_table), "features x", ncol(peak_table) - ncol(feat_info), "samples\n")
  invisible(list(raw = raw, peaks = peaks, refined = refined, aligned = aligned, 
                 grouped = filled, peak_table = peak_table, peak_simp = peak_simp))
}

#' Phase 2: Masscleaner QC for ONE ion mode
run_masscleaner_for_mode <- function(ion_mode, spinfo, rsd_threshold = RSD_THRESHOLD) {
  cat("\n========== Masscleaner QC:", toupper(ion_mode), "Mode ==========\n")
  xcms_file <- file.path(XCMS_OUT, paste0(ion_mode, ".Peak_table.simplify.csv"))
  peak <- read.csv(xcms_file, stringsAsFactors = FALSE, row.names = 1)
  var_info <- data.frame(
  variable_id = rownames(peak),
  mz = peak$mz,
  rt = peak$rt,
  stringsAsFactors = FALSE
  exp_data <- peak[, !(colnames(peak) %in% c("mz", "rt"))]

  rownames(exp_data) <- var_info$variable_id

  # Map sample names
  if (any(colnames(exp_data) %in% spinfo$result_name)) {
    name_map <- setNames(spinfo$sample.name, spinfo$result_name)
    colnames(exp_data) <- name_map[colnames(exp_data)]
  }

  # Align
  common_samples <- intersect(colnames(exp_data), spinfo$sample_id)
  if (length(common_samples) == 0) {
    common_samples <- intersect(colnames(exp_data), spinfo$sample.name)
    if (length(common_samples) == 0) stop("No matching samples for ", ion_mode)
    spinfo <- spinfo[spinfo$sample.name %in% common_samples, ]
    spinfo$sample_id <- spinfo$sample.name
  } else {
    spinfo <- spinfo[spinfo$sample_id %in% common_samples, ]
  }
  exp_data <- exp_data[, common_samples, drop = FALSE]

  # RT unit check
  var_info$rt <- as.numeric(var_info$rt)
  if (max(var_info$rt, na.rm = TRUE) < 1000) var_info$rt <- var_info$rt * 60

  # Build mass_dataset
  object <- create_mass_dataset(
    expression_data = exp_data,
    sample_info = spinfo,
    variable_info = var_info[, c("variable_id", "mz", "rt")]
  )
  cat("  Features:", nrow(var_info), "| Samples:", ncol(exp_data), "\n")

  # MV plot
  mv_sp <- data.frame(
    sample = names(get_mv_number(object, by = "sample")),
    feature_number = get_mv_number(object, by = "sample"),
    stringsAsFactors = FALSE
  )
  mv_sp$group <- case_when(
    grepl("QC", mv_sp$sample) ~ "QC",
    grepl("SCZ", mv_sp$sample) ~ "SCZ",
    grepl("NOR", mv_sp$sample) ~ "NOR"
  )
  mv_sp$group <- factor(mv_sp$group, levels = c("QC", "NOR", "SCZ"))

  p_mv <- ggplot(mv_sp, aes(x = reorder(sample, feature_number), y = feature_number, fill = group)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_manual(values = c("#63AA83", "#33CCD0", "#FBABA5")) +
    theme_classic() + labs(title = paste(toupper(ion_mode), "Missing Values"))
  ggsave(file.path(MASS_DIR, paste0('massclean.', ion_mode, ".mv_plot.pdf")), p_mv, width = 8, height = 6)

  # Filter MV by group
  qc_id      <- spinfo$sample_id[spinfo$group == "QC"]
  control_id <- spinfo$sample_id[spinfo$group == "control"]
  case_id    <- spinfo$sample_id[spinfo$group == "case"]

  object <- object %>%
    mutate_variable_na_freq(according_to_samples = qc_id) %>%
    mutate_variable_na_freq(according_to_samples = control_id) %>%
    mutate_variable_na_freq(according_to_samples = case_id) %>%
    activate_mass_dataset(what = "variable_info") %>%
    filter(na_freq <= 0.2, na_freq.1 <= 0.5, na_freq.2 <= 0.5)
  cat("  After MV filter:", nrow(extract_variable_info(object)), "features\n")

  # Outliers
  outlier_samples <- object %>% `+`(1) %>% log(2) %>% scale() %>% detect_outlier()
  outlier_table <- extract_outlier_table(outlier_samples)
  if (nrow(outlier_table) > 0) {
    cat("  Outliers:\n"); print(outlier_table)
  }

  # KNN imputation
  object <- impute_mv(object, method = "knn")
  cat("  MV after KNN:", get_mv_number(object), "\n")

  # Normalization
  object <- object %>% activate_mass_dataset(what = "sample_info") %>% mutate(batch = as.character(batch))
  object_norm <- normalize_data(object, method = "pqn")

  # RSD filter
  object_norm <- object_norm %>% mutate_rsd(according_to_samples = qc_id)
  object_norm <- object_norm %>% activate_mass_dataset(what = "variable_info") %>% filter(rsd <= rsd_threshold)
  cat("  After RSD filter:", nrow(extract_variable_info(object_norm)), "features\n")

  # Save
  save(object_norm, file = file.path(MASS_DIR, paste0('massclean.', ion_mode, ".normalized.Rdata")))
  write.csv(extract_expression_data(object_norm),
            file = file.path(MASS_DIR, paste0('massclean.', ion_mode, ".normalized-intensity.csv")),
            row.names = TRUE, quote = FALSE)

  invisible(object_norm)
}

#' Phase 3: Statistics for ONE ion mode
run_statistics_for_mode <- function(ion_mode, obj, annt_file, spinfo_stat) {
  cat("\n========== Statistics:", toupper(ion_mode), "Mode ==========\n")

  sp_scz_all <- subset(spinfo_stat, group == "case")
  sp_nor     <- subset(spinfo_stat, group == "control")
  sp_scz_sub <- subset(sp_scz_all, age <= 25)

  # Univariate t-test
  run_univariate_t <- function(obj, sample_info) {
    scz <- subset(sample_info, group == "case")$sample_id
    nor <- subset(sample_info, group == "control")$sample_id
    obj_stat <- obj %>%
      mutate_fc(control_sample_id = nor, case_sample_id = scz, mean_median = "mean") %>%
      mutate_p_value(control_sample_id = nor, case_sample_id = scz, method = "t.test", p_adjust_methods = "BH")
    extract_variable_info(obj_stat) %>% select(variable_id, fc, p_value, p_value_adjust)
  }

  t_all <- run_univariate_t(obj, rbind(sp_scz_all, sp_nor))
  t_sub <- run_univariate_t(obj, rbind(sp_scz_sub, sp_nor))

  # LM all ages
  all_age_lm <- function(obj, sample_info) {
    cat("  LM all ages...\n")
    expr_mat <- obj@expression_data
    lm_res <- lapply(rownames(expr_mat), function(met_id) {
      y <- log2(as.numeric(expr_mat[met_id, sample_info$sample_id]))
      df <- data.frame(
        y = y,
        group = factor(sample_info$group, levels = c("control", "case")),
        age = as.numeric(sample_info$age),
        gender = sample_info$gender,
        disease.duration = as.numeric(sample_info$disease.duration),
        drug.dose = as.numeric(sample_info$drug.dose)
      )
      fit_main <- lm(y ~ group + age + gender, data = df)
      coef_main <- summary(fit_main)$coefficients
      df_case <- subset(df, group == "case")
      fit_case <- lm(y ~ age + gender + disease.duration + drug.dose, data = df_case)
      coef_case <- summary(fit_case)$coefficients
      fit_reduced <- lm(y ~ age + gender, data = df_case)
      anova_res <- anova(fit_reduced, fit_case)
      data.frame(
        variable_id = met_id,
        group_beta = coef_main["groupcase", "Estimate"],
        group_p = coef_main["groupcase", "Pr(>|t|)"],
        duration_beta = coef_case["disease.duration", "Estimate"],
        duration_p = coef_case["disease.duration", "Pr(>|t|)"],
        dose_beta = coef_case["drug.dose", "Estimate"],
        dose_p = coef_case["drug.dose", "Pr(>|t|)"],
        duration_dose_joint_p = anova_res$`Pr(>F)`[2],
        r_squared_main = summary(fit_main)$r.squared,
        r_squared_case = summary(fit_case)$r.squared,
        stringsAsFactors = FALSE
      )
    })
    res <- do.call(rbind, lm_res)
    res$group_p_adj <- p.adjust(res$group_p, method = "BH")
    res$duration_p_adj <- p.adjust(res$duration_p, method = "BH")
    res$dose_p_adj <- p.adjust(res$dose_p, method = "BH")
    res$duration_dose_joint_p_adj <- p.adjust(res$duration_dose_joint_p, method = "BH")
    res
  }

  lm_all <- all_age_lm(obj, spinfo_stat)

  # LM subgroup
  sub_age_lm <- function(obj, sample_info) {
    cat("  LM subgroup (<=25)...\n")
    sp_sub <- subset(sample_info, age <= 25)
    expr_mat <- obj@expression_data[, sp_sub$sample_id, drop = FALSE]
    lm_res <- lapply(rownames(expr_mat), function(met_id) {
      y <- log2(as.numeric(expr_mat[met_id, sp_sub$sample_id]))
      df <- data.frame(
        y = y,
        group = factor(sp_sub$group, levels = c("control", "case")),
        age = sp_sub$age,
        gender = sp_sub$gender,
        disease.duration = sp_sub$disease.duration,
        drug.dose = sp_sub$drug.dose
      )
      fit_main <- lm(y ~ group + age + gender, data = df)
      coef_main <- summary(fit_main)$coefficients
      df_case <- subset(df, group == "case")
      fit_case <- lm(y ~ age + gender + disease.duration + drug.dose, data = df_case)
      coef_case <- summary(fit_case)$coefficients
      data.frame(
        variable_id = met_id,
        subgroup_beta = coef_main["groupcase", "Estimate"],
        subgroup_p = coef_main["groupcase", "Pr(>|t|)"],
        subgroup_duration_beta = coef_case["disease.duration", "Estimate"],
        subgroup_duration_p = coef_case["disease.duration", "Pr(>|t|)"],
        subgroup_dose_beta = coef_case["drug.dose", "Estimate"],
        subgroup_dose_p = coef_case["drug.dose", "Pr(>|t|)"],
        subgroup_r_squared_main = summary(fit_main)$r.squared,
        subgroup_r_squared_case = summary(fit_case)$r.squared,
        stringsAsFactors = FALSE
      )
    })
    res <- do.call(rbind, lm_res)
    res$subgroup_p_adj <- p.adjust(res$subgroup_p, method = "BH")
    res$subgroup_duration_p_adj <- p.adjust(res$subgroup_duration_p, method = "BH")
    res$subgroup_dose_p_adj <- p.adjust(res$subgroup_dose_p, method = "BH")
    res
  }

  lm_sub <- sub_age_lm(obj, spinfo_stat)

  # PLS-DA
  cat("  PLS-DA...\n")
  dataMatrix <- t(obj@expression_data[, spinfo_stat$sample_id])
  dataMatrix <- dataMatrix[order(rownames(dataMatrix)), , drop = FALSE]
  sampleMetadata <- spinfo_stat[order(spinfo_stat$sample_id), ]
  stopifnot(all(rownames(dataMatrix) == sampleMetadata$sample_id))
  plsda <- opls(dataMatrix, sampleMetadata$group, permI = 1000)
  vip <- getVipVn(plsda)
  vip_df <- data.frame(XCMS.id = names(vip), VIP_group = vip, stringsAsFactors = FALSE)

  # Merge results
  result <- extract_variable_info(obj, with_expression_data = FALSE)[, 1:3]
  colnames(result)[1] <- "XCMS.id"

  result <- result %>%
    left_join(t_all, by = c("XCMS.id" = "variable_id")) %>%
    rename(SCZ_allvsNOR.FC = fc, SCZ_allvsNOR.ttest.p = p_value, SCZ_allvsNOR.ttest.padj = p_value_adjust) %>%
    left_join(t_sub, by = c("XCMS.id" = "variable_id")) %>%
    rename(SCZ_AvsNOR.FC = fc, SCZ_AvsNOR.ttest.p = p_value, SCZ_AvsNOR.ttest.padj = p_value_adjust)

  result <- result %>%
    left_join(lm_all %>% select(variable_id, group_beta, group_p, group_p_adj,
                                duration_beta, duration_p, duration_p_adj,
                                dose_beta, dose_p, dose_p_adj),
              by = c("XCMS.id" = "variable_id")) %>%
    rename_with(~ paste0("SCZ_allvsNOR.", .), .cols = c(group_beta, group_p, group_p_adj,
                                                         duration_beta, duration_p, duration_p_adj,
                                                         dose_beta, dose_p, dose_p_adj)) %>%
    left_join(lm_sub %>% select(variable_id, subgroup_beta, subgroup_p, subgroup_p_adj,
                                subgroup_duration_beta, subgroup_duration_p, subgroup_duration_p_adj,
                                subgroup_dose_beta, subgroup_dose_p, subgroup_dose_p_adj),
              by = c("XCMS.id" = "variable_id")) %>%
    rename_with(~ paste0("SCZ_AvsNOR.", .), .cols = starts_with("subgroup_"))

  result <- result %>% left_join(vip_df, by = "XCMS.id") %>%
    mutate(
      is.sig.ttest.SCZ_allvsNOR = ifelse(SCZ_allvsNOR.ttest.padj < 0.05, 1, 0),
      is.sig.ttest.SCZ_AvsNOR = ifelse(SCZ_AvsNOR.ttest.padj < 0.05, 1, 0),
      is.sig.lm.SCZ_allvsNOR = ifelse(SCZ_allvsNOR.group_p_adj < 0.05, 1, 0),
      is.sig.lm.SCZ_AvsNOR = ifelse(SCZ_AvsNOR.subgroup_p_adj < 0.05, 1, 0)
    )

  # Add annotation
  if (file.exists(annt_file)) {
    annt <- read.csv(annt_file, stringsAsFactors = FALSE)
    colnames(annt)[2] <- "XCMS.id"
    annt_cols <- intersect(c("XCMS.id", "MS2Metabolite", "MS2superclass", "NumberMS1hmdb", "NumberMS1kegg",
                             "MS1hmdbName", "MS1keggName", "MS2hmd", "MS2kegg",
                             "MS1hmdbID", "MS1hmdbTokegg", "MS1keggID"),
                           colnames(annt))
    result <- result %>% left_join(annt[, annt_cols], by = "XCMS.id") %>% mutate(model = ion_mode)
  } else {
    warning("Annotation file not found: ", annt_file)
    result$model <- ion_mode
  }

  exp <- extract_expression_data(obj)
  exp$XCMS.id <- rownames(exp)
  result <- merge(result, exp, by = "XCMS.id")

  write.table(result, file = file.path(MASS_DIR, paste0(ion_mode, ".results.txt")),
              row.names = FALSE, sep = "\t", quote = FALSE)

  # Venn
  venn.diagram(
    list(
      t_all = result$XCMS.id[result$is.sig.ttest.SCZ_allvsNOR == 1],
      t_sub = result$XCMS.id[result$is.sig.ttest.SCZ_AvsNOR == 1],
      lm_all = result$XCMS.id[result$is.sig.lm.SCZ_allvsNOR == 1],
      lm_sub = result$XCMS.id[result$is.sig.lm.SCZ_AvsNOR == 1]
    ),
    filename = file.path(STAT_DIR, paste0('LM.', ion_mode, ".venn.tiff")),
    imagetype = "tiff"
  )

  invisible(result)
}

#-------------------------------------------------------------------------------
# 4. Main Execution
#-------------------------------------------------------------------------------

cat("\n========== Loading Data ==========\n")

file <- read.table(SAMPLE_LIST, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file$file_name <- file.path(SEQDATA_DIR, file$file_name)
file.input <- file[file$sample_name != "Not_used", ]

spinfo <- read.csv(SPINFO_FILE, stringsAsFactors = FALSE)
colnames(spinfo)[1:2] <- c("sample_id", "subject_id")
spinfo$gender <- factor(spinfo$gender, levels = c("M", "F"))

sp_scz <- subset(spinfo, group == "case")
cat("Correlation (duration vs dose):", cor(sp_scz$disease.duration, sp_scz$drug.dose, use = "complete.obs"), "\n")
cat("Correlation (age vs dose):", cor(sp_scz$age, sp_scz$drug.dose, use = "complete.obs"), "\n")

spinfo_stat <- subset(spinfo, class != "QC" & !(sample.name %in% c("SCZ_04", "SCZ_08", "SCZ_23")))

# Determine which modes to run
modes_to_run <- if (RUN_MODE == "both") c("positive", "negative") else RUN_MODE

# Phase 1 & 2 parameters
cent_params  <- CentWaveParam(snthresh = 6, ppm = 45, peakwidth = c(5, 25), mzdiff = 0.01)
merge_params <- MergeNeighboringPeaksParam(ppm = 45)

# Run for each mode
all_results <- list()

for (mode in modes_to_run) {

  # ---- Phase 1: XCMS ----
  cat("\n>>>>>>>>>> Processing", toupper(mode), "Mode <<<<<<<<<<\n")
  register(bpstart(MulticoreParam(N_CORES)))

  mode_files <- file.input[file.input$type == ifelse(mode == "positive", "Pos", "Neg"), ]
  xcms_res <- run_xcms_for_mode(mode, mode_files, cent_params, merge_params)

  register(bpstop(MulticoreParam(N_CORES)))

  # ---- Phase 2: Masscleaner ----
  qc_obj <- run_masscleaner_for_mode(mode, spinfo)

  # ---- Phase 3: Statistics ----
  annt_file <- file.path(SEQDATA_DIR, paste0(mode, "-all-identification.csv"))
  stat_res <- run_statistics_for_mode(mode, qc_obj, annt_file, spinfo_stat)

  all_results[[mode]] <- stat_res
}

# ---- Cross-mode annotation check (only if both modes run) ----
if (RUN_MODE == "both" && length(all_results) == 2) {
  cat("\n========== Cross-Mode Annotation Check ==========\n")

  mg_res <- bind_rows(all_results[["pos"]], all_results[["neg"]])

  ms1_only <- mg_res %>% filter((NumberMS1hmdb > 0 | NumberMS1kegg > 0) & MS2Metabolite == "-")
  ms2_only <- mg_res %>% filter(MS2Metabolite != "-" & NumberMS1hmdb == 0 & NumberMS1kegg == 0)
  ms1_ms2  <- mg_res %>% filter(MS2Metabolite != "-" & (NumberMS1hmdb > 0 | NumberMS1kegg > 0))

  cat("MS1 only:", nrow(ms1_only), "\n")
  cat("MS2 only:", nrow(ms2_only), "\n")
  cat("MS1 + MS2:", nrow(ms1_ms2), "\n")

  ms1_ms2_match <- ms1_ms2 %>%
    filter(
      grepl(MS2Metabolite, MS1hmdbName, fixed = TRUE) |
      grepl(MS2Metabolite, MS1keggName, fixed = TRUE) |
      (MS2hmd != "-" & !is.na(MS2hmd) & grepl(MS2hmd, MS1hmdbID, fixed = TRUE)) |
      (MS2kegg != "-" & !is.na(MS2kegg) & grepl(MS2kegg, MS1keggID, fixed = TRUE))
    )

  cat("MS1-MS2 matched:", nrow(ms1_ms2_match), "\n")

  high_conf <- bind_rows(ms1_ms2_match, ms2_only)
  write.table(high_conf, file = file.path(STAT_DIR, "LM.high_confidence_annotations.txt"),
              row.names = FALSE, sep = "\t", quote = FALSE)
}

cat("\n========== Pipeline Complete ==========\n")
cat("Run mode:", RUN_MODE, "\n")
cat("Outputs:\n")
cat("  XCMS:", XCMS_OUT, "\n")
cat("  QC:", MASS_DIR, "\n")
cat("  Stats:", STAT_DIR, "\n")
