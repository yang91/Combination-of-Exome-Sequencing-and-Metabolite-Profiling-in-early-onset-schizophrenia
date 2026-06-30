# =============================================================================
# utils.R — Shared Utility Functions for EOSCZ Analysis Pipeline
# =============================================================================
# Source this file at the top of every analysis script (after config.R):
#   source("config.R")
#   source("utils.R")
#
# These functions eliminate code duplication across:
#   3.2.Metabolite_symptom_corr.R
#   3.3.Variant_metabolite_assoc.R
#   3.4.Metabolite_metabolite_correlation.R
#   3.5.Protein_metabolite_network.R
# =============================================================================

library(magrittr)
library(dplyr)

#-------------------------------------------------------------------------------
# 1. Metabolite Data Loading & Preprocessing
#-------------------------------------------------------------------------------

#' Load and combine positive and negative mode metabolite data
#'
#' @param pos_file Path to positive mode statistics file (default: POS_STAT_FILE from config)
#' @param neg_file Path to negative mode statistics file (default: NEG_STAT_FILE from config)
#' @param select_cols Column names or indices to retain. If NULL, returns all columns.
#' @return Combined data frame with both ion modes.
#' @examples
#' all_metabolites <- load_metabolite_data()
#' all_metabolites <- load_metabolite_data(select_cols = c("XCMS.id", "model", "MS2Metabolite", "All_sample.LM.padj"))
load_metabolite_data <- function(pos_file = POS_STAT_FILE,
                                  neg_file = NEG_STAT_FILE,
                                  select_cols = NULL) {
  if (!file.exists(pos_file)) {
    warning("Positive mode file not found: ", pos_file)
    pos_metabolites <- data.frame()
  } else {
    pos_metabolites <- read.csv(pos_file, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
  }

  if (!file.exists(neg_file)) {
    warning("Negative mode file not found: ", neg_file)
    neg_metabolites <- data.frame()
  } else {
    neg_metabolites <- read.csv(neg_file, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
  }

  all_metabolites <- rbind(pos_metabolites, neg_metabolites)

  if (!is.null(select_cols)) {
    available_cols <- intersect(select_cols, colnames(all_metabolites))
    missing_cols <- setdiff(select_cols, colnames(all_metabolites))
    if (length(missing_cols) > 0) {
      warning("Missing columns in metabolite data: ", paste(missing_cols, collapse = ", "))
    }
    all_metabolites <- all_metabolites[, available_cols, drop = FALSE]
  }

  return(all_metabolites)
}


#' Extract sample intensity column names from metabolite data frame
#'
#' Automatically identifies sample columns by excluding known metadata columns.
#' This replaces hardcoded column indices (e.g., 39:106, 25:72).
#'
#' @param metabolite_df Data frame containing metabolite annotations and intensities
#' @param meta_patterns Regex patterns for metadata column names (customizable)
#' @return Character vector of sample column names.
#' @examples
#' sample_names <- extract_sample_names(dem_metabolites)
extract_sample_names <- function(metabolite_df,
                                 meta_patterns = c("^XCMS\\.id$", "^model$", "^mz$", "^rt$",
                                                   "^MS2", "^MS1", "^annotation$", "^Number",
                                                   "^p_value", "^padj", "^FC", "^fc$",
                                                   "^group_", "^duration_", "^dose_")) {
  is_meta <- sapply(colnames(metabolite_df), function(col) {
    any(sapply(meta_patterns, function(p) grepl(p, col, ignore.case = TRUE)))
  })
  sample_names <- colnames(metabolite_df)[!is_meta]
  return(sample_names)
}


#' Filter differentially expressed metabolites (DEMs)
#'
#' @param metabolite_df Combined metabolite data frame
#' @param padj_col Column name for adjusted p-value (default: "All_sample.LM.padj")
#' @param threshold Significance threshold (default: PADJ_THRESHOLD from config)
#' @return Data frame of DEMs.
filter_dem <- function(metabolite_df,
                       padj_col = "All_sample.LM.padj",
                       threshold = PADJ_THRESHOLD) {
  if (!padj_col %in% colnames(metabolite_df)) {
    stop("p-value column not found: ", padj_col,
         "\nAvailable columns: ", paste(head(colnames(metabolite_df), 10), collapse = ", "), "...")
  }
  metabolite_df %>%
    dplyr::filter(!!rlang::sym(padj_col) < threshold)
}


#' Merge metabolites with duplicate MS2 annotations
#'
#' For metabolites sharing the same MS2 annotation, concatenates metadata
#' columns and averages sample intensity columns.
#'
#' @param metabolite_df Data frame with metabolite annotations and intensities.
#' @param id_var Name of the annotation column to check for duplicates (default: "MS2Metabolite")
#' @param concat_cols Column names/indices to concatenate for duplicates (default: first 4 columns + annotation column)
#' @param first_cols Column names/indices to take the first value (default: all remaining metadata)
#' @return Data frame with duplicate annotations merged.
#'
#' @note Replaces the fragile hardcoded column indices (25:72, etc.) in the
#'   original scripts with dynamic detection.
merge_duplicate_metabolites <- function(metabolite_df,
                                         id_var = "MS2Metabolite",
                                         concat_cols = NULL,
                                         first_cols = NULL) {
  if (!id_var %in% colnames(metabolite_df)) {
    warning("Annotation column '", id_var, "' not found. Returning data unchanged.")
    return(metabolite_df)
  }

  # Find duplicate annotation names (excluding NA, empty, and "-")
  dup_names <- metabolite_df[[id_var]][duplicated(metabolite_df[[id_var]])]
  dup_names <- unique(dup_names)
  dup_names <- dup_names[!is.na(dup_names) & dup_names != "-" & dup_names != ""]

  if (length(dup_names) == 0) {
    return(metabolite_df)
  }

  dup_rows <- metabolite_df %>% dplyr::filter(!!rlang::sym(id_var) %in% dup_names)
  unique_rows <- metabolite_df %>% dplyr::filter(!(!!rlang::sym(id_var) %in% dup_names))

  # Auto-detect column roles if not specified
  if (is.null(concat_cols)) {
    # Default: concatenate first 4 columns and the id_var column
    concat_cols <- c(seq_len(min(4, ncol(metabolite_df))), which(colnames(metabolite_df) == id_var))
    concat_cols <- unique(concat_cols)
  }
  if (is.null(first_cols)) {
    # All other non-intensity columns
    sample_cols <- extract_sample_names(metabolite_df)
    first_cols <- setdiff(seq_len(ncol(metabolite_df)), c(concat_cols, which(colnames(metabolite_df) %in% sample_cols)))
  }

  # Intensity columns are all remaining columns
  sample_cols <- extract_sample_names(metabolite_df)
  intensity_cols <- which(colnames(metabolite_df) %in% sample_cols)

  merged_list <- vector("list", length(dup_names))

  for (i in seq_along(dup_names)) {
    met_name <- dup_names[i]
    met_group <- dplyr::filter(dup_rows, !!rlang::sym(id_var) == met_name)

    new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(metabolite_df)))
    colnames(new_row) <- colnames(metabolite_df)

    # Concatenate specified columns
    for (j in intersect(concat_cols, seq_len(ncol(metabolite_df)))) {
      vals <- met_group[, j]
      vals <- vals[!is.na(vals)]
      new_row[1, j] <- paste(unique(vals), collapse = ";")
    }

    # Take first value for specified columns
    for (j in intersect(first_cols, seq_len(ncol(metabolite_df)))) {
      new_row[1, j] <- met_group[1, j]
    }

    # Average intensity columns
    if (length(intensity_cols) > 0) {
      new_row[1, intensity_cols] <- apply(met_group[, intensity_cols, drop = FALSE], 2, function(x) {
        mean(as.numeric(x), na.rm = TRUE)
      })
    }

    merged_list[[i]] <- new_row
  }

  merged_df <- do.call(rbind, merged_list)
  result <- rbind(unique_rows, merged_df)

  # Restore row names if original had them
  if (!is.null(rownames(metabolite_df))) {
    rownames(result) <- NULL  # Reset to avoid conflicts
  }

  return(result)
}


#' Prepare intensity data matrix for downstream analysis
#'
#' Extracts sample intensity columns, creates a composite XCMS.id, and reorders
#' columns for analysis. Replaces the hardcoded column indexing in original scripts.
#'
#' @param metabolite_df Data frame with metabolite data
#' @param sample_names Character vector of sample column names. If NULL, auto-detected.
#' @param id_col Name of the metabolite ID column (default: "XCMS.id")
#' @param mode_col Name of the ion mode column (default: "model")
#' @return Data frame with columns: XCMS.id, sample intensities..., model
prepare_intensity_data <- function(metabolite_df,
                                   sample_names = NULL,
                                   id_col = "XCMS.id",
                                   mode_col = "model") {
  if (is.null(sample_names)) {
    sample_names <- extract_sample_names(metabolite_df)
  }

  if (length(sample_names) == 0) {
    stop("No sample columns detected in metabolite data.")
  }

  # Validate required columns exist
  if (!id_col %in% colnames(metabolite_df)) {
    stop("ID column not found: ", id_col)
  }
  if (!mode_col %in% colnames(metabolite_df)) {
    stop("Mode column not found: ", mode_col)
  }

  # Extract intensities and create composite ID
  intensity_data <- metabolite_df[, sample_names, drop = FALSE]
  intensity_data$XCMS.id <- paste0(metabolite_df[[id_col]], "_", metabolite_df[[mode_col]])
  intensity_data$model <- metabolite_df[[mode_col]]

  # Reorder: XCMS.id, samples, model
  intensity_data <- intensity_data[, c("XCMS.id", sample_names, "model")]

  return(intensity_data)
}


#' Select metabolites with reliable MS2 or matched MS1+MS2 annotations
#'
#' Implements the annotation filtering logic from the original pipeline:
#' 1. Keep MS2-only annotations (no MS1 match)
#' 2. Keep MS1+MS2 annotations where names or IDs match
#'
#' @param metabolite_df Data frame with MS1 and MS2 annotation columns
#' @return Filtered data frame with only reliably annotated metabolites.
#'
#' @note Fixes the original bug where 'tmp$MS2Metabolite' was passed as a
#'   string literal inside grepl() instead of the variable value.
select_annt <- function(metabolite_df) {
  # MS1 annotations from HMDB or KEGG
  ms1_ant <- dplyr::bind_rows(
    metabolite_df %>% dplyr::filter(NumberMS1hmdb > 0),
    metabolite_df %>% dplyr::filter(NumberMS1kegg > 0)
  ) %>% dplyr::distinct()

  # MS2 annotations (non-empty)
  ms2_ant <- metabolite_df %>%
    dplyr::filter(MS2Metabolite != "-" & MS2Metabolite != "" & !is.na(MS2Metabolite))

  # MS2-only annotations (no MS1)
  ms2_ant_only <- ms2_ant %>%
    dplyr::filter(NumberMS1hmdb == 0, NumberMS1kegg == 0)

  # Overlap: MS1 + MS2
  ms1_2_ant <- dplyr::bind_rows(
    ms1_ant %>% dplyr::filter(MS2Metabolite != "-" & !is.na(MS2Metabolite)),
    ms2_ant %>% dplyr::filter(NumberMS1hmdb > 0),
    ms2_ant %>% dplyr::filter(NumberMS1kegg > 0)
  ) %>% dplyr::distinct()

  # Check MS1-MS2 matches (name or ID)
  has_match <- sapply(seq_len(nrow(ms1_2_ant)), function(i) {
    tmp <- ms1_2_ant[i, ]
    j <- 0

    # Name match in HMDB
    if (!is.na(tmp$MS2Metabolite) && !is.na(tmp$MS1hmdbName) &&
        grepl(tmp$MS2Metabolite, tmp$MS1hmdbName, perl = TRUE)) {
      j <- j + 1
    }
    # Name match in KEGG
    if (!is.na(tmp$MS2Metabolite) && !is.na(tmp$MS1keggName) &&
        grepl(tmp$MS2Metabolite, tmp$MS1keggName, perl = TRUE)) {
      j <- j + 1
    }
    # ID match in HMDB
    if (!is.na(tmp$MS2hmd) && tmp$MS2hmd != "-" && tmp$MS2hmd != "" &&
        !is.na(tmp$MS1hmdbID) && grepl(tmp$MS2hmd, tmp$MS1hmdbID)) {
      j <- j + 1
    }
    # ID match via KEGG bridge
    if (!is.na(tmp$MS2kegg) && tmp$MS2kegg != "-" && tmp$MS2kegg != "" &&
        !is.na(tmp$MS1hmdbTokegg) && grepl(tmp$MS2kegg, tmp$MS1hmdbTokegg)) {
      j <- j + 1
    }
    # ID match in KEGG
    if (!is.na(tmp$MS2kegg) && tmp$MS2kegg != "-" && tmp$MS2kegg != "" &&
        !is.na(tmp$MS1keggID) && grepl(tmp$MS2kegg, tmp$MS1keggID)) {
      j <- j + 1
    }

    j > 0
  })

  ms1_2_ant_match <- ms1_2_ant[has_match, ]

  dplyr::bind_rows(ms2_ant_only, ms1_2_ant_match)
}


#-------------------------------------------------------------------------------
# 2. Clinical Data Loading
#-------------------------------------------------------------------------------

#' Load and clean sample/clinical information
#'
#' @param spinfo_file Path to sample information CSV (default: SPINFO_FILE from config)
#' @param remove_qc Whether to remove QC samples (default: TRUE)
#' @param with_panss_only Whether to retain only samples with PANSS data (default: FALSE)
#' @return Cleaned data frame with sample information.
load_sample_info <- function(spinfo_file = SPINFO_FILE,
                             remove_qc = TRUE,
                             with_panss_only = FALSE) {
  if (!file.exists(spinfo_file)) {
    stop("Sample info file not found: ", spinfo_file)
  }

  spinfo <- read.csv(spinfo_file, stringsAsFactors = FALSE)

  # Standardize first two column names
  if (ncol(spinfo) >= 2) {
    colnames(spinfo)[1:2] <- c("sample_id", "subject_id")
  }

  # Factorize gender
  if ("gender" %in% colnames(spinfo)) {
    spinfo$gender <- factor(spinfo$gender, levels = c("M", "F"))
  }

  # Remove QC samples
  if (remove_qc && "group" %in% colnames(spinfo)) {
    spinfo <- spinfo %>% dplyr::filter(group != "QC")
  }

  # Retain only PANSS samples
  if (with_panss_only && "with_panss" %in% colnames(spinfo)) {
    spinfo <- spinfo %>% dplyr::filter(with_panss == "Y")
  }

  # Sort by group if available
  if ("group" %in% colnames(spinfo)) {
    spinfo <- spinfo %>% dplyr::arrange(group)
  }

  # Set row names if sample name column exists
  if ("sample_name" %in% colnames(spinfo)) {
    rownames(spinfo) <- spinfo$sample_name
  } else if ("sample_id" %in% colnames(spinfo)) {
    rownames(spinfo) <- spinfo$sample_id
  }

  return(spinfo)
}


#' Extract PANSS score vectors from sample info
#'
#' @param spinfo Sample info data frame (output from load_sample_info)
#' @return List of named numeric vectors: total, positive, negative.
extract_panss_scores <- function(spinfo) {
  required_cols <- c("PANSS.total", "PANSS.positive", "PANSS.negative")
  missing <- setdiff(required_cols, colnames(spinfo))
  if (length(missing) > 0) {
    stop("Missing PANSS columns: ", paste(missing, collapse = ", "))
  }

  # Determine sample ID vector (fix: was returning the string "Sample" as a single name)
  id_vec <- if ("Sample" %in% colnames(spinfo)) spinfo$Sample else rownames(spinfo)

  list(
    total    = setNames(as.numeric(spinfo$PANSS.total), id_vec),
    positive = setNames(as.numeric(spinfo$PANSS.positive), id_vec),
    negative = setNames(as.numeric(spinfo$PANSS.negative), id_vec)
  )
}


#-------------------------------------------------------------------------------
# 3. Statistical Utilities
#-------------------------------------------------------------------------------

#' Flatten correlation matrix to long-format data frame
#'
#' Converts a symmetric correlation matrix (and its p-value matrix) to a tidy
#' data frame with one row per upper-triangle pair.
#'
#' @param cormat Correlation matrix (square, symmetric)
#' @param pmat Matrix of p-values (same dimensions as cormat)
#' @param method Correction method for q-values (default: "BH")
#' @return Data frame with columns: row, column, cor, p, q
#'
#' @note Replaces the typo-ridden 'flatternCorrMatrix' in the original scripts.
flatten_corr_matrix <- function(cormat, pmat, method = "BH") {
  ut <- upper.tri(cormat)
  result <- data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = cormat[ut],
    p = pmat[ut],
    stringsAsFactors = FALSE
  )
  result$q <- p.adjust(result$p, method = method)
  return(result)
}


#' Safe t-test wrapper with error handling
#'
#' Runs t.test and returns a list with estimate, p.value, and conf.int.
#' If the test fails (e.g., identical groups), returns NA instead of crashing.
#'
#' @param x Numeric vector (group 1)
#' @param y Numeric vector (group 2)
#' @param ... Additional arguments passed to t.test()
#' @return List with elements: estimate, p.value, conf.low, conf.high, error (NULL if success)
safe_t_test <- function(x, y, ...) {
  res <- tryCatch({
    tt <- t.test(x, y, ...)
    list(
      estimate = as.numeric(tt$estimate[1] - tt$estimate[2]),
      p.value = tt$p.value,
      conf.low = tt$conf.int[1],
      conf.high = tt$conf.int[2],
      error = NULL
    )
  }, error = function(e) {
    warning("t.test failed: ", conditionMessage(e))
    list(
      estimate = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      error = conditionMessage(e)
    )
  })
  return(res)
}


#' Fisher's exact test wrapper with edge case handling
#'
#' Handles cases where all/none have variant in a group, which causes
#' fisher.test to fail with "workspace" errors.
#'
#' @param mat 2x2 contingency table
#' @return List with OR estimate and p-value (NA if test fails or is degenerate)
safe_fisher_test <- function(mat) {
  if (any(mat < 0)) {
    return(list(OR = NA_real_, p_value = NA_real_))
  }

  # Degenerate cases: all or none in a group
  if (sum(mat[1, ]) == 0 || sum(mat[2, ]) == 0 ||
      mat[1, 1] == sum(mat[1, ]) || mat[2, 1] == sum(mat[2, ])) {
    OR <- if (mat[1, 1] > 0 && mat[2, 1] == 0) Inf
          else if (mat[1, 1] == 0 && mat[2, 1] > 0) 0
          else 1
    return(list(OR = OR, p_value = NA_real_))
  }

  ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
  if (is.null(ft)) {
    return(list(OR = NA_real_, p_value = NA_real_))
  }

  list(OR = as.numeric(ft$estimate), p_value = ft$p.value)
}


#-------------------------------------------------------------------------------
# 4. Variant Data Processing
#-------------------------------------------------------------------------------

#' Extract sample IDs and impact levels from VEP records
#'
#' Parses the VEP output format to extract per-sample variant information.
#' Replaces the fragile apply()+strsplit() logic in the original scripts.
#'
#' @param vep_records Data frame subset of VEP results for a single variant
#' @param sample_col Column name or index containing the sample/IND field (default: 14)
#' @return Data frame with columns: sample, impact, impact_num
extract_vep_samples <- function(vep_records, sample_col = 14) {
  if (ncol(vep_records) < sample_col) {
    stop("VEP records have fewer columns than expected (", sample_col, ")")
  }

  variant_info <- apply(vep_records, 1, function(row) {
    fields <- unlist(strsplit(as.character(row[sample_col]), ";"))

    sample_field <- fields[grepl("^IND=", fields)]
    sample_id <- gsub("^IND=", "", sample_field)

    impact_field <- fields[grepl("^IMPACT=", fields)]
    impact_level <- gsub("^IMPACT=", "", impact_field)

    data.frame(
      sample = sample_id,
      impact = impact_level,
      stringsAsFactors = FALSE
    )
  })

  variant_df <- do.call(rbind, variant_info)
  variant_df$impact_num <- EFFECT_LEVELS[variant_df$impact]
  variant_df$impact_num[is.na(variant_df$impact_num)] <- 1  # MODIFIER as default

  return(variant_df)
}


#' Build variant presence and effect matrices from VEP data
#'
#' Creates binary presence and impact-level matrices for downstream analysis.
#'
#' @param filtered_vep VEP data frame filtered to target genes/variants
#' @param sample_names Character vector of sample names
#' @param id_col Column name for variant ID (default: "X.Uploaded_variation")
#' @return List with matrices: presence (0/1), effect (1-4)
build_variant_matrices <- function(filtered_vep, sample_names,
                                   id_col = "X.Uploaded_variation") {
  variant_list <- unique(filtered_vep[[id_col]])
  n_variants <- length(variant_list)
  n_samples <- length(sample_names)

  var_presence <- matrix(0, nrow = n_variants, ncol = n_samples)
  var_effect <- matrix(0, nrow = n_variants, ncol = n_samples)
  colnames(var_presence) <- colnames(var_effect) <- sample_names
  rownames(var_presence) <- rownames(var_effect) <- variant_list

  for (i in seq_len(n_variants)) {
    variant_id <- variant_list[i]
    variant_records <- dplyr::filter(filtered_vep, !!rlang::sym(id_col) == variant_id)

    variant_df <- extract_vep_samples(variant_records)
    carrier_samples <- unique(variant_df$sample)
    carrier_samples <- carrier_samples[carrier_samples %in% sample_names]

    for (sample_id in carrier_samples) {
      sample_effects <- variant_df$impact_num[variant_df$sample == sample_id]
      var_presence[i, sample_id] <- 1
      var_effect[i, sample_id] <- max(sample_effects, na.rm = TRUE)
    }
  }

  list(presence = var_presence, effect = var_effect)
}


#-------------------------------------------------------------------------------
# 5. Network Data Preparation
#-------------------------------------------------------------------------------

#' Extract significant variant-metabolite association pairs
#'
#' Filters association results for pairs significant in all four tests
#' (intensity and rank, SCZ vs NOR and SCZ vs SCZ).
#'
#' @param assoc_df Association results data frame (from 3.3 output)
#' @param threshold Significance threshold (default: PADJ_THRESHOLD)
#' @param prefix Column name prefix indicating test type (e.g., "intensity", "rank")
#' @return Filtered data frame with only significant associations.
filter_significant_assoc <- function(assoc_df,
                                     threshold = PADJ_THRESHOLD,
                                     padj_cols = NULL) {
  if (is.null(padj_cols)) {
    # Auto-detect all columns ending in .padj
    padj_cols <- grep("\\.padj$", colnames(assoc_df), value = TRUE)
  }

  if (length(padj_cols) < 4) {
    warning("Expected at least 4 padj columns but found ", length(padj_cols), ": ",
            paste(padj_cols, collapse = ", "))
  }

  mask <- apply(assoc_df[, padj_cols, drop = FALSE], 1, function(x) {
    all(x < threshold, na.rm = TRUE)
  })

  assoc_df[mask, ]
}


#' Determine metabolite node type (DEM or NonDEM)
#'
#' @param node_ids Character vector of metabolite node IDs
#' @param dem_ids Character vector of DEM XCMS.id_model identifiers
#' @param nondem_ids Character vector of non-DEM XCMS.id_model identifiers
#' @return Character vector of node types ("DEM", "NonDEM", or "-")
classify_metabolite_node <- function(node_ids, dem_ids, nondem_ids) {
  ifelse(node_ids %in% dem_ids, "DEM",
         ifelse(node_ids %in% nondem_ids, "NonDEM", "-"))
}


#-------------------------------------------------------------------------------
# 6. Validation & Debugging
#-------------------------------------------------------------------------------

#' Validate that all expected columns exist in a data frame
#'
#' @param df Data frame to check
#' @param required_cols Character vector of required column names
#' @param context Description of the data frame (for error messages)
#' @return TRUE if all columns exist, otherwise stops with informative message.
validate_columns <- function(df, required_cols, context = "data frame") {
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(
      "Missing required columns in ", context, ":\n",
      "  Required: ", paste(required_cols, collapse = ", "), "\n",
      "  Missing:  ", paste(missing, collapse = ", "), "\n",
      "  Present:  ", paste(head(colnames(df), 10), collapse = ", "),
      if (ncol(df) > 10) " ..." else ""
    )
  }
  invisible(TRUE)
}


#' Print summary of a data frame (dimensions, key columns, NA counts)
#'
#' @param df Data frame to summarize
#' @param name Optional name for the summary header
#' @param key_cols Columns to show NA summary for (if NULL, auto-selected)
print_df_summary <- function(df, name = NULL, key_cols = NULL) {
  header <- if (!is.null(name)) paste0("Summary: ", name) else "Data Frame Summary"
  cat("\n========== ", header, " ==========\n", sep = "")
  cat("Dimensions: ", nrow(df), " rows x ", ncol(df), " columns\n", sep = "")

  if (!is.null(key_cols)) {
    na_counts <- sapply(df[, key_cols, drop = FALSE], function(x) sum(is.na(x)))
    cat("NA counts:\n")
    print(na_counts)
  }
  cat("========================================\n\n")
}


#-------------------------------------------------------------------------------
# 7. Backward Compatibility Aliases
#-------------------------------------------------------------------------------
# These allow old scripts to call the old function names while using the
# new, corrected implementations.

#' @rdname flatten_corr_matrix
flatternCorrMatrix <- flatten_corr_matrix  # legacy alias (fixes typo)
