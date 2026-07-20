# ============================================================
# SDY2583 CP22
# STEP 5 SAFE: Targeted threshold sensitivity
# B-cell / humoral remodeling panel
#
# Threshold configurations:
#   main
#   permissive: positive-marker thresholds -0.2; DUMP_LOW +0.2
#   stringent : positive-marker thresholds +0.2; DUMP_LOW -0.2
#
# Models:
#   variable ~ disease_group + age_for_model + sex
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "purrr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install("flowCore", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(broom)
  library(purrr)
})

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP22")

step1_rdata <- file.path(
  analysis_dir,
  "11_RData",
  "SDY2583_CP22_STEP1_fcs_inventory_marker_QC.RData"
)

step3a_rdata <- file.path(
  analysis_dir,
  "11_RData",
  "SDY2583_CP22_STEP3A_metadata_merge_age_QC.RData"
)

step4_rdata <- file.path(
  analysis_dir,
  "11_RData",
  "SDY2583_CP22_STEP4_composite_scores.RData"
)

if (!file.exists(step1_rdata)) stop("Step 1 RData bulunamadı: ", step1_rdata)
if (!file.exists(step3a_rdata)) stop("Step 3A RData bulunamadı: ", step3a_rdata)
if (!file.exists(step4_rdata)) stop("Step 4 RData bulunamadı: ", step4_rdata)

out_dir <- file.path(analysis_dir, "07_threshold_sensitivity")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Load required objects
# ------------------------------------------------------------

load(step1_rdata)
if (!exists("fcs_files")) stop("Step 1 RData içinde fcs_files bulunamadı.")

load(step3a_rdata)
if (!exists("cp22_analysis_data")) stop("Step 3A RData içinde cp22_analysis_data bulunamadı.")

metadata_model <- cp22_analysis_data %>%
  dplyr::select(
    subject_id,
    disease_group,
    age_for_model,
    sex,
    any_of(c(
      "cancer_subgroup",
      "therapy_status_4level",
      "chemotherapy",
      "targeted_therapy",
      "any_immunotherapy",
      "ici_immunotherapy",
      "therapy_line_number",
      "time_from_start_days"
    ))
  ) %>%
  dplyr::distinct(subject_id, .keep_all = TRUE)

load(step4_rdata)
if (!exists("composite_definitions")) stop("Step 4 RData içinde composite_definitions bulunamadı.")

# ------------------------------------------------------------
# 4. Threshold definitions
# ------------------------------------------------------------

main_thresholds <- list(
  DUMP_LOW = 1.5,
  CD19 = 2.0,
  CD27 = 1.5,
  IgD = 1.5,
  IgM = 1.5,
  IgA = 1.5,
  IgG = 1.5,
  CD38 = 1.5,
  CD38_HIGH = 2.2,
  CD138 = 1.5,
  CD24 = 1.5,
  CD10 = 1.5,
  CD39 = 1.5
)

make_threshold_config <- function(config_name) {

  th <- main_thresholds

  if (config_name == "main") {
    return(th)
  }

  positive_markers <- setdiff(names(th), "DUMP_LOW")

  if (config_name == "permissive") {
    for (m in positive_markers) th[[m]] <- th[[m]] - 0.2
    th$DUMP_LOW <- th$DUMP_LOW + 0.2
    return(th)
  }

  if (config_name == "stringent") {
    for (m in positive_markers) th[[m]] <- th[[m]] + 0.2
    th$DUMP_LOW <- th$DUMP_LOW - 0.2
    return(th)
  }

  stop("Unknown threshold config: ", config_name)
}

threshold_configs <- c("main", "permissive", "stringent")

threshold_config_table <- dplyr::bind_rows(
  lapply(threshold_configs, function(cfg) {
    th <- make_threshold_config(cfg)
    tibble::tibble(
      threshold_set = cfg,
      marker_or_gate = names(th),
      threshold = unlist(th)
    )
  })
)

readr::write_csv(
  threshold_config_table,
  file.path(out_dir, "SDY2583_CP22_threshold_configurations_STEP5.csv")
)

# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

extract_subject_id <- function(x) {
  b <- basename(x)
  id <- stringr::str_extract(b, "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(b, "\\.fcs$"), id)
}

get_marker_map <- function(ff) {
  pp <- Biobase::pData(flowCore::parameters(ff))
  channels <- as.character(pp$name)
  markers <- as.character(pp$desc)
  idx <- is.na(markers) | markers == ""
  markers[idx] <- channels[idx]

  tibble::tibble(
    parameter_index = seq_along(channels),
    channel_name = channels,
    marker_desc = markers
  )
}

find_channel_exact_marker <- function(map, marker) {
  idx <- which(map$marker_desc == marker)
  if (length(idx) >= 1) return(map$channel_name[idx[1]])
  NA_character_
}

find_channel_regex_marker <- function(map, pattern) {
  idx <- which(grepl(pattern, map$marker_desc, ignore.case = TRUE))
  if (length(idx) >= 1) return(map$channel_name[idx[1]])
  NA_character_
}

find_channel_by_channel <- function(map, channel) {
  idx <- which(map$channel_name == channel)
  if (length(idx) >= 1) return(map$channel_name[idx[1]])
  NA_character_
}

build_channel_map <- function(map) {

  ch <- list()

  ch$CD19 <- find_channel_exact_marker(map, "CD19")
  ch$CD27 <- find_channel_exact_marker(map, "CD27")
  ch$IgD  <- find_channel_exact_marker(map, "IgD")
  ch$IgM  <- find_channel_exact_marker(map, "IgM")
  ch$IgA  <- find_channel_exact_marker(map, "IgA")
  ch$IgG  <- find_channel_exact_marker(map, "IgG")
  ch$CD38 <- find_channel_exact_marker(map, "CD38")
  ch$CD138 <- find_channel_exact_marker(map, "CD138")
  ch$CD24 <- find_channel_exact_marker(map, "CD24")
  ch$CD10 <- find_channel_exact_marker(map, "CD10")
  ch$CD39 <- find_channel_exact_marker(map, "CD39")

  ch$DUMP <- find_channel_regex_marker(map, "^Viability")
  if (is.na(ch$DUMP)) {
    ch$DUMP <- find_channel_by_channel(map, "BV510-A")
  }

  dump_marker <- NA_character_
  if (!is.na(ch$DUMP)) {
    dump_marker <- map$marker_desc[which(map$channel_name == ch$DUMP)[1]]
  }

  list(
    channels = ch,
    dump_marker_variant = dump_marker
  )
}

safe_pct <- function(num, den) {
  ifelse(is.na(den) | den <= 0, NA_real_, 100 * num / den)
}

safe_median <- function(x) {
  if (length(x) == 0) return(NA_real_)
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}

get_spill_matrix <- function(ff) {
  kk <- flowCore::keyword(ff)
  for (nm in c("SPILL", "$SPILLOVER", "SPILLOVER")) {
    if (nm %in% names(kk)) {
      sp <- kk[[nm]]
      if (is.matrix(sp)) return(sp)
    }
  }
  NULL
}

apply_compensation_safely <- function(ff) {
  sp <- get_spill_matrix(ff)

  if (is.null(sp)) {
    return(list(ff = ff, compensation_applied = FALSE, compensation_error = NA_character_))
  }

  tryCatch({
    ff2 <- flowCore::compensate(ff, sp)
    list(ff = ff2, compensation_applied = TRUE, compensation_error = NA_character_)
  }, error = function(e) {
    list(ff = ff, compensation_applied = FALSE, compensation_error = as.character(e$message))
  })
}

transform_safely <- function(ff) {

  ex_names <- colnames(flowCore::exprs(ff))
  fluor_ch <- ex_names[!grepl("FSC|SSC|Time", ex_names, ignore.case = TRUE)]

  if (length(fluor_ch) == 0) {
    return(list(ff = ff, transform_applied = FALSE, transform_error = "No fluorescence channels detected"))
  }

  trans <- flowCore::logicleTransform(
    transformationId = "fixed_logicle",
    w = 0.5,
    t = 262144,
    m = 4.5,
    a = 0
  )

  tryCatch({
    tf <- flowCore::transformList(fluor_ch, trans)
    ff2 <- flowCore::transform(ff, tf)
    list(ff = ff2, transform_applied = TRUE, transform_error = NA_character_)
  }, error = function(e) {
    list(ff = ff, transform_applied = FALSE, transform_error = as.character(e$message))
  })
}

get_vec <- function(ex, ch) {
  if (is.na(ch) || !(ch %in% colnames(ex))) {
    return(rep(NA_real_, nrow(ex)))
  }
  as.numeric(ex[, ch])
}

# ------------------------------------------------------------
# 6. Feature extraction for one file under one threshold set
# ------------------------------------------------------------

extract_one_cp22_threshold <- function(fp, threshold_set, thresholds) {

  tryCatch({

    ff_raw <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE
    )

    marker_map_file <- get_marker_map(ff_raw)
    channel_info <- build_channel_map(marker_map_file)
    ch <- channel_info$channels

    required_markers <- c(
      "DUMP", "CD19", "CD27", "IgD", "IgM", "IgA", "IgG",
      "CD38", "CD138", "CD24", "CD10", "CD39"
    )

    missing_markers <- required_markers[
      vapply(ch[required_markers], function(x) is.na(x), logical(1))
    ]

    if (length(missing_markers) > 0) {
      stop("Missing required channels/markers: ", paste(missing_markers, collapse = ", "))
    }

    comp <- apply_compensation_safely(ff_raw)
    trans <- transform_safely(comp$ff)

    ff <- trans$ff
    ex <- flowCore::exprs(ff)

    total_events <- nrow(ex)

    dump <- get_vec(ex, ch$DUMP)
    cd19 <- get_vec(ex, ch$CD19)
    cd27 <- get_vec(ex, ch$CD27)
    igd  <- get_vec(ex, ch$IgD)
    igm  <- get_vec(ex, ch$IgM)
    iga  <- get_vec(ex, ch$IgA)
    igg  <- get_vec(ex, ch$IgG)
    cd38 <- get_vec(ex, ch$CD38)
    cd138 <- get_vec(ex, ch$CD138)
    cd24 <- get_vec(ex, ch$CD24)
    cd10 <- get_vec(ex, ch$CD10)
    cd39 <- get_vec(ex, ch$CD39)

    dump_low <- dump <= thresholds$DUMP_LOW
    b_cell <- dump_low & cd19 > thresholds$CD19

    n_dump_low <- sum(dump_low, na.rm = TRUE)
    n_b <- sum(b_cell, na.rm = TRUE)

    b_idx <- if (n_b == 0) rep(FALSE, total_events) else b_cell

    cd27_pos <- b_idx & cd27 > thresholds$CD27
    igd_pos <- b_idx & igd > thresholds$IgD
    igm_pos <- b_idx & igm > thresholds$IgM
    iga_pos <- b_idx & iga > thresholds$IgA
    igg_pos <- b_idx & igg > thresholds$IgG
    cd38_pos <- b_idx & cd38 > thresholds$CD38
    cd38_high <- b_idx & cd38 > thresholds$CD38_HIGH
    cd138_pos <- b_idx & cd138 > thresholds$CD138
    cd24_pos <- b_idx & cd24 > thresholds$CD24
    cd10_pos <- b_idx & cd10 > thresholds$CD10
    cd39_pos <- b_idx & cd39 > thresholds$CD39

    naive_like <- b_idx & igd > thresholds$IgD & cd27 <= thresholds$CD27
    unswitched_memory_like <- b_idx & igd > thresholds$IgD & cd27 > thresholds$CD27
    switched_memory_like <- b_idx & igd <= thresholds$IgD & cd27 > thresholds$CD27
    double_negative_like <- b_idx & igd <= thresholds$IgD & cd27 <= thresholds$CD27

    igm_unswitched_memory_like <- unswitched_memory_like & igm > thresholds$IgM
    iga_switched_memory_like <- switched_memory_like & iga > thresholds$IgA
    igg_switched_memory_like <- switched_memory_like & igg > thresholds$IgG
    iga_igg_double_negative_switched_like <- switched_memory_like &
      iga <= thresholds$IgA &
      igg <= thresholds$IgG

    plasmablast_like <- b_idx & cd38 > thresholds$CD38_HIGH & cd27 > thresholds$CD27
    class_switched_plasmablast_like <- plasmablast_like & igd <= thresholds$IgD
    iga_plasmablast_like <- plasmablast_like & iga > thresholds$IgA
    igg_plasmablast_like <- plasmablast_like & igg > thresholds$IgG

    plasma_cell_like <- b_idx & cd138 > thresholds$CD138
    cd38_cd138_plasma_cell_like <- b_idx & cd38 > thresholds$CD38 & cd138 > thresholds$CD138

    transitional_like <- b_idx & cd24 > thresholds$CD24 & cd38 > thresholds$CD38
    immature_transitional_like <- b_idx & cd10 > thresholds$CD10
    cd10_cd24_cd38_transitional_like <- b_idx &
      cd10 > thresholds$CD10 &
      cd24 > thresholds$CD24 &
      cd38 > thresholds$CD38

    cd39_regulatory_like <- b_idx & cd39 > thresholds$CD39
    cd39_cd24_regulatory_like <- b_idx & cd39 > thresholds$CD39 & cd24 > thresholds$CD24
    cd39_cd24_cd38_regulatory_transitional_like <- b_idx &
      cd39 > thresholds$CD39 &
      cd24 > thresholds$CD24 &
      cd38 > thresholds$CD38

    n_naive_like <- sum(naive_like, na.rm = TRUE)
    n_unswitched_memory_like <- sum(unswitched_memory_like, na.rm = TRUE)
    n_switched_memory_like <- sum(switched_memory_like, na.rm = TRUE)
    n_double_negative_like <- sum(double_negative_like, na.rm = TRUE)

    n_cd27_pos <- sum(cd27_pos, na.rm = TRUE)
    n_igd_pos <- sum(igd_pos, na.rm = TRUE)
    n_igm_pos <- sum(igm_pos, na.rm = TRUE)
    n_iga_pos <- sum(iga_pos, na.rm = TRUE)
    n_igg_pos <- sum(igg_pos, na.rm = TRUE)
    n_cd38_pos <- sum(cd38_pos, na.rm = TRUE)
    n_cd38_high <- sum(cd38_high, na.rm = TRUE)
    n_cd138_pos <- sum(cd138_pos, na.rm = TRUE)
    n_cd24_pos <- sum(cd24_pos, na.rm = TRUE)
    n_cd10_pos <- sum(cd10_pos, na.rm = TRUE)
    n_cd39_pos <- sum(cd39_pos, na.rm = TRUE)

    n_igm_unswitched_memory_like <- sum(igm_unswitched_memory_like, na.rm = TRUE)
    n_iga_switched_memory_like <- sum(iga_switched_memory_like, na.rm = TRUE)
    n_igg_switched_memory_like <- sum(igg_switched_memory_like, na.rm = TRUE)
    n_iga_igg_double_negative_switched_like <- sum(iga_igg_double_negative_switched_like, na.rm = TRUE)

    n_plasmablast_like <- sum(plasmablast_like, na.rm = TRUE)
    n_class_switched_plasmablast_like <- sum(class_switched_plasmablast_like, na.rm = TRUE)
    n_iga_plasmablast_like <- sum(iga_plasmablast_like, na.rm = TRUE)
    n_igg_plasmablast_like <- sum(igg_plasmablast_like, na.rm = TRUE)
    n_plasma_cell_like <- sum(plasma_cell_like, na.rm = TRUE)
    n_cd38_cd138_plasma_cell_like <- sum(cd38_cd138_plasma_cell_like, na.rm = TRUE)

    n_transitional_like <- sum(transitional_like, na.rm = TRUE)
    n_immature_transitional_like <- sum(immature_transitional_like, na.rm = TRUE)
    n_cd10_cd24_cd38_transitional_like <- sum(cd10_cd24_cd38_transitional_like, na.rm = TRUE)

    n_cd39_regulatory_like <- sum(cd39_regulatory_like, na.rm = TRUE)
    n_cd39_cd24_regulatory_like <- sum(cd39_cd24_regulatory_like, na.rm = TRUE)
    n_cd39_cd24_cd38_regulatory_transitional_like <- sum(cd39_cd24_cd38_regulatory_transitional_like, na.rm = TRUE)

    tibble::tibble(
      threshold_set = threshold_set,
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      feature_ok = TRUE,
      feature_error = NA_character_,
      transform_method = ifelse(trans$transform_applied, "compensated_fixed_logicle", "compensated_untransformed_fallback"),
      compensation_applied = comp$compensation_applied,
      compensation_error = comp$compensation_error,
      transform_applied = trans$transform_applied,
      transform_error = trans$transform_error,
      dump_channel = ch$DUMP,
      dump_marker_variant = channel_info$dump_marker_variant,
      total_events = total_events,
      n_dump_low = n_dump_low,
      n_cd19_b = n_b,
      pct_dump_low_total = safe_pct(n_dump_low, total_events),
      pct_cd19_b_total = safe_pct(n_b, total_events),
      pct_cd19_b_within_dump_low = safe_pct(n_b, n_dump_low),

      pct_naive_like_within_b = safe_pct(n_naive_like, n_b),
      pct_unswitched_memory_like_within_b = safe_pct(n_unswitched_memory_like, n_b),
      pct_switched_memory_like_within_b = safe_pct(n_switched_memory_like, n_b),
      pct_double_negative_like_within_b = safe_pct(n_double_negative_like, n_b),

      pct_cd27_pos_within_b = safe_pct(n_cd27_pos, n_b),
      pct_igd_pos_within_b = safe_pct(n_igd_pos, n_b),
      pct_igm_pos_within_b = safe_pct(n_igm_pos, n_b),
      pct_iga_pos_within_b = safe_pct(n_iga_pos, n_b),
      pct_igg_pos_within_b = safe_pct(n_igg_pos, n_b),
      pct_cd38_pos_within_b = safe_pct(n_cd38_pos, n_b),
      pct_cd38high_within_b = safe_pct(n_cd38_high, n_b),
      pct_cd138_pos_within_b = safe_pct(n_cd138_pos, n_b),
      pct_cd24_pos_within_b = safe_pct(n_cd24_pos, n_b),
      pct_cd10_pos_within_b = safe_pct(n_cd10_pos, n_b),
      pct_cd39_pos_within_b = safe_pct(n_cd39_pos, n_b),

      pct_igm_unswitched_memory_like_within_b = safe_pct(n_igm_unswitched_memory_like, n_b),
      pct_iga_switched_memory_like_within_b = safe_pct(n_iga_switched_memory_like, n_b),
      pct_igg_switched_memory_like_within_b = safe_pct(n_igg_switched_memory_like, n_b),
      pct_iga_igg_double_negative_switched_like_within_b = safe_pct(n_iga_igg_double_negative_switched_like, n_b),

      pct_iga_within_switched_memory_like = safe_pct(n_iga_switched_memory_like, n_switched_memory_like),
      pct_igg_within_switched_memory_like = safe_pct(n_igg_switched_memory_like, n_switched_memory_like),
      pct_iga_igg_double_negative_within_switched_memory_like = safe_pct(n_iga_igg_double_negative_switched_like, n_switched_memory_like),
      pct_igm_within_unswitched_memory_like = safe_pct(n_igm_unswitched_memory_like, n_unswitched_memory_like),

      pct_plasmablast_like_within_b = safe_pct(n_plasmablast_like, n_b),
      pct_class_switched_plasmablast_like_within_b = safe_pct(n_class_switched_plasmablast_like, n_b),
      pct_iga_plasmablast_like_within_b = safe_pct(n_iga_plasmablast_like, n_b),
      pct_igg_plasmablast_like_within_b = safe_pct(n_igg_plasmablast_like, n_b),
      pct_plasma_cell_like_within_b = safe_pct(n_plasma_cell_like, n_b),
      pct_cd38_cd138_plasma_cell_like_within_b = safe_pct(n_cd38_cd138_plasma_cell_like, n_b),

      pct_transitional_like_within_b = safe_pct(n_transitional_like, n_b),
      pct_immature_transitional_like_within_b = safe_pct(n_immature_transitional_like, n_b),
      pct_cd10_cd24_cd38_transitional_like_within_b = safe_pct(n_cd10_cd24_cd38_transitional_like, n_b),

      pct_cd39_regulatory_like_within_b = safe_pct(n_cd39_regulatory_like, n_b),
      pct_cd39_cd24_regulatory_like_within_b = safe_pct(n_cd39_cd24_regulatory_like, n_b),
      pct_cd39_cd24_cd38_regulatory_transitional_like_within_b = safe_pct(n_cd39_cd24_cd38_regulatory_transitional_like, n_b),

      median_CD19_in_B = safe_median(cd19[b_idx]),
      median_CD27_in_B = safe_median(cd27[b_idx]),
      median_IgD_in_B = safe_median(igd[b_idx]),
      median_IgM_in_B = safe_median(igm[b_idx]),
      median_IgA_in_B = safe_median(iga[b_idx]),
      median_IgG_in_B = safe_median(igg[b_idx]),
      median_CD38_in_B = safe_median(cd38[b_idx]),
      median_CD138_in_B = safe_median(cd138[b_idx]),
      median_CD24_in_B = safe_median(cd24[b_idx]),
      median_CD10_in_B = safe_median(cd10[b_idx]),
      median_CD39_in_B = safe_median(cd39[b_idx])
    )

  }, error = function(e) {

    tibble::tibble(
      threshold_set = threshold_set,
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      feature_ok = FALSE,
      feature_error = as.character(e$message)
    )
  })
}

# ------------------------------------------------------------
# 7. Run threshold extraction
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 5 THRESHOLD SENSITIVITY STARTED\n")
cat("============================================================\n")

threshold_feature_tables <- list()

for (cfg in threshold_configs) {

  cat("\nThreshold set:", cfg, "\n")

  th <- make_threshold_config(cfg)
  rows <- vector("list", length(fcs_files))

  for (i in seq_along(fcs_files)) {
    if (i %% 50 == 0) {
      cat("  processed", i, "of", length(fcs_files), "files\n")
    }
    rows[[i]] <- extract_one_cp22_threshold(fcs_files[i], cfg, th)
  }

  threshold_feature_tables[[cfg]] <- dplyr::bind_rows(rows)
}

threshold_feature_table <- dplyr::bind_rows(threshold_feature_tables)

# ------------------------------------------------------------
# 8. Merge metadata
# ------------------------------------------------------------

threshold_analysis_data <- threshold_feature_table %>%
  dplyr::left_join(metadata_model, by = "subject_id") %>%
  dplyr::mutate(
    disease_group = factor(
      disease_group,
      levels = c("Healthy control", "Cancer patient")
    ),
    sex = factor(sex),
    model_ready_age_sex = feature_ok == TRUE &
      !is.na(disease_group) &
      !is.na(age_for_model) &
      !is.na(sex)
  )

# ------------------------------------------------------------
# 9. Composite scores under each threshold set
# ------------------------------------------------------------

zscore <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (is.na(sig) || sig == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mu) / sig
}

make_composite_for_df <- function(df, comp_name, defs) {

  defs_comp <- defs %>%
    dplyr::filter(composite_score == comp_name, feature %in% names(df))

  if (nrow(defs_comp) == 0) {
    return(rep(NA_real_, nrow(df)))
  }

  component_mat <- matrix(NA_real_, nrow = nrow(df), ncol = nrow(defs_comp))

  for (i in seq_len(nrow(defs_comp))) {
    f <- defs_comp$feature[i]
    direction <- defs_comp$direction[i]
    component_mat[, i] <- zscore(df[[f]]) * direction
  }

  out <- rowMeans(component_mat, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

composite_names <- unique(composite_definitions$composite_score)

threshold_analysis_data_with_scores <- threshold_analysis_data %>%
  dplyr::group_by(threshold_set) %>%
  dplyr::group_modify(function(.x, .y) {
    temp <- .x
    for (cs in composite_names) {
      temp[[cs]] <- make_composite_for_df(temp, cs, composite_definitions)
    }
    temp
  }) %>%
  dplyr::ungroup()

# ------------------------------------------------------------
# 10. Variable dictionary
# ------------------------------------------------------------

feature_dictionary <- tibble::tribble(
  ~variable, ~variable_type, ~module,
  "pct_cd19_b_total", "feature", "B_cell_composition",
  "pct_cd19_b_within_dump_low", "feature", "B_cell_composition",
  "pct_naive_like_within_b", "feature", "core_B_cell_states",
  "pct_unswitched_memory_like_within_b", "feature", "core_B_cell_states",
  "pct_switched_memory_like_within_b", "feature", "core_B_cell_states",
  "pct_double_negative_like_within_b", "feature", "core_B_cell_states",
  "pct_cd27_pos_within_b", "feature", "core_B_cell_states",
  "pct_igd_pos_within_b", "feature", "core_B_cell_states",
  "pct_igm_pos_within_b", "feature", "immunoglobulin_isotype_B_cells",
  "pct_iga_pos_within_b", "feature", "immunoglobulin_isotype_B_cells",
  "pct_igg_pos_within_b", "feature", "immunoglobulin_isotype_B_cells",
  "pct_igm_unswitched_memory_like_within_b", "feature", "immunoglobulin_isotype_B_cells",
  "pct_iga_switched_memory_like_within_b", "feature", "immunoglobulin_isotype_B_cells",
  "pct_igg_switched_memory_like_within_b", "feature", "immunoglobulin_isotype_B_cells",
  "pct_iga_igg_double_negative_switched_like_within_b", "feature", "immunoglobulin_isotype_B_cells",
  "pct_iga_within_switched_memory_like", "feature", "immunoglobulin_isotype_B_cells",
  "pct_igg_within_switched_memory_like", "feature", "immunoglobulin_isotype_B_cells",
  "pct_iga_igg_double_negative_within_switched_memory_like", "feature", "immunoglobulin_isotype_B_cells",
  "pct_igm_within_unswitched_memory_like", "feature", "immunoglobulin_isotype_B_cells",
  "pct_cd38_pos_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_cd38high_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_cd138_pos_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_plasmablast_like_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_class_switched_plasmablast_like_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_iga_plasmablast_like_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_igg_plasmablast_like_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_plasma_cell_like_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_cd38_cd138_plasma_cell_like_within_b", "feature", "plasmablast_plasma_cell_axis",
  "pct_cd24_pos_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_cd10_pos_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_transitional_like_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_immature_transitional_like_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_cd10_cd24_cd38_transitional_like_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_cd39_pos_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_cd39_regulatory_like_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_cd39_cd24_regulatory_like_within_b", "feature", "transitional_regulatory_like_B_cells",
  "pct_cd39_cd24_cd38_regulatory_transitional_like_within_b", "feature", "transitional_regulatory_like_B_cells",
  "median_CD19_in_B", "feature", "B_cell_marker_medians",
  "median_CD27_in_B", "feature", "B_cell_marker_medians",
  "median_IgD_in_B", "feature", "B_cell_marker_medians",
  "median_IgM_in_B", "feature", "B_cell_marker_medians",
  "median_IgA_in_B", "feature", "B_cell_marker_medians",
  "median_IgG_in_B", "feature", "B_cell_marker_medians",
  "median_CD38_in_B", "feature", "B_cell_marker_medians",
  "median_CD138_in_B", "feature", "B_cell_marker_medians",
  "median_CD24_in_B", "feature", "B_cell_marker_medians",
  "median_CD10_in_B", "feature", "B_cell_marker_medians",
  "median_CD39_in_B", "feature", "B_cell_marker_medians"
)

composite_dictionary <- tibble::tibble(
  variable = composite_names,
  variable_type = "composite",
  module = "composite_scores"
)

variable_dictionary <- dplyr::bind_rows(feature_dictionary, composite_dictionary) %>%
  dplyr::filter(variable %in% names(threshold_analysis_data_with_scores))

target_variables <- variable_dictionary$variable

# ------------------------------------------------------------
# 11. Modeling
# ------------------------------------------------------------

run_lm_variable <- function(df, variable_name, threshold_set_name) {

  model_df <- df %>%
    dplyr::filter(threshold_set == threshold_set_name) %>%
    dplyr::transmute(
      value = suppressWarnings(as.numeric(.data[[variable_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = sex
    ) %>%
    dplyr::filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex)
    ) %>%
    dplyr::mutate(
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex = droplevels(factor(as.character(sex)))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")

  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) {
    return(NULL)
  }

  fit <- tryCatch(
    stats::lm(value ~ disease_group + age_for_model + sex, data = model_df),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  tt <- tryCatch(
    broom::tidy(fit, conf.int = TRUE),
    error = function(e) NULL
  )

  if (is.null(tt)) return(NULL)

  disease_term <- "disease_groupCancer patient"
  if (!(disease_term %in% tt$term)) return(NULL)

  out <- tt %>% dplyr::filter(term == disease_term)

  desc <- model_df %>%
    dplyr::group_by(disease_group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      .groups = "drop"
    )

  healthy_mean <- desc$mean[desc$disease_group == "Healthy control"]
  cancer_mean <- desc$mean[desc$disease_group == "Cancer patient"]
  healthy_median <- desc$median[desc$disease_group == "Healthy control"]
  cancer_median <- desc$median[desc$disease_group == "Cancer patient"]

  if (length(healthy_mean) == 0) healthy_mean <- NA_real_
  if (length(cancer_mean) == 0) cancer_mean <- NA_real_
  if (length(healthy_median) == 0) healthy_median <- NA_real_
  if (length(cancer_median) == 0) cancer_median <- NA_real_

  tibble::tibble(
    threshold_set = threshold_set_name,
    variable = variable_name,
    n_model = n_model,
    n_healthy = n_healthy,
    n_cancer = n_cancer,
    healthy_mean = healthy_mean,
    cancer_mean = cancer_mean,
    healthy_median = healthy_median,
    cancer_median = cancer_median,
    beta_cancer_vs_healthy = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = dplyr::case_when(
      out$estimate[1] > 0 ~ "higher_in_cancer",
      out$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

threshold_model_results <- dplyr::bind_rows(
  lapply(threshold_configs, function(cfg) {
    dplyr::bind_rows(
      lapply(target_variables, function(vv) {
        run_lm_variable(threshold_analysis_data_with_scores, vv, cfg)
      })
    )
  })
) %>%
  dplyr::left_join(variable_dictionary, by = "variable") %>%
  dplyr::group_by(threshold_set, module) %>%
  dplyr::mutate(FDR_within_module = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(threshold_set) %>%
  dplyr::mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(threshold_set, FDR_global, p_value)

# ------------------------------------------------------------
# 12. Robustness summary
# ------------------------------------------------------------

get_one <- function(df, cfg, col) {
  x <- df[df$threshold_set == cfg, col, drop = TRUE]
  if (length(x) == 0) return(NA)
  x[1]
}

make_robustness_one <- function(vv) {

  rr <- threshold_model_results %>%
    dplyr::filter(variable == vv)

  if (nrow(rr) == 0) return(NULL)

  main_direction <- get_one(rr, "main", "direction")
  permissive_direction <- get_one(rr, "permissive", "direction")
  stringent_direction <- get_one(rr, "stringent", "direction")

  main_FDR_global <- as.numeric(get_one(rr, "main", "FDR_global"))
  permissive_FDR_global <- as.numeric(get_one(rr, "permissive", "FDR_global"))
  stringent_FDR_global <- as.numeric(get_one(rr, "stringent", "FDR_global"))

  main_FDR_module <- as.numeric(get_one(rr, "main", "FDR_within_module"))
  permissive_FDR_module <- as.numeric(get_one(rr, "permissive", "FDR_within_module"))
  stringent_FDR_module <- as.numeric(get_one(rr, "stringent", "FDR_within_module"))

  directions <- c(main_direction, permissive_direction, stringent_direction)
  direction_preserved <- all(!is.na(directions)) && length(unique(directions)) == 1

  tibble::tibble(
    variable = vv,
    variable_type = get_one(rr, "main", "variable_type"),
    module = get_one(rr, "main", "module"),

    main_beta = as.numeric(get_one(rr, "main", "beta_cancer_vs_healthy")),
    permissive_beta = as.numeric(get_one(rr, "permissive", "beta_cancer_vs_healthy")),
    stringent_beta = as.numeric(get_one(rr, "stringent", "beta_cancer_vs_healthy")),

    main_direction = main_direction,
    permissive_direction = permissive_direction,
    stringent_direction = stringent_direction,

    main_FDR_global = main_FDR_global,
    permissive_FDR_global = permissive_FDR_global,
    stringent_FDR_global = stringent_FDR_global,

    main_FDR_within_module = main_FDR_module,
    permissive_FDR_within_module = permissive_FDR_module,
    stringent_FDR_within_module = stringent_FDR_module,

    direction_preserved_all_sets = direction_preserved,
    global_FDR_lt_0p05_all_sets =
      direction_preserved &&
      main_FDR_global < 0.05 &&
      permissive_FDR_global < 0.05 &&
      stringent_FDR_global < 0.05,
    module_FDR_lt_0p05_all_sets =
      direction_preserved &&
      main_FDR_module < 0.05 &&
      permissive_FDR_module < 0.05 &&
      stringent_FDR_module < 0.05
  )
}

threshold_robustness_summary <- dplyr::bind_rows(
  lapply(target_variables, make_robustness_one)
) %>%
  dplyr::mutate(
    robustness_class = dplyr::case_when(
      global_FDR_lt_0p05_all_sets == TRUE ~ "direction_and_global_FDR_preserved_all_sets",
      module_FDR_lt_0p05_all_sets == TRUE ~ "direction_and_module_FDR_preserved_all_sets",
      direction_preserved_all_sets == TRUE ~ "direction_preserved_FDR_not_all_sets",
      TRUE ~ "direction_not_preserved"
    )
  ) %>%
  dplyr::arrange(
    factor(
      robustness_class,
      levels = c(
        "direction_and_global_FDR_preserved_all_sets",
        "direction_and_module_FDR_preserved_all_sets",
        "direction_preserved_FDR_not_all_sets",
        "direction_not_preserved"
      )
    ),
    main_FDR_global
  )

threshold_robustness_counts <- threshold_robustness_summary %>%
  dplyr::count(variable_type, robustness_class, name = "n") %>%
  dplyr::arrange(variable_type, robustness_class)

threshold_extraction_summary <- threshold_analysis_data_with_scores %>%
  dplyr::group_by(threshold_set) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_unique_subjects = dplyr::n_distinct(subject_id),
    n_feature_ok = sum(feature_ok == TRUE, na.rm = TRUE),
    n_feature_failed = sum(feature_ok == FALSE, na.rm = TRUE),
    median_total_events = median(total_events, na.rm = TRUE),
    median_dump_low_events = median(n_dump_low, na.rm = TRUE),
    median_cd19_b_events = median(n_cd19_b, na.rm = TRUE),
    min_cd19_b_events = min(n_cd19_b, na.rm = TRUE),
    max_cd19_b_events = max(n_cd19_b, na.rm = TRUE),
    median_pct_cd19_b_total = median(pct_cd19_b_total, na.rm = TRUE),
    median_pct_cd19_b_within_dump_low = median(pct_cd19_b_within_dump_low, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 13. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  threshold_extraction_summary,
  file.path(out_dir, "SDY2583_CP22_threshold_extraction_summary_STEP5.csv")
)

readr::write_csv(
  threshold_analysis_data_with_scores,
  file.path(out_dir, "SDY2583_CP22_threshold_analysis_data_with_scores_STEP5.csv")
)

readr::write_csv(
  variable_dictionary,
  file.path(out_dir, "SDY2583_CP22_threshold_variable_dictionary_STEP5.csv")
)

readr::write_csv(
  threshold_model_results,
  file.path(out_dir, "SDY2583_CP22_threshold_model_results_STEP5.csv")
)

readr::write_csv(
  threshold_robustness_summary,
  file.path(out_dir, "SDY2583_CP22_threshold_robustness_summary_STEP5.csv")
)

readr::write_csv(
  threshold_robustness_counts,
  file.path(out_dir, "SDY2583_CP22_threshold_robustness_counts_STEP5.csv")
)

save(
  threshold_configs,
  threshold_config_table,
  threshold_feature_table,
  threshold_analysis_data_with_scores,
  variable_dictionary,
  threshold_model_results,
  threshold_robustness_summary,
  threshold_robustness_counts,
  threshold_extraction_summary,
  composite_definitions,
  file = file.path(
    rdata_dir,
    "SDY2583_CP22_STEP5_threshold_sensitivity_TARGETED.RData"
  )
)

# ------------------------------------------------------------
# 14. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 5 COMPLETE: THRESHOLD SENSITIVITY\n")
cat("============================================================\n")

cat("\nThreshold extraction summary:\n")
print(as.data.frame(threshold_extraction_summary), row.names = FALSE)

cat("\nThreshold robustness counts:\n")
print(as.data.frame(threshold_robustness_counts), row.names = FALSE)

cat("\nTop threshold-robust variables:\n")
print(
  as.data.frame(
    threshold_robustness_summary %>%
      dplyr::filter(robustness_class == "direction_and_global_FDR_preserved_all_sets") %>%
      dplyr::slice_head(n = 40)
  ),
  row.names = FALSE
)

cat("\nDirection preserved but FDR not all sets:\n")
print(
  as.data.frame(
    threshold_robustness_summary %>%
      dplyr::filter(robustness_class == "direction_preserved_FDR_not_all_sets") %>%
      dplyr::slice_head(n = 40)
  ),
  row.names = FALSE
)

cat("\nDirection not preserved:\n")
print(
  as.data.frame(
    threshold_robustness_summary %>%
      dplyr::filter(robustness_class == "direction_not_preserved")
  ),
  row.names = FALSE
)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
