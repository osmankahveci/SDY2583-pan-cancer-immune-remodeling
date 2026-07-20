# ============================================================
# SDY2583 CP10
# STEP 5 SAFE: Targeted threshold sensitivity
#
# Threshold sets:
#   main
#   permissive: positive-marker thresholds -0.2; viability-low threshold +0.2
#   stringent:  positive-marker thresholds +0.2; viability-low threshold -0.2
#
# Re-extracts CP10 targeted features and rebuilds CP10 composite scores.
#
# Interpretation:
# - Phenotype-like terms only.
# - No definitive neutrophil/eosinophil/DC/MDSC claim.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "tidyr")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
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
  library(tidyr)
})

# Namespace safety
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
count <- dplyr::count
distinct <- dplyr::distinct
case_when <- dplyr::case_when
bind_rows <- dplyr::bind_rows
left_join <- dplyr::left_join
n_distinct <- dplyr::n_distinct
n <- dplyr::n

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP10")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "07_threshold_sensitivity")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP1_fcs_inventory_marker_QC.RData")
step3a_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP3A_metadata_merge_age_QC.RData")
step3b_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP3B_age_sex_adjusted_statistics.RData")
step4_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP4_composite_scores.RData")

if (!file.exists(step1_rdata)) stop("Step 1 RData bulunamadı: ", step1_rdata)
if (!file.exists(step3a_rdata)) stop("Step 3A RData bulunamadı: ", step3a_rdata)
if (!file.exists(step3b_rdata)) stop("Step 3B RData bulunamadı: ", step3b_rdata)
if (!file.exists(step4_rdata)) stop("Step 4 RData bulunamadı: ", step4_rdata)

load(step1_rdata)
load(step3a_rdata)
load(step3b_rdata)
load(step4_rdata)

if (!exists("fcs_files")) stop("fcs_files bulunamadı.")
if (!exists("cp10_analysis_data")) stop("cp10_analysis_data bulunamadı.")
if (!exists("composite_definitions")) stop("composite_definitions bulunamadı.")
if (!exists("primary_results")) stop("primary_results bulunamadı.")
if (!exists("primary_composite_results")) stop("primary_composite_results bulunamadı.")

# ------------------------------------------------------------
# 3. Threshold sets
# ------------------------------------------------------------

base_thresholds <- list(
  VIABILITY_LOW = 1.5,
  CD45 = 2.0,
  CD3 = 1.5,
  CD19 = 1.5,
  CD56 = 1.5,
  CD14 = 1.5,
  HLA_DR = 1.5,
  CD11c = 1.5,
  CD13 = 1.5,
  CD66b = 1.5,
  CCR3 = 1.5,
  CD123 = 1.5
)

make_threshold_set <- function(set_name) {
  th <- base_thresholds

  if (set_name == "main") {
    return(th)
  }

  if (set_name == "permissive") {
    # low gate becomes more permissive by increasing allowed low threshold
    th$VIABILITY_LOW <- base_thresholds$VIABILITY_LOW + 0.2

    # positive gates become more permissive by lowering threshold
    for (nm in setdiff(names(th), "VIABILITY_LOW")) {
      th[[nm]] <- base_thresholds[[nm]] - 0.2
    }
    return(th)
  }

  if (set_name == "stringent") {
    # low gate becomes more stringent by lowering allowed low threshold
    th$VIABILITY_LOW <- base_thresholds$VIABILITY_LOW - 0.2

    # positive gates become more stringent by raising threshold
    for (nm in setdiff(names(th), "VIABILITY_LOW")) {
      th[[nm]] <- base_thresholds[[nm]] + 0.2
    }
    return(th)
  }

  stop("Unknown threshold set: ", set_name)
}

threshold_sets <- c("main", "permissive", "stringent")

threshold_set_table <- bind_rows(lapply(threshold_sets, function(ss) {
  th <- make_threshold_set(ss)
  tibble(
    threshold_set = ss,
    marker_or_gate = names(th),
    threshold = unlist(th)
  )
}))

write_csv(threshold_set_table, file.path(out_dir, "SDY2583_CP10_threshold_set_definitions_STEP5.csv"))

# ------------------------------------------------------------
# 4. Helper functions for FCS extraction
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

  tibble(
    parameter_index = seq_along(channels),
    channel_name = channels,
    marker_desc = markers
  )
}

find_marker_channel_exact <- function(map, marker) {
  idx <- which(map$marker_desc == marker)
  if (length(idx) >= 1) return(map$channel_name[idx[1]])
  NA_character_
}

find_marker_channel_regex <- function(map, pattern) {
  idx <- which(grepl(pattern, map$marker_desc, ignore.case = TRUE))
  if (length(idx) >= 1) return(map$channel_name[idx[1]])
  NA_character_
}

build_channel_map <- function(map) {
  ch <- list()
  ch$CD45 <- find_marker_channel_exact(map, "CD45")
  ch$CD3 <- find_marker_channel_exact(map, "CD3")
  ch$CD19 <- find_marker_channel_exact(map, "CD19")
  ch$CD56 <- find_marker_channel_exact(map, "CD56")
  ch$CD14 <- find_marker_channel_exact(map, "CD14")
  ch$HLA_DR <- find_marker_channel_regex(map, "^HLA-DR$|HLA_DR|HLADR")
  ch$CD11c <- find_marker_channel_exact(map, "CD11c")
  ch$CD13 <- find_marker_channel_exact(map, "CD13")
  ch$CD66b <- find_marker_channel_exact(map, "CD66b")
  ch$CD123 <- find_marker_channel_exact(map, "CD123")
  ch$CCR3 <- find_marker_channel_regex(map, "CCR3|CD193")
  ch$VIABILITY <- find_marker_channel_regex(map, "^Viability$")
  if (is.na(ch$VIABILITY)) ch$VIABILITY <- find_marker_channel_regex(map, "Viability|Live|Dead")
  ch
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
  if (is.null(sp)) return(ff)
  tryCatch(flowCore::compensate(ff, sp), error = function(e) ff)
}

transform_safely <- function(ff) {
  ex_names <- colnames(flowCore::exprs(ff))
  fluor_ch <- ex_names[!grepl("FSC|SSC|Time", ex_names, ignore.case = TRUE)]
  if (length(fluor_ch) == 0) return(ff)
  trans <- flowCore::logicleTransform(
    transformationId = "fixed_logicle",
    w = 0.5,
    t = 262144,
    m = 4.5,
    a = 0
  )
  tryCatch({
    tf <- flowCore::transformList(fluor_ch, trans)
    flowCore::transform(ff, tf)
  }, error = function(e) ff)
}

get_vec <- function(ex, ch) {
  if (is.na(ch) || !(ch %in% colnames(ex))) return(rep(NA_real_, nrow(ex)))
  as.numeric(ex[, ch])
}

safe_pct <- function(num, den) {
  ifelse(is.na(den) | den <= 0, NA_real_, 100 * num / den)
}

safe_ratio <- function(num, den) {
  ifelse(is.na(den) | den <= 0, NA_real_, num / den)
}

safe_median <- function(x) {
  if (length(x) == 0) return(NA_real_)
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}

safe_sum <- function(x) sum(x, na.rm = TRUE)

# ------------------------------------------------------------
# 5. Extract targeted CP10 features from one FCS file
# ------------------------------------------------------------

extract_one_cp10_threshold <- function(fp, threshold_set_name) {

  th <- make_threshold_set(threshold_set_name)

  tryCatch({

    ff_raw <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE
    )

    marker_map_file <- get_marker_map(ff_raw)
    ch <- build_channel_map(marker_map_file)

    required <- c("VIABILITY", "CD45", "CD3", "CD19", "CD56", "CD14", "HLA_DR", "CD11c", "CD13", "CD66b", "CCR3", "CD123")
    missing <- required[vapply(ch[required], function(x) is.na(x), logical(1))]
    if (length(missing) > 0) stop("Missing required markers: ", paste(missing, collapse = ", "))

    ff <- ff_raw %>%
      apply_compensation_safely() %>%
      transform_safely()

    ex <- flowCore::exprs(ff)

    total_events <- nrow(ex)

    viability <- get_vec(ex, ch$VIABILITY)
    cd45 <- get_vec(ex, ch$CD45)
    cd3 <- get_vec(ex, ch$CD3)
    cd19 <- get_vec(ex, ch$CD19)
    cd56 <- get_vec(ex, ch$CD56)
    cd14 <- get_vec(ex, ch$CD14)
    hla_dr <- get_vec(ex, ch$HLA_DR)
    cd11c <- get_vec(ex, ch$CD11c)
    cd13 <- get_vec(ex, ch$CD13)
    cd66b <- get_vec(ex, ch$CD66b)
    ccr3 <- get_vec(ex, ch$CCR3)
    cd123 <- get_vec(ex, ch$CD123)

    viable <- viability <= th$VIABILITY_LOW
    primary <- viable & cd45 > th$CD45

    n_viable <- safe_sum(viable)
    n_cd45 <- safe_sum(primary)

    t_like <- primary & cd3 > th$CD3
    b_like <- primary & cd19 > th$CD19
    nk_like <- primary & cd3 <= th$CD3 & cd56 > th$CD56
    lymphoid_like_any <- primary & (cd3 > th$CD3 | cd19 > th$CD19 | cd56 > th$CD56)

    cd14_mono_like <- primary & cd14 > th$CD14
    cd14_hladr_pos <- cd14_mono_like & hla_dr > th$HLA_DR
    cd14_hladr_low <- cd14_mono_like & hla_dr <= th$HLA_DR
    cd14_cd11c_hladr_pos <- cd14_mono_like & cd11c > th$CD11c & hla_dr > th$HLA_DR

    cd11c_hladr_apc_like <- primary & cd11c > th$CD11c & hla_dr > th$HLA_DR

    cd66b_gran_like <- primary & cd66b > th$CD66b
    cd13_cd66b_gran_like <- primary & cd13 > th$CD13 & cd66b > th$CD66b
    cd13_pos_within_cd66b <- cd66b_gran_like & cd13 > th$CD13

    ccr3_pos <- primary & ccr3 > th$CCR3
    ccr3_cd66b_eosinophil_like <- primary & ccr3 > th$CCR3 & cd66b > th$CD66b
    ccr3_cd13_cd66b_gran_like <- primary & ccr3 > th$CCR3 & cd13 > th$CD13 & cd66b > th$CD66b

    cd123_pos <- primary & cd123 > th$CD123
    cd123_hladr_dc_like <- primary & cd123 > th$CD123 & hla_dr > th$HLA_DR
    cd123_hladr_cd11c_neg_pdc_like <- primary &
      cd123 > th$CD123 & hla_dr > th$HLA_DR &
      cd11c <= th$CD11c & cd14 <= th$CD14 &
      cd3 <= th$CD3 & cd19 <= th$CD19
    cd123_cd11c_hladr_mixed_apc_like <- primary & cd123 > th$CD123 & cd11c > th$CD11c & hla_dr > th$HLA_DR

    myeloid_granulocytic_like_any <- primary & (
      cd14 > th$CD14 | cd11c > th$CD11c | cd13 > th$CD13 |
        cd66b > th$CD66b | ccr3 > th$CCR3 | cd123 > th$CD123
    )

    n_t_like <- safe_sum(t_like)
    n_b_like <- safe_sum(b_like)
    n_lymphoid_like_any <- safe_sum(lymphoid_like_any)
    n_cd14_mono_like <- safe_sum(cd14_mono_like)
    n_cd14_hladr_pos <- safe_sum(cd14_hladr_pos)
    n_cd14_hladr_low <- safe_sum(cd14_hladr_low)
    n_cd14_cd11c_hladr_pos <- safe_sum(cd14_cd11c_hladr_pos)
    n_cd11c_hladr_apc_like <- safe_sum(cd11c_hladr_apc_like)
    n_cd66b_gran_like <- safe_sum(cd66b_gran_like)
    n_cd13_cd66b_gran_like <- safe_sum(cd13_cd66b_gran_like)
    n_ccr3_pos <- safe_sum(ccr3_pos)
    n_ccr3_cd66b_eosinophil_like <- safe_sum(ccr3_cd66b_eosinophil_like)
    n_ccr3_cd13_cd66b_gran_like <- safe_sum(ccr3_cd13_cd66b_gran_like)
    n_cd123_pos <- safe_sum(cd123_pos)
    n_cd123_hladr_dc_like <- safe_sum(cd123_hladr_dc_like)
    n_cd123_hladr_cd11c_neg_pdc_like <- safe_sum(cd123_hladr_cd11c_neg_pdc_like)
    n_cd123_cd11c_hladr_mixed_apc_like <- safe_sum(cd123_cd11c_hladr_mixed_apc_like)
    n_myeloid_granulocytic_like_any <- safe_sum(myeloid_granulocytic_like_any)

    tibble(
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      threshold_set = threshold_set_name,
      feature_ok = TRUE,
      feature_error = NA_character_,
      total_events = total_events,
      n_viable = n_viable,
      n_cd45_viable = n_cd45,

      pct_t_like_within_cd45 = safe_pct(n_t_like, n_cd45),
      pct_b_like_within_cd45 = safe_pct(n_b_like, n_cd45),
      pct_lymphoid_like_any_within_cd45 = safe_pct(n_lymphoid_like_any, n_cd45),

      pct_cd14_mono_like_within_cd45 = safe_pct(n_cd14_mono_like, n_cd45),
      pct_cd14_hladr_pos_within_cd45 = safe_pct(n_cd14_hladr_pos, n_cd45),
      pct_cd14_hladr_low_within_cd45 = safe_pct(n_cd14_hladr_low, n_cd45),
      pct_cd14_cd11c_hladr_pos_within_cd45 = safe_pct(n_cd14_cd11c_hladr_pos, n_cd45),
      pct_cd11c_hladr_apc_like_within_cd45 = safe_pct(n_cd11c_hladr_apc_like, n_cd45),

      pct_cd66b_gran_like_within_cd45 = safe_pct(n_cd66b_gran_like, n_cd45),
      pct_cd13_cd66b_gran_like_within_cd45 = safe_pct(n_cd13_cd66b_gran_like, n_cd45),
      pct_cd13_pos_within_cd66b_gran_like = safe_pct(n_cd13_cd66b_gran_like, n_cd66b_gran_like),

      pct_ccr3_pos_within_cd45 = safe_pct(n_ccr3_pos, n_cd45),
      pct_ccr3_cd66b_eosinophil_like_within_cd45 = safe_pct(n_ccr3_cd66b_eosinophil_like, n_cd45),
      pct_ccr3_cd13_cd66b_gran_like_within_cd45 = safe_pct(n_ccr3_cd13_cd66b_gran_like, n_cd45),
      pct_ccr3_pos_within_cd66b_gran_like = safe_pct(n_ccr3_cd66b_eosinophil_like, n_cd66b_gran_like),

      pct_cd123_pos_within_cd45 = safe_pct(n_cd123_pos, n_cd45),
      pct_cd123_hladr_dc_like_within_cd45 = safe_pct(n_cd123_hladr_dc_like, n_cd45),
      pct_cd123_hladr_cd11c_neg_pdc_like_within_cd45 = safe_pct(n_cd123_hladr_cd11c_neg_pdc_like, n_cd45),
      pct_cd123_cd11c_hladr_mixed_apc_like_within_cd45 = safe_pct(n_cd123_cd11c_hladr_mixed_apc_like, n_cd45),

      pct_myeloid_granulocytic_like_any_within_cd45 = safe_pct(n_myeloid_granulocytic_like_any, n_cd45),
      ratio_myeloid_granulocytic_to_lymphoid_like = safe_ratio(n_myeloid_granulocytic_like_any, n_lymphoid_like_any),
      ratio_cd66b_gran_like_to_t_like = safe_ratio(n_cd66b_gran_like, n_t_like),
      ratio_cd14_mono_like_to_t_like = safe_ratio(n_cd14_mono_like, n_t_like),

      median_CD13_in_CD45 = safe_median(cd13[primary]),
      median_CD66b_in_CD45 = safe_median(cd66b[primary]),
      median_CD14_in_CD45 = safe_median(cd14[primary]),
      median_CD11c_in_CD45 = safe_median(cd11c[primary]),
      median_CCR3_in_CD45 = safe_median(ccr3[primary]),

      median_CD13_in_CD66b_gran_like = safe_median(cd13[cd66b_gran_like]),
      median_CD66b_in_CD66b_gran_like = safe_median(cd66b[cd66b_gran_like]),
      median_CCR3_in_CD66b_gran_like = safe_median(ccr3[cd66b_gran_like])
    )

  }, error = function(e) {
    tibble(
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      threshold_set = threshold_set_name,
      feature_ok = FALSE,
      feature_error = as.character(e$message)
    )
  })
}

# ------------------------------------------------------------
# 6. Extract all threshold sets
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 5 STARTED: TARGETED THRESHOLD SENSITIVITY\n")
cat("============================================================\n")

threshold_feature_list <- list()

for (ss in threshold_sets) {

  cat("\nThreshold set:", ss, "\n")
  rows <- vector("list", length(fcs_files))

  for (i in seq_along(fcs_files)) {
    if (i %% 50 == 0) cat("Processed", i, "of", length(fcs_files), "files for", ss, "\n")
    rows[[i]] <- extract_one_cp10_threshold(fcs_files[i], ss)
  }

  threshold_feature_list[[ss]] <- bind_rows(rows)
}

threshold_features_long <- bind_rows(threshold_feature_list)

# ------------------------------------------------------------
# 7. Add metadata and rebuild composites for each threshold set
# ------------------------------------------------------------

metadata_cols <- cp10_analysis_data %>%
  select(
    subject_id,
    disease_group,
    age_for_model,
    sex,
    sex_binary,
    model_ready_age_sex,
    channel_order_mismatch_file,
    marker_mismatch_file
  ) %>%
  distinct(subject_id, .keep_all = TRUE)

threshold_features_long <- threshold_features_long %>%
  left_join(metadata_cols, by = "subject_id") %>%
  mutate(
    model_ready_threshold = feature_ok == TRUE &
      !is.na(disease_group) &
      !is.na(age_for_model) &
      !is.na(sex)
  )

# Only keep composite components present in threshold_features_long.
composite_definitions_threshold <- composite_definitions %>%
  filter(feature %in% names(threshold_features_long))

zscore_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

build_scores_for_set <- function(df_set, set_name) {

  out <- df_set

  for (score_name in unique(composite_definitions_threshold$composite_score)) {

    defs <- composite_definitions_threshold %>%
      filter(composite_score == score_name)

    component_mat <- matrix(NA_real_, nrow = nrow(out), ncol = nrow(defs))

    for (j in seq_len(nrow(defs))) {
      ff <- defs$feature[j]
      dd <- defs$component_direction[j]
      component_mat[, j] <- zscore_safe(out[[ff]]) * dd
    }

    out[[score_name]] <- rowMeans(component_mat, na.rm = TRUE)
  }

  out
}

threshold_scored_long <- bind_rows(lapply(threshold_sets, function(ss) {
  df_set <- threshold_features_long %>% filter(threshold_set == ss)
  build_scores_for_set(df_set, ss)
}))

# ------------------------------------------------------------
# 8. Target variables and modeling
# ------------------------------------------------------------

target_features <- unique(c(
  composite_definitions_threshold$feature,
  "pct_cd13_cd66b_gran_like_within_cd45",
  "pct_cd66b_gran_like_within_cd45",
  "pct_cd14_mono_like_within_cd45",
  "pct_myeloid_granulocytic_like_any_within_cd45",
  "ratio_myeloid_granulocytic_to_lymphoid_like",
  "ratio_cd66b_gran_like_to_t_like",
  "pct_lymphoid_like_any_within_cd45",
  "pct_t_like_within_cd45",
  "pct_cd11c_hladr_apc_like_within_cd45",
  "pct_ccr3_cd66b_eosinophil_like_within_cd45"
))

target_features <- target_features[target_features %in% names(threshold_scored_long)]
target_composites <- unique(composite_definitions_threshold$composite_score)

target_variables <- unique(c(target_composites, target_features))

run_lm_variable <- function(df, variable_name, threshold_set_name) {

  model_df <- df %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[variable_name]])),
      disease_group = disease_group,
      age_for_model = suppressWarnings(as.numeric(age_for_model)),
      sex = sex
    ) %>%
    filter(!is.na(value), !is.na(disease_group), !is.na(age_for_model), !is.na(sex)) %>%
    mutate(
      disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
      sex = droplevels(factor(as.character(sex)))
    )

  n_model <- nrow(model_df)
  n_healthy <- sum(model_df$disease_group == "Healthy control")
  n_cancer <- sum(model_df$disease_group == "Cancer patient")

  if (n_model < 50 || n_healthy < 20 || n_cancer < 20) return(NULL)

  fit <- tryCatch(lm(value ~ disease_group + age_for_model + sex, data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(NULL)

  term <- "disease_groupCancer patient"
  if (!(term %in% tt$term)) return(NULL)
  out <- tt %>% filter(term == !!term)

  desc <- model_df %>%
    group_by(disease_group) %>%
    summarise(
      n = n(),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      .groups = "drop"
    )

  tibble(
    threshold_set = threshold_set_name,
    variable = variable_name,
    variable_type = ifelse(variable_name %in% target_composites, "composite", "feature"),
    n_model = n_model,
    n_healthy = n_healthy,
    n_cancer = n_cancer,
    healthy_mean = desc$mean[desc$disease_group == "Healthy control"][1],
    cancer_mean = desc$mean[desc$disease_group == "Cancer patient"][1],
    healthy_median = desc$median[desc$disease_group == "Healthy control"][1],
    cancer_median = desc$median[desc$disease_group == "Cancer patient"][1],
    beta_cancer_vs_healthy = out$estimate[1],
    conf_low = out$conf.low[1],
    conf_high = out$conf.high[1],
    p_value = out$p.value[1],
    direction = case_when(
      out$estimate[1] > 0 ~ "higher_in_cancer",
      out$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

threshold_model_results <- bind_rows(lapply(threshold_sets, function(ss) {
  df_set <- threshold_scored_long %>%
    filter(threshold_set == ss, model_ready_threshold == TRUE)

  bind_rows(lapply(target_variables, function(vv) {
    run_lm_variable(df_set, vv, ss)
  })) %>%
    mutate(
      FDR_global = p.adjust(p_value, method = "BH")
    )
})) %>%
  arrange(variable, threshold_set)

# ------------------------------------------------------------
# 9. Robustness classification
# ------------------------------------------------------------

wide_results <- threshold_model_results %>%
  select(
    threshold_set,
    variable,
    variable_type,
    beta_cancer_vs_healthy,
    direction,
    p_value,
    FDR_global
  ) %>%
  pivot_wider(
    names_from = threshold_set,
    values_from = c(beta_cancer_vs_healthy, direction, p_value, FDR_global),
    names_sep = "__"
  )

threshold_robustness_summary <- wide_results %>%
  mutate(
    direction_preserved_all_sets =
      !is.na(direction__main) &
      direction__main == direction__permissive &
      direction__main == direction__stringent,
    global_FDR_lt_0p05_all_sets =
      FDR_global__main < 0.05 &
      FDR_global__permissive < 0.05 &
      FDR_global__stringent < 0.05,
    robustness_class = case_when(
      direction_preserved_all_sets == TRUE & global_FDR_lt_0p05_all_sets == TRUE ~ "direction_and_global_FDR_preserved_all_sets",
      direction_preserved_all_sets == TRUE & global_FDR_lt_0p05_all_sets == FALSE ~ "direction_preserved_FDR_not_all_sets",
      TRUE ~ "direction_not_preserved"
    )
  ) %>%
  arrange(
    variable_type,
    factor(
      robustness_class,
      levels = c(
        "direction_and_global_FDR_preserved_all_sets",
        "direction_preserved_FDR_not_all_sets",
        "direction_not_preserved"
      )
    ),
    FDR_global__main
  )

threshold_extraction_summary <- threshold_features_long %>%
  group_by(threshold_set) %>%
  summarise(
    n_rows = n(),
    n_unique_subjects = n_distinct(subject_id),
    n_feature_ok = sum(feature_ok == TRUE, na.rm = TRUE),
    n_feature_failed = sum(feature_ok == FALSE, na.rm = TRUE),
    median_total_events = median(total_events, na.rm = TRUE),
    median_viable_events = median(n_viable, na.rm = TRUE),
    median_cd45_viable_events = median(n_cd45_viable, na.rm = TRUE),
    median_pct_cd13_cd66b_gran_like_within_cd45 = median(pct_cd13_cd66b_gran_like_within_cd45, na.rm = TRUE),
    median_pct_cd66b_gran_like_within_cd45 = median(pct_cd66b_gran_like_within_cd45, na.rm = TRUE),
    median_pct_cd14_mono_like_within_cd45 = median(pct_cd14_mono_like_within_cd45, na.rm = TRUE),
    median_ratio_myeloid_granulocytic_to_lymphoid_like = median(ratio_myeloid_granulocytic_to_lymphoid_like, na.rm = TRUE),
    .groups = "drop"
  )

robustness_counts <- threshold_robustness_summary %>%
  count(variable_type, robustness_class, name = "n") %>%
  arrange(variable_type, robustness_class)

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

write_csv(threshold_features_long, file.path(out_dir, "SDY2583_CP10_threshold_features_long_STEP5.csv"))
write_csv(threshold_scored_long, file.path(out_dir, "SDY2583_CP10_threshold_scored_long_STEP5.csv"))
write_csv(threshold_model_results, file.path(out_dir, "SDY2583_CP10_threshold_model_results_STEP5.csv"))
write_csv(threshold_robustness_summary, file.path(out_dir, "SDY2583_CP10_threshold_robustness_summary_STEP5.csv"))
write_csv(threshold_extraction_summary, file.path(out_dir, "SDY2583_CP10_threshold_extraction_summary_STEP5.csv"))
write_csv(robustness_counts, file.path(out_dir, "SDY2583_CP10_threshold_robustness_counts_STEP5.csv"))
write_csv(composite_definitions_threshold, file.path(out_dir, "SDY2583_CP10_composite_definitions_threshold_components_STEP5.csv"))

save(
  threshold_features_long,
  threshold_scored_long,
  threshold_model_results,
  threshold_robustness_summary,
  threshold_extraction_summary,
  robustness_counts,
  threshold_set_table,
  composite_definitions_threshold,
  target_variables,
  target_composites,
  target_features,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP5_threshold_sensitivity_TARGETED.RData")
)

# ------------------------------------------------------------
# 11. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 5 COMPLETE: TARGETED THRESHOLD SENSITIVITY\n")
cat("============================================================\n")

cat("\nThreshold extraction summary:\n")
print(as.data.frame(threshold_extraction_summary), row.names = FALSE)

cat("\nThreshold robustness counts:\n")
print(as.data.frame(robustness_counts), row.names = FALSE)

cat("\nTop threshold-robust variables:\n")
print(
  as.data.frame(
    threshold_robustness_summary %>%
      filter(robustness_class == "direction_and_global_FDR_preserved_all_sets") %>%
      arrange(variable_type, FDR_global__main)
  ),
  row.names = FALSE
)

cat("\nDirection preserved but FDR not all sets:\n")
print(
  as.data.frame(
    threshold_robustness_summary %>%
      filter(robustness_class == "direction_preserved_FDR_not_all_sets") %>%
      arrange(variable_type, FDR_global__main)
  ),
  row.names = FALSE
)

cat("\nDirection not preserved:\n")
print(
  as.data.frame(
    threshold_robustness_summary %>%
      filter(robustness_class == "direction_not_preserved") %>%
      arrange(variable_type, FDR_global__main)
  ),
  row.names = FALSE
)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
