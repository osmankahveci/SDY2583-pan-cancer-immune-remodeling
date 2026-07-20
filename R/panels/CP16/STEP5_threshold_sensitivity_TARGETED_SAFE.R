# ============================================================
# SDY2583 CP16
# STEP 5 SAFE: Targeted threshold sensitivity
#
# Purpose:
#   Re-extract CP16 APC/DC-like features under three threshold
#   configurations and test whether the primary cancer-associated
#   directions and FDR significance are preserved.
#
# Threshold sets:
#   main:
#     original Step 2 thresholds
#   permissive:
#     positive-marker thresholds -0.2 transformed units
#     dump-low threshold +0.2 transformed units
#   stringent:
#     positive-marker thresholds +0.2 transformed units
#     dump-low threshold -0.2 transformed units
#
# Variables tested:
#   1) all CP16 composite scores
#   2) all CP16 feature-dictionary features available after extraction
#
# Output:
#   outputs/CP16/07_threshold_sensitivity
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "broom", "purrr", "tidyr")

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
  library(purrr)
  library(tidyr)
})

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

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "07_threshold_sensitivity")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP1_fcs_inventory_marker_QC.RData")
step3a_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP3A_metadata_merge_age_QC.RData")
step3b_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP3B_age_sex_adjusted_statistics.RData")
step4_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP4_composite_scores.RData")

required_rdata <- c(step1_rdata, step3a_rdata, step3b_rdata, step4_rdata)

missing_rdata <- required_rdata[!file.exists(required_rdata)]
if (length(missing_rdata) > 0) {
  stop("Eksik RData dosyaları:\n", paste(missing_rdata, collapse = "\n"))
}

load(step1_rdata)
load(step3a_rdata)
load(step3b_rdata)
load(step4_rdata)

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "07_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("fcs_files")) stop("fcs_files bulunamadı.")
if (!exists("cp16_analysis_data")) stop("cp16_analysis_data bulunamadı.")
if (!exists("feature_dictionary_full")) stop("feature_dictionary_full bulunamadı.")
if (!exists("composite_definitions")) stop("composite_definitions bulunamadı.")
if (!exists("composite_score_cols")) stop("composite_score_cols bulunamadı.")

main_thresholds <- list(
  DUMP_LOW = 1.5,
  CD45 = 2.0,
  HLA_DR = 1.5,
  CD11c = 1.5,
  CD14 = 1.5,
  CD16 = 1.5,
  CD1c = 1.5,
  CD141 = 1.5,
  CLEC9A = 1.5,
  CD123 = 1.5,
  FceRI = 1.5,
  CD13 = 1.5
)

make_threshold_set <- function(kind) {
  th <- main_thresholds
  positive_markers <- setdiff(names(th), "DUMP_LOW")
  if (kind == "main") return(th)
  if (kind == "permissive") {
    for (mm in positive_markers) th[[mm]] <- th[[mm]] - 0.2
    th$DUMP_LOW <- th$DUMP_LOW + 0.2
    return(th)
  }
  if (kind == "stringent") {
    for (mm in positive_markers) th[[mm]] <- th[[mm]] + 0.2
    th$DUMP_LOW <- th$DUMP_LOW - 0.2
    return(th)
  }
  stop("Unknown threshold kind: ", kind)
}

threshold_sets <- list(
  main = make_threshold_set("main"),
  permissive = make_threshold_set("permissive"),
  stringent = make_threshold_set("stringent")
)

threshold_table <- bind_rows(lapply(names(threshold_sets), function(ss) {
  tibble(
    threshold_set = ss,
    marker_or_gate = names(threshold_sets[[ss]]),
    threshold = unlist(threshold_sets[[ss]])
  )
}))

write_csv(threshold_table, file.path(out_dir, "SDY2583_CP16_threshold_sets_STEP5.csv"))

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
  tibble(parameter_index = seq_along(channels), channel_name = channels, marker_desc = markers)
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
  ch$DUMP <- find_marker_channel_regex(map, "Viability.*CD15.*CD3.*CD19.*CCR3.*CD7|Viability|Dump")
  ch$HLA_DR <- find_marker_channel_regex(map, "^HLA-DR$|HLA_DR|HLADR")
  ch$CD11c <- find_marker_channel_exact(map, "CD11c")
  ch$CD14 <- find_marker_channel_exact(map, "CD14")
  ch$CD16 <- find_marker_channel_exact(map, "CD16")
  ch$CD1c <- find_marker_channel_exact(map, "CD1c")
  ch$CD141 <- find_marker_channel_exact(map, "CD141")
  ch$CLEC9A <- find_marker_channel_exact(map, "CLEC9A")
  ch$CD123 <- find_marker_channel_exact(map, "CD123")
  ch$FceRI <- find_marker_channel_regex(map, "FceRI|FcER1|FcERI|FCER1")
  ch$CD13 <- find_marker_channel_exact(map, "CD13")
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

safe_sum <- function(x) sum(x, na.rm = TRUE)
safe_pct <- function(num, den) ifelse(is.na(den) | den <= 0, NA_real_, 100 * num / den)
safe_ratio <- function(num, den) ifelse(is.na(den) | den <= 0, NA_real_, num / den)
safe_median <- function(x) {
  if (length(x) == 0) return(NA_real_)
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}

extract_one_cp16_threshold <- function(fp, threshold_set_name, th) {
  tryCatch({
    ff_raw <- flowCore::read.FCS(fp, transformation = FALSE, truncate_max_range = FALSE)
    marker_map_file <- get_marker_map(ff_raw)
    ch <- build_channel_map(marker_map_file)
    required <- c("DUMP", "CD45", "HLA_DR", "CD11c", "CD14", "CD16", "CD1c", "CD141", "CLEC9A", "CD123", "FceRI", "CD13")
    missing <- required[vapply(ch[required], function(x) is.na(x), logical(1))]
    if (length(missing) > 0) stop("Missing required markers: ", paste(missing, collapse = ", "))
    ff <- ff_raw %>% apply_compensation_safely() %>% transform_safely()
    ex <- flowCore::exprs(ff)
    total_events <- nrow(ex)
    dump <- get_vec(ex, ch$DUMP)
    cd45 <- get_vec(ex, ch$CD45)
    hla_dr <- get_vec(ex, ch$HLA_DR)
    cd11c <- get_vec(ex, ch$CD11c)
    cd14 <- get_vec(ex, ch$CD14)
    cd16 <- get_vec(ex, ch$CD16)
    cd1c <- get_vec(ex, ch$CD1c)
    cd141 <- get_vec(ex, ch$CD141)
    clec9a <- get_vec(ex, ch$CLEC9A)
    cd123 <- get_vec(ex, ch$CD123)
    fceri <- get_vec(ex, ch$FceRI)
    cd13 <- get_vec(ex, ch$CD13)
    dump_low <- dump <= th$DUMP_LOW
    primary <- dump_low & cd45 > th$CD45
    hladr_core <- primary & hla_dr > th$HLA_DR
    cd11c_hladr_apc_like <- hladr_core & cd11c > th$CD11c
    cd14_mono_like <- primary & cd14 > th$CD14
    cd14_hladr_pos <- cd14_mono_like & hla_dr > th$HLA_DR
    cd14_hladr_low <- cd14_mono_like & hla_dr <= th$HLA_DR
    cd14pos_cd16neg_mono_like <- cd14_mono_like & cd16 <= th$CD16
    cd14pos_cd16pos_mono_like <- cd14_mono_like & cd16 > th$CD16
    cd14lowneg_cd16pos_hladr_like <- primary & cd14 <= th$CD14 & cd16 > th$CD16 & hla_dr > th$HLA_DR
    cd1c_apc_like <- hladr_core & cd11c > th$CD11c & cd1c > th$CD1c
    cd1c_cd14low_cdc2_like <- cd1c_apc_like & cd14 <= th$CD14
    cd1c_fceri_apc_like <- cd1c_apc_like & fceri > th$FceRI
    fceri_pos_hladr <- hladr_core & fceri > th$FceRI
    cd141_pos_hladr <- hladr_core & cd141 > th$CD141
    clec9a_pos_hladr <- hladr_core & clec9a > th$CLEC9A
    cd141_clec9a_cdc1_like <- hladr_core & cd11c > th$CD11c & cd141 > th$CD141 & clec9a > th$CLEC9A & cd14 <= th$CD14
    cd123_hladr_dc_like <- hladr_core & cd123 > th$CD123
    cd123_pdc_like <- hladr_core & cd123 > th$CD123 & cd11c <= th$CD11c & cd14 <= th$CD14 & cd16 <= th$CD16
    cd123_cd11c_mixed_apc_like <- hladr_core & cd123 > th$CD123 & cd11c > th$CD11c
    cd13_pos_primary <- primary & cd13 > th$CD13
    cd13_pos_hladr <- hladr_core & cd13 > th$CD13
    cd13_cd11c_hladr_apc_like <- hladr_core & cd13 > th$CD13 & cd11c > th$CD11c
    n_dump_low <- safe_sum(dump_low)
    n_primary <- safe_sum(primary)
    n_hladr_core <- safe_sum(hladr_core)
    n_cd11c_hladr <- safe_sum(cd11c_hladr_apc_like)
    n_cd14_mono <- safe_sum(cd14_mono_like)
    n_cd14_hladr_pos <- safe_sum(cd14_hladr_pos)
    n_cd14_hladr_low <- safe_sum(cd14_hladr_low)
    n_cd14pos_cd16neg <- safe_sum(cd14pos_cd16neg_mono_like)
    n_cd14pos_cd16pos <- safe_sum(cd14pos_cd16pos_mono_like)
    n_cd14lowneg_cd16pos <- safe_sum(cd14lowneg_cd16pos_hladr_like)
    n_cd1c_apc <- safe_sum(cd1c_apc_like)
    n_cd1c_cd14low_cdc2 <- safe_sum(cd1c_cd14low_cdc2_like)
    n_cd1c_fceri_apc <- safe_sum(cd1c_fceri_apc_like)
    n_fceri_pos_hladr <- safe_sum(fceri_pos_hladr)
    n_cd141_pos_hladr <- safe_sum(cd141_pos_hladr)
    n_clec9a_pos_hladr <- safe_sum(clec9a_pos_hladr)
    n_cd141_clec9a_cdc1 <- safe_sum(cd141_clec9a_cdc1_like)
    n_cd123_hladr <- safe_sum(cd123_hladr_dc_like)
    n_cd123_pdc <- safe_sum(cd123_pdc_like)
    n_cd123_cd11c_mixed <- safe_sum(cd123_cd11c_mixed_apc_like)
    n_cd13_primary <- safe_sum(cd13_pos_primary)
    n_cd13_hladr <- safe_sum(cd13_pos_hladr)
    n_cd13_cd11c_hladr <- safe_sum(cd13_cd11c_hladr_apc_like)
    tibble(
      threshold_set = threshold_set_name,
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      feature_ok = TRUE,
      feature_error = NA_character_,
      total_events = total_events,
      n_dump_low = n_dump_low,
      n_cd45_dump_low = n_primary,
      n_hladr_apc_core = n_hladr_core,
      n_cd11c_hladr_apc_like = n_cd11c_hladr,
      pct_dump_low_within_total = safe_pct(n_dump_low, total_events),
      pct_cd45_dump_low_within_total = safe_pct(n_primary, total_events),
      pct_hladr_apc_core_within_total = safe_pct(n_hladr_core, total_events),
      pct_hladr_apc_core_within_cd45_dump_low = safe_pct(n_hladr_core, n_primary),
      pct_cd11c_hladr_apc_like_within_cd45_dump_low = safe_pct(n_cd11c_hladr, n_primary),
      pct_cd11c_hladr_apc_like_within_hladr_core = safe_pct(n_cd11c_hladr, n_hladr_core),
      pct_cd14_mono_like_within_cd45_dump_low = safe_pct(n_cd14_mono, n_primary),
      pct_cd14_hladr_pos_within_cd45_dump_low = safe_pct(n_cd14_hladr_pos, n_primary),
      pct_cd14_hladr_low_within_cd45_dump_low = safe_pct(n_cd14_hladr_low, n_primary),
      pct_hladr_pos_within_cd14_mono_like = safe_pct(n_cd14_hladr_pos, n_cd14_mono),
      pct_hladr_low_within_cd14_mono_like = safe_pct(n_cd14_hladr_low, n_cd14_mono),
      pct_cd14pos_cd16neg_mono_like_within_cd45_dump_low = safe_pct(n_cd14pos_cd16neg, n_primary),
      pct_cd14pos_cd16pos_mono_like_within_cd45_dump_low = safe_pct(n_cd14pos_cd16pos, n_primary),
      pct_cd14lowneg_cd16pos_hladr_like_within_cd45_dump_low = safe_pct(n_cd14lowneg_cd16pos, n_primary),
      pct_cd16_pos_within_cd14_mono_like = safe_pct(n_cd14pos_cd16pos, n_cd14_mono),
      pct_cd1c_apc_like_within_cd45_dump_low = safe_pct(n_cd1c_apc, n_primary),
      pct_cd1c_apc_like_within_hladr_core = safe_pct(n_cd1c_apc, n_hladr_core),
      pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low = safe_pct(n_cd1c_cd14low_cdc2, n_primary),
      pct_cd1c_cd14low_cdc2_like_within_hladr_core = safe_pct(n_cd1c_cd14low_cdc2, n_hladr_core),
      pct_cd1c_fceri_apc_like_within_cd45_dump_low = safe_pct(n_cd1c_fceri_apc, n_primary),
      pct_cd1c_fceri_apc_like_within_cd1c_apc_like = safe_pct(n_cd1c_fceri_apc, n_cd1c_apc),
      pct_fceri_pos_hladr_within_cd45_dump_low = safe_pct(n_fceri_pos_hladr, n_primary),
      pct_cd141_pos_hladr_within_cd45_dump_low = safe_pct(n_cd141_pos_hladr, n_primary),
      pct_clec9a_pos_hladr_within_cd45_dump_low = safe_pct(n_clec9a_pos_hladr, n_primary),
      pct_cd141_clec9a_cdc1_like_within_cd45_dump_low = safe_pct(n_cd141_clec9a_cdc1, n_primary),
      pct_cd141_clec9a_cdc1_like_within_hladr_core = safe_pct(n_cd141_clec9a_cdc1, n_hladr_core),
      pct_cd123_hladr_dc_like_within_cd45_dump_low = safe_pct(n_cd123_hladr, n_primary),
      pct_cd123_hladr_dc_like_within_hladr_core = safe_pct(n_cd123_hladr, n_hladr_core),
      pct_cd123_pdc_like_within_cd45_dump_low = safe_pct(n_cd123_pdc, n_primary),
      pct_cd123_pdc_like_within_hladr_core = safe_pct(n_cd123_pdc, n_hladr_core),
      pct_cd123_cd11c_mixed_apc_like_within_cd45_dump_low = safe_pct(n_cd123_cd11c_mixed, n_primary),
      pct_cd13_pos_within_cd45_dump_low = safe_pct(n_cd13_primary, n_primary),
      pct_cd13_pos_within_hladr_core = safe_pct(n_cd13_hladr, n_hladr_core),
      pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low = safe_pct(n_cd13_cd11c_hladr, n_primary),
      ratio_cd11c_hladr_apc_to_cd14_mono_like = safe_ratio(n_cd11c_hladr, n_cd14_mono),
      ratio_cd1c_cdc2_like_to_cd123_pdc_like = safe_ratio(n_cd1c_cd14low_cdc2, n_cd123_pdc),
      ratio_cd141_clec9a_cdc1_like_to_cd123_pdc_like = safe_ratio(n_cd141_clec9a_cdc1, n_cd123_pdc),
      ratio_cd11c_hladr_apc_to_cd123_pdc_like = safe_ratio(n_cd11c_hladr, n_cd123_pdc),
      ratio_cd14_mono_like_to_cd123_pdc_like = safe_ratio(n_cd14_mono, n_cd123_pdc),
      median_HLA_DR_in_cd45_dump_low = safe_median(hla_dr[primary]),
      median_CD11c_in_cd45_dump_low = safe_median(cd11c[primary]),
      median_CD14_in_cd45_dump_low = safe_median(cd14[primary]),
      median_CD16_in_cd45_dump_low = safe_median(cd16[primary]),
      median_CD1c_in_cd45_dump_low = safe_median(cd1c[primary]),
      median_CD123_in_cd45_dump_low = safe_median(cd123[primary]),
      median_CD141_in_cd45_dump_low = safe_median(cd141[primary]),
      median_CLEC9A_in_cd45_dump_low = safe_median(clec9a[primary]),
      median_FceRI_in_cd45_dump_low = safe_median(fceri[primary]),
      median_CD13_in_cd45_dump_low = safe_median(cd13[primary]),
      median_HLA_DR_in_cd14_mono_like = safe_median(hla_dr[cd14_mono_like]),
      median_CD16_in_cd14_mono_like = safe_median(cd16[cd14_mono_like]),
      median_CD1c_in_hladr_core = safe_median(cd1c[hladr_core]),
      median_FceRI_in_cd1c_apc_like = safe_median(fceri[cd1c_apc_like]),
      median_CLEC9A_in_cd141_pos_hladr = safe_median(clec9a[cd141_pos_hladr]),
      median_CD123_in_hladr_core = safe_median(cd123[hladr_core]),
      median_CD13_in_hladr_core = safe_median(cd13[hladr_core])
    )
  }, error = function(e) {
    tibble(
      threshold_set = threshold_set_name,
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      feature_ok = FALSE,
      feature_error = as.character(e$message)
    )
  })
}

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 5 STARTED: TARGETED THRESHOLD SENSITIVITY\n")
cat("============================================================\n")

all_threshold_features <- list()

for (ss in names(threshold_sets)) {
  cat("\nProcessing threshold set: ", ss, "\n", sep = "")
  th <- threshold_sets[[ss]]
  rows <- vector("list", length(fcs_files))
  for (i in seq_along(fcs_files)) {
    if (i %% 50 == 0) cat("  ", ss, ": processed ", i, " of ", length(fcs_files), " FCS files\n", sep = "")
    rows[[i]] <- extract_one_cp16_threshold(fcs_files[i], ss, th)
  }
  all_threshold_features[[ss]] <- bind_rows(rows)
}

threshold_features <- bind_rows(all_threshold_features)

metadata_keep <- cp16_analysis_data %>%
  select(subject_id, disease_group, age_for_model, sex, sex_binary, channel_order_mismatch_file, marker_order_mismatch_file, marker_set_mismatch_file) %>%
  distinct(subject_id, .keep_all = TRUE)

threshold_features <- threshold_features %>%
  left_join(metadata_keep, by = "subject_id") %>%
  mutate(
    model_ready_age_sex = feature_ok == TRUE & !is.na(disease_group) & !is.na(age_for_model) & !is.na(sex)
  )

zscore_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

build_scores_one_threshold <- function(df_one) {
  df_out <- df_one
  for (sc in composite_score_cols) {
    defs <- composite_definitions %>% filter(composite_score == sc, feature %in% names(df_out))
    if (nrow(defs) == 0) {
      df_out[[sc]] <- NA_real_
      df_out[[paste0(sc, "_n_components_available")]] <- 0
      next
    }
    component_mat <- matrix(NA_real_, nrow = nrow(df_out), ncol = nrow(defs))
    colnames(component_mat) <- defs$feature
    for (j in seq_len(nrow(defs))) {
      ff <- defs$feature[j]
      dd <- defs$component_direction[j]
      component_mat[, j] <- zscore_safe(df_out[[ff]]) * dd
    }
    n_available_components <- rowSums(!is.na(component_mat))
    min_required <- ceiling(nrow(defs) / 2)
    score_value <- rowMeans(component_mat, na.rm = TRUE)
    score_value[n_available_components < min_required] <- NA_real_
    df_out[[sc]] <- score_value
    df_out[[paste0(sc, "_n_components_available")]] <- n_available_components
  }
  df_out
}

threshold_scores_data <- threshold_features %>%
  group_split(threshold_set) %>%
  lapply(build_scores_one_threshold) %>%
  bind_rows()

target_feature_cols <- feature_dictionary_full$feature
target_feature_cols <- target_feature_cols[target_feature_cols %in% names(threshold_scores_data)]
target_variables <- unique(c(composite_score_cols, target_feature_cols))

variable_dictionary <- bind_rows(
  tibble(variable = composite_score_cols, variable_type = "composite", module = "composite", interpretation = composite_score_cols),
  feature_dictionary_full %>% transmute(variable = feature, variable_type = "feature", module = module, interpretation = interpretation)
) %>%
  filter(variable %in% target_variables) %>%
  distinct(variable, .keep_all = TRUE)

run_variable_model <- function(df, variable_name, threshold_set_name) {
  model_df <- df %>%
    filter(threshold_set == threshold_set_name) %>%
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
    summarise(n = n(), mean = mean(value, na.rm = TRUE), median = median(value, na.rm = TRUE), .groups = "drop")
  tibble(
    threshold_set = threshold_set_name,
    variable = variable_name,
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

threshold_model_results <- bind_rows(lapply(names(threshold_sets), function(ss) {
  res <- bind_rows(lapply(target_variables, function(vv) run_variable_model(threshold_scores_data, vv, ss)))
  res %>%
    left_join(variable_dictionary, by = "variable") %>%
    group_by(variable_type) %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    ungroup()
})) %>%
  arrange(variable_type, variable, threshold_set)

wide_results <- threshold_model_results %>%
  select(variable, variable_type, module, interpretation, threshold_set, beta_cancer_vs_healthy, direction, p_value, FDR_global) %>%
  pivot_wider(
    names_from = threshold_set,
    values_from = c(beta_cancer_vs_healthy, direction, p_value, FDR_global),
    names_sep = "__"
  )

threshold_robustness_summary <- wide_results %>%
  mutate(
    direction_preserved_all_sets =
      !is.na(direction__main) &
      !is.na(direction__permissive) &
      !is.na(direction__stringent) &
      direction__main == direction__permissive &
      direction__main == direction__stringent,
    global_FDR_lt_0p05_all_sets =
      !is.na(FDR_global__main) &
      !is.na(FDR_global__permissive) &
      !is.na(FDR_global__stringent) &
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
    factor(robustness_class, levels = c("direction_and_global_FDR_preserved_all_sets", "direction_preserved_FDR_not_all_sets", "direction_not_preserved")),
    FDR_global__main
  )

threshold_robustness_counts <- threshold_robustness_summary %>%
  count(variable_type, robustness_class, name = "n") %>%
  arrange(variable_type, robustness_class)

top_threshold_robust_variables <- threshold_robustness_summary %>%
  filter(robustness_class == "direction_and_global_FDR_preserved_all_sets") %>%
  arrange(FDR_global__main) %>%
  slice_head(n = 40)

direction_preserved_FDR_not_all_sets <- threshold_robustness_summary %>%
  filter(robustness_class == "direction_preserved_FDR_not_all_sets") %>%
  arrange(FDR_global__main)

direction_not_preserved_threshold <- threshold_robustness_summary %>%
  filter(robustness_class == "direction_not_preserved") %>%
  arrange(FDR_global__main)

selected_summary_features <- c(
  "pct_cd45_dump_low_within_total",
  "pct_hladr_apc_core_within_cd45_dump_low",
  "pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low",
  "pct_cd123_pdc_like_within_cd45_dump_low",
  "pct_cd141_clec9a_cdc1_like_within_cd45_dump_low",
  "pct_cd14_hladr_low_within_cd45_dump_low",
  "pct_hladr_low_within_cd14_mono_like"
)
selected_summary_features <- selected_summary_features[selected_summary_features %in% names(threshold_scores_data)]

threshold_extraction_summary <- threshold_scores_data %>%
  group_by(threshold_set) %>%
  summarise(
    n_rows = n(),
    n_unique_subjects = n_distinct(subject_id),
    n_feature_ok = sum(feature_ok == TRUE, na.rm = TRUE),
    n_feature_failed = sum(feature_ok == FALSE, na.rm = TRUE),
    median_total_events = median(total_events, na.rm = TRUE),
    median_cd45_dump_low_events = median(n_cd45_dump_low, na.rm = TRUE),
    median_hladr_apc_core_events = median(n_hladr_apc_core, na.rm = TRUE),
    across(all_of(selected_summary_features), ~ median(.x, na.rm = TRUE), .names = "median_{.col}"),
    .groups = "drop"
  ) %>%
  arrange(factor(threshold_set, levels = c("main", "permissive", "stringent")))

failed_threshold_files <- threshold_scores_data %>%
  filter(feature_ok == FALSE) %>%
  select(threshold_set, subject_id, file_name, file_path, feature_error) %>%
  arrange(threshold_set, subject_id)

write_csv(threshold_features, file.path(out_dir, "SDY2583_CP16_threshold_features_STEP5.csv"))
write_csv(threshold_scores_data, file.path(out_dir, "SDY2583_CP16_threshold_scores_data_STEP5.csv"))
write_csv(threshold_extraction_summary, file.path(out_dir, "SDY2583_CP16_threshold_extraction_summary_STEP5.csv"))
write_csv(threshold_model_results, file.path(out_dir, "SDY2583_CP16_threshold_model_results_STEP5.csv"))
write_csv(threshold_robustness_summary, file.path(out_dir, "SDY2583_CP16_threshold_robustness_summary_STEP5.csv"))
write_csv(threshold_robustness_counts, file.path(out_dir, "SDY2583_CP16_threshold_robustness_counts_STEP5.csv"))
write_csv(top_threshold_robust_variables, file.path(out_dir, "SDY2583_CP16_top_threshold_robust_variables_STEP5.csv"))
write_csv(direction_preserved_FDR_not_all_sets, file.path(out_dir, "SDY2583_CP16_direction_preserved_FDR_not_all_sets_STEP5.csv"))
write_csv(direction_not_preserved_threshold, file.path(out_dir, "SDY2583_CP16_direction_not_preserved_threshold_STEP5.csv"))
write_csv(failed_threshold_files, file.path(out_dir, "SDY2583_CP16_failed_threshold_files_STEP5.csv"))

save(
  threshold_sets,
  threshold_table,
  threshold_features,
  threshold_scores_data,
  threshold_extraction_summary,
  threshold_model_results,
  threshold_robustness_summary,
  threshold_robustness_counts,
  top_threshold_robust_variables,
  direction_preserved_FDR_not_all_sets,
  direction_not_preserved_threshold,
  failed_threshold_files,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP5_threshold_sensitivity_TARGETED.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 5 COMPLETE: TARGETED THRESHOLD SENSITIVITY\n")
cat("============================================================\n")

cat("\nThreshold extraction summary:\n")
print(as.data.frame(threshold_extraction_summary), row.names = FALSE)

cat("\nThreshold robustness counts:\n")
print(as.data.frame(threshold_robustness_counts), row.names = FALSE)

cat("\nTop threshold-robust variables:\n")
print(as.data.frame(top_threshold_robust_variables %>% slice_head(n = 40)), row.names = FALSE)

cat("\nDirection preserved but FDR not all sets:\n")
print(as.data.frame(direction_preserved_FDR_not_all_sets), row.names = FALSE)

cat("\nDirection not preserved:\n")
print(as.data.frame(direction_not_preserved_threshold), row.names = FALSE)

cat("\nFailed threshold files:\n")
print(as.data.frame(failed_threshold_files), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
