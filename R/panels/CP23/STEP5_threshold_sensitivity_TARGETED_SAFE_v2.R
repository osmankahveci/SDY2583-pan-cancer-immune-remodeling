# ============================================================
# SDY2583 CP23
# STEP 5 SAFE: Targeted threshold sensitivity
#
# Re-extracts CP23 targeted composite components under:
#   main       : original thresholds
#   permissive : positive-marker thresholds -0.2; dump-low threshold +0.2
#   stringent  : positive-marker thresholds +0.2; dump-low threshold -0.2
#
# Output:
#   outputs/CP23/07_threshold_sensitivity
#
# Important interpretation:
#   CP23 is a monocyte/macrophage-like myeloid remodeling panel.
#   Phenotypes are threshold-defined "like" populations.
#   Do not claim definitive macrophages, MDSCs, basophils, dendritic cells,
#   or tissue macrophage polarization from this panel alone.
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

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "07_threshold_sensitivity")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP1_fcs_inventory_marker_QC.RData")
step3a_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP3A_metadata_merge_age_QC.RData")
step4_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP4_composite_scores.RData")

if (!file.exists(step1_rdata)) stop("Step 1 RData bulunamadı: ", step1_rdata)
if (!file.exists(step3a_rdata)) stop("Step 3A RData bulunamadı: ", step3a_rdata)
if (!file.exists(step4_rdata)) stop("Step 4 RData bulunamadı: ", step4_rdata)

load(step1_rdata)
if (!exists("fcs_files")) stop("fcs_files bulunamadı.")
if (!exists("fcs_inventory")) stop("fcs_inventory bulunamadı.")
fcs_files_STEP1 <- fcs_files
fcs_inventory_STEP1 <- fcs_inventory

load(step3a_rdata)
if (!exists("cp23_analysis_data")) stop("cp23_analysis_data bulunamadı.")
metadata_for_models <- cp23_analysis_data %>%
  select(
    subject_id,
    disease_group,
    age_for_model,
    sex,
    sex_binary,
    model_ready_age_sex,
    model_ready_binary_sex
  ) %>%
  distinct(subject_id, .keep_all = TRUE)

load(step4_rdata)
if (!exists("composite_definitions_available")) stop("composite_definitions_available bulunamadı.")
if (!exists("composite_score_dictionary")) stop("composite_score_dictionary bulunamadı.")

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "07_threshold_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fcs_files <- fcs_files_STEP1
fcs_inventory <- fcs_inventory_STEP1

# ------------------------------------------------------------
# Threshold sets
# ------------------------------------------------------------

main_thresholds <- list(
  DUMP_LOW = 1.5,
  CD45 = 2.0,
  HLA_DR = 1.5,
  CD14 = 1.5,
  CD16 = 1.5,
  CD33 = 1.5,
  CD15 = 1.5,
  FceRI = 1.5,
  CD9 = 1.5,
  CD84 = 1.5,
  CD206 = 1.5,
  CD169 = 1.5
)

make_threshold_set <- function(base, positive_shift = 0, dump_low_shift = 0) {
  out <- base
  out$DUMP_LOW <- base$DUMP_LOW + dump_low_shift
  positive_names <- setdiff(names(base), "DUMP_LOW")
  for (nm in positive_names) out[[nm]] <- base[[nm]] + positive_shift
  out
}

threshold_sets <- list(
  main = main_thresholds,
  permissive = make_threshold_set(main_thresholds, positive_shift = -0.2, dump_low_shift = +0.2),
  stringent = make_threshold_set(main_thresholds, positive_shift = +0.2, dump_low_shift = -0.2)
)

threshold_table <- bind_rows(lapply(names(threshold_sets), function(ts) {
  tibble(
    threshold_set = ts,
    marker_or_gate = names(threshold_sets[[ts]]),
    threshold = unlist(threshold_sets[[ts]])
  )
}))

# ------------------------------------------------------------
# Helper functions copied from CP23 Step 2
# ------------------------------------------------------------

extract_subject_id <- function(x) {
  b <- basename(x)
  id <- stringr::str_extract(b, "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(b, "\\.fcs$"), id)
}

get_marker_map_safe <- function(ff) {
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
  ch$DUMP <- find_marker_channel_regex(map, "Viability.*CD3.*CCR3.*CD19.*CD7|Viability|Dump")
  ch$HLA_DR <- find_marker_channel_regex(map, "^HLA-DR$|HLA_DR|HLADR")
  ch$CD14 <- find_marker_channel_exact(map, "CD14")
  ch$CD16 <- find_marker_channel_exact(map, "CD16")
  ch$CD33 <- find_marker_channel_exact(map, "CD33")
  ch$CD15 <- find_marker_channel_exact(map, "CD15")
  ch$FceRI <- find_marker_channel_regex(map, "FceRI|FcER1|FcERI|FCER1")
  ch$CD9 <- find_marker_channel_exact(map, "CD9")
  ch$CD84 <- find_marker_channel_exact(map, "CD84")
  ch$CD206 <- find_marker_channel_exact(map, "CD206")
  ch$CD169 <- find_marker_channel_exact(map, "CD169")

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

event_qc_bin <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x < 50 ~ "<50",
    x < 100 ~ "50-99",
    x < 300 ~ "100-299",
    x < 1000 ~ "300-999",
    TRUE ~ ">=1000"
  )
}

# ------------------------------------------------------------
# Feature extraction for one file and one threshold set
# ------------------------------------------------------------

extract_one_cp23_threshold <- function(fp, thresholds, threshold_set) {

  tryCatch({

    ff_raw <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE,
      emptyValue = FALSE
    )

    marker_map_file <- get_marker_map_safe(ff_raw)
    ch <- build_channel_map(marker_map_file)

    required <- c(
      "DUMP", "CD45", "HLA_DR", "CD14", "CD16", "CD33",
      "CD15", "FceRI", "CD9", "CD84", "CD206", "CD169"
    )

    missing <- required[vapply(ch[required], function(x) is.na(x), logical(1))]
    if (length(missing) > 0) {
      stop("Missing required markers: ", paste(missing, collapse = ", "))
    }

    ff <- ff_raw %>%
      apply_compensation_safely() %>%
      transform_safely()

    ex <- flowCore::exprs(ff)

    total_events <- nrow(ex)

    dump <- get_vec(ex, ch$DUMP)
    cd45 <- get_vec(ex, ch$CD45)
    hla_dr <- get_vec(ex, ch$HLA_DR)
    cd14 <- get_vec(ex, ch$CD14)
    cd16 <- get_vec(ex, ch$CD16)
    cd33 <- get_vec(ex, ch$CD33)
    cd15 <- get_vec(ex, ch$CD15)
    fceri <- get_vec(ex, ch$FceRI)
    cd9 <- get_vec(ex, ch$CD9)
    cd84 <- get_vec(ex, ch$CD84)
    cd206 <- get_vec(ex, ch$CD206)
    cd169 <- get_vec(ex, ch$CD169)

    # Parent gates
    dump_low <- dump <= thresholds$DUMP_LOW
    cd45_dump_low <- dump_low & cd45 > thresholds$CD45

    # Broad myeloid/APC-like gates
    hladr_pos <- cd45_dump_low & hla_dr > thresholds$HLA_DR
    cd33_pos <- cd45_dump_low & cd33 > thresholds$CD33
    cd33_hladr_myeloid_like <- cd45_dump_low & cd33 > thresholds$CD33 & hla_dr > thresholds$HLA_DR
    cd15_gran_like <- cd45_dump_low & cd15 > thresholds$CD15
    cd33_cd15_myeloid_gran_like <- cd45_dump_low & cd33 > thresholds$CD33 & cd15 > thresholds$CD15

    # Monocyte-like phenotypes
    cd14_mono_like <- cd45_dump_low & cd14 > thresholds$CD14
    cd14_hladr_pos <- cd14_mono_like & hla_dr > thresholds$HLA_DR
    cd14_hladr_low <- cd14_mono_like & hla_dr <= thresholds$HLA_DR
    cd14_cd16neg_mono_like <- cd14_mono_like & cd16 <= thresholds$CD16
    cd14_cd16pos_mono_like <- cd14_mono_like & cd16 > thresholds$CD16
    cd14lowneg_cd16pos_hladr_like <- cd45_dump_low & cd14 <= thresholds$CD14 & cd16 > thresholds$CD16 & hla_dr > thresholds$HLA_DR

    # FcERI-associated myeloid/APC-like phenotypes
    fceri_pos <- cd45_dump_low & fceri > thresholds$FceRI
    fceri_hladr_apc_like <- cd45_dump_low & fceri > thresholds$FceRI & hla_dr > thresholds$HLA_DR
    cd14_fceri_mono_like <- cd14_mono_like & fceri > thresholds$FceRI
    cd33_fceri_hladr_myeloid_like <- cd45_dump_low & cd33 > thresholds$CD33 & fceri > thresholds$FceRI & hla_dr > thresholds$HLA_DR

    # Macrophage-like activation marker phenotypes
    cd14_cd206_pos <- cd14_mono_like & cd206 > thresholds$CD206
    cd14_cd169_pos <- cd14_mono_like & cd169 > thresholds$CD169
    cd14_cd206_cd169_pos <- cd14_mono_like & cd206 > thresholds$CD206 & cd169 > thresholds$CD169
    cd14_cd9_pos <- cd14_mono_like & cd9 > thresholds$CD9
    cd14_cd84_pos <- cd14_mono_like & cd84 > thresholds$CD84
    cd14_cd9_cd84_pos <- cd14_mono_like & cd9 > thresholds$CD9 & cd84 > thresholds$CD84

    hladr_cd206_pos <- hladr_pos & cd206 > thresholds$CD206
    hladr_cd169_pos <- hladr_pos & cd169 > thresholds$CD169
    hladr_cd206_cd169_pos <- hladr_pos & cd206 > thresholds$CD206 & cd169 > thresholds$CD169
    cd33_hladr_cd206_pos <- cd33_hladr_myeloid_like & cd206 > thresholds$CD206
    cd33_hladr_cd169_pos <- cd33_hladr_myeloid_like & cd169 > thresholds$CD169

    # Count denominators
    n_dump_low <- safe_sum(dump_low)
    n_cd45_dump_low <- safe_sum(cd45_dump_low)
    n_hladr_pos <- safe_sum(hladr_pos)
    n_cd33_pos <- safe_sum(cd33_pos)
    n_cd33_hladr <- safe_sum(cd33_hladr_myeloid_like)
    n_cd15 <- safe_sum(cd15_gran_like)
    n_cd33_cd15 <- safe_sum(cd33_cd15_myeloid_gran_like)

    n_cd14 <- safe_sum(cd14_mono_like)
    n_cd14_hladr_pos <- safe_sum(cd14_hladr_pos)
    n_cd14_hladr_low <- safe_sum(cd14_hladr_low)
    n_cd14_cd16neg <- safe_sum(cd14_cd16neg_mono_like)
    n_cd14_cd16pos <- safe_sum(cd14_cd16pos_mono_like)
    n_cd14lowneg_cd16pos <- safe_sum(cd14lowneg_cd16pos_hladr_like)

    n_fceri <- safe_sum(fceri_pos)
    n_fceri_hladr <- safe_sum(fceri_hladr_apc_like)
    n_cd14_fceri <- safe_sum(cd14_fceri_mono_like)
    n_cd33_fceri_hladr <- safe_sum(cd33_fceri_hladr_myeloid_like)

    n_cd14_cd206 <- safe_sum(cd14_cd206_pos)
    n_cd14_cd169 <- safe_sum(cd14_cd169_pos)
    n_cd14_cd206_cd169 <- safe_sum(cd14_cd206_cd169_pos)
    n_cd14_cd9 <- safe_sum(cd14_cd9_pos)
    n_cd14_cd84 <- safe_sum(cd14_cd84_pos)
    n_cd14_cd9_cd84 <- safe_sum(cd14_cd9_cd84_pos)

    n_hladr_cd206 <- safe_sum(hladr_cd206_pos)
    n_hladr_cd169 <- safe_sum(hladr_cd169_pos)
    n_hladr_cd206_cd169 <- safe_sum(hladr_cd206_cd169_pos)
    n_cd33_hladr_cd206 <- safe_sum(cd33_hladr_cd206_pos)
    n_cd33_hladr_cd169 <- safe_sum(cd33_hladr_cd169_pos)

    tibble(
      threshold_set = threshold_set,
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      feature_ok = TRUE,
      feature_error = NA_character_,
      total_events = total_events,
      n_dump_low = n_dump_low,
      n_cd45_dump_low = n_cd45_dump_low,
      n_hladr_pos = n_hladr_pos,
      n_cd33_hladr_myeloid_like = n_cd33_hladr,
      n_cd14_mono_like = n_cd14,

      pct_dump_low_within_total = safe_pct(n_dump_low, total_events),
      pct_cd45_dump_low_within_total = safe_pct(n_cd45_dump_low, total_events),

      pct_hladr_pos_within_cd45_dump_low = safe_pct(n_hladr_pos, n_cd45_dump_low),
      pct_cd33_pos_within_cd45_dump_low = safe_pct(n_cd33_pos, n_cd45_dump_low),
      pct_cd33_hladr_myeloid_like_within_cd45_dump_low = safe_pct(n_cd33_hladr, n_cd45_dump_low),
      pct_cd33_hladr_myeloid_like_within_hladr_pos = safe_pct(n_cd33_hladr, n_hladr_pos),
      pct_cd15_gran_like_within_cd45_dump_low = safe_pct(n_cd15, n_cd45_dump_low),
      pct_cd33_cd15_myeloid_gran_like_within_cd45_dump_low = safe_pct(n_cd33_cd15, n_cd45_dump_low),
      pct_cd15_pos_within_cd33_pos = safe_pct(n_cd33_cd15, n_cd33_pos),

      pct_cd14_mono_like_within_cd45_dump_low = safe_pct(n_cd14, n_cd45_dump_low),
      pct_cd14_hladr_pos_within_cd45_dump_low = safe_pct(n_cd14_hladr_pos, n_cd45_dump_low),
      pct_cd14_hladr_low_within_cd45_dump_low = safe_pct(n_cd14_hladr_low, n_cd45_dump_low),
      pct_hladr_pos_within_cd14_mono_like = safe_pct(n_cd14_hladr_pos, n_cd14),
      pct_hladr_low_within_cd14_mono_like = safe_pct(n_cd14_hladr_low, n_cd14),
      pct_cd14_cd16neg_mono_like_within_cd45_dump_low = safe_pct(n_cd14_cd16neg, n_cd45_dump_low),
      pct_cd14_cd16pos_mono_like_within_cd45_dump_low = safe_pct(n_cd14_cd16pos, n_cd45_dump_low),
      pct_cd16_pos_within_cd14_mono_like = safe_pct(n_cd14_cd16pos, n_cd14),
      pct_cd14lowneg_cd16pos_hladr_like_within_cd45_dump_low = safe_pct(n_cd14lowneg_cd16pos, n_cd45_dump_low),

      pct_fceri_pos_within_cd45_dump_low = safe_pct(n_fceri, n_cd45_dump_low),
      pct_fceri_hladr_apc_like_within_cd45_dump_low = safe_pct(n_fceri_hladr, n_cd45_dump_low),
      pct_fceri_hladr_apc_like_within_hladr_pos = safe_pct(n_fceri_hladr, n_hladr_pos),
      pct_cd14_fceri_mono_like_within_cd45_dump_low = safe_pct(n_cd14_fceri, n_cd45_dump_low),
      pct_fceri_pos_within_cd14_mono_like = safe_pct(n_cd14_fceri, n_cd14),
      pct_cd33_fceri_hladr_myeloid_like_within_cd45_dump_low = safe_pct(n_cd33_fceri_hladr, n_cd45_dump_low),

      pct_cd14_cd206_pos_within_cd45_dump_low = safe_pct(n_cd14_cd206, n_cd45_dump_low),
      pct_cd206_pos_within_cd14_mono_like = safe_pct(n_cd14_cd206, n_cd14),
      pct_cd14_cd169_pos_within_cd45_dump_low = safe_pct(n_cd14_cd169, n_cd45_dump_low),
      pct_cd169_pos_within_cd14_mono_like = safe_pct(n_cd14_cd169, n_cd14),
      pct_cd14_cd206_cd169_pos_within_cd45_dump_low = safe_pct(n_cd14_cd206_cd169, n_cd45_dump_low),
      pct_cd206_cd169_pos_within_cd14_mono_like = safe_pct(n_cd14_cd206_cd169, n_cd14),
      pct_cd14_cd9_pos_within_cd45_dump_low = safe_pct(n_cd14_cd9, n_cd45_dump_low),
      pct_cd9_pos_within_cd14_mono_like = safe_pct(n_cd14_cd9, n_cd14),
      pct_cd14_cd84_pos_within_cd45_dump_low = safe_pct(n_cd14_cd84, n_cd45_dump_low),
      pct_cd84_pos_within_cd14_mono_like = safe_pct(n_cd14_cd84, n_cd14),
      pct_cd14_cd9_cd84_pos_within_cd45_dump_low = safe_pct(n_cd14_cd9_cd84, n_cd45_dump_low),
      pct_cd9_cd84_pos_within_cd14_mono_like = safe_pct(n_cd14_cd9_cd84, n_cd14),

      pct_hladr_cd206_pos_within_cd45_dump_low = safe_pct(n_hladr_cd206, n_cd45_dump_low),
      pct_hladr_cd169_pos_within_cd45_dump_low = safe_pct(n_hladr_cd169, n_cd45_dump_low),
      pct_hladr_cd206_cd169_pos_within_cd45_dump_low = safe_pct(n_hladr_cd206_cd169, n_cd45_dump_low),
      pct_cd33_hladr_cd206_pos_within_cd45_dump_low = safe_pct(n_cd33_hladr_cd206, n_cd45_dump_low),
      pct_cd33_hladr_cd169_pos_within_cd45_dump_low = safe_pct(n_cd33_hladr_cd169, n_cd45_dump_low),

      ratio_cd33_hladr_myeloid_to_cd14_mono_like = safe_ratio(n_cd33_hladr, n_cd14),
      ratio_cd15_gran_like_to_cd14_mono_like = safe_ratio(n_cd15, n_cd14),
      ratio_cd33_cd15_to_cd14_mono_like = safe_ratio(n_cd33_cd15, n_cd14),
      ratio_cd14_hladr_low_to_hladr_pos = safe_ratio(n_cd14_hladr_low, n_hladr_pos),
      ratio_cd206_cd169_cd14_to_cd15_gran_like = safe_ratio(n_cd14_cd206_cd169, n_cd15),

      median_HLA_DR_in_cd45_dump_low = safe_median(hla_dr[cd45_dump_low]),
      median_CD14_in_cd45_dump_low = safe_median(cd14[cd45_dump_low]),
      median_CD16_in_cd45_dump_low = safe_median(cd16[cd45_dump_low]),
      median_CD33_in_cd45_dump_low = safe_median(cd33[cd45_dump_low]),
      median_CD15_in_cd45_dump_low = safe_median(cd15[cd45_dump_low]),
      median_FceRI_in_cd45_dump_low = safe_median(fceri[cd45_dump_low]),
      median_CD9_in_cd45_dump_low = safe_median(cd9[cd45_dump_low]),
      median_CD84_in_cd45_dump_low = safe_median(cd84[cd45_dump_low]),
      median_CD206_in_cd45_dump_low = safe_median(cd206[cd45_dump_low]),
      median_CD169_in_cd45_dump_low = safe_median(cd169[cd45_dump_low]),

      median_HLA_DR_in_cd14_mono_like = safe_median(hla_dr[cd14_mono_like]),
      median_CD16_in_cd14_mono_like = safe_median(cd16[cd14_mono_like]),
      median_CD33_in_cd14_mono_like = safe_median(cd33[cd14_mono_like]),
      median_CD15_in_cd14_mono_like = safe_median(cd15[cd14_mono_like]),
      median_FceRI_in_cd14_mono_like = safe_median(fceri[cd14_mono_like]),
      median_CD9_in_cd14_mono_like = safe_median(cd9[cd14_mono_like]),
      median_CD84_in_cd14_mono_like = safe_median(cd84[cd14_mono_like]),
      median_CD206_in_cd14_mono_like = safe_median(cd206[cd14_mono_like]),
      median_CD169_in_cd14_mono_like = safe_median(cd169[cd14_mono_like]),

      median_CD33_in_hladr_pos = safe_median(cd33[hladr_pos]),
      median_CD14_in_hladr_pos = safe_median(cd14[hladr_pos]),
      median_CD15_in_cd33_pos = safe_median(cd15[cd33_pos]),
      median_CD206_in_cd33_hladr_myeloid_like = safe_median(cd206[cd33_hladr_myeloid_like]),
      median_CD169_in_cd33_hladr_myeloid_like = safe_median(cd169[cd33_hladr_myeloid_like])
    )

  }, error = function(e) {
    tibble(
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
# Re-extraction across threshold sets
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 5 STARTED: TARGETED THRESHOLD SENSITIVITY\n")
cat("============================================================\n")

all_feature_rows <- list()

for (ts in names(threshold_sets)) {
  cat("\nThreshold set: ", ts, "\n", sep = "")

  th <- threshold_sets[[ts]]
  rows <- vector("list", length(fcs_files))

  for (i in seq_along(fcs_files)) {
    if (i %% 50 == 0) cat("  processed ", i, " / ", length(fcs_files), "\n", sep = "")
    rows[[i]] <- extract_one_cp23_threshold(fcs_files[i], th, ts)
  }

  all_feature_rows[[ts]] <- bind_rows(rows)
}

threshold_features_raw <- bind_rows(all_feature_rows)

technical_flags <- fcs_inventory %>%
  select(
    subject_id,
    file_name,
    channel_order_mismatch_file,
    marker_order_mismatch_file,
    marker_set_mismatch_file,
    has_spillover_keyword,
    spillover_keyword_name
  ) %>%
  distinct(subject_id, file_name, .keep_all = TRUE)

threshold_features <- threshold_features_raw %>%
  left_join(technical_flags, by = c("subject_id", "file_name")) %>%
  left_join(metadata_for_models, by = "subject_id") %>%
  mutate(
    cd45_dump_low_event_qc_bin = event_qc_bin(n_cd45_dump_low),
    hladr_myeloid_event_qc_bin = event_qc_bin(n_cd33_hladr_myeloid_like)
  )

# ------------------------------------------------------------
# Build threshold-specific composite scores
# ------------------------------------------------------------

z_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)
  if (is.na(sig) || sig == 0) return(rep(NA_real_, length(x)))
  (x - mu) / sig
}

score_names <- unique(composite_definitions_available$composite_score)

build_scores_for_threshold <- function(df_one_set) {
  out <- df_one_set

  for (i in seq_len(nrow(composite_definitions_available))) {
    ff <- composite_definitions_available$feature[i]
    direction <- composite_definitions_available$component_direction[i]
    zname <- paste0("z_oriented__", ff)

    if (ff %in% names(out) && !(zname %in% names(out))) {
      out[[zname]] <- z_safe(out[[ff]]) * direction
    }
  }

  for (ss in score_names) {
    comp_features <- composite_definitions_available %>%
      filter(composite_score == ss) %>%
      pull(feature)

    zcols <- paste0("z_oriented__", comp_features)
    zcols <- zcols[zcols %in% names(out)]

    if (length(zcols) == 0) {
      out[[ss]] <- NA_real_
    } else {
      out[[ss]] <- rowMeans(out[, zcols, drop = FALSE], na.rm = TRUE)
      all_na <- apply(out[, zcols, drop = FALSE], 1, function(x) all(is.na(x)))
      out[[ss]][all_na] <- NA_real_
    }
  }

  out
}

threshold_score_data <- bind_rows(lapply(names(threshold_sets), function(ts) {
  threshold_features %>%
    filter(threshold_set == ts) %>%
    build_scores_for_threshold()
}))

# ------------------------------------------------------------
# Model helpers
# ------------------------------------------------------------

run_model_one <- function(df, variable_name, variable_type, threshold_set, sex_var = "sex") {

  if (!(variable_name %in% names(df))) return(NULL)

  work <- df %>%
    filter(threshold_set == !!threshold_set, model_ready_age_sex == TRUE) %>%
    transmute(
      value = suppressWarnings(as.numeric(.data[[variable_name]])),
      disease_group = disease_group,
      age_for_model = age_for_model,
      sex_model = .data[[sex_var]]
    ) %>%
    filter(
      !is.na(value),
      !is.na(disease_group),
      !is.na(age_for_model),
      !is.na(sex_model)
    ) %>%
    mutate(
      disease_group = factor(
        as.character(disease_group),
        levels = c("Healthy control", "Cancer patient")
      ),
      sex_model = factor(sex_model)
    )

  if (nrow(work) < 30) return(NULL)
  if (length(unique(work$disease_group)) < 2) return(NULL)

  n_healthy <- sum(work$disease_group == "Healthy control")
  n_cancer <- sum(work$disease_group == "Cancer patient")
  if (n_healthy < 10 || n_cancer < 10) return(NULL)

  fit <- tryCatch(
    lm(value ~ disease_group + age_for_model + sex_model, data = work),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  tt <- tryCatch(broom::tidy(fit, conf.int = TRUE), error = function(e) NULL)
  if (is.null(tt)) return(NULL)

  disease_term <- "disease_groupCancer patient"
  if (!(disease_term %in% tt$term)) return(NULL)
  term_row <- tt %>% filter(term == disease_term)

  tibble(
    threshold_set = threshold_set,
    variable = variable_name,
    variable_type = variable_type,
    n_model = nrow(work),
    n_healthy = n_healthy,
    n_cancer = n_cancer,
    healthy_mean = mean(work$value[work$disease_group == "Healthy control"], na.rm = TRUE),
    cancer_mean = mean(work$value[work$disease_group == "Cancer patient"], na.rm = TRUE),
    healthy_median = median(work$value[work$disease_group == "Healthy control"], na.rm = TRUE),
    cancer_median = median(work$value[work$disease_group == "Cancer patient"], na.rm = TRUE),
    beta_cancer_vs_healthy = term_row$estimate[1],
    conf_low = term_row$conf.low[1],
    conf_high = term_row$conf.high[1],
    p_value = term_row$p.value[1],
    direction = case_when(
      term_row$estimate[1] > 0 ~ "higher_in_cancer",
      term_row$estimate[1] < 0 ~ "lower_in_cancer",
      TRUE ~ "no_direction"
    )
  )
}

target_feature_names <- unique(composite_definitions_available$feature)
target_feature_names <- target_feature_names[target_feature_names %in% names(threshold_score_data)]

target_variables <- bind_rows(
  tibble(variable = score_names, variable_type = "composite"),
  tibble(variable = target_feature_names, variable_type = "feature")
) %>%
  distinct(variable, .keep_all = TRUE)

threshold_model_results <- bind_rows(lapply(names(threshold_sets), function(ts) {
  bind_rows(lapply(seq_len(nrow(target_variables)), function(i) {
    run_model_one(
      threshold_score_data,
      variable_name = target_variables$variable[i],
      variable_type = target_variables$variable_type[i],
      threshold_set = ts,
      sex_var = "sex"
    )
  })) %>%
    group_by(threshold_set) %>%
    mutate(FDR_global = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}))

# ------------------------------------------------------------
# Annotate and reshape robustness
# ------------------------------------------------------------

# Robust annotation builder.
# Some earlier Step 4 objects do not contain an `interpretation` column.
# Therefore all annotation columns are handled with any_of() and safe fallbacks.

feature_annotation <- composite_definitions_available %>%
  select(any_of(c("feature", "component_axis", "interpretation"))) %>%
  distinct(feature, .keep_all = TRUE) %>%
  rename(variable = feature)

if (!("module" %in% names(feature_annotation))) {
  if ("component_axis" %in% names(feature_annotation)) {
    feature_annotation <- feature_annotation %>% rename(module = component_axis)
  } else {
    feature_annotation$module <- "feature"
  }
}

if (!("interpretation" %in% names(feature_annotation))) {
  feature_annotation$interpretation <- feature_annotation$variable
}

composite_annotation <- composite_score_dictionary %>%
  select(any_of(c("composite_score", "score_label", "interpretation"))) %>%
  distinct(composite_score, .keep_all = TRUE) %>%
  rename(variable = composite_score)

if (!("score_label" %in% names(composite_annotation))) {
  composite_annotation$score_label <- composite_annotation$variable
}

if (!("interpretation" %in% names(composite_annotation))) {
  composite_annotation$interpretation <- composite_annotation$score_label
}

composite_annotation <- composite_annotation %>%
  transmute(
    variable = variable,
    module = "composite",
    interpretation = interpretation
  )

variable_annotation <- bind_rows(
  feature_annotation %>%
    transmute(
      variable = variable,
      module = module,
      interpretation = interpretation
    ),
  composite_annotation
) %>%
  distinct(variable, .keep_all = TRUE)

threshold_model_results_annotated <- threshold_model_results %>%
  left_join(variable_annotation, by = "variable") %>%
  arrange(variable_type, FDR_global, p_value)

wide_robustness <- threshold_model_results_annotated %>%
  select(
    variable, variable_type, module, interpretation,
    threshold_set, beta_cancer_vs_healthy, direction, p_value, FDR_global
  ) %>%
  pivot_wider(
    names_from = threshold_set,
    values_from = c(beta_cancer_vs_healthy, direction, p_value, FDR_global),
    names_sep = "__"
  ) %>%
  mutate(
    direction_preserved_all_sets =
      direction__main == direction__permissive &
      direction__main == direction__stringent,
    global_FDR_lt_0p05_all_sets =
      FDR_global__main < 0.05 &
      FDR_global__permissive < 0.05 &
      FDR_global__stringent < 0.05,
    robustness_class = case_when(
      direction_preserved_all_sets == TRUE & global_FDR_lt_0p05_all_sets == TRUE ~
        "direction_and_global_FDR_preserved_all_sets",
      direction_preserved_all_sets == TRUE & global_FDR_lt_0p05_all_sets == FALSE ~
        "direction_preserved_FDR_not_all_sets",
      TRUE ~ "direction_not_preserved"
    )
  ) %>%
  arrange(
    factor(variable_type, levels = c("composite", "feature")),
    desc(robustness_class == "direction_and_global_FDR_preserved_all_sets"),
    FDR_global__main
  )

threshold_robustness_counts <- wide_robustness %>%
  count(variable_type, robustness_class, name = "n") %>%
  arrange(variable_type, robustness_class)

top_threshold_robust_variables <- wide_robustness %>%
  filter(robustness_class == "direction_and_global_FDR_preserved_all_sets") %>%
  arrange(
    factor(variable_type, levels = c("composite", "feature")),
    FDR_global__main
  )

direction_preserved_FDR_not_all_sets <- wide_robustness %>%
  filter(robustness_class == "direction_preserved_FDR_not_all_sets") %>%
  arrange(variable_type, FDR_global__main)

direction_not_preserved <- wide_robustness %>%
  filter(robustness_class == "direction_not_preserved") %>%
  arrange(variable_type, FDR_global__main)

threshold_extraction_summary <- threshold_score_data %>%
  group_by(threshold_set) %>%
  summarise(
    n_rows = n(),
    n_unique_subjects = n_distinct(subject_id),
    n_feature_ok = sum(feature_ok == TRUE, na.rm = TRUE),
    n_feature_failed = sum(feature_ok == FALSE, na.rm = TRUE),
    median_total_events = median(total_events, na.rm = TRUE),
    median_cd45_dump_low_events = median(n_cd45_dump_low, na.rm = TRUE),
    median_cd33_hladr_myeloid_like_events = median(n_cd33_hladr_myeloid_like, na.rm = TRUE),
    median_pct_cd45_dump_low_within_total = median(pct_cd45_dump_low_within_total, na.rm = TRUE),
    median_pct_cd14_mono_like_within_cd45_dump_low =
      median(pct_cd14_mono_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd33_hladr_myeloid_like_within_cd45_dump_low =
      median(pct_cd33_hladr_myeloid_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_cd9_cd84_pos_within_cd45_dump_low =
      median(pct_cd14_cd9_cd84_pos_within_cd45_dump_low, na.rm = TRUE),
    median_pct_fceri_hladr_apc_like_within_cd45_dump_low =
      median(pct_fceri_hladr_apc_like_within_cd45_dump_low, na.rm = TRUE),
    .groups = "drop"
  )

failed_threshold_files <- threshold_score_data %>%
  filter(feature_ok == FALSE) %>%
  select(threshold_set, subject_id, file_name, file_path, feature_error)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write_csv(threshold_table, file.path(out_dir, "SDY2583_CP23_threshold_table_STEP5.csv"))
write_csv(threshold_score_data, file.path(out_dir, "SDY2583_CP23_threshold_score_data_STEP5.csv"))
write_csv(threshold_extraction_summary, file.path(out_dir, "SDY2583_CP23_threshold_extraction_summary_STEP5.csv"))
write_csv(target_variables, file.path(out_dir, "SDY2583_CP23_target_variables_STEP5.csv"))
write_csv(threshold_model_results_annotated, file.path(out_dir, "SDY2583_CP23_threshold_model_results_STEP5.csv"))
write_csv(wide_robustness, file.path(out_dir, "SDY2583_CP23_threshold_wide_robustness_STEP5.csv"))
write_csv(threshold_robustness_counts, file.path(out_dir, "SDY2583_CP23_threshold_robustness_counts_STEP5.csv"))
write_csv(top_threshold_robust_variables, file.path(out_dir, "SDY2583_CP23_top_threshold_robust_variables_STEP5.csv"))
write_csv(direction_preserved_FDR_not_all_sets, file.path(out_dir, "SDY2583_CP23_direction_preserved_FDR_not_all_sets_STEP5.csv"))
write_csv(direction_not_preserved, file.path(out_dir, "SDY2583_CP23_direction_not_preserved_STEP5.csv"))
write_csv(failed_threshold_files, file.path(out_dir, "SDY2583_CP23_failed_threshold_files_STEP5.csv"))

save(
  threshold_sets,
  threshold_table,
  threshold_score_data,
  threshold_extraction_summary,
  target_variables,
  threshold_model_results_annotated,
  wide_robustness,
  threshold_robustness_counts,
  top_threshold_robust_variables,
  direction_preserved_FDR_not_all_sets,
  direction_not_preserved,
  failed_threshold_files,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP5_threshold_sensitivity.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 5 COMPLETE: TARGETED THRESHOLD SENSITIVITY\n")
cat("============================================================\n")

cat("\nThreshold extraction summary:\n")
print(as.data.frame(threshold_extraction_summary), row.names = FALSE)

cat("\nThreshold robustness counts:\n")
print(as.data.frame(threshold_robustness_counts), row.names = FALSE)

cat("\nTop threshold-robust variables:\n")
print(as.data.frame(top_threshold_robust_variables), row.names = FALSE)

cat("\nDirection preserved but FDR not all sets:\n")
print(as.data.frame(direction_preserved_FDR_not_all_sets), row.names = FALSE)

cat("\nDirection not preserved:\n")
print(as.data.frame(direction_not_preserved), row.names = FALSE)

cat("\nFailed threshold files:\n")
print(as.data.frame(failed_threshold_files), row.names = FALSE)

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
