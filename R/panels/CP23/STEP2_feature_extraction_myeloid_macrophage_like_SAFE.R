# ============================================================
# SDY2583 CP23
# STEP 2 SAFE: Feature extraction
#
# Panel interpretation:
#   CP23 is treated as a circulating myeloid / monocyte /
#   macrophage-like remodeling panel.
#
# Markers:
#   CD45, dump/viability_CD3_CCR3_CD19_CD7, HLA-DR,
#   CD14, CD16, CD33, CD15, FcERI, CD9, CD84, CD206, CD169
#
# Key caution:
#   Phenotypes are threshold-defined "like" populations.
#   Do not claim definitive tissue macrophages, MDSCs, M1/M2,
#   basophils, or dendritic cells from this panel alone.
#
# Output:
#   outputs/CP23/02_feature_extraction
#   outputs/CP23/11_RData
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr")

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
  library(purrr)
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
# 2. Paths and Step 1 load
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "02_feature_extraction")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP23_STEP1_fcs_inventory_marker_QC.RData")

if (!file.exists(step1_rdata)) {
  stop("Step 1 RData bulunamadı: ", step1_rdata)
}

load(step1_rdata)

# Reset paths after RData load.
analysis_dir <- sd_analysis_dir("CP23")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "02_feature_extraction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("fcs_files")) stop("fcs_files bulunamadı.")
if (!exists("fcs_inventory")) stop("fcs_inventory bulunamadı.")

# ------------------------------------------------------------
# 3. Thresholds
# ------------------------------------------------------------

thresholds <- list(
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

threshold_table <- tibble(
  marker_or_gate = names(thresholds),
  threshold = unlist(thresholds)
)

# ------------------------------------------------------------
# 4. Helper functions
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
# 5. Feature extraction for one file
# ------------------------------------------------------------

extract_one_cp23 <- function(fp) {

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

    # Macrophage-like activation/polarization marker phenotypes
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

      # Composition
      pct_dump_low_within_total = safe_pct(n_dump_low, total_events),
      pct_cd45_dump_low_within_total = safe_pct(n_cd45_dump_low, total_events),

      # HLA-DR / CD33 / CD15 myeloid architecture
      pct_hladr_pos_within_cd45_dump_low = safe_pct(n_hladr_pos, n_cd45_dump_low),
      pct_cd33_pos_within_cd45_dump_low = safe_pct(n_cd33_pos, n_cd45_dump_low),
      pct_cd33_hladr_myeloid_like_within_cd45_dump_low = safe_pct(n_cd33_hladr, n_cd45_dump_low),
      pct_cd33_hladr_myeloid_like_within_hladr_pos = safe_pct(n_cd33_hladr, n_hladr_pos),
      pct_cd15_gran_like_within_cd45_dump_low = safe_pct(n_cd15, n_cd45_dump_low),
      pct_cd33_cd15_myeloid_gran_like_within_cd45_dump_low = safe_pct(n_cd33_cd15, n_cd45_dump_low),
      pct_cd15_pos_within_cd33_pos = safe_pct(n_cd33_cd15, n_cd33_pos),

      # Monocyte-like architecture
      pct_cd14_mono_like_within_cd45_dump_low = safe_pct(n_cd14, n_cd45_dump_low),
      pct_cd14_hladr_pos_within_cd45_dump_low = safe_pct(n_cd14_hladr_pos, n_cd45_dump_low),
      pct_cd14_hladr_low_within_cd45_dump_low = safe_pct(n_cd14_hladr_low, n_cd45_dump_low),
      pct_hladr_pos_within_cd14_mono_like = safe_pct(n_cd14_hladr_pos, n_cd14),
      pct_hladr_low_within_cd14_mono_like = safe_pct(n_cd14_hladr_low, n_cd14),
      pct_cd14_cd16neg_mono_like_within_cd45_dump_low = safe_pct(n_cd14_cd16neg, n_cd45_dump_low),
      pct_cd14_cd16pos_mono_like_within_cd45_dump_low = safe_pct(n_cd14_cd16pos, n_cd45_dump_low),
      pct_cd16_pos_within_cd14_mono_like = safe_pct(n_cd14_cd16pos, n_cd14),
      pct_cd14lowneg_cd16pos_hladr_like_within_cd45_dump_low = safe_pct(n_cd14lowneg_cd16pos, n_cd45_dump_low),

      # FcERI-associated myeloid/APC-like axis
      pct_fceri_pos_within_cd45_dump_low = safe_pct(n_fceri, n_cd45_dump_low),
      pct_fceri_hladr_apc_like_within_cd45_dump_low = safe_pct(n_fceri_hladr, n_cd45_dump_low),
      pct_fceri_hladr_apc_like_within_hladr_pos = safe_pct(n_fceri_hladr, n_hladr_pos),
      pct_cd14_fceri_mono_like_within_cd45_dump_low = safe_pct(n_cd14_fceri, n_cd45_dump_low),
      pct_fceri_pos_within_cd14_mono_like = safe_pct(n_cd14_fceri, n_cd14),
      pct_cd33_fceri_hladr_myeloid_like_within_cd45_dump_low = safe_pct(n_cd33_fceri_hladr, n_cd45_dump_low),

      # CD206/CD169/CD9/CD84 macrophage-like marker phenotypes
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

      # Ratios
      ratio_cd33_hladr_myeloid_to_cd14_mono_like = safe_ratio(n_cd33_hladr, n_cd14),
      ratio_cd15_gran_like_to_cd14_mono_like = safe_ratio(n_cd15, n_cd14),
      ratio_cd33_cd15_to_cd14_mono_like = safe_ratio(n_cd33_cd15, n_cd14),
      ratio_cd14_hladr_low_to_hladr_pos = safe_ratio(n_cd14_hladr_low, n_hladr_pos),
      ratio_cd206_cd169_cd14_to_cd15_gran_like = safe_ratio(n_cd14_cd206_cd169, n_cd15),

      # Marker medians in CD45+ dump-low parent
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

      # Marker medians in CD14+ monocyte-like parent
      median_HLA_DR_in_cd14_mono_like = safe_median(hla_dr[cd14_mono_like]),
      median_CD16_in_cd14_mono_like = safe_median(cd16[cd14_mono_like]),
      median_CD33_in_cd14_mono_like = safe_median(cd33[cd14_mono_like]),
      median_CD15_in_cd14_mono_like = safe_median(cd15[cd14_mono_like]),
      median_FceRI_in_cd14_mono_like = safe_median(fceri[cd14_mono_like]),
      median_CD9_in_cd14_mono_like = safe_median(cd9[cd14_mono_like]),
      median_CD84_in_cd14_mono_like = safe_median(cd84[cd14_mono_like]),
      median_CD206_in_cd14_mono_like = safe_median(cd206[cd14_mono_like]),
      median_CD169_in_cd14_mono_like = safe_median(cd169[cd14_mono_like]),

      # Marker medians in HLA-DR+ / CD33+ parents
      median_CD33_in_hladr_pos = safe_median(cd33[hladr_pos]),
      median_CD14_in_hladr_pos = safe_median(cd14[hladr_pos]),
      median_CD15_in_cd33_pos = safe_median(cd15[cd33_pos]),
      median_CD206_in_cd33_hladr_myeloid_like = safe_median(cd206[cd33_hladr_myeloid_like]),
      median_CD169_in_cd33_hladr_myeloid_like = safe_median(cd169[cd33_hladr_myeloid_like])
    )

  }, error = function(e) {
    tibble(
      subject_id = extract_subject_id(fp),
      file_name = basename(fp),
      file_path = fp,
      feature_ok = FALSE,
      feature_error = as.character(e$message)
    )
  })
}

# ------------------------------------------------------------
# 6. Run extraction
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 2 STARTED: FEATURE EXTRACTION\n")
cat("============================================================\n")

feature_rows <- vector("list", length(fcs_files))

for (i in seq_along(fcs_files)) {
  if (i %% 50 == 0) cat("  processed ", i, " / ", length(fcs_files), "\n", sep = "")
  feature_rows[[i]] <- extract_one_cp23(fcs_files[i])
}

cp23_features_raw <- bind_rows(feature_rows)

# ------------------------------------------------------------
# 7. Add technical mismatch flags from Step 1
# ------------------------------------------------------------

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

cp23_features <- cp23_features_raw %>%
  left_join(technical_flags, by = c("subject_id", "file_name")) %>%
  mutate(
    cd45_dump_low_event_qc_bin = event_qc_bin(n_cd45_dump_low),
    hladr_myeloid_event_qc_bin = event_qc_bin(n_cd33_hladr_myeloid_like)
  )

# ------------------------------------------------------------
# 8. Feature dictionary
# ------------------------------------------------------------

feature_dictionary_full <- tribble(
  ~feature, ~module, ~interpretation,

  "pct_cd45_dump_low_within_total", "composition", "CD45+ dump-low myeloid-enriched fraction among total events",

  "pct_hladr_pos_within_cd45_dump_low", "HLA_DR_CD33_myeloid_axis", "HLA-DR+ fraction within CD45+ dump-low events",
  "pct_cd33_pos_within_cd45_dump_low", "HLA_DR_CD33_myeloid_axis", "CD33+ myeloid-like fraction within CD45+ dump-low events",
  "pct_cd33_hladr_myeloid_like_within_cd45_dump_low", "HLA_DR_CD33_myeloid_axis", "CD33+HLA-DR+ myeloid/APC-like fraction within CD45+ dump-low events",
  "pct_cd33_hladr_myeloid_like_within_hladr_pos", "HLA_DR_CD33_myeloid_axis", "CD33+ fraction within HLA-DR+ CD45+ dump-low events",
  "pct_cd15_gran_like_within_cd45_dump_low", "CD15_granulocytic_like_axis", "CD15+ granulocytic-like fraction within CD45+ dump-low events",
  "pct_cd33_cd15_myeloid_gran_like_within_cd45_dump_low", "CD15_granulocytic_like_axis", "CD33+CD15+ myeloid/granulocytic-like fraction within CD45+ dump-low events",
  "pct_cd15_pos_within_cd33_pos", "CD15_granulocytic_like_axis", "CD15+ fraction within CD33+ CD45+ dump-low events",

  "pct_cd14_mono_like_within_cd45_dump_low", "monocyte_like_axis", "CD14+ monocyte-like fraction within CD45+ dump-low events",
  "pct_cd14_hladr_pos_within_cd45_dump_low", "monocyte_like_axis", "CD14+HLA-DR+ monocyte-like fraction within CD45+ dump-low events",
  "pct_cd14_hladr_low_within_cd45_dump_low", "monocyte_like_axis", "CD14+HLA-DR-low monocyte-like fraction within CD45+ dump-low events",
  "pct_hladr_pos_within_cd14_mono_like", "monocyte_like_axis", "HLA-DR+ fraction within CD14+ monocyte-like cells",
  "pct_hladr_low_within_cd14_mono_like", "monocyte_like_axis", "HLA-DR-low fraction within CD14+ monocyte-like cells",
  "pct_cd14_cd16neg_mono_like_within_cd45_dump_low", "monocyte_like_axis", "CD14+CD16- monocyte-like fraction within CD45+ dump-low events",
  "pct_cd14_cd16pos_mono_like_within_cd45_dump_low", "monocyte_like_axis", "CD14+CD16+ monocyte-like fraction within CD45+ dump-low events",
  "pct_cd16_pos_within_cd14_mono_like", "monocyte_like_axis", "CD16+ fraction within CD14+ monocyte-like cells",
  "pct_cd14lowneg_cd16pos_hladr_like_within_cd45_dump_low", "monocyte_like_axis", "CD14-low/negative CD16+HLA-DR+ monocyte-like fraction within CD45+ dump-low events",

  "pct_fceri_pos_within_cd45_dump_low", "FcERI_myeloid_axis", "FcERI+ fraction within CD45+ dump-low events",
  "pct_fceri_hladr_apc_like_within_cd45_dump_low", "FcERI_myeloid_axis", "FcERI+HLA-DR+ APC-like fraction within CD45+ dump-low events",
  "pct_fceri_hladr_apc_like_within_hladr_pos", "FcERI_myeloid_axis", "FcERI+ fraction within HLA-DR+ CD45+ dump-low events",
  "pct_cd14_fceri_mono_like_within_cd45_dump_low", "FcERI_myeloid_axis", "CD14+FcERI+ monocyte-like fraction within CD45+ dump-low events",
  "pct_fceri_pos_within_cd14_mono_like", "FcERI_myeloid_axis", "FcERI+ fraction within CD14+ monocyte-like cells",
  "pct_cd33_fceri_hladr_myeloid_like_within_cd45_dump_low", "FcERI_myeloid_axis", "CD33+FcERI+HLA-DR+ myeloid/APC-like fraction within CD45+ dump-low events",

  "pct_cd14_cd206_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "CD14+CD206+ macrophage-like marker fraction within CD45+ dump-low events",
  "pct_cd206_pos_within_cd14_mono_like", "CD206_CD169_macrophage_like_axis", "CD206+ fraction within CD14+ monocyte-like cells",
  "pct_cd14_cd169_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "CD14+CD169+ macrophage-like marker fraction within CD45+ dump-low events",
  "pct_cd169_pos_within_cd14_mono_like", "CD206_CD169_macrophage_like_axis", "CD169+ fraction within CD14+ monocyte-like cells",
  "pct_cd14_cd206_cd169_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "CD14+CD206+CD169+ macrophage-like marker fraction within CD45+ dump-low events",
  "pct_cd206_cd169_pos_within_cd14_mono_like", "CD206_CD169_macrophage_like_axis", "CD206+CD169+ fraction within CD14+ monocyte-like cells",
  "pct_hladr_cd206_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "HLA-DR+CD206+ macrophage-like marker fraction within CD45+ dump-low events",
  "pct_hladr_cd169_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "HLA-DR+CD169+ macrophage-like marker fraction within CD45+ dump-low events",
  "pct_hladr_cd206_cd169_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "HLA-DR+CD206+CD169+ macrophage-like marker fraction within CD45+ dump-low events",
  "pct_cd33_hladr_cd206_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "CD33+HLA-DR+CD206+ myeloid/APC-like fraction within CD45+ dump-low events",
  "pct_cd33_hladr_cd169_pos_within_cd45_dump_low", "CD206_CD169_macrophage_like_axis", "CD33+HLA-DR+CD169+ myeloid/APC-like fraction within CD45+ dump-low events",

  "pct_cd14_cd9_pos_within_cd45_dump_low", "CD9_CD84_activation_axis", "CD14+CD9+ activation-marker fraction within CD45+ dump-low events",
  "pct_cd9_pos_within_cd14_mono_like", "CD9_CD84_activation_axis", "CD9+ fraction within CD14+ monocyte-like cells",
  "pct_cd14_cd84_pos_within_cd45_dump_low", "CD9_CD84_activation_axis", "CD14+CD84+ activation-marker fraction within CD45+ dump-low events",
  "pct_cd84_pos_within_cd14_mono_like", "CD9_CD84_activation_axis", "CD84+ fraction within CD14+ monocyte-like cells",
  "pct_cd14_cd9_cd84_pos_within_cd45_dump_low", "CD9_CD84_activation_axis", "CD14+CD9+CD84+ activation-marker fraction within CD45+ dump-low events",
  "pct_cd9_cd84_pos_within_cd14_mono_like", "CD9_CD84_activation_axis", "CD9+CD84+ fraction within CD14+ monocyte-like cells",

  "ratio_cd33_hladr_myeloid_to_cd14_mono_like", "ratio_balance", "Ratio of CD33+HLA-DR+ myeloid/APC-like to CD14+ monocyte-like events",
  "ratio_cd15_gran_like_to_cd14_mono_like", "ratio_balance", "Ratio of CD15+ granulocytic-like to CD14+ monocyte-like events",
  "ratio_cd33_cd15_to_cd14_mono_like", "ratio_balance", "Ratio of CD33+CD15+ myeloid/granulocytic-like to CD14+ monocyte-like events",
  "ratio_cd14_hladr_low_to_hladr_pos", "ratio_balance", "Ratio of CD14+HLA-DR-low monocyte-like to HLA-DR+ myeloid/APC-like events",
  "ratio_cd206_cd169_cd14_to_cd15_gran_like", "ratio_balance", "Ratio of CD14+CD206+CD169+ macrophage-like marker events to CD15+ granulocytic-like events"
)

marker_median_features <- tibble(
  feature = c(
    "median_HLA_DR_in_cd45_dump_low", "median_CD14_in_cd45_dump_low",
    "median_CD16_in_cd45_dump_low", "median_CD33_in_cd45_dump_low",
    "median_CD15_in_cd45_dump_low", "median_FceRI_in_cd45_dump_low",
    "median_CD9_in_cd45_dump_low", "median_CD84_in_cd45_dump_low",
    "median_CD206_in_cd45_dump_low", "median_CD169_in_cd45_dump_low",
    "median_HLA_DR_in_cd14_mono_like", "median_CD16_in_cd14_mono_like",
    "median_CD33_in_cd14_mono_like", "median_CD15_in_cd14_mono_like",
    "median_FceRI_in_cd14_mono_like", "median_CD9_in_cd14_mono_like",
    "median_CD84_in_cd14_mono_like", "median_CD206_in_cd14_mono_like",
    "median_CD169_in_cd14_mono_like", "median_CD33_in_hladr_pos",
    "median_CD14_in_hladr_pos", "median_CD15_in_cd33_pos",
    "median_CD206_in_cd33_hladr_myeloid_like",
    "median_CD169_in_cd33_hladr_myeloid_like"
  ),
  module = "marker_medians",
  interpretation = c(
    "HLA-DR median intensity in CD45+ dump-low events",
    "CD14 median intensity in CD45+ dump-low events",
    "CD16 median intensity in CD45+ dump-low events",
    "CD33 median intensity in CD45+ dump-low events",
    "CD15 median intensity in CD45+ dump-low events",
    "FcERI median intensity in CD45+ dump-low events",
    "CD9 median intensity in CD45+ dump-low events",
    "CD84 median intensity in CD45+ dump-low events",
    "CD206 median intensity in CD45+ dump-low events",
    "CD169 median intensity in CD45+ dump-low events",
    "HLA-DR median intensity in CD14+ monocyte-like cells",
    "CD16 median intensity in CD14+ monocyte-like cells",
    "CD33 median intensity in CD14+ monocyte-like cells",
    "CD15 median intensity in CD14+ monocyte-like cells",
    "FcERI median intensity in CD14+ monocyte-like cells",
    "CD9 median intensity in CD14+ monocyte-like cells",
    "CD84 median intensity in CD14+ monocyte-like cells",
    "CD206 median intensity in CD14+ monocyte-like cells",
    "CD169 median intensity in CD14+ monocyte-like cells",
    "CD33 median intensity in HLA-DR+ CD45+ dump-low events",
    "CD14 median intensity in HLA-DR+ CD45+ dump-low events",
    "CD15 median intensity in CD33+ CD45+ dump-low events",
    "CD206 median intensity in CD33+HLA-DR+ myeloid/APC-like events",
    "CD169 median intensity in CD33+HLA-DR+ myeloid/APC-like events"
  )
)

feature_dictionary_full <- bind_rows(feature_dictionary_full, marker_median_features) %>%
  filter(feature %in% names(cp23_features)) %>%
  distinct(feature, .keep_all = TRUE)

# ------------------------------------------------------------
# 9. Summaries
# ------------------------------------------------------------

feature_extraction_summary <- cp23_features %>%
  summarise(
    n_fcs_files = n(),
    n_unique_subjects = n_distinct(subject_id),
    n_feature_ok = sum(feature_ok == TRUE, na.rm = TRUE),
    n_feature_failed = sum(feature_ok == FALSE, na.rm = TRUE),
    n_channel_order_mismatch_files = sum(channel_order_mismatch_file == TRUE, na.rm = TRUE),
    n_marker_order_mismatch_files = sum(marker_order_mismatch_file == TRUE, na.rm = TRUE),
    n_marker_set_mismatch_files = sum(marker_set_mismatch_file == TRUE, na.rm = TRUE),
    median_total_events = median(total_events, na.rm = TRUE),
    median_cd45_dump_low_events = median(n_cd45_dump_low, na.rm = TRUE),
    min_cd45_dump_low_events = min(n_cd45_dump_low, na.rm = TRUE),
    max_cd45_dump_low_events = max(n_cd45_dump_low, na.rm = TRUE),
    median_cd33_hladr_myeloid_like_events = median(n_cd33_hladr_myeloid_like, na.rm = TRUE),
    min_cd33_hladr_myeloid_like_events = min(n_cd33_hladr_myeloid_like, na.rm = TRUE),
    max_cd33_hladr_myeloid_like_events = max(n_cd33_hladr_myeloid_like, na.rm = TRUE)
  )

cd45_event_qc_summary <- cp23_features %>%
  count(cd45_dump_low_event_qc_bin, name = "n") %>%
  arrange(factor(cd45_dump_low_event_qc_bin, levels = c("<50", "50-99", "100-299", "300-999", ">=1000")))

hladr_myeloid_event_qc_summary <- cp23_features %>%
  count(hladr_myeloid_event_qc_bin, name = "n") %>%
  arrange(factor(hladr_myeloid_event_qc_bin, levels = c("<50", "50-99", "100-299", "300-999", ">=1000")))

selected_feature_medians <- cp23_features %>%
  filter(feature_ok == TRUE) %>%
  summarise(
    median_pct_cd45_dump_low_within_total = median(pct_cd45_dump_low_within_total, na.rm = TRUE),
    median_pct_cd14_mono_like_within_cd45_dump_low = median(pct_cd14_mono_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_hladr_low_within_cd45_dump_low = median(pct_cd14_hladr_low_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd33_hladr_myeloid_like_within_cd45_dump_low = median(pct_cd33_hladr_myeloid_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd15_gran_like_within_cd45_dump_low = median(pct_cd15_gran_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_cd206_pos_within_cd45_dump_low = median(pct_cd14_cd206_pos_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_cd169_pos_within_cd45_dump_low = median(pct_cd14_cd169_pos_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_cd9_pos_within_cd45_dump_low = median(pct_cd14_cd9_pos_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_cd84_pos_within_cd45_dump_low = median(pct_cd14_cd84_pos_within_cd45_dump_low, na.rm = TRUE),
    median_pct_fceri_hladr_apc_like_within_cd45_dump_low = median(pct_fceri_hladr_apc_like_within_cd45_dump_low, na.rm = TRUE)
  )

failed_feature_files <- cp23_features %>%
  filter(feature_ok == FALSE) %>%
  select(subject_id, file_name, file_path, feature_error)

# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

write_csv(threshold_table, file.path(out_dir, "SDY2583_CP23_threshold_table_STEP2.csv"))
write_csv(cp23_features, file.path(out_dir, "SDY2583_CP23_features_STEP2.csv"))
write_csv(feature_dictionary_full, file.path(out_dir, "SDY2583_CP23_feature_dictionary_STEP2.csv"))
write_csv(feature_extraction_summary, file.path(out_dir, "SDY2583_CP23_feature_extraction_summary_STEP2.csv"))
write_csv(cd45_event_qc_summary, file.path(out_dir, "SDY2583_CP23_cd45_dump_low_event_QC_summary_STEP2.csv"))
write_csv(hladr_myeloid_event_qc_summary, file.path(out_dir, "SDY2583_CP23_hladr_myeloid_event_QC_summary_STEP2.csv"))
write_csv(selected_feature_medians, file.path(out_dir, "SDY2583_CP23_selected_feature_medians_STEP2.csv"))
write_csv(failed_feature_files, file.path(out_dir, "SDY2583_CP23_failed_feature_files_STEP2.csv"))

save(
  thresholds,
  threshold_table,
  cp23_features,
  feature_dictionary_full,
  feature_extraction_summary,
  cd45_event_qc_summary,
  hladr_myeloid_event_qc_summary,
  selected_feature_medians,
  failed_feature_files,
  file = file.path(rdata_dir, "SDY2583_CP23_STEP2_feature_extraction_myeloid_macrophage_like.RData")
)

# ------------------------------------------------------------
# 11. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP23 STEP 2 COMPLETE: FEATURE EXTRACTION\n")
cat("============================================================\n")

cat("\nFeature extraction summary:\n")
print(as.data.frame(feature_extraction_summary), row.names = FALSE)

cat("\nEvent QC summary - CD45 dump-low denominator:\n")
print(as.data.frame(cd45_event_qc_summary), row.names = FALSE)

cat("\nEvent QC summary - CD33+HLA-DR+ myeloid/APC-like denominator:\n")
print(as.data.frame(hladr_myeloid_event_qc_summary), row.names = FALSE)

cat("\nSelected feature medians:\n")
print(as.data.frame(selected_feature_medians), row.names = FALSE)

cat("\nFeature dictionary:\n")
print(as.data.frame(feature_dictionary_full), row.names = FALSE)

cat("\nFailed feature files:\n")
print(as.data.frame(failed_feature_files), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
