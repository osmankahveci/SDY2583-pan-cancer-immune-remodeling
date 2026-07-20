# ============================================================
# SDY2583 CP16
# STEP 2 SAFE: Feature extraction for dendritic/APC-like myeloid panel
#
# Panel markers:
#   CD45
#   Dump/viability: Viability_CD15_CD3_CD19_CCR3_CD7
#   HLA-DR, CD11c, CD14, CD16, CD1c, CD141, CLEC9A,
#   CD123, FceRI, CD13
#
# Critical:
#   - Marker-name tolerant extraction is used.
#   - Do NOT use fixed channel indices.
#   - DBG161/DBG177 have Time at parameter 1 and shifted channels;
#     they are retained because marker set is identical.
#
# Interpretation:
#   - This is a phenotype-like APC/DC/myeloid panel.
#   - Do not claim definitive cDC1, cDC2, pDC, monocytes, or MDSCs.
#   - Use "cDC1-like", "cDC2-like", "pDC-like",
#     "monocyte-like", "APC-like myeloid".
#
# Output:
#   outputs/CP16/02_feature_extraction
#   outputs/CP16/11_RData
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
# 2. Paths and load Step 1
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "02_feature_extraction")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP1_fcs_inventory_marker_QC.RData")

if (!file.exists(step1_rdata)) {
  stop("Step 1 RData bulunamadı: ", step1_rdata)
}

load(step1_rdata)

if (!exists("fcs_files")) stop("fcs_files bulunamadı.")
if (!exists("fcs_inventory")) stop("fcs_inventory bulunamadı.")

# ------------------------------------------------------------
# 3. Main thresholds
# ------------------------------------------------------------
# These are initial transformed-intensity gates.
# Threshold sensitivity will be performed later.

thresholds <- list(
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

threshold_table <- tibble(
  marker_or_gate = names(thresholds),
  threshold = unlist(thresholds)
)

write_csv(threshold_table, file.path(out_dir, "SDY2583_CP16_thresholds_STEP2.csv"))

# ------------------------------------------------------------
# 4. Helper functions
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
  ch$DUMP <- find_marker_channel_regex(map, "Viability.*CD15.*CD3.*CD19.*CCR3.*CD7|Viability|Dump")
  ch$HLA_DR <- find_marker_channel_regex(map, "^HLA-DR$|HLA_DR|HLADR")
  ch$CD11c <- find_marker_channel_exact(map, "CD11c")
  ch$CD14 <- find_marker_channel_exact(map, "CD14")
  ch$CD16 <- find_marker_channel_exact(map, "CD16")
  ch$CD1c <- find_marker_channel_exact(map, "CD1c")
  ch$CD141 <- find_marker_channel_exact(map, "CD141")
  ch$CLEC9A <- find_marker_channel_exact(map, "CLEC9A")
  ch$CD123 <- find_marker_channel_exact(map, "CD123")
  ch$FceRI <- find_marker_channel_regex(map, "FceRI|FcER1|FcERI|FCER1|FcERI")
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

# ------------------------------------------------------------
# 5. Feature extraction from one FCS
# ------------------------------------------------------------

extract_one_cp16 <- function(fp) {

  tryCatch({

    ff_raw <- flowCore::read.FCS(
      fp,
      transformation = FALSE,
      truncate_max_range = FALSE
    )

    marker_map_file <- get_marker_map(ff_raw)
    ch <- build_channel_map(marker_map_file)

    required <- c(
      "DUMP", "CD45", "HLA_DR", "CD11c", "CD14", "CD16",
      "CD1c", "CD141", "CLEC9A", "CD123", "FceRI", "CD13"
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
    cd11c <- get_vec(ex, ch$CD11c)
    cd14 <- get_vec(ex, ch$CD14)
    cd16 <- get_vec(ex, ch$CD16)
    cd1c <- get_vec(ex, ch$CD1c)
    cd141 <- get_vec(ex, ch$CD141)
    clec9a <- get_vec(ex, ch$CLEC9A)
    cd123 <- get_vec(ex, ch$CD123)
    fceri <- get_vec(ex, ch$FceRI)
    cd13 <- get_vec(ex, ch$CD13)

    dump_low <- dump <= thresholds$DUMP_LOW
    primary <- dump_low & cd45 > thresholds$CD45

    # APC-like core
    hladr_core <- primary & hla_dr > thresholds$HLA_DR
    cd11c_hladr_apc_like <- hladr_core & cd11c > thresholds$CD11c

    # Monocyte-like phenotypes
    cd14_mono_like <- primary & cd14 > thresholds$CD14
    cd14_hladr_pos <- cd14_mono_like & hla_dr > thresholds$HLA_DR
    cd14_hladr_low <- cd14_mono_like & hla_dr <= thresholds$HLA_DR
    cd14pos_cd16neg_mono_like <- cd14_mono_like & cd16 <= thresholds$CD16
    cd14pos_cd16pos_mono_like <- cd14_mono_like & cd16 > thresholds$CD16
    cd14lowneg_cd16pos_hladr_like <- primary & cd14 <= thresholds$CD14 & cd16 > thresholds$CD16 & hla_dr > thresholds$HLA_DR

    # cDC2-like / CD1c-associated APC phenotypes
    cd1c_apc_like <- hladr_core & cd11c > thresholds$CD11c & cd1c > thresholds$CD1c
    cd1c_cd14low_cdc2_like <- cd1c_apc_like & cd14 <= thresholds$CD14
    cd1c_fceri_apc_like <- cd1c_apc_like & fceri > thresholds$FceRI
    fceri_pos_hladr <- hladr_core & fceri > thresholds$FceRI

    # cDC1-like / CD141-CLEC9A-associated APC phenotypes
    cd141_pos_hladr <- hladr_core & cd141 > thresholds$CD141
    clec9a_pos_hladr <- hladr_core & clec9a > thresholds$CLEC9A
    cd141_clec9a_cdc1_like <- hladr_core &
      cd11c > thresholds$CD11c &
      cd141 > thresholds$CD141 &
      clec9a > thresholds$CLEC9A &
      cd14 <= thresholds$CD14

    # pDC-like / CD123-associated phenotypes
    cd123_hladr_dc_like <- hladr_core & cd123 > thresholds$CD123
    cd123_pdc_like <- hladr_core &
      cd123 > thresholds$CD123 &
      cd11c <= thresholds$CD11c &
      cd14 <= thresholds$CD14 &
      cd16 <= thresholds$CD16
    cd123_cd11c_mixed_apc_like <- hladr_core &
      cd123 > thresholds$CD123 &
      cd11c > thresholds$CD11c

    # CD13-associated myeloid/APC phenotype
    cd13_pos_primary <- primary & cd13 > thresholds$CD13
    cd13_pos_hladr <- hladr_core & cd13 > thresholds$CD13
    cd13_cd11c_hladr_apc_like <- hladr_core & cd13 > thresholds$CD13 & cd11c > thresholds$CD11c

    # Counts
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

      # primary composition
      pct_dump_low_within_total = safe_pct(n_dump_low, total_events),
      pct_cd45_dump_low_within_total = safe_pct(n_primary, total_events),
      pct_hladr_apc_core_within_total = safe_pct(n_hladr_core, total_events),
      pct_hladr_apc_core_within_cd45_dump_low = safe_pct(n_hladr_core, n_primary),
      pct_cd11c_hladr_apc_like_within_cd45_dump_low = safe_pct(n_cd11c_hladr, n_primary),
      pct_cd11c_hladr_apc_like_within_hladr_core = safe_pct(n_cd11c_hladr, n_hladr_core),

      # monocyte-like phenotypes
      pct_cd14_mono_like_within_cd45_dump_low = safe_pct(n_cd14_mono, n_primary),
      pct_cd14_hladr_pos_within_cd45_dump_low = safe_pct(n_cd14_hladr_pos, n_primary),
      pct_cd14_hladr_low_within_cd45_dump_low = safe_pct(n_cd14_hladr_low, n_primary),
      pct_hladr_pos_within_cd14_mono_like = safe_pct(n_cd14_hladr_pos, n_cd14_mono),
      pct_hladr_low_within_cd14_mono_like = safe_pct(n_cd14_hladr_low, n_cd14_mono),
      pct_cd14pos_cd16neg_mono_like_within_cd45_dump_low = safe_pct(n_cd14pos_cd16neg, n_primary),
      pct_cd14pos_cd16pos_mono_like_within_cd45_dump_low = safe_pct(n_cd14pos_cd16pos, n_primary),
      pct_cd14lowneg_cd16pos_hladr_like_within_cd45_dump_low = safe_pct(n_cd14lowneg_cd16pos, n_primary),
      pct_cd16_pos_within_cd14_mono_like = safe_pct(n_cd14pos_cd16pos, n_cd14_mono),

      # cDC2-like / CD1c-associated APC phenotypes
      pct_cd1c_apc_like_within_cd45_dump_low = safe_pct(n_cd1c_apc, n_primary),
      pct_cd1c_apc_like_within_hladr_core = safe_pct(n_cd1c_apc, n_hladr_core),
      pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low = safe_pct(n_cd1c_cd14low_cdc2, n_primary),
      pct_cd1c_cd14low_cdc2_like_within_hladr_core = safe_pct(n_cd1c_cd14low_cdc2, n_hladr_core),
      pct_cd1c_fceri_apc_like_within_cd45_dump_low = safe_pct(n_cd1c_fceri_apc, n_primary),
      pct_cd1c_fceri_apc_like_within_cd1c_apc_like = safe_pct(n_cd1c_fceri_apc, n_cd1c_apc),
      pct_fceri_pos_hladr_within_cd45_dump_low = safe_pct(n_fceri_pos_hladr, n_primary),

      # cDC1-like / CD141-CLEC9A phenotypes
      pct_cd141_pos_hladr_within_cd45_dump_low = safe_pct(n_cd141_pos_hladr, n_primary),
      pct_clec9a_pos_hladr_within_cd45_dump_low = safe_pct(n_clec9a_pos_hladr, n_primary),
      pct_cd141_clec9a_cdc1_like_within_cd45_dump_low = safe_pct(n_cd141_clec9a_cdc1, n_primary),
      pct_cd141_clec9a_cdc1_like_within_hladr_core = safe_pct(n_cd141_clec9a_cdc1, n_hladr_core),

      # pDC-like / CD123-associated phenotypes
      pct_cd123_hladr_dc_like_within_cd45_dump_low = safe_pct(n_cd123_hladr, n_primary),
      pct_cd123_hladr_dc_like_within_hladr_core = safe_pct(n_cd123_hladr, n_hladr_core),
      pct_cd123_pdc_like_within_cd45_dump_low = safe_pct(n_cd123_pdc, n_primary),
      pct_cd123_pdc_like_within_hladr_core = safe_pct(n_cd123_pdc, n_hladr_core),
      pct_cd123_cd11c_mixed_apc_like_within_cd45_dump_low = safe_pct(n_cd123_cd11c_mixed, n_primary),

      # CD13-associated myeloid/APC phenotypes
      pct_cd13_pos_within_cd45_dump_low = safe_pct(n_cd13_primary, n_primary),
      pct_cd13_pos_within_hladr_core = safe_pct(n_cd13_hladr, n_hladr_core),
      pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low = safe_pct(n_cd13_cd11c_hladr, n_primary),

      # ratios
      ratio_cd11c_hladr_apc_to_cd14_mono_like = safe_ratio(n_cd11c_hladr, n_cd14_mono),
      ratio_cd1c_cdc2_like_to_cd123_pdc_like = safe_ratio(n_cd1c_cd14low_cdc2, n_cd123_pdc),
      ratio_cd141_clec9a_cdc1_like_to_cd123_pdc_like = safe_ratio(n_cd141_clec9a_cdc1, n_cd123_pdc),
      ratio_cd11c_hladr_apc_to_cd123_pdc_like = safe_ratio(n_cd11c_hladr, n_cd123_pdc),
      ratio_cd14_mono_like_to_cd123_pdc_like = safe_ratio(n_cd14_mono, n_cd123_pdc),

      # medians in denominators
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
cat("SDY2583 CP16 STEP 2 STARTED: FEATURE EXTRACTION\n")
cat("============================================================\n")

feature_rows <- vector("list", length(fcs_files))

for (i in seq_along(fcs_files)) {
  if (i %% 50 == 0) cat("Processed", i, "of", length(fcs_files), "FCS files\n")
  feature_rows[[i]] <- extract_one_cp16(fcs_files[i])
}

cp16_features <- bind_rows(feature_rows)

# Add technical mismatch flags from Step 1.
technical_flags <- fcs_inventory %>%
  select(
    subject_id,
    channel_order_mismatch_file,
    marker_order_mismatch_file,
    marker_set_mismatch_file,
    has_spillover_keyword,
    spillover_keyword_name
  ) %>%
  distinct(subject_id, .keep_all = TRUE)

cp16_features <- cp16_features %>%
  left_join(technical_flags, by = "subject_id")

# Event QC bins
cp16_features <- cp16_features %>%
  mutate(
    cd45_dump_low_event_qc_bin = case_when(
      is.na(n_cd45_dump_low) ~ "missing",
      n_cd45_dump_low < 50 ~ "<50",
      n_cd45_dump_low < 100 ~ "50-99",
      n_cd45_dump_low < 300 ~ "100-299",
      n_cd45_dump_low < 1000 ~ "300-999",
      TRUE ~ ">=1000"
    ),
    hladr_apc_core_event_qc_bin = case_when(
      is.na(n_hladr_apc_core) ~ "missing",
      n_hladr_apc_core < 50 ~ "<50",
      n_hladr_apc_core < 100 ~ "50-99",
      n_hladr_apc_core < 300 ~ "100-299",
      n_hladr_apc_core < 1000 ~ "300-999",
      TRUE ~ ">=1000"
    )
  )

# ------------------------------------------------------------
# 7. Feature dictionary
# ------------------------------------------------------------

feature_dictionary <- tibble::tribble(
  ~feature, ~module, ~interpretation,

  "pct_cd45_dump_low_within_total", "composition", "CD45+ dump-low APC-enriched fraction among total events",
  "pct_hladr_apc_core_within_cd45_dump_low", "APC_core", "HLA-DR+ APC-like core fraction within CD45+ dump-low events",
  "pct_cd11c_hladr_apc_like_within_cd45_dump_low", "APC_core", "CD11c+HLA-DR+ APC-like myeloid fraction within CD45+ dump-low events",
  "pct_cd11c_hladr_apc_like_within_hladr_core", "APC_core", "CD11c+ fraction within HLA-DR+ APC-like core",

  "pct_cd14_mono_like_within_cd45_dump_low", "monocyte_like_axis", "CD14+ monocyte-like fraction within CD45+ dump-low events",
  "pct_cd14_hladr_pos_within_cd45_dump_low", "monocyte_like_axis", "CD14+HLA-DR+ monocyte-like fraction within CD45+ dump-low events",
  "pct_cd14_hladr_low_within_cd45_dump_low", "monocyte_like_axis", "CD14+HLA-DR-low monocyte-like fraction within CD45+ dump-low events",
  "pct_hladr_pos_within_cd14_mono_like", "monocyte_like_axis", "HLA-DR+ fraction within CD14+ monocyte-like cells",
  "pct_hladr_low_within_cd14_mono_like", "monocyte_like_axis", "HLA-DR-low fraction within CD14+ monocyte-like cells",
  "pct_cd14pos_cd16neg_mono_like_within_cd45_dump_low", "monocyte_like_axis", "CD14+CD16- monocyte-like fraction within CD45+ dump-low events",
  "pct_cd14pos_cd16pos_mono_like_within_cd45_dump_low", "monocyte_like_axis", "CD14+CD16+ monocyte-like fraction within CD45+ dump-low events",
  "pct_cd14lowneg_cd16pos_hladr_like_within_cd45_dump_low", "monocyte_like_axis", "CD14-low/negative CD16+HLA-DR+ monocyte-like fraction within CD45+ dump-low events",
  "pct_cd16_pos_within_cd14_mono_like", "monocyte_like_axis", "CD16+ fraction within CD14+ monocyte-like cells",

  "pct_cd1c_apc_like_within_cd45_dump_low", "CD1c_cDC2_like_axis", "CD1c+CD11c+HLA-DR+ APC-like fraction within CD45+ dump-low events",
  "pct_cd1c_apc_like_within_hladr_core", "CD1c_cDC2_like_axis", "CD1c+CD11c+ fraction within HLA-DR+ APC-like core",
  "pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low", "CD1c_cDC2_like_axis", "CD1c+CD11c+HLA-DR+CD14-low cDC2-like fraction within CD45+ dump-low events",
  "pct_cd1c_cd14low_cdc2_like_within_hladr_core", "CD1c_cDC2_like_axis", "CD1c+CD11c+CD14-low cDC2-like fraction within HLA-DR+ APC-like core",
  "pct_cd1c_fceri_apc_like_within_cd45_dump_low", "CD1c_cDC2_like_axis", "CD1c+FcERI+ APC-like fraction within CD45+ dump-low events",
  "pct_cd1c_fceri_apc_like_within_cd1c_apc_like", "CD1c_cDC2_like_axis", "FcERI+ fraction within CD1c+ APC-like cells",
  "pct_fceri_pos_hladr_within_cd45_dump_low", "CD1c_cDC2_like_axis", "FcERI+HLA-DR+ fraction within CD45+ dump-low events",

  "pct_cd141_pos_hladr_within_cd45_dump_low", "CD141_CLEC9A_cDC1_like_axis", "CD141+HLA-DR+ fraction within CD45+ dump-low events",
  "pct_clec9a_pos_hladr_within_cd45_dump_low", "CD141_CLEC9A_cDC1_like_axis", "CLEC9A+HLA-DR+ fraction within CD45+ dump-low events",
  "pct_cd141_clec9a_cdc1_like_within_cd45_dump_low", "CD141_CLEC9A_cDC1_like_axis", "CD141+CLEC9A+CD11c+HLA-DR+CD14-low cDC1-like fraction within CD45+ dump-low events",
  "pct_cd141_clec9a_cdc1_like_within_hladr_core", "CD141_CLEC9A_cDC1_like_axis", "CD141+CLEC9A+ cDC1-like fraction within HLA-DR+ APC-like core",

  "pct_cd123_hladr_dc_like_within_cd45_dump_low", "CD123_pDC_like_axis", "CD123+HLA-DR+ DC-like fraction within CD45+ dump-low events",
  "pct_cd123_hladr_dc_like_within_hladr_core", "CD123_pDC_like_axis", "CD123+ fraction within HLA-DR+ APC-like core",
  "pct_cd123_pdc_like_within_cd45_dump_low", "CD123_pDC_like_axis", "CD123+HLA-DR+CD11c-lowCD14-lowCD16-low pDC-like fraction within CD45+ dump-low events",
  "pct_cd123_pdc_like_within_hladr_core", "CD123_pDC_like_axis", "pDC-like fraction within HLA-DR+ APC-like core",
  "pct_cd123_cd11c_mixed_apc_like_within_cd45_dump_low", "CD123_pDC_like_axis", "CD123+CD11c+HLA-DR+ mixed APC-like fraction within CD45+ dump-low events",

  "pct_cd13_pos_within_cd45_dump_low", "CD13_myeloid_axis", "CD13+ fraction within CD45+ dump-low events",
  "pct_cd13_pos_within_hladr_core", "CD13_myeloid_axis", "CD13+ fraction within HLA-DR+ APC-like core",
  "pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low", "CD13_myeloid_axis", "CD13+CD11c+HLA-DR+ APC-like myeloid fraction within CD45+ dump-low events",

  "ratio_cd11c_hladr_apc_to_cd14_mono_like", "ratio_balance", "Ratio of CD11c+HLA-DR+ APC-like to CD14+ monocyte-like events",
  "ratio_cd1c_cdc2_like_to_cd123_pdc_like", "ratio_balance", "Ratio of CD1c+ cDC2-like to CD123+ pDC-like events",
  "ratio_cd141_clec9a_cdc1_like_to_cd123_pdc_like", "ratio_balance", "Ratio of CD141+CLEC9A+ cDC1-like to CD123+ pDC-like events",
  "ratio_cd11c_hladr_apc_to_cd123_pdc_like", "ratio_balance", "Ratio of CD11c+HLA-DR+ APC-like to CD123+ pDC-like events",
  "ratio_cd14_mono_like_to_cd123_pdc_like", "ratio_balance", "Ratio of CD14+ monocyte-like to CD123+ pDC-like events"
)

# Keep only features that exist in extraction output.
feature_dictionary <- feature_dictionary %>%
  filter(feature %in% names(cp16_features))

# ------------------------------------------------------------
# 8. Summaries
# ------------------------------------------------------------

feature_extraction_summary <- tibble(
  n_fcs_files = length(fcs_files),
  n_unique_subjects = n_distinct(cp16_features$subject_id),
  n_feature_ok = sum(cp16_features$feature_ok == TRUE, na.rm = TRUE),
  n_feature_failed = sum(cp16_features$feature_ok == FALSE, na.rm = TRUE),
  n_channel_order_mismatch_files = sum(cp16_features$channel_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_order_mismatch_files = sum(cp16_features$marker_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_set_mismatch_files = sum(cp16_features$marker_set_mismatch_file == TRUE, na.rm = TRUE),
  median_total_events = median(cp16_features$total_events, na.rm = TRUE),
  median_cd45_dump_low_events = median(cp16_features$n_cd45_dump_low, na.rm = TRUE),
  median_hladr_apc_core_events = median(cp16_features$n_hladr_apc_core, na.rm = TRUE),
  min_hladr_apc_core_events = min(cp16_features$n_hladr_apc_core, na.rm = TRUE),
  max_hladr_apc_core_events = max(cp16_features$n_hladr_apc_core, na.rm = TRUE)
)

event_qc_summary_cd45_dump_low <- cp16_features %>%
  count(cd45_dump_low_event_qc_bin, name = "n") %>%
  arrange(factor(cd45_dump_low_event_qc_bin, levels = c("<50", "50-99", "100-299", "300-999", ">=1000", "missing")))

event_qc_summary_hladr_core <- cp16_features %>%
  count(hladr_apc_core_event_qc_bin, name = "n") %>%
  arrange(factor(hladr_apc_core_event_qc_bin, levels = c("<50", "50-99", "100-299", "300-999", ">=1000", "missing")))

selected_feature_medians <- cp16_features %>%
  filter(feature_ok == TRUE) %>%
  summarise(
    median_pct_cd45_dump_low_within_total = median(pct_cd45_dump_low_within_total, na.rm = TRUE),
    median_pct_hladr_apc_core_within_cd45_dump_low = median(pct_hladr_apc_core_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd11c_hladr_apc_like_within_cd45_dump_low = median(pct_cd11c_hladr_apc_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd14_mono_like_within_cd45_dump_low = median(pct_cd14_mono_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low = median(pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd141_clec9a_cdc1_like_within_cd45_dump_low = median(pct_cd141_clec9a_cdc1_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd123_pdc_like_within_cd45_dump_low = median(pct_cd123_pdc_like_within_cd45_dump_low, na.rm = TRUE),
    median_pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low = median(pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low, na.rm = TRUE)
  )

failed_feature_files <- cp16_features %>%
  filter(feature_ok == FALSE) %>%
  select(subject_id, file_name, file_path, feature_error)

# ------------------------------------------------------------
# 9. Save outputs
# ------------------------------------------------------------

write_csv(cp16_features, file.path(out_dir, "SDY2583_CP16_features_STEP2.csv"))
write_csv(feature_dictionary, file.path(out_dir, "SDY2583_CP16_feature_dictionary_STEP2.csv"))
write_csv(feature_extraction_summary, file.path(out_dir, "SDY2583_CP16_feature_extraction_summary_STEP2.csv"))
write_csv(event_qc_summary_cd45_dump_low, file.path(out_dir, "SDY2583_CP16_event_QC_summary_CD45_dump_low_STEP2.csv"))
write_csv(event_qc_summary_hladr_core, file.path(out_dir, "SDY2583_CP16_event_QC_summary_HLADR_core_STEP2.csv"))
write_csv(selected_feature_medians, file.path(out_dir, "SDY2583_CP16_selected_feature_medians_STEP2.csv"))
write_csv(failed_feature_files, file.path(out_dir, "SDY2583_CP16_failed_feature_files_STEP2.csv"))

save(
  cp16_features,
  feature_dictionary,
  thresholds,
  threshold_table,
  feature_extraction_summary,
  event_qc_summary_cd45_dump_low,
  event_qc_summary_hladr_core,
  selected_feature_medians,
  failed_feature_files,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP2_feature_extraction_APC_DC_myeloid.RData")
)

# ------------------------------------------------------------
# 10. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 2 COMPLETE: FEATURE EXTRACTION\n")
cat("============================================================\n")

cat("\nFeature extraction summary:\n")
print(as.data.frame(feature_extraction_summary), row.names = FALSE)

cat("\nEvent QC summary - CD45 dump-low denominator:\n")
print(as.data.frame(event_qc_summary_cd45_dump_low), row.names = FALSE)

cat("\nEvent QC summary - HLA-DR APC-like core denominator:\n")
print(as.data.frame(event_qc_summary_hladr_core), row.names = FALSE)

cat("\nSelected feature medians:\n")
print(as.data.frame(selected_feature_medians), row.names = FALSE)

cat("\nFeature dictionary:\n")
print(as.data.frame(feature_dictionary), row.names = FALSE)

cat("\nFailed feature files:\n")
print(as.data.frame(failed_feature_files), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)

cat("============================================================\n")
