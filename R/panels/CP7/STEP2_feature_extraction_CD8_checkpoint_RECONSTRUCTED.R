# ============================================================
# SDY2583 CP7
# STEP 2 RECONSTRUCTED: CD8 differentiation/checkpoint features
#
# Reconstruction basis:
#   - archived CP7 threshold table;
#   - archived 850-subject feature-table schema;
#   - archived marker-channel map and Methods record;
#   - shared compensation/logicle conventions in recovered scripts.
#
# This is reconstructed source, not the original archived script.
# Validate all generated columns and summary values against the archived
# SDY2583_CP7_FULL_850_feature_table_STEP2.csv before release.
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("flowCore", quietly = TRUE)) BiocManager::install("flowCore", ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(flowCore)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
})

analysis_dir <- sd_analysis_dir("CP7")
out_dir <- file.path(analysis_dir, "02_feature_extraction")
threshold_dir <- file.path(analysis_dir, "03_thresholds")
rdata_dir <- file.path(analysis_dir, "11_RData")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step1_file <- file.path(rdata_dir, "SDY2583_CP7_STEP1_fcs_inventory_marker_QC_RECONSTRUCTED.RData")
if (!file.exists(step1_file)) stop("Run reconstructed CP7 Step 1 first: ", step1_file)
load(step1_file)
if (!exists("fcs_files")) stop("Step 1 RData does not contain fcs_files.")

thresholds <- c(
  CD3 = 2.0,
  CD8 = 2.2,
  CD45RA = 2.0,
  CD62L = 2.2,
  CD27 = 2.0,
  PD1_pos = 1.5,
  PD1_high = 2.0,
  TIM3_pos = 1.5,
  LAG3_pos = 1.5,
  TIGIT_pos = 1.5,
  ICOS_pos = 1.5,
  CD39_pos = 1.5
)

threshold_table <- tibble(
  feature = names(thresholds),
  channel = c(
    "BV605-A", "PerCP-Cy5-5-A", "PE-Cy5-A", "BV650-A", "BV711-A",
    "PE-CF594-A", "PE-CF594-A", "BV421-A", "BV786-A", "PE-A", "BB515-A", "PE-Cy7-A"
  ),
  threshold = as.numeric(thresholds),
  interpretation = c(
    "CD3 positive", "CD8 positive", "CD45RA positive", "CD62L positive", "CD27 positive",
    "PD-1 positive", "PD-1 high", "TIM-3 positive", "LAG-3 positive", "TIGIT positive",
    "ICOS positive", "CD39 positive"
  )
)
readr::write_csv(threshold_table, file.path(threshold_dir, "SDY2583_CP7_threshold_definitions_STEP2.csv"))

extract_subject_id <- function(path) {
  id <- stringr::str_extract(basename(path), "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(basename(path), "\\.fcs$"), id)
}

get_marker_map <- function(ff) {
  pp <- Biobase::pData(flowCore::parameters(ff))
  channel <- as.character(pp$name)
  marker <- as.character(pp$desc)
  marker[is.na(marker) | marker == ""] <- channel[is.na(marker) | marker == ""]
  tibble(channel = channel, marker = marker)
}

find_marker <- function(map, pattern) {
  idx <- which(grepl(pattern, map$marker, ignore.case = TRUE))
  if (length(idx) == 0L) NA_character_ else map$channel[idx[[1L]]]
}

build_channel_map <- function(map) {
  c(
    CD3 = find_marker(map, "^CD3$"),
    CD8 = find_marker(map, "^CD8$"),
    CD45RA = find_marker(map, "CD45RA"),
    CD62L = find_marker(map, "CD62L"),
    CD27 = find_marker(map, "^CD27$"),
    PD1 = find_marker(map, "PD-?1"),
    TIM3 = find_marker(map, "TIM-?3"),
    LAG3 = find_marker(map, "LAG-?3"),
    TIGIT = find_marker(map, "TIGIT"),
    ICOS = find_marker(map, "ICOS"),
    CD39 = find_marker(map, "CD39")
  )
}

get_spill_matrix <- function(ff) {
  kk <- flowCore::keyword(ff)
  for (nm in c("SPILL", "$SPILLOVER", "SPILLOVER")) {
    if (nm %in% names(kk) && is.matrix(kk[[nm]])) return(kk[[nm]])
  }
  NULL
}

apply_compensation <- function(ff) {
  sp <- get_spill_matrix(ff)
  if (is.null(sp)) return(ff)
  tryCatch(flowCore::compensate(ff, sp), error = function(e) ff)
}

apply_fixed_logicle <- function(ff) {
  channels <- colnames(flowCore::exprs(ff))
  fluorescence <- channels[!grepl("FSC|SSC|Time", channels, ignore.case = TRUE)]
  transform <- flowCore::logicleTransform(
    transformationId = "fixed_logicle", w = 0.5, t = 262144, m = 4.5, a = 0
  )
  flowCore::transform(ff, flowCore::transformList(fluorescence, transform))
}

safe_pct <- function(mask, denominator_mask = NULL) {
  if (is.null(denominator_mask)) denominator_mask <- rep(TRUE, length(mask))
  den <- sum(denominator_mask, na.rm = TRUE)
  if (den <= 0L) return(NA_real_)
  100 * sum(mask & denominator_mask, na.rm = TRUE) / den
}

safe_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

extract_one <- function(path) {
  tryCatch({
    ff_raw <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
    map <- get_marker_map(ff_raw)
    ch <- build_channel_map(map)
    if (any(is.na(ch))) stop("Missing markers: ", paste(names(ch)[is.na(ch)], collapse = ", "))

    ff <- apply_fixed_logicle(apply_compensation(ff_raw))
    ex <- flowCore::exprs(ff)
    v <- lapply(ch, function(channel) as.numeric(ex[, channel]))

    cd3 <- v$CD3 > thresholds[["CD3"]]
    cd8 <- v$CD8 > thresholds[["CD8"]]
    cd3cd8 <- cd3 & cd8
    cd45ra <- v$CD45RA > thresholds[["CD45RA"]]
    cd62l <- v$CD62L > thresholds[["CD62L"]]
    cd27 <- v$CD27 > thresholds[["CD27"]]
    pd1 <- v$PD1 > thresholds[["PD1_pos"]]
    pd1_high <- v$PD1 > thresholds[["PD1_high"]]
    tim3 <- v$TIM3 > thresholds[["TIM3_pos"]]
    lag3 <- v$LAG3 > thresholds[["LAG3_pos"]]
    tigit <- v$TIGIT > thresholds[["TIGIT_pos"]]
    icos <- v$ICOS > thresholds[["ICOS_pos"]]
    cd39 <- v$CD39 > thresholds[["CD39_pos"]]

    naive <- cd45ra & cd62l
    tcm <- !cd45ra & cd62l
    tem <- !cd45ra & !cd62l
    temra <- cd45ra & !cd62l
    temra_cd27neg <- temra & !cd27
    temra_cd27pos <- temra & cd27
    cd62lneg_cd27neg <- !cd62l & !cd27
    cd45ra_cd62lneg_cd27neg <- cd45ra & !cd62l & !cd27

    n_cd3cd8 <- sum(cd3cd8)
    row <- tibble(
      file_name = basename(path),
      file_path = normalizePath(path, mustWork = FALSE),
      subject_id = extract_subject_id(path),
      panel = "CP7",
      feature_ok = TRUE,
      error_message = NA_character_,
      total_events = nrow(ex),
      n_channels = ncol(ex),
      n_cd3_pos = sum(cd3),
      n_cd8_pos = sum(cd8),
      n_cd3_cd8_pos = n_cd3cd8,
      pct_cd3_pos_total = safe_pct(cd3),
      pct_cd8_pos_total = safe_pct(cd8),
      pct_cd3_cd8_pos_total = safe_pct(cd3cd8),
      pct_cd8_within_cd3 = safe_pct(cd8, cd3),
      median_CD62L_in_CD3CD8 = safe_median(v$CD62L[cd3cd8]),
      median_CD27_in_CD3CD8 = safe_median(v$CD27[cd3cd8]),
      median_PD1_in_CD3CD8 = safe_median(v$PD1[cd3cd8]),
      median_TIM3_in_CD3CD8 = safe_median(v$TIM3[cd3cd8]),
      median_LAG3_in_CD3CD8 = safe_median(v$LAG3[cd3cd8]),
      median_TIGIT_in_CD3CD8 = safe_median(v$TIGIT[cd3cd8]),
      median_ICOS_in_CD3CD8 = safe_median(v$ICOS[cd3cd8]),
      median_CD39_in_CD3CD8 = safe_median(v$CD39[cd3cd8]),
      median_CD45RA_in_CD3CD8 = safe_median(v$CD45RA[cd3cd8]),
      pct_naive_like = safe_pct(naive, cd3cd8),
      pct_tcm_like = safe_pct(tcm, cd3cd8),
      pct_tem_like = safe_pct(tem, cd3cd8),
      pct_temra_like = safe_pct(temra, cd3cd8),
      pct_cd62l_pos_within_cd3cd8 = safe_pct(cd62l, cd3cd8),
      pct_cd62l_neg_within_cd3cd8 = safe_pct(!cd62l, cd3cd8),
      pct_cd45ra_pos_within_cd3cd8 = safe_pct(cd45ra, cd3cd8),
      pct_cd45ra_neg_within_cd3cd8 = safe_pct(!cd45ra, cd3cd8),
      pct_cd27_pos_within_cd3cd8 = safe_pct(cd27, cd3cd8),
      pct_cd27_neg_within_cd3cd8 = safe_pct(!cd27, cd3cd8),
      pct_temra_cd27neg = safe_pct(temra_cd27neg, cd3cd8),
      pct_temra_cd27pos = safe_pct(temra_cd27pos, cd3cd8),
      pct_cd62lneg_cd27neg = safe_pct(cd62lneg_cd27neg, cd3cd8),
      pct_cd45ra_pos_cd62lneg_cd27neg = safe_pct(cd45ra_cd62lneg_cd27neg, cd3cd8),
      pct_naive_cd27pos = safe_pct(naive & cd27, cd3cd8),
      pct_naive_cd27neg = safe_pct(naive & !cd27, cd3cd8),
      pct_pd1_pos = safe_pct(pd1, cd3cd8),
      pct_pd1_high = safe_pct(pd1_high, cd3cd8),
      pct_tim3_pos = safe_pct(tim3, cd3cd8),
      pct_lag3_pos = safe_pct(lag3, cd3cd8),
      pct_tigit_pos = safe_pct(tigit, cd3cd8),
      pct_icos_pos = safe_pct(icos, cd3cd8),
      pct_cd39_pos = safe_pct(cd39, cd3cd8),
      pct_pd1_tim3_pos = safe_pct(pd1 & tim3, cd3cd8),
      pct_pd1_lag3_pos = safe_pct(pd1 & lag3, cd3cd8),
      pct_pd1_tigit_pos = safe_pct(pd1 & tigit, cd3cd8),
      pct_pd1_cd39_pos = safe_pct(pd1 & cd39, cd3cd8),
      pct_tim3_lag3_pos = safe_pct(tim3 & lag3, cd3cd8),
      pct_tigit_cd39_pos = safe_pct(tigit & cd39, cd3cd8),
      pct_icos_cd39_pos = safe_pct(icos & cd39, cd3cd8),
      pct_pd1_tigit_cd39_pos = safe_pct(pd1 & tigit & cd39, cd3cd8),
      pct_pd1_tim3_lag3_pos = safe_pct(pd1 & tim3 & lag3, cd3cd8),
      pct_pd1_tigit_icos_cd39_pos = safe_pct(pd1 & tigit & icos & cd39, cd3cd8),
      pct_pd1_pos_within_temra_like = safe_pct(pd1, cd3cd8 & temra),
      pct_tim3_pos_within_temra_like = safe_pct(tim3, cd3cd8 & temra),
      pct_lag3_pos_within_temra_like = safe_pct(lag3, cd3cd8 & temra),
      pct_tigit_pos_within_temra_like = safe_pct(tigit, cd3cd8 & temra),
      pct_cd39_pos_within_temra_like = safe_pct(cd39, cd3cd8 & temra),
      pct_pd1_pos_within_temra_cd27neg = safe_pct(pd1, cd3cd8 & temra_cd27neg),
      pct_tigit_pos_within_temra_cd27neg = safe_pct(tigit, cd3cd8 & temra_cd27neg),
      pct_cd39_pos_within_temra_cd27neg = safe_pct(cd39, cd3cd8 & temra_cd27neg),
      pct_pd1_temra = safe_pct(pd1 & temra, cd3cd8),
      pct_tigit_temra = safe_pct(tigit & temra, cd3cd8),
      pct_cd39_temra = safe_pct(cd39 & temra, cd3cd8),
      pct_pd1_temra_cd27neg = safe_pct(pd1 & temra_cd27neg, cd3cd8),
      pct_tigit_temra_cd27neg = safe_pct(tigit & temra_cd27neg, cd3cd8),
      pct_cd39_temra_cd27neg = safe_pct(cd39 & temra_cd27neg, cd3cd8),
      pct_pd1_tigit_cd39_temra_cd27neg = safe_pct(pd1 & tigit & cd39 & temra_cd27neg, cd3cd8)
    )
    row
  }, error = function(e) {
    tibble(
      file_name = basename(path), file_path = normalizePath(path, mustWork = FALSE),
      subject_id = extract_subject_id(path), panel = "CP7", feature_ok = FALSE,
      error_message = conditionMessage(e), total_events = NA_integer_, n_channels = NA_integer_,
      n_cd3_pos = NA_integer_, n_cd8_pos = NA_integer_, n_cd3_cd8_pos = NA_integer_
    )
  })
}

message("Extracting reconstructed CP7 features from ", length(fcs_files), " files...")
feature_table <- purrr::map_dfr(fcs_files, extract_one)

qc_summary <- tibble(
  n_files = nrow(feature_table),
  n_unique_subjects = n_distinct(feature_table$subject_id),
  feature_ok = sum(feature_table$feature_ok, na.rm = TRUE),
  feature_failed = sum(!feature_table$feature_ok, na.rm = TRUE),
  median_total_events = median(feature_table$total_events, na.rm = TRUE),
  min_total_events = min(feature_table$total_events, na.rm = TRUE),
  max_total_events = max(feature_table$total_events, na.rm = TRUE),
  median_cd3cd8_events = median(feature_table$n_cd3_cd8_pos, na.rm = TRUE),
  min_cd3cd8_events = min(feature_table$n_cd3_cd8_pos, na.rm = TRUE),
  max_cd3cd8_events = max(feature_table$n_cd3_cd8_pos, na.rm = TRUE)
)

readr::write_csv(feature_table, file.path(out_dir, "SDY2583_CP7_FULL_850_feature_table_STEP2_RECONSTRUCTED.csv"))
readr::write_csv(qc_summary, file.path(out_dir, "SDY2583_CP7_FULL_850_STEP2_QC_summary_RECONSTRUCTED.csv"))
save(feature_table, qc_summary, thresholds, threshold_table,
     file = file.path(rdata_dir, "SDY2583_CP7_STEP2_feature_extraction_RECONSTRUCTED.RData"))

print(qc_summary)
message("Expected archived benchmark: 850/850 successful; median total events 301989.5; median CD3+CD8+ events 12411.5.")
