# ============================================================
# SDY2583 CP10
# STEP 2 SAFE: Feature extraction
# Myeloid / granulocyte-like / leukocyte-lineage composition panel
#
# Marker-name tolerant extraction is REQUIRED because DBG161 and DBG177
# have Time as parameter 1 and shifted fluorescence/scatter positions.
#
# Primary gate:
#   Viability-low CD45+ leukocytes
#
# Main axes:
#   - leukocyte-lineage composition: CD3, CD19, CD56
#   - monocyte-like / HLA-DR axis: CD14, HLA-DR, CD11c
#   - granulocyte-like axis: CD13, CD66b
#   - CCR3/CD193 eosinophil-like axis
#   - CD123 HLA-DR DC-like axis
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("flowCore", quietly = TRUE)) BiocManager::install("flowCore", ask = FALSE, update = FALSE)

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

analysis_dir <- sd_analysis_dir("CP10")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "02_feature_extraction_myeloid_granulocyte")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP1_fcs_inventory_marker_QC.RData")
step1b_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP1B_mismatch_diagnostic.RData")
if (!file.exists(step1_rdata)) stop("Step 1 RData bulunamadı: ", step1_rdata)
load(step1_rdata)
if (!exists("fcs_files")) stop("Step 1 RData içinde fcs_files bulunamadı.")
step1b_available <- file.exists(step1b_rdata)
if (step1b_available) load(step1b_rdata)

# Main thresholds on transformed intensities. Step 5 will stress-test them.
thresholds <- list(
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
threshold_table <- tibble(marker_or_gate = names(thresholds), threshold = unlist(thresholds))
readr::write_csv(threshold_table, file.path(out_dir, "SDY2583_CP10_STEP2_thresholds_main.csv"))

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
  ch$VIABILITY <- find_marker_channel_regex(map, "^Viability$|Viability|Live|Dead")
  list(channels = ch)
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
  if (is.null(sp)) return(list(ff = ff, compensation_applied = FALSE, compensation_error = NA_character_))
  tryCatch({
    list(ff = flowCore::compensate(ff, sp), compensation_applied = TRUE, compensation_error = NA_character_)
  }, error = function(e) {
    list(ff = ff, compensation_applied = FALSE, compensation_error = as.character(e$message))
  })
}
transform_safely <- function(ff) {
  ex_names <- colnames(flowCore::exprs(ff))
  fluor_ch <- ex_names[!grepl("FSC|SSC|Time", ex_names, ignore.case = TRUE)]
  if (length(fluor_ch) == 0) return(list(ff = ff, transform_applied = FALSE, transform_error = "No fluorescence channels detected"))
  trans <- flowCore::logicleTransform(transformationId = "fixed_logicle", w = 0.5, t = 262144, m = 4.5, a = 0)
  tryCatch({
    tf <- flowCore::transformList(fluor_ch, trans)
    list(ff = flowCore::transform(ff, tf), transform_applied = TRUE, transform_error = NA_character_)
  }, error = function(e) {
    list(ff = ff, transform_applied = FALSE, transform_error = as.character(e$message))
  })
}
get_vec <- function(ex, ch) {
  if (is.na(ch) || !(ch %in% colnames(ex))) return(rep(NA_real_, nrow(ex)))
  as.numeric(ex[, ch])
}
safe_pct <- function(num, den) ifelse(is.na(den) | den <= 0, NA_real_, 100 * num / den)
safe_ratio <- function(num, den) ifelse(is.na(den) | den <= 0, NA_real_, num / den)
safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}
safe_sum <- function(x) sum(x, na.rm = TRUE)

extract_one_cp10 <- function(fp) {
  tryCatch({
    ff_raw <- flowCore::read.FCS(fp, transformation = FALSE, truncate_max_range = FALSE)
    marker_map_file <- get_marker_map(ff_raw)
    ch <- build_channel_map(marker_map_file)$channels
    required <- c("VIABILITY", "CD45", "CD3", "CD19", "CD56", "CD14", "HLA_DR", "CD11c", "CD13", "CD66b", "CCR3", "CD123")
    missing_markers <- required[vapply(ch[required], function(x) is.na(x), logical(1))]
    if (length(missing_markers) > 0) stop("Missing required markers/channels: ", paste(missing_markers, collapse = ", "))

    comp <- apply_compensation_safely(ff_raw)
    trans <- transform_safely(comp$ff)
    ff <- trans$ff
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

    viable <- viability <= thresholds$VIABILITY_LOW
    primary <- viable & cd45 > thresholds$CD45
    n_viable <- safe_sum(viable)
    n_cd45 <- safe_sum(primary)

    t_like <- primary & cd3 > thresholds$CD3
    b_like <- primary & cd19 > thresholds$CD19
    nk_like <- primary & cd3 <= thresholds$CD3 & cd56 > thresholds$CD56
    lymphoid_like_any <- primary & (cd3 > thresholds$CD3 | cd19 > thresholds$CD19 | cd56 > thresholds$CD56)

    cd14_mono_like <- primary & cd14 > thresholds$CD14
    cd14_hladr_pos <- cd14_mono_like & hla_dr > thresholds$HLA_DR
    cd14_hladr_low <- cd14_mono_like & hla_dr <= thresholds$HLA_DR
    cd14_cd11c_pos <- cd14_mono_like & cd11c > thresholds$CD11c
    cd14_cd11c_hladr_pos <- cd14_mono_like & cd11c > thresholds$CD11c & hla_dr > thresholds$HLA_DR
    cd11c_pos <- primary & cd11c > thresholds$CD11c
    cd11c_hladr_apc_like <- primary & cd11c > thresholds$CD11c & hla_dr > thresholds$HLA_DR
    hladr_pos <- primary & hla_dr > thresholds$HLA_DR

    cd66b_gran_like <- primary & cd66b > thresholds$CD66b
    cd13_pos <- primary & cd13 > thresholds$CD13
    cd13_cd66b_gran_like <- primary & cd13 > thresholds$CD13 & cd66b > thresholds$CD66b
    ccr3_pos <- primary & ccr3 > thresholds$CCR3
    ccr3_cd66b_eosinophil_like <- primary & ccr3 > thresholds$CCR3 & cd66b > thresholds$CD66b
    ccr3_cd13_cd66b_gran_like <- primary & ccr3 > thresholds$CCR3 & cd13 > thresholds$CD13 & cd66b > thresholds$CD66b

    cd123_pos <- primary & cd123 > thresholds$CD123
    cd123_hladr_dc_like <- primary & cd123 > thresholds$CD123 & hla_dr > thresholds$HLA_DR
    cd123_hladr_cd11c_neg_pdc_like <- primary & cd123 > thresholds$CD123 & hla_dr > thresholds$HLA_DR & cd11c <= thresholds$CD11c & cd14 <= thresholds$CD14 & cd3 <= thresholds$CD3 & cd19 <= thresholds$CD19
    cd123_cd11c_hladr_mixed_apc_like <- primary & cd123 > thresholds$CD123 & cd11c > thresholds$CD11c & hla_dr > thresholds$HLA_DR

    myeloid_granulocytic_like_any <- primary & (cd14 > thresholds$CD14 | cd11c > thresholds$CD11c | cd13 > thresholds$CD13 | cd66b > thresholds$CD66b | ccr3 > thresholds$CCR3 | cd123 > thresholds$CD123)

    n_t_like <- safe_sum(t_like); n_b_like <- safe_sum(b_like); n_nk_like <- safe_sum(nk_like); n_lymphoid_like_any <- safe_sum(lymphoid_like_any)
    n_cd14_mono_like <- safe_sum(cd14_mono_like); n_cd14_hladr_pos <- safe_sum(cd14_hladr_pos); n_cd14_hladr_low <- safe_sum(cd14_hladr_low); n_cd14_cd11c_pos <- safe_sum(cd14_cd11c_pos); n_cd14_cd11c_hladr_pos <- safe_sum(cd14_cd11c_hladr_pos); n_cd11c_pos <- safe_sum(cd11c_pos); n_cd11c_hladr_apc_like <- safe_sum(cd11c_hladr_apc_like); n_hladr_pos <- safe_sum(hladr_pos)
    n_cd66b_gran_like <- safe_sum(cd66b_gran_like); n_cd13_pos <- safe_sum(cd13_pos); n_cd13_cd66b_gran_like <- safe_sum(cd13_cd66b_gran_like); n_ccr3_pos <- safe_sum(ccr3_pos); n_ccr3_cd66b_eosinophil_like <- safe_sum(ccr3_cd66b_eosinophil_like); n_ccr3_cd13_cd66b_gran_like <- safe_sum(ccr3_cd13_cd66b_gran_like)
    n_cd123_pos <- safe_sum(cd123_pos); n_cd123_hladr_dc_like <- safe_sum(cd123_hladr_dc_like); n_cd123_hladr_cd11c_neg_pdc_like <- safe_sum(cd123_hladr_cd11c_neg_pdc_like); n_cd123_cd11c_hladr_mixed_apc_like <- safe_sum(cd123_cd11c_hladr_mixed_apc_like)
    n_myeloid_granulocytic_like_any <- safe_sum(myeloid_granulocytic_like_any)

    channel_order_mismatch_file <- FALSE; marker_mismatch_file <- FALSE
    if (exists("file_mismatch_summary")) {
      mm_row <- file_mismatch_summary %>% filter(file_path == fp)
      if (nrow(mm_row) > 0) {
        channel_order_mismatch_file <- any(mm_row$any_channel_mismatch, na.rm = TRUE)
        marker_mismatch_file <- any(mm_row$any_marker_mismatch, na.rm = TRUE)
      }
    }

    tibble(
      subject_id = extract_subject_id(fp), file_name = basename(fp), file_path = fp,
      feature_ok = TRUE, feature_error = NA_character_,
      transform_method = ifelse(trans$transform_applied, "fixed_logicle", "untransformed_fallback"),
      compensation_applied = comp$compensation_applied, compensation_error = comp$compensation_error,
      transform_applied = trans$transform_applied, transform_error = trans$transform_error,
      channel_order_mismatch_file = channel_order_mismatch_file, marker_mismatch_file = marker_mismatch_file,
      channel_Viability = ch$VIABILITY, channel_CD45 = ch$CD45, channel_CD3 = ch$CD3, channel_CD19 = ch$CD19, channel_CD56 = ch$CD56, channel_CD14 = ch$CD14, channel_HLA_DR = ch$HLA_DR, channel_CD11c = ch$CD11c, channel_CD13 = ch$CD13, channel_CD66b = ch$CD66b, channel_CCR3 = ch$CCR3, channel_CD123 = ch$CD123,
      total_events = total_events, n_viable = n_viable, n_cd45_viable = n_cd45,
      pct_viable_total = safe_pct(n_viable, total_events), pct_cd45_viable_total = safe_pct(n_cd45, total_events), pct_cd45_within_viable = safe_pct(n_cd45, n_viable),

      n_t_like = n_t_like, n_b_like = n_b_like, n_nk_like = n_nk_like, n_lymphoid_like_any = n_lymphoid_like_any,
      pct_t_like_within_cd45 = safe_pct(n_t_like, n_cd45), pct_b_like_within_cd45 = safe_pct(n_b_like, n_cd45), pct_nk_like_within_cd45 = safe_pct(n_nk_like, n_cd45), pct_lymphoid_like_any_within_cd45 = safe_pct(n_lymphoid_like_any, n_cd45),

      n_cd14_mono_like = n_cd14_mono_like, n_cd14_hladr_pos = n_cd14_hladr_pos, n_cd14_hladr_low = n_cd14_hladr_low, n_cd14_cd11c_pos = n_cd14_cd11c_pos, n_cd14_cd11c_hladr_pos = n_cd14_cd11c_hladr_pos, n_cd11c_pos = n_cd11c_pos, n_cd11c_hladr_apc_like = n_cd11c_hladr_apc_like, n_hladr_pos = n_hladr_pos,
      pct_cd14_mono_like_within_cd45 = safe_pct(n_cd14_mono_like, n_cd45), pct_cd14_hladr_pos_within_cd45 = safe_pct(n_cd14_hladr_pos, n_cd45), pct_cd14_hladr_low_within_cd45 = safe_pct(n_cd14_hladr_low, n_cd45), pct_cd14_cd11c_pos_within_cd45 = safe_pct(n_cd14_cd11c_pos, n_cd45), pct_cd14_cd11c_hladr_pos_within_cd45 = safe_pct(n_cd14_cd11c_hladr_pos, n_cd45), pct_cd11c_pos_within_cd45 = safe_pct(n_cd11c_pos, n_cd45), pct_cd11c_hladr_apc_like_within_cd45 = safe_pct(n_cd11c_hladr_apc_like, n_cd45), pct_hladr_pos_within_cd45 = safe_pct(n_hladr_pos, n_cd45),
      pct_hladr_pos_within_cd14_mono_like = safe_pct(n_cd14_hladr_pos, n_cd14_mono_like), pct_hladr_low_within_cd14_mono_like = safe_pct(n_cd14_hladr_low, n_cd14_mono_like), pct_cd11c_pos_within_cd14_mono_like = safe_pct(n_cd14_cd11c_pos, n_cd14_mono_like), pct_cd11c_hladr_pos_within_cd14_mono_like = safe_pct(n_cd14_cd11c_hladr_pos, n_cd14_mono_like),

      n_cd66b_gran_like = n_cd66b_gran_like, n_cd13_pos = n_cd13_pos, n_cd13_cd66b_gran_like = n_cd13_cd66b_gran_like, n_ccr3_pos = n_ccr3_pos, n_ccr3_cd66b_eosinophil_like = n_ccr3_cd66b_eosinophil_like, n_ccr3_cd13_cd66b_gran_like = n_ccr3_cd13_cd66b_gran_like,
      pct_cd66b_gran_like_within_cd45 = safe_pct(n_cd66b_gran_like, n_cd45), pct_cd13_pos_within_cd45 = safe_pct(n_cd13_pos, n_cd45), pct_cd13_cd66b_gran_like_within_cd45 = safe_pct(n_cd13_cd66b_gran_like, n_cd45), pct_ccr3_pos_within_cd45 = safe_pct(n_ccr3_pos, n_cd45), pct_ccr3_cd66b_eosinophil_like_within_cd45 = safe_pct(n_ccr3_cd66b_eosinophil_like, n_cd45), pct_ccr3_cd13_cd66b_gran_like_within_cd45 = safe_pct(n_ccr3_cd13_cd66b_gran_like, n_cd45),
      pct_ccr3_pos_within_cd66b_gran_like = safe_pct(n_ccr3_cd66b_eosinophil_like, n_cd66b_gran_like), pct_cd13_pos_within_cd66b_gran_like = safe_pct(n_cd13_cd66b_gran_like, n_cd66b_gran_like),

      n_cd123_pos = n_cd123_pos, n_cd123_hladr_dc_like = n_cd123_hladr_dc_like, n_cd123_hladr_cd11c_neg_pdc_like = n_cd123_hladr_cd11c_neg_pdc_like, n_cd123_cd11c_hladr_mixed_apc_like = n_cd123_cd11c_hladr_mixed_apc_like,
      pct_cd123_pos_within_cd45 = safe_pct(n_cd123_pos, n_cd45), pct_cd123_hladr_dc_like_within_cd45 = safe_pct(n_cd123_hladr_dc_like, n_cd45), pct_cd123_hladr_cd11c_neg_pdc_like_within_cd45 = safe_pct(n_cd123_hladr_cd11c_neg_pdc_like, n_cd45), pct_cd123_cd11c_hladr_mixed_apc_like_within_cd45 = safe_pct(n_cd123_cd11c_hladr_mixed_apc_like, n_cd45), pct_hladr_pos_within_cd123_pos = safe_pct(n_cd123_hladr_dc_like, n_cd123_pos),

      n_myeloid_granulocytic_like_any = n_myeloid_granulocytic_like_any,
      pct_myeloid_granulocytic_like_any_within_cd45 = safe_pct(n_myeloid_granulocytic_like_any, n_cd45),
      ratio_myeloid_granulocytic_to_lymphoid_like = safe_ratio(n_myeloid_granulocytic_like_any, n_lymphoid_like_any),
      ratio_cd66b_gran_like_to_t_like = safe_ratio(n_cd66b_gran_like, n_t_like),
      ratio_cd14_mono_like_to_t_like = safe_ratio(n_cd14_mono_like, n_t_like),

      median_CD45_in_CD45 = safe_median(cd45[primary]), median_CD3_in_CD45 = safe_median(cd3[primary]), median_CD19_in_CD45 = safe_median(cd19[primary]), median_CD56_in_CD45 = safe_median(cd56[primary]), median_CD14_in_CD45 = safe_median(cd14[primary]), median_HLA_DR_in_CD45 = safe_median(hla_dr[primary]), median_CD11c_in_CD45 = safe_median(cd11c[primary]), median_CD13_in_CD45 = safe_median(cd13[primary]), median_CD66b_in_CD45 = safe_median(cd66b[primary]), median_CCR3_in_CD45 = safe_median(ccr3[primary]), median_CD123_in_CD45 = safe_median(cd123[primary]),
      median_HLA_DR_in_CD14_mono_like = safe_median(hla_dr[cd14_mono_like]), median_CD11c_in_CD14_mono_like = safe_median(cd11c[cd14_mono_like]), median_CD14_in_CD14_mono_like = safe_median(cd14[cd14_mono_like]),
      median_CD66b_in_CD66b_gran_like = safe_median(cd66b[cd66b_gran_like]), median_CD13_in_CD66b_gran_like = safe_median(cd13[cd66b_gran_like]), median_CCR3_in_CD66b_gran_like = safe_median(ccr3[cd66b_gran_like]),
      median_CCR3_in_CCR3_CD66b_eosinophil_like = safe_median(ccr3[ccr3_cd66b_eosinophil_like]), median_CD123_in_CD123_HLADR_dc_like = safe_median(cd123[cd123_hladr_dc_like]), median_HLA_DR_in_CD123_HLADR_dc_like = safe_median(hla_dr[cd123_hladr_dc_like])
    )
  }, error = function(e) {
    tibble(subject_id = extract_subject_id(fp), file_name = basename(fp), file_path = fp, feature_ok = FALSE, feature_error = as.character(e$message))
  })
}

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 2 STARTED: FEATURE EXTRACTION\n")
cat("============================================================\n")
cat("\nNumber of FCS files:", length(fcs_files), "\n")

feature_rows <- vector("list", length(fcs_files))
for (i in seq_along(fcs_files)) {
  if (i %% 50 == 0) cat("Processed", i, "of", length(fcs_files), "files\n")
  feature_rows[[i]] <- extract_one_cp10(fcs_files[i])
}
cp10_features <- bind_rows(feature_rows)

feature_summary <- tibble(
  n_fcs_files = length(fcs_files),
  n_unique_subjects = n_distinct(cp10_features$subject_id),
  n_feature_ok = sum(cp10_features$feature_ok == TRUE, na.rm = TRUE),
  n_feature_failed = sum(cp10_features$feature_ok == FALSE, na.rm = TRUE),
  n_channel_order_mismatch_files = sum(cp10_features$channel_order_mismatch_file == TRUE, na.rm = TRUE),
  n_marker_mismatch_files = sum(cp10_features$marker_mismatch_file == TRUE, na.rm = TRUE),
  median_total_events = median(cp10_features$total_events, na.rm = TRUE),
  median_viable_events = median(cp10_features$n_viable, na.rm = TRUE),
  median_cd45_viable_events = median(cp10_features$n_cd45_viable, na.rm = TRUE),
  min_cd45_viable_events = min(cp10_features$n_cd45_viable, na.rm = TRUE),
  max_cd45_viable_events = max(cp10_features$n_cd45_viable, na.rm = TRUE),
  median_pct_cd45_viable_total = median(cp10_features$pct_cd45_viable_total, na.rm = TRUE),
  median_pct_cd45_within_viable = median(cp10_features$pct_cd45_within_viable, na.rm = TRUE)
)

event_qc_summary <- cp10_features %>%
  mutate(cd45_event_qc_bin = case_when(
    is.na(n_cd45_viable) ~ "missing",
    n_cd45_viable < 50 ~ "<50",
    n_cd45_viable < 100 ~ "50-99",
    n_cd45_viable < 300 ~ "100-299",
    n_cd45_viable < 1000 ~ "300-999",
    TRUE ~ ">=1000"
  )) %>%
  count(cd45_event_qc_bin, name = "n_files") %>%
  arrange(cd45_event_qc_bin)

feature_dictionary <- tibble::tribble(
  ~feature, ~module, ~interpretation,
  "pct_t_like_within_cd45", "lineage_composition", "CD3+ T-cell-like fraction within viable CD45+ leukocytes",
  "pct_b_like_within_cd45", "lineage_composition", "CD19+ B-cell-like fraction within viable CD45+ leukocytes",
  "pct_nk_like_within_cd45", "lineage_composition", "CD3- CD56+ NK-like fraction within viable CD45+ leukocytes",
  "pct_lymphoid_like_any_within_cd45", "lineage_composition", "CD3/CD19/CD56-defined lymphoid-like fraction within viable CD45+ leukocytes",
  "pct_cd14_mono_like_within_cd45", "monocyte_HLA_DR_axis", "CD14+ monocyte-like fraction within viable CD45+ leukocytes",
  "pct_cd14_hladr_pos_within_cd45", "monocyte_HLA_DR_axis", "CD14+HLA-DR+ monocyte-like fraction within viable CD45+ leukocytes",
  "pct_cd14_hladr_low_within_cd45", "monocyte_HLA_DR_axis", "CD14+HLA-DR-low monocyte-like fraction within viable CD45+ leukocytes",
  "pct_hladr_pos_within_cd14_mono_like", "monocyte_HLA_DR_axis", "HLA-DR+ fraction within CD14+ monocyte-like cells",
  "pct_hladr_low_within_cd14_mono_like", "monocyte_HLA_DR_axis", "HLA-DR-low fraction within CD14+ monocyte-like cells",
  "median_HLA_DR_in_CD14_mono_like", "monocyte_HLA_DR_axis", "HLA-DR median intensity in CD14+ monocyte-like cells",
  "pct_cd11c_hladr_apc_like_within_cd45", "APC_like_myeloid_axis", "CD11c+HLA-DR+ APC-like myeloid fraction within viable CD45+ leukocytes",
  "pct_cd14_cd11c_hladr_pos_within_cd45", "APC_like_myeloid_axis", "CD14+CD11c+HLA-DR+ myeloid APC-like fraction within viable CD45+ leukocytes",
  "pct_cd66b_gran_like_within_cd45", "granulocyte_like_axis", "CD66b+ granulocyte-like fraction within viable CD45+ leukocytes",
  "pct_cd13_cd66b_gran_like_within_cd45", "granulocyte_like_axis", "CD13+CD66b+ granulocyte-like fraction within viable CD45+ leukocytes",
  "pct_cd13_pos_within_cd66b_gran_like", "granulocyte_like_axis", "CD13+ fraction within CD66b+ granulocyte-like cells",
  "median_CD66b_in_CD66b_gran_like", "granulocyte_like_axis", "CD66b median intensity in CD66b+ granulocyte-like cells",
  "median_CD13_in_CD66b_gran_like", "granulocyte_like_axis", "CD13 median intensity in CD66b+ granulocyte-like cells",
  "pct_ccr3_pos_within_cd45", "CCR3_eosinophil_like_axis", "CCR3/CD193+ fraction within viable CD45+ leukocytes",
  "pct_ccr3_cd66b_eosinophil_like_within_cd45", "CCR3_eosinophil_like_axis", "CCR3+CD66b+ eosinophil-like granulocytic fraction within viable CD45+ leukocytes",
  "pct_ccr3_cd13_cd66b_gran_like_within_cd45", "CCR3_eosinophil_like_axis", "CCR3+CD13+CD66b+ granulocyte-like fraction within viable CD45+ leukocytes",
  "pct_ccr3_pos_within_cd66b_gran_like", "CCR3_eosinophil_like_axis", "CCR3+ fraction within CD66b+ granulocyte-like cells",
  "median_CCR3_in_CD66b_gran_like", "CCR3_eosinophil_like_axis", "CCR3 median intensity in CD66b+ granulocyte-like cells",
  "pct_cd123_pos_within_cd45", "CD123_DC_like_axis", "CD123+ fraction within viable CD45+ leukocytes",
  "pct_cd123_hladr_dc_like_within_cd45", "CD123_DC_like_axis", "CD123+HLA-DR+ DC-like fraction within viable CD45+ leukocytes",
  "pct_cd123_hladr_cd11c_neg_pdc_like_within_cd45", "CD123_DC_like_axis", "CD123+HLA-DR+CD11c-CD14- lineage-negative pDC-like fraction within viable CD45+ leukocytes",
  "pct_cd123_cd11c_hladr_mixed_apc_like_within_cd45", "CD123_DC_like_axis", "CD123+CD11c+HLA-DR+ mixed APC-like fraction within viable CD45+ leukocytes",
  "pct_myeloid_granulocytic_like_any_within_cd45", "integrated_myeloid_granulocytic_balance", "Union myeloid/granulocytic-like fraction within viable CD45+ leukocytes",
  "ratio_myeloid_granulocytic_to_lymphoid_like", "integrated_myeloid_granulocytic_balance", "Ratio of myeloid/granulocytic-like to lymphoid-like events",
  "ratio_cd66b_gran_like_to_t_like", "integrated_myeloid_granulocytic_balance", "Ratio of CD66b+ granulocyte-like to CD3+ T-like events",
  "ratio_cd14_mono_like_to_t_like", "integrated_myeloid_granulocytic_balance", "Ratio of CD14+ monocyte-like to CD3+ T-like events"
)

readr::write_csv(cp10_features, file.path(out_dir, "SDY2583_CP10_features_STEP2.csv"))
readr::write_csv(feature_summary, file.path(out_dir, "SDY2583_CP10_feature_extraction_summary_STEP2.csv"))
readr::write_csv(event_qc_summary, file.path(out_dir, "SDY2583_CP10_event_QC_summary_STEP2.csv"))
readr::write_csv(feature_dictionary, file.path(out_dir, "SDY2583_CP10_feature_dictionary_STEP2.csv"))

save(cp10_features, feature_summary, event_qc_summary, feature_dictionary, thresholds, threshold_table, step1b_available,
     file = file.path(rdata_dir, "SDY2583_CP10_STEP2_feature_extraction_myeloid_granulocyte.RData"))

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 2 COMPLETE: FEATURE EXTRACTION\n")
cat("============================================================\n")

cat("\nFeature extraction summary:\n")
print(as.data.frame(feature_summary), row.names = FALSE)

cat("\nCD45 event QC summary:\n")
print(as.data.frame(event_qc_summary), row.names = FALSE)

cat("\nSelected CP10 feature medians:\n")
selected_medians <- cp10_features %>%
  summarise(
    median_pct_cd45_viable_total = median(pct_cd45_viable_total, na.rm = TRUE),
    median_pct_t_like_within_cd45 = median(pct_t_like_within_cd45, na.rm = TRUE),
    median_pct_b_like_within_cd45 = median(pct_b_like_within_cd45, na.rm = TRUE),
    median_pct_nk_like_within_cd45 = median(pct_nk_like_within_cd45, na.rm = TRUE),
    median_pct_cd14_mono_like_within_cd45 = median(pct_cd14_mono_like_within_cd45, na.rm = TRUE),
    median_pct_hladr_low_within_cd14_mono_like = median(pct_hladr_low_within_cd14_mono_like, na.rm = TRUE),
    median_pct_cd66b_gran_like_within_cd45 = median(pct_cd66b_gran_like_within_cd45, na.rm = TRUE),
    median_pct_ccr3_cd66b_eosinophil_like_within_cd45 = median(pct_ccr3_cd66b_eosinophil_like_within_cd45, na.rm = TRUE),
    median_pct_cd123_hladr_dc_like_within_cd45 = median(pct_cd123_hladr_dc_like_within_cd45, na.rm = TRUE),
    median_ratio_myeloid_granulocytic_to_lymphoid_like = median(ratio_myeloid_granulocytic_to_lymphoid_like, na.rm = TRUE)
  )
print(as.data.frame(selected_medians), row.names = FALSE)

cat("\nFeature dictionary:\n")
print(as.data.frame(feature_dictionary), row.names = FALSE)

cat("\nFailed files, if any:\n")
print(as.data.frame(cp10_features %>% filter(feature_ok == FALSE) %>% select(subject_id, file_name, feature_error)), row.names = FALSE)

cat("\nOutputs saved in:\n")
print(out_dir)
cat("============================================================\n")
