# ============================================================
# SDY2583 CP22
# STEP 2 SAFE: B-cell / humoral remodeling feature extraction
# Marker-name tolerant extraction for CP22 mismatch files
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "purrr")

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

out_dir <- file.path(analysis_dir, "02_feature_extraction")
rdata_dir <- file.path(analysis_dir, "11_RData")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(step1_rdata)) {
  stop("Step 1 RData bulunamadı: ", step1_rdata)
}

load(step1_rdata)

if (!exists("fcs_files")) {
  stop("Step 1 RData içinde fcs_files bulunamadı.")
}

# ------------------------------------------------------------
# 3. Threshold definitions
# ------------------------------------------------------------

thresholds <- list(
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

threshold_table <- tibble::tibble(
  marker_or_gate = names(thresholds),
  threshold = unlist(thresholds)
)

readr::write_csv(
  threshold_table,
  file.path(out_dir, "SDY2583_CP22_thresholds_STEP2.csv")
)

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

  # Dump channel: allow both CP22 marker variants
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

  out <- tryCatch({
    ff2 <- flowCore::compensate(ff, sp)
    list(ff = ff2, compensation_applied = TRUE, compensation_error = NA_character_)
  }, error = function(e) {
    list(ff = ff, compensation_applied = FALSE, compensation_error = as.character(e$message))
  })

  out
}

transform_safely <- function(ff) {
  ex_names <- colnames(flowCore::exprs(ff))

  fluor_ch <- ex_names[
    !grepl("FSC|SSC|Time", ex_names, ignore.case = TRUE)
  ]

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

  out <- tryCatch({
    tf <- flowCore::transformList(fluor_ch, trans)
    ff2 <- flowCore::transform(ff, tf)
    list(ff = ff2, transform_applied = TRUE, transform_error = NA_character_)
  }, error = function(e) {
    list(ff = ff, transform_applied = FALSE, transform_error = as.character(e$message))
  })

  out
}

get_vec <- function(ex, ch) {
  if (is.na(ch) || !(ch %in% colnames(ex))) {
    return(rep(NA_real_, nrow(ex)))
  }
  as.numeric(ex[, ch])
}

# ------------------------------------------------------------
# 5. Feature extraction for one FCS file
# ------------------------------------------------------------

extract_one_cp22 <- function(fp) {

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

    if (n_b == 0) {
      b_idx <- rep(FALSE, total_events)
    } else {
      b_idx <- b_cell
    }

    # B-cell marker positivity
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

    # Core B-cell states
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
    cd38_cd138_plasma_cell_like <- b_idx &
      cd38 > thresholds$CD38 &
      cd138 > thresholds$CD138

    transitional_like <- b_idx &
      cd24 > thresholds$CD24 &
      cd38 > thresholds$CD38

    immature_transitional_like <- b_idx &
      cd10 > thresholds$CD10

    cd10_cd24_cd38_transitional_like <- b_idx &
      cd10 > thresholds$CD10 &
      cd24 > thresholds$CD24 &
      cd38 > thresholds$CD38

    cd39_regulatory_like <- b_idx & cd39 > thresholds$CD39
    cd39_cd24_regulatory_like <- b_idx &
      cd39 > thresholds$CD39 &
      cd24 > thresholds$CD24

    cd39_cd24_cd38_regulatory_transitional_like <- b_idx &
      cd39 > thresholds$CD39 &
      cd24 > thresholds$CD24 &
      cd38 > thresholds$CD38

    # Counts
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

    # Output row
    tibble::tibble(
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

      n_naive_like = n_naive_like,
      n_unswitched_memory_like = n_unswitched_memory_like,
      n_switched_memory_like = n_switched_memory_like,
      n_double_negative_like = n_double_negative_like,

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
cat("SDY2583 CP22 STEP 2 FEATURE EXTRACTION STARTED\n")
cat("============================================================\n")

feature_list <- vector("list", length(fcs_files))

for (i in seq_along(fcs_files)) {
  if (i %% 25 == 0) {
    cat("  processed", i, "of", length(fcs_files), "files\n")
  }
  feature_list[[i]] <- extract_one_cp22(fcs_files[i])
}

cp22_feature_table <- dplyr::bind_rows(feature_list)

# ------------------------------------------------------------
# 7. Summaries
# ------------------------------------------------------------

feature_summary <- tibble::tibble(
  n_rows = nrow(cp22_feature_table),
  n_unique_subjects = dplyr::n_distinct(cp22_feature_table$subject_id),
  n_feature_ok = sum(cp22_feature_table$feature_ok == TRUE, na.rm = TRUE),
  n_feature_failed = sum(cp22_feature_table$feature_ok == FALSE, na.rm = TRUE),
  n_compensation_applied = sum(cp22_feature_table$compensation_applied == TRUE, na.rm = TRUE),
  n_transform_applied = sum(cp22_feature_table$transform_applied == TRUE, na.rm = TRUE),
  median_total_events = median(cp22_feature_table$total_events, na.rm = TRUE),
  median_dump_low_events = median(cp22_feature_table$n_dump_low, na.rm = TRUE),
  median_cd19_b_events = median(cp22_feature_table$n_cd19_b, na.rm = TRUE),
  min_cd19_b_events = min(cp22_feature_table$n_cd19_b, na.rm = TRUE),
  max_cd19_b_events = max(cp22_feature_table$n_cd19_b, na.rm = TRUE),
  median_pct_cd19_b_total = median(cp22_feature_table$pct_cd19_b_total, na.rm = TRUE),
  median_pct_cd19_b_within_dump_low = median(cp22_feature_table$pct_cd19_b_within_dump_low, na.rm = TRUE)
)

dump_variant_summary <- cp22_feature_table %>%
  dplyr::count(dump_marker_variant, sort = TRUE)

event_qc_by_subject <- cp22_feature_table %>%
  dplyr::mutate(
    b_event_qc_bin = dplyr::case_when(
      is.na(n_cd19_b) ~ "missing",
      n_cd19_b < 50 ~ "<50",
      n_cd19_b < 100 ~ "50-99",
      n_cd19_b < 300 ~ "100-299",
      n_cd19_b < 1000 ~ "300-999",
      TRUE ~ ">=1000"
    )
  ) %>%
  dplyr::count(b_event_qc_bin, sort = FALSE)

selected_feature_medians <- cp22_feature_table %>%
  dplyr::summarise(
    median_pct_cd19_b_total = median(pct_cd19_b_total, na.rm = TRUE),
    median_pct_cd19_b_within_dump_low = median(pct_cd19_b_within_dump_low, na.rm = TRUE),
    median_pct_naive_like_within_b = median(pct_naive_like_within_b, na.rm = TRUE),
    median_pct_unswitched_memory_like_within_b = median(pct_unswitched_memory_like_within_b, na.rm = TRUE),
    median_pct_switched_memory_like_within_b = median(pct_switched_memory_like_within_b, na.rm = TRUE),
    median_pct_double_negative_like_within_b = median(pct_double_negative_like_within_b, na.rm = TRUE),
    median_pct_iga_switched_memory_like_within_b = median(pct_iga_switched_memory_like_within_b, na.rm = TRUE),
    median_pct_igg_switched_memory_like_within_b = median(pct_igg_switched_memory_like_within_b, na.rm = TRUE),
    median_pct_plasmablast_like_within_b = median(pct_plasmablast_like_within_b, na.rm = TRUE),
    median_pct_plasma_cell_like_within_b = median(pct_plasma_cell_like_within_b, na.rm = TRUE),
    median_pct_transitional_like_within_b = median(pct_transitional_like_within_b, na.rm = TRUE),
    median_pct_cd39_regulatory_like_within_b = median(pct_cd39_regulatory_like_within_b, na.rm = TRUE)
  )

failed_files <- cp22_feature_table %>%
  dplyr::filter(feature_ok == FALSE)

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  cp22_feature_table,
  file.path(out_dir, "SDY2583_CP22_FULL_850_feature_table_STEP2_Bcell_humoral.csv")
)

readr::write_csv(
  feature_summary,
  file.path(out_dir, "SDY2583_CP22_feature_extraction_summary_STEP2_Bcell_humoral.csv")
)

readr::write_csv(
  dump_variant_summary,
  file.path(out_dir, "SDY2583_CP22_dump_marker_variant_summary_STEP2.csv")
)

readr::write_csv(
  event_qc_by_subject,
  file.path(out_dir, "SDY2583_CP22_Bcell_event_QC_summary_STEP2.csv")
)

readr::write_csv(
  selected_feature_medians,
  file.path(out_dir, "SDY2583_CP22_selected_feature_medians_STEP2.csv")
)

readr::write_csv(
  failed_files,
  file.path(out_dir, "SDY2583_CP22_failed_files_STEP2.csv")
)

save(
  cp22_feature_table,
  feature_summary,
  dump_variant_summary,
  event_qc_by_subject,
  selected_feature_medians,
  failed_files,
  thresholds,
  threshold_table,
  file = file.path(
    rdata_dir,
    "SDY2583_CP22_STEP2_feature_extraction_Bcell_humoral.RData"
  )
)

# Backward/simple name for later scripts
save(
  cp22_feature_table,
  feature_summary,
  dump_variant_summary,
  event_qc_by_subject,
  selected_feature_medians,
  failed_files,
  thresholds,
  threshold_table,
  file = file.path(
    rdata_dir,
    "SDY2583_CP22_STEP2_feature_extraction.RData"
  )
)

# ------------------------------------------------------------
# 9. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 2 COMPLETE: B-CELL FEATURE EXTRACTION\n")
cat("============================================================\n")

cat("\nFeature extraction summary:\n")
print(as.data.frame(feature_summary), row.names = FALSE)

cat("\nDump marker variant summary:\n")
print(as.data.frame(dump_variant_summary), row.names = FALSE)

cat("\nB-cell event QC summary:\n")
print(as.data.frame(event_qc_by_subject), row.names = FALSE)

cat("\nSelected feature medians:\n")
print(as.data.frame(selected_feature_medians), row.names = FALSE)

cat("\nFailed files:\n")
if (nrow(failed_files) == 0) {
  cat("No failed files.\n")
} else {
  print(as.data.frame(failed_files), row.names = FALSE)
}

cat("\nFiles saved in:\n")
print(out_dir)

cat("============================================================\n")
