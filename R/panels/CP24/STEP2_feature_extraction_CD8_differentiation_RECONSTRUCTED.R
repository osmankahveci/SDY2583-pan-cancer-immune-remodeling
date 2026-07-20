# CP24 full-cohort Step 2 reconstructed CD8 differentiation feature extraction.
# The output includes the complete official 55-row model family, including the
# full 19-endpoint PD-1 module.

rm(list = ls())
source(file.path(
  Sys.getenv("SDY2583_REPO_ROOT", unset = "."),
  "R", "shared", "reconstructed_panel_framework.R"
))
rp_install_and_load(
  c("dplyr", "readr", "stringr", "tibble", "purrr"),
  "flowCore"
)
source(file.path(
  sd_repo_root(), "R", "panels", "CP24", "MANIFEST_RECONSTRUCTED.R"
))

analysis_dir <- sd_analysis_dir("CP24")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "02_feature_extraction")
threshold_dir <- file.path(analysis_dir, "03_thresholds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(threshold_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(
  rdata_dir,
  "SDY2583_CP24_STEP1_inventory_RECONSTRUCTED.RData"
))

thresholds_main <- CP24_MANIFEST$thresholds
readr::write_csv(
  tibble::tibble(
    marker = names(thresholds_main),
    threshold = as.numeric(thresholds_main)
  ),
  file.path(threshold_dir, "SDY2583_CP24_thresholds_RECONSTRUCTED.csv")
)

extract_one <- function(path, thresholds = thresholds_main) {
  tryCatch({
    object <- rp_read_transform(path)
    marker_map <- object$marker_map
    expr <- flowCore::exprs(object$ff)

    channels <- c(
      CD3 = rp_find_channel(marker_map, c("^CD3$"), "BV605-A"),
      CD8 = rp_find_channel(marker_map, c("^CD8$"), "PerCP-Cy5-5-A"),
      CD62L = rp_find_channel(marker_map, c("CD62L"), "BV650-A"),
      CD95 = rp_find_channel(marker_map, c("CD95"), "BV711-A"),
      CD57 = rp_find_channel(marker_map, c("^CD57$"), "BV786-A"),
      CD27 = rp_find_channel(marker_map, c("^CD27$"), "BB515-A"),
      CX3CR1 = rp_find_channel(marker_map, c("CX3CR1"), "PE-A"),
      PD1 = rp_find_channel(marker_map, c("PD-?1"), "PE-CF594-A"),
      CD45RA = rp_find_channel(marker_map, c("CD45RA"), "PE-Cy5-A"),
      CXCR3 = rp_find_channel(marker_map, c("CXCR3"), "BV421-A"),
      CXCR5 = rp_find_channel(marker_map, c("CXCR5"), "PE-Cy7-A")
    )

    if (any(is.na(channels))) {
      stop(
        "Missing CP24 markers: ",
        paste(names(channels)[is.na(channels)], collapse = ", ")
      )
    }

    values <- lapply(channels, function(channel) rp_vec(expr, channel))

    cd3 <- values$CD3 > thresholds[["CD3"]]
    cd8 <- values$CD8 > thresholds[["CD8"]]
    cd3_cd8 <- cd3 & cd8

    cd45ra <- values$CD45RA > thresholds[["CD45RA"]]
    cd62l <- values$CD62L > thresholds[["CD62L"]]
    cd27 <- values$CD27 > thresholds[["CD27"]]
    cd57 <- values$CD57 > thresholds[["CD57"]]
    cx3cr1 <- values$CX3CR1 > thresholds[["CX3CR1"]]
    cd95 <- values$CD95 > thresholds[["CD95"]]
    pd1 <- values$PD1 > thresholds[["PD1"]]
    pd1_high <- values$PD1 > thresholds[["PD1_HIGH"]]

    naive_like <- cd45ra & cd62l
    tcm_like <- !cd45ra & cd62l
    tem_like <- !cd45ra & !cd62l
    temra_like <- cd45ra & !cd62l
    temra_cd27neg <- temra_like & !cd27

    tibble::tibble(
      subject_id = rp_extract_subject_id(path),
      result_file_name = basename(path),
      feature_ok = TRUE,
      error_message = NA_character_,
      total_events = nrow(expr),
      n_cd3_pos = rp_n(cd3),
      n_cd8_pos = rp_n(cd8),
      n_cd3_cd8_pos = rp_n(cd3_cd8),

      pct_cd3_pos_total = rp_pct(cd3),
      pct_cd8_pos_total = rp_pct(cd8),
      pct_cd3_cd8_pos_total = rp_pct(cd3_cd8),
      pct_cd8_within_cd3 = rp_pct(cd8, cd3),

      median_CD62L_in_CD3CD8 = rp_median(values$CD62L[cd3_cd8]),
      median_CD27_in_CD3CD8 = rp_median(values$CD27[cd3_cd8]),
      median_PD1_in_CD3CD8 = rp_median(values$PD1[cd3_cd8]),
      median_CD57_in_CD3CD8 = rp_median(values$CD57[cd3_cd8]),
      median_CX3CR1_in_CD3CD8 = rp_median(values$CX3CR1[cd3_cd8]),
      median_CD95_in_CD3CD8 = rp_median(values$CD95[cd3_cd8]),
      median_CD45RA_in_CD3CD8 = rp_median(values$CD45RA[cd3_cd8]),
      median_CXCR3_in_CD3CD8 = rp_median(values$CXCR3[cd3_cd8]),
      median_CXCR5_in_CD3CD8 = rp_median(values$CXCR5[cd3_cd8]),

      pct_naive_like = rp_pct(naive_like, cd3_cd8),
      pct_tcm_like = rp_pct(tcm_like, cd3_cd8),
      pct_tem_like = rp_pct(tem_like, cd3_cd8),
      pct_temra_like = rp_pct(temra_like, cd3_cd8),
      pct_cd62l_pos_within_cd3cd8 = rp_pct(cd62l, cd3_cd8),
      pct_cd62l_neg_within_cd3cd8 = rp_pct(!cd62l, cd3_cd8),
      pct_cd45ra_pos_within_cd3cd8 = rp_pct(cd45ra, cd3_cd8),
      pct_cd45ra_neg_within_cd3cd8 = rp_pct(!cd45ra, cd3_cd8),

      pct_cd27_pos_within_cd3cd8 = rp_pct(cd27, cd3_cd8),
      pct_cd27_neg_within_cd3cd8 = rp_pct(!cd27, cd3_cd8),
      pct_temra_cd27neg = rp_pct(temra_cd27neg, cd3_cd8),
      pct_temra_cd27pos = rp_pct(temra_like & cd27, cd3_cd8),
      pct_cd62lneg_cd27neg = rp_pct(!cd62l & !cd27, cd3_cd8),
      pct_cd45ra_pos_cd62lneg_cd27neg = rp_pct(
        cd45ra & !cd62l & !cd27,
        cd3_cd8
      ),
      pct_naive_cd27pos = rp_pct(naive_like & cd27, cd3_cd8),
      pct_naive_cd27neg = rp_pct(naive_like & !cd27, cd3_cd8),

      pct_cd57_pos = rp_pct(cd57, cd3_cd8),
      pct_cx3cr1_pos = rp_pct(cx3cr1, cd3_cd8),
      pct_cd95_pos = rp_pct(cd95, cd3_cd8),
      pct_cd57_cx3cr1_pos = rp_pct(cd57 & cx3cr1, cd3_cd8),
      pct_cd57_cd95_pos = rp_pct(cd57 & cd95, cd3_cd8),
      pct_cx3cr1_cd95_pos = rp_pct(cx3cr1 & cd95, cd3_cd8),
      pct_cd57_cx3cr1_cd95_pos = rp_pct(cd57 & cx3cr1 & cd95, cd3_cd8),

      pct_pd1_pos = rp_pct(pd1, cd3_cd8),
      pct_pd1_high = rp_pct(pd1_high, cd3_cd8),
      pct_pd1_pos_within_naive_like = rp_pct(pd1, cd3_cd8 & naive_like),
      pct_pd1_pos_within_tcm_like = rp_pct(pd1, cd3_cd8 & tcm_like),
      pct_pd1_pos_within_tem_like = rp_pct(pd1, cd3_cd8 & tem_like),
      pct_pd1_pos_within_temra_like = rp_pct(pd1, cd3_cd8 & temra_like),
      pct_pd1_pos_within_temra_cd27neg = rp_pct(
        pd1,
        cd3_cd8 & temra_cd27neg
      ),
      pct_pd1_temra = rp_pct(pd1 & temra_like, cd3_cd8),
      pct_pd1_temra_cd27neg = rp_pct(pd1 & temra_cd27neg, cd3_cd8),
      pct_pd1high_temra = rp_pct(pd1_high & temra_like, cd3_cd8),
      pct_pd1high_temra_cd27neg = rp_pct(
        pd1_high & temra_cd27neg,
        cd3_cd8
      ),
      pct_pd1_cd57_pos = rp_pct(pd1 & cd57, cd3_cd8),
      pct_pd1_cx3cr1_pos = rp_pct(pd1 & cx3cr1, cd3_cd8),
      pct_pd1_cd95_pos = rp_pct(pd1 & cd95, cd3_cd8),
      pct_pd1_cd57_cx3cr1_pos = rp_pct(pd1 & cd57 & cx3cr1, cd3_cd8),
      pct_pd1_cd57_cd95_pos = rp_pct(pd1 & cd57 & cd95, cd3_cd8),
      pct_pd1_cx3cr1_cd95_pos = rp_pct(pd1 & cx3cr1 & cd95, cd3_cd8),
      pct_pd1_cd57_cx3cr1_cd95_pos = rp_pct(
        pd1 & cd57 & cx3cr1 & cd95,
        cd3_cd8
      )
    )
  }, error = function(error) {
    tibble::tibble(
      subject_id = rp_extract_subject_id(path),
      result_file_name = basename(path),
      feature_ok = FALSE,
      error_message = conditionMessage(error),
      total_events = NA_integer_,
      n_cd3_pos = NA_integer_,
      n_cd8_pos = NA_integer_,
      n_cd3_cd8_pos = NA_integer_
    )
  })
}

feature_table <- purrr::map_dfr(fcs_files, extract_one)
summary <- tibble::tibble(
  n_files = nrow(feature_table),
  n_subjects = dplyr::n_distinct(feature_table$subject_id),
  n_success = sum(feature_table$feature_ok, na.rm = TRUE),
  n_failed = sum(!feature_table$feature_ok, na.rm = TRUE),
  median_total_events = median(feature_table$total_events, na.rm = TRUE),
  median_cd3_cd8_events = median(feature_table$n_cd3_cd8_pos, na.rm = TRUE),
  n_ge500 = sum(feature_table$n_cd3_cd8_pos >= 500, na.rm = TRUE),
  n_ge1000 = sum(feature_table$n_cd3_cd8_pos >= 1000, na.rm = TRUE)
)

readr::write_csv(
  feature_table,
  file.path(out_dir, "SDY2583_CP24_FULL_850_feature_table_RECONSTRUCTED.csv")
)
readr::write_csv(
  summary,
  file.path(out_dir, "SDY2583_CP24_feature_extraction_summary_RECONSTRUCTED.csv")
)
save(
  feature_table,
  summary,
  thresholds_main,
  extract_one,
  file = file.path(
    rdata_dir,
    "SDY2583_CP24_STEP2_feature_extraction_RECONSTRUCTED.RData"
  )
)
print(summary)
