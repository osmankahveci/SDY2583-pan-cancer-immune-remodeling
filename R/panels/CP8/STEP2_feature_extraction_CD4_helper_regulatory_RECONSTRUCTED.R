# CP8 Step 2 reconstructed CD4 helper/regulatory feature extraction.
# Rebuilt from the archived threshold table, 46-feature schema, and final
# Methods record. Validate all columns against the archived Step 2 CSV.
rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "reconstructed_panel_framework.R"))
rp_install_and_load(c("dplyr", "readr", "stringr", "tibble", "purrr"), "flowCore")
source(file.path(sd_repo_root(), "R", "panels", "CP8", "MANIFEST_RECONSTRUCTED.R"))

panel <- "CP8"; analysis_dir <- sd_analysis_dir(panel)
rdata_dir <- file.path(analysis_dir, "11_RData"); out_dir <- file.path(analysis_dir, "02_feature_extraction")
threshold_dir <- file.path(analysis_dir, "03_thresholds")
dir.create(rdata_dir, recursive=TRUE, showWarnings=FALSE); dir.create(out_dir, recursive=TRUE, showWarnings=FALSE); dir.create(threshold_dir, recursive=TRUE, showWarnings=FALSE)
step1 <- file.path(rdata_dir, "SDY2583_CP8_STEP1_inventory_RECONSTRUCTED.RData")
if (!file.exists(step1)) stop("Run CP8 Step 1 first.")
load(step1)
th <- CP8_MANIFEST$thresholds

threshold_table <- tibble::tibble(marker=names(th), threshold_transformed_scale=as.numeric(th))
readr::write_csv(threshold_table, file.path(threshold_dir, "SDY2583_CP8_feature_extraction_thresholds_RECONSTRUCTED.csv"))

extract_one <- function(path, thresholds=th) {
  tryCatch({
    obj <- rp_read_transform(path); map <- obj$marker_map; ex <- flowCore::exprs(obj$ff)
    ch <- c(
      CD3=rp_find_channel(map, c("^CD3$")), CD4=rp_find_channel(map, c("^CD4$"), "BB515-A"),
      CD45RA=rp_find_channel(map, c("CD45RA")), CD62L=rp_find_channel(map, c("CD62L")),
      CD27=rp_find_channel(map, c("^CD27$")), CCR6=rp_find_channel(map, c("CCR6", "CD196")),
      IL7RA=rp_find_channel(map, c("IL7RA", "CD127")), CCR4=rp_find_channel(map, c("CCR4", "CD194")),
      CD25=rp_find_channel(map, c("^CD25$")), CXCR5=rp_find_channel(map, c("CXCR5", "CD185"))
    )
    if (any(is.na(ch))) stop("Missing markers: ", paste(names(ch)[is.na(ch)], collapse=", "))
    v <- lapply(ch, function(x) rp_vec(ex, x))
    cd3 <- v$CD3 > thresholds["CD3"]; cd4 <- v$CD4 > thresholds["CD4"]; primary <- cd3 & cd4
    ra <- v$CD45RA > thresholds["CD45RA"]; l62 <- v$CD62L > thresholds["CD62L"]; cd27 <- v$CD27 > thresholds["CD27"]
    ccr6 <- v$CCR6 > thresholds["CCR6"]; il7ra <- v$IL7RA > thresholds["IL7RA"]; il7ra_low <- v$IL7RA <= thresholds["IL7RA"]
    ccr4 <- v$CCR4 > thresholds["CCR4"]; cd25 <- v$CD25 > thresholds["CD25"]; cd25h <- v$CD25 > thresholds["CD25_HIGH"]
    cxcr5 <- v$CXCR5 > thresholds["CXCR5"]
    naive <- ra & l62; tcm <- !ra & l62; tem <- !ra & !l62; temra <- ra & !l62; memory <- !ra
    treg <- cd25 & il7ra_low; treg_enriched <- cd25h & il7ra_low
    tibble::tibble(
      subject_id=rp_extract_subject_id(path), result_file_name=basename(path), feature_ok=TRUE, error_message=NA_character_,
      total_events=nrow(ex), n_cd3_pos=rp_n(cd3), n_cd4_pos=rp_n(cd4), n_cd3_cd4_pos=rp_n(primary),
      pct_cd3_pos_total=rp_pct(cd3), pct_cd4_pos_total=rp_pct(cd4), pct_cd3_cd4_pos_total=rp_pct(primary), pct_cd4_within_cd3=rp_pct(cd4, cd3),
      median_CD45RA_in_CD3CD4=rp_median(v$CD45RA[primary]), median_CD62L_in_CD3CD4=rp_median(v$CD62L[primary]),
      median_CD27_in_CD3CD4=rp_median(v$CD27[primary]), median_CCR6_in_CD3CD4=rp_median(v$CCR6[primary]),
      median_IL7RA_in_CD3CD4=rp_median(v$IL7RA[primary]), median_CCR4_in_CD3CD4=rp_median(v$CCR4[primary]),
      median_CD25_in_CD3CD4=rp_median(v$CD25[primary]), median_CXCR5_in_CD3CD4=rp_median(v$CXCR5[primary]),
      pct_naive_like=rp_pct(naive, primary), pct_tcm_like=rp_pct(tcm, primary), pct_tem_like=rp_pct(tem, primary), pct_temra_like=rp_pct(temra, primary),
      pct_memory_like=rp_pct(memory, primary), pct_cd27_pos=rp_pct(cd27, primary), pct_cd27_neg=rp_pct(!cd27, primary),
      pct_naive_cd27pos=rp_pct(naive & cd27, primary), pct_memory_cd27pos=rp_pct(memory & cd27, primary),
      pct_memory_cd27neg=rp_pct(memory & !cd27, primary), pct_cd62lneg_cd27neg=rp_pct(!l62 & !cd27, primary),
      pct_ccr6_pos=rp_pct(ccr6, primary), pct_il7ra_pos=rp_pct(il7ra, primary), pct_il7ra_low=rp_pct(il7ra_low, primary),
      pct_ccr4_pos=rp_pct(ccr4, primary), pct_cd25_pos=rp_pct(cd25, primary), pct_cd25_high=rp_pct(cd25h, primary), pct_cxcr5_pos=rp_pct(cxcr5, primary),
      pct_cd25pos_il7ralow_treg_like=rp_pct(treg, primary), pct_cd25high_il7ralow_treg_enriched=rp_pct(treg_enriched, primary),
      pct_ccr4_cd25pos_il7ralow=rp_pct(ccr4 & treg, primary), pct_ccr4_cd25high_il7ralow=rp_pct(ccr4 & treg_enriched, primary),
      pct_cxcr5pos_tfh_like=rp_pct(cxcr5, primary), pct_ccr6pos_th17_like=rp_pct(ccr6, primary),
      pct_ccr6pos_cxcr5pos=rp_pct(ccr6 & cxcr5, primary), pct_ccr4pos_ccr6pos=rp_pct(ccr4 & ccr6, primary),
      pct_ccr4pos_cxcr5pos=rp_pct(ccr4 & cxcr5, primary), pct_cd25pos_within_cxcr5pos=rp_pct(cd25, primary & cxcr5),
      pct_ccr6pos_within_cxcr5pos=rp_pct(ccr6, primary & cxcr5), pct_cxcr5pos_within_ccr6pos=rp_pct(cxcr5, primary & ccr6),
      pct_ccr4pos_within_treg_like=rp_pct(ccr4, primary & treg), pct_cd25pos_il7ralow_within_memory=rp_pct(treg, primary & memory)
    )
  }, error=function(e) tibble::tibble(subject_id=rp_extract_subject_id(path), result_file_name=basename(path), feature_ok=FALSE, error_message=conditionMessage(e), total_events=NA_integer_, n_cd3_pos=NA_integer_, n_cd4_pos=NA_integer_, n_cd3_cd4_pos=NA_integer_))
}

feature_table <- purrr::map_dfr(fcs_files, extract_one)
summary <- tibble::tibble(
  n_files=nrow(feature_table), n_unique_subjects=dplyr::n_distinct(feature_table$subject_id), n_success=sum(feature_table$feature_ok), n_failed=sum(!feature_table$feature_ok),
  median_total_events=median(feature_table$total_events, na.rm=TRUE), min_total_events=min(feature_table$total_events, na.rm=TRUE), max_total_events=max(feature_table$total_events, na.rm=TRUE),
  median_cd3_cd4_events=median(feature_table$n_cd3_cd4_pos, na.rm=TRUE), min_cd3_cd4_events=min(feature_table$n_cd3_cd4_pos, na.rm=TRUE), max_cd3_cd4_events=max(feature_table$n_cd3_cd4_pos, na.rm=TRUE)
)
readr::write_csv(feature_table, file.path(out_dir, "SDY2583_CP8_FULL_850_feature_table_STEP2_RECONSTRUCTED.csv"))
readr::write_csv(summary, file.path(out_dir, "SDY2583_CP8_feature_extraction_summary_STEP2_RECONSTRUCTED.csv"))
save(feature_table, summary, th, extract_one, file=file.path(rdata_dir, "SDY2583_CP8_STEP2_feature_extraction_RECONSTRUCTED.RData"))
print(summary)
message("Archived benchmarks: 850/850; median total events about 299655; median CD3+CD4+ events about 26156.")
