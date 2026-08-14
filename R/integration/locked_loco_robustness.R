# SDY2583 locked leave-one-cancer-type-out (LOCO) robustness analysis
# Public aggregate-output implementation.
#
# Expected objects in the calling environment:
#   S     - participant-level locked principal-score matrix created by
#           R/integration/locked_internal_validation.R
#   integ - local integrated ALL10 matrix containing cancer_subgroup_model
#
# Critical rule: the locked scores are used unchanged. No score membership,
# direction, centering, scaling, or weighting is re-estimated after any cancer
# subgroup is omitted.

req_loco <- c("dplyr", "tidyr", "purrr", "readr", "tibble", "lmtest", "sandwich")
miss_loco <- req_loco[!vapply(req_loco, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss_loco)) stop("Missing packages for LOCO: ", paste(miss_loco, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(lmtest)
  library(sandwich)
})

if (!exists("S", inherits = TRUE)) {
  stop("Object `S` is missing. Run/source locked_internal_validation.R first.")
}
if (!exists("integ", inherits = TRUE)) {
  stop("Object `integ` is missing. Run/source locked_internal_validation.R first.")
}
if (!"cancer_subgroup_model" %in% names(integ)) {
  stop("Integrated matrix is missing `cancer_subgroup_model`.")
}

repo_root_loco <- Sys.getenv("SDY2583_PUBLIC_REPO_ROOT", "")
if (!nzchar(repo_root_loco)) repo_root_loco <- getwd()
out_dir_loco <- Sys.getenv(
  "SDY2583_LOCKED_LOCO_OUTPUT_DIR",
  file.path(repo_root_loco, "results", "local", "locked_loco_robustness")
)
dir.create(out_dir_loco, recursive = TRUE, showWarnings = FALSE)

meta_loco <- integ %>%
  select(subject_id, cancer_subgroup_model) %>%
  distinct(subject_id, .keep_all = TRUE)

dat_loco <- S %>% left_join(meta_loco, by = "subject_id")

principal_scores_loco <- if (exists("prin", inherits = TRUE)) {
  as.character(prin$score)
} else {
  grep("^CP[0-9]+_integrated_.*_score$", names(dat_loco), value = TRUE)
}

if (length(principal_scores_loco) != 10L) {
  stop("Expected 10 locked principal scores; found ", length(principal_scores_loco))
}

default_groups_loco <- c(
  "Sarkom", "Diğer kanserler", "Meme", "Kolorektal",
  "Pankreas", "Akciğer", "Deri", "Prostat"
)
groups_env_loco <- Sys.getenv("SDY2583_LOCO_CANCER_GROUPS", "")
cancer_groups_loco <- if (nzchar(groups_env_loco)) {
  trimws(strsplit(groups_env_loco, ",", fixed = TRUE)[[1]])
} else {
  default_groups_loco
}

missing_groups_loco <- setdiff(
  cancer_groups_loco,
  unique(as.character(dat_loco$cancer_subgroup_model))
)
if (length(missing_groups_loco)) {
  stop("Requested cancer subgroup labels not found: ", paste(missing_groups_loco, collapse = ", "))
}

fit_hc3_loco <- function(data, score) {
  d <- data %>%
    transmute(
      y = suppressWarnings(as.numeric(.data[[score]])),
      disease = as.character(disease_group_model),
      age = suppressWarnings(as.numeric(age_for_clinical_model)),
      sex = factor(as.character(sex_for_clinical_model))
    ) %>%
    filter(
      disease %in% c("Cancer patient", "Healthy control"),
      !is.na(y), !is.na(age), !is.na(sex)
    ) %>%
    mutate(cancer = ifelse(disease == "Cancer patient", 1, 0))

  if (nrow(d) < 20L || length(unique(d$cancer)) < 2L) {
    return(tibble(
      n = nrow(d), n_cancer = sum(d$cancer == 1), n_healthy = sum(d$cancer == 0),
      beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_
    ))
  }

  fit <- lm(y ~ cancer + age + sex, data = d)
  ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  b <- unname(ct["cancer", 1])
  se <- unname(ct["cancer", 2])
  p <- unname(ct["cancer", 4])
  z <- qnorm(0.975)

  tibble(
    n = nrow(d),
    n_cancer = sum(d$cancer == 1),
    n_healthy = sum(d$cancer == 0),
    beta = b,
    se = se,
    ci_low = b - z * se,
    ci_high = b + z * se,
    p = p
  )
}

baseline_loco <- map_dfr(principal_scores_loco, function(s) {
  fit_hc3_loco(dat_loco, s) %>%
    mutate(panel = sub("_.*$", "", s), score = s, .before = 1)
}) %>%
  mutate(q_bh_10 = p.adjust(p, method = "BH"))

loco_batches <- map(cancer_groups_loco, function(g) {
  omit_n <- sum(
    dat_loco$disease_group_model == "Cancer patient" &
      dat_loco$cancer_subgroup_model == g,
    na.rm = TRUE
  )

  d <- dat_loco %>%
    filter(!(disease_group_model == "Cancer patient" & cancer_subgroup_model == g))

  ans <- map_dfr(principal_scores_loco, function(s) {
    fit_hc3_loco(d, s) %>%
      mutate(
        panel = sub("_.*$", "", s),
        score = s,
        omitted_cancer_group = g,
        omitted_n = omit_n,
        .before = 1
      )
  })
  ans$q_bh_10_within_omission <- p.adjust(ans$p, method = "BH")
  ans
})

loco <- bind_rows(loco_batches)
loco$q_bh_80_global <- p.adjust(loco$p, method = "BH")

baseline_map_loco <- baseline_loco %>% select(score, baseline_beta = beta)
loco <- loco %>%
  left_join(baseline_map_loco, by = "score") %>%
  mutate(
    beta_ratio_to_baseline = beta / baseline_beta,
    relative_beta_change_pct = (beta_ratio_to_baseline - 1) * 100,
    direction_preserved = sign(beta) == sign(baseline_beta),
    positive_beta = beta > 0,
    ci_excludes_zero = ci_low > 0 | ci_high < 0,
    within_omission_FDR_sig = q_bh_10_within_omission < 0.05,
    global_80_FDR_sig = q_bh_80_global < 0.05
  )

score_summary_loco <- loco %>%
  group_by(panel, score) %>%
  summarise(
    baseline_beta = first(baseline_beta),
    min_LOCO_beta = min(beta, na.rm = TRUE),
    max_LOCO_beta = max(beta, na.rm = TRUE),
    all_8_directions_preserved = all(direction_preserved, na.rm = TRUE),
    all_8_beta_positive = all(positive_beta, na.rm = TRUE),
    n_of_8_within_omission_FDR_sig = sum(within_omission_FDR_sig, na.rm = TRUE),
    n_of_8_global_80_FDR_sig = sum(global_80_FDR_sig, na.rm = TRUE),
    worst_within_omission_q = max(q_bh_10_within_omission, na.rm = TRUE),
    omission_with_min_beta = omitted_cancer_group[which.min(beta)],
    largest_abs_relative_change_pct = relative_beta_change_pct[which.max(abs(relative_beta_change_pct))],
    omission_with_largest_abs_change = omitted_cancer_group[which.max(abs(relative_beta_change_pct))],
    .groups = "drop"
  )

omission_summary_loco <- loco %>%
  group_by(omitted_cancer_group, omitted_n) %>%
  summarise(
    n_principal_scores = n(),
    n_direction_preserved = sum(direction_preserved, na.rm = TRUE),
    n_positive_beta = sum(positive_beta, na.rm = TRUE),
    n_within_omission_FDR_sig = sum(within_omission_FDR_sig, na.rm = TRUE),
    n_global_80_FDR_sig = sum(global_80_FDR_sig, na.rm = TRUE),
    min_beta = min(beta, na.rm = TRUE),
    max_beta = max(beta, na.rm = TRUE),
    max_within_omission_q = max(q_bh_10_within_omission, na.rm = TRUE),
    .groups = "drop"
  )

counts_loco <- dat_loco %>%
  filter(disease_group_model == "Cancer patient") %>%
  count(cancer_subgroup_model, name = "n") %>%
  filter(cancer_subgroup_model %in% cancer_groups_loco) %>%
  arrange(match(cancer_subgroup_model, cancer_groups_loco))

beta_matrix_loco <- loco %>%
  select(panel, omitted_cancer_group, beta) %>%
  pivot_wider(names_from = omitted_cancer_group, values_from = beta)

write_csv(baseline_loco, file.path(out_dir_loco, "locked_loco_baseline_results.csv"))
write_csv(loco, file.path(out_dir_loco, "locked_loco_principal_results.csv"))
write_csv(score_summary_loco, file.path(out_dir_loco, "locked_loco_score_summary.csv"))
write_csv(omission_summary_loco, file.path(out_dir_loco, "locked_loco_omission_summary.csv"))
write_csv(counts_loco, file.path(out_dir_loco, "locked_loco_cancer_group_counts.csv"))
write_csv(beta_matrix_loco, file.path(out_dir_loco, "locked_loco_beta_matrix.csv"))

cat("Locked LOCO robustness analysis complete. Aggregate outputs: ", out_dir_loco, "\n", sep = "")
cat("Direction preserved: ", sum(loco$direction_preserved, na.rm = TRUE), "/", nrow(loco), " models\n", sep = "")
cat("Within-omission BH-FDR significant: ", sum(loco$within_omission_FDR_sig, na.rm = TRUE), "/", nrow(loco), "\n", sep = "")
cat("Global 80-test BH-FDR significant: ", sum(loco$global_80_FDR_sig, na.rm = TRUE), "/", nrow(loco), "\n", sep = "")
