# Shared framework for provenance-tracked reconstructed SDY2583 panel workflows.
# This file contains generic mechanics only. Panel biology, marker mapping,
# gating expressions, feature modules, composite directions, and benchmarks
# must remain in panel-specific manifests/scripts.

source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

rp_install_and_load <- function(cran = character(), bioc = character()) {
  for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  if (length(bioc) > 0L) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
  }
  invisible(lapply(c(cran, bioc), function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

rp_extract_subject_id <- function(path) {
  id <- stringr::str_extract(basename(path), "DBG[0-9]+")
  ifelse(is.na(id), stringr::str_remove(basename(path), "\\.fcs$"), id)
}

rp_marker_map <- function(ff) {
  pp <- Biobase::pData(flowCore::parameters(ff))
  channel <- as.character(pp$name)
  marker <- as.character(pp$desc)
  marker[is.na(marker) | marker == ""] <- channel[is.na(marker) | marker == ""]
  tibble::tibble(
    parameter_index = seq_along(channel),
    channel = channel,
    marker = marker,
    range = suppressWarnings(as.numeric(pp$range)),
    minRange = suppressWarnings(as.numeric(pp$minRange)),
    maxRange = suppressWarnings(as.numeric(pp$maxRange))
  )
}

rp_find_channel <- function(map, patterns, fallback = NA_character_) {
  for (pattern in patterns) {
    idx <- which(grepl(pattern, map$marker, ignore.case = TRUE))
    if (length(idx) > 0L) return(map$channel[idx[[1L]]])
  }
  if (!is.na(fallback) && fallback %in% map$channel) return(fallback)
  NA_character_
}

rp_spill_matrix <- function(ff) {
  keys <- flowCore::keyword(ff)
  for (name in c("SPILL", "$SPILLOVER", "SPILLOVER")) {
    if (name %in% names(keys) && is.matrix(keys[[name]])) return(keys[[name]])
  }
  NULL
}

rp_compensate <- function(ff) {
  spill <- rp_spill_matrix(ff)
  if (is.null(spill)) return(list(ff = ff, applied = FALSE, error = NA_character_))
  tryCatch(
    list(ff = flowCore::compensate(ff, spill), applied = TRUE, error = NA_character_),
    error = function(e) list(ff = ff, applied = FALSE, error = conditionMessage(e))
  )
}

rp_logicle <- function(ff, channels = NULL) {
  all_channels <- colnames(flowCore::exprs(ff))
  if (is.null(channels)) channels <- all_channels[!grepl("FSC|SSC|Time", all_channels, ignore.case = TRUE)]
  transform <- flowCore::logicleTransform(
    transformationId = "fixed_logicle", w = 0.5, t = 262144, m = 4.5, a = 0
  )
  tryCatch(
    list(ff = flowCore::transform(ff, flowCore::transformList(channels, transform)), applied = TRUE, error = NA_character_),
    error = function(e) list(ff = ff, applied = FALSE, error = conditionMessage(e))
  )
}

rp_read_transform <- function(path) {
  ff_raw <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
  compensation <- rp_compensate(ff_raw)
  transformed <- rp_logicle(compensation$ff)
  list(
    raw = ff_raw,
    ff = transformed$ff,
    marker_map = rp_marker_map(ff_raw),
    compensation_applied = compensation$applied,
    compensation_error = compensation$error,
    transform_applied = transformed$applied,
    transform_error = transformed$error
  )
}

rp_vec <- function(exprs, channel) {
  if (is.na(channel) || !(channel %in% colnames(exprs))) return(rep(NA_real_, nrow(exprs)))
  as.numeric(exprs[, channel])
}

rp_n <- function(mask) sum(mask, na.rm = TRUE)
rp_pct <- function(mask, denominator = NULL) {
  if (is.null(denominator)) denominator <- rep(TRUE, length(mask))
  n_denominator <- sum(denominator, na.rm = TRUE)
  if (n_denominator <= 0L) return(NA_real_)
  100 * sum(mask & denominator, na.rm = TRUE) / n_denominator
}
rp_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}
rp_ratio <- function(numerator, denominator) {
  if (is.na(denominator) || denominator <= 0) return(NA_real_)
  numerator / denominator
}
rp_z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

rp_inventory <- function(panel, expected_markers = character(), fallback_channels = list()) {
  rp_install_and_load(c("dplyr", "readr", "stringr", "tibble", "purrr", "tidyr"), "flowCore")
  panel <- toupper(panel)
  analysis_dir <- sd_analysis_dir(panel)
  out_dir <- file.path(analysis_dir, "01_inventory_marker_QC")
  rdata_dir <- file.path(analysis_dir, "11_RData")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)

  fcs_files <- list.files(sd_fcs_dir(panel), "\\.fcs$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  fcs_files <- sort(unique(fcs_files[grepl(panel, basename(fcs_files), ignore.case = TRUE)]))
  if (length(fcs_files) == 0L) stop("No ", panel, " FCS files found under ", sd_fcs_dir(panel))

  read_one <- function(path) {
    tryCatch({
      ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
      map <- rp_marker_map(ff)
      keys <- names(flowCore::keyword(ff))
      tibble::tibble(
        file_name = basename(path), file_path = normalizePath(path, mustWork = FALSE),
        subject_id = rp_extract_subject_id(path), panel = panel, read_ok = TRUE,
        read_error = NA_character_, n_events = nrow(flowCore::exprs(ff)),
        n_channels = ncol(flowCore::exprs(ff)),
        has_spillover = any(c("SPILL", "$SPILLOVER", "SPILLOVER") %in% keys),
        channel_signature = paste(map$channel, collapse = "|"),
        marker_signature = paste(map$marker, collapse = "|"),
        marker_set_signature = paste(sort(unique(map$marker)), collapse = "|")
      )
    }, error = function(e) tibble::tibble(
      file_name = basename(path), file_path = normalizePath(path, mustWork = FALSE),
      subject_id = rp_extract_subject_id(path), panel = panel, read_ok = FALSE,
      read_error = conditionMessage(e), n_events = NA_integer_, n_channels = NA_integer_,
      has_spillover = NA, channel_signature = NA_character_, marker_signature = NA_character_,
      marker_set_signature = NA_character_
    ))
  }

  inventory <- purrr::map_dfr(fcs_files, read_one)
  if (!any(inventory$read_ok)) stop("No ", panel, " file could be read.")
  modal_channels <- inventory |> dplyr::filter(read_ok) |> dplyr::count(n_channels) |>
    dplyr::arrange(dplyr::desc(n), n_channels) |> dplyr::slice(1) |> dplyr::pull(n_channels)
  reference_file <- inventory |> dplyr::filter(read_ok, n_channels == modal_channels) |>
    dplyr::arrange(file_name) |> dplyr::slice(1) |> dplyr::pull(file_path)
  reference_map <- rp_marker_map(flowCore::read.FCS(reference_file, transformation = FALSE, truncate_max_range = FALSE))

  reference_channel_signature <- paste(reference_map$channel, collapse = "|")
  reference_marker_signature <- paste(reference_map$marker, collapse = "|")
  reference_marker_set_signature <- paste(sort(unique(reference_map$marker)), collapse = "|")
  inventory <- inventory |> dplyr::mutate(
    channel_order_mismatch = read_ok & channel_signature != reference_channel_signature,
    marker_order_mismatch = read_ok & marker_signature != reference_marker_signature,
    marker_set_mismatch = read_ok & marker_set_signature != reference_marker_set_signature
  )

  marker_presence <- tibble::tibble(marker = expected_markers) |> dplyr::rowwise() |> dplyr::mutate(
    reference_channel = rp_find_channel(reference_map, marker, fallback = fallback_channels[[marker]] %||% NA_character_),
    present = !is.na(reference_channel)
  ) |> dplyr::ungroup()

  summary <- tibble::tibble(
    panel = panel, n_files = nrow(inventory), n_subjects = dplyr::n_distinct(inventory$subject_id),
    n_read_ok = sum(inventory$read_ok, na.rm = TRUE), n_read_failed = sum(!inventory$read_ok, na.rm = TRUE),
    n_with_spillover = sum(inventory$has_spillover, na.rm = TRUE), modal_n_channels = modal_channels,
    median_events = stats::median(inventory$n_events, na.rm = TRUE),
    min_events = min(inventory$n_events, na.rm = TRUE), max_events = max(inventory$n_events, na.rm = TRUE),
    n_channel_patterns = dplyr::n_distinct(inventory$channel_signature, na.rm = TRUE),
    n_marker_patterns = dplyr::n_distinct(inventory$marker_signature, na.rm = TRUE),
    n_channel_order_mismatch = sum(inventory$channel_order_mismatch, na.rm = TRUE),
    n_marker_order_mismatch = sum(inventory$marker_order_mismatch, na.rm = TRUE),
    n_marker_set_mismatch = sum(inventory$marker_set_mismatch, na.rm = TRUE)
  )

  readr::write_csv(inventory, file.path(out_dir, paste0("SDY2583_", panel, "_FCS_inventory_RECONSTRUCTED.csv")))
  readr::write_csv(reference_map, file.path(out_dir, paste0("SDY2583_", panel, "_reference_marker_map_RECONSTRUCTED.csv")))
  readr::write_csv(marker_presence, file.path(out_dir, paste0("SDY2583_", panel, "_expected_marker_presence_RECONSTRUCTED.csv")))
  readr::write_csv(summary, file.path(out_dir, paste0("SDY2583_", panel, "_inventory_QC_summary_RECONSTRUCTED.csv")))
  save(fcs_files, inventory, reference_map, marker_presence, summary,
       file = file.path(rdata_dir, paste0("SDY2583_", panel, "_STEP1_inventory_RECONSTRUCTED.RData")))
  invisible(list(fcs_files = fcs_files, inventory = inventory, reference_map = reference_map, summary = summary))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

rp_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

rp_standardize_disease <- function(x) {
  y <- trimws(as.character(x))
  dplyr::case_when(
    stringr::str_detect(y, stringr::regex("healthy|control", ignore_case = TRUE)) ~ "Healthy control",
    stringr::str_detect(y, stringr::regex("cancer|patient|tumou?r", ignore_case = TRUE)) ~ "Cancer patient",
    TRUE ~ y
  )
}

rp_standardize_sex <- function(x) {
  y <- trimws(as.character(x))
  dplyr::case_when(
    tolower(y) %in% c("female", "f", "woman") ~ "Female",
    tolower(y) %in% c("male", "m", "man") ~ "Male",
    is.na(y) | y == "" | tolower(y) %in% c("na", "unknown") ~ NA_character_,
    TRUE ~ y
  )
}

rp_find_metadata_file <- function() {
  explicit <- path.expand(Sys.getenv("SDY2583_METADATA_MATRIX_FILE", unset = ""))
  if (nzchar(explicit) && file.exists(explicit)) return(explicit)
  dirs <- unique(c(sd_integrated_dir(), sd_metadata_dir()))
  dirs <- dirs[dir.exists(dirs)]
  files <- unique(unlist(lapply(dirs, function(d) list.files(d, "\\.csv$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE))))
  if (length(files) == 0L) stop("No metadata CSV found. Set SDY2583_METADATA_MATRIX_FILE.")
  inspect <- function(path) {
    x <- tryCatch(readr::read_csv(path, n_max = 5, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(x)) return(tibble::tibble(path = path, score = -Inf))
    subject <- rp_first_col(x, c("subject_id", "dbg_id", "subject"))
    disease <- rp_first_col(x, c("disease_group", "disease_group_model", "disease_group_clinical", "group", "condition", "arm_name"))
    age <- rp_first_col(x, c("age_for_model", "age_raw", "age_clinical", "age_years", "age"))
    sex <- rp_first_col(x, c("sex", "sex_for_clinical_model", "sex_clinical", "gender"))
    score <- 1000 * !is.na(subject) + 300 * !is.na(disease) + 300 * !is.na(age) + 300 * !is.na(sex)
    score <- score + 100 * grepl("integrated|clinical|metadata|matrix", basename(path), ignore.case = TRUE)
    score <- score - 500 * grepl("summary|statistics|result|figure|dictionary|manifest", basename(path), ignore.case = TRUE)
    tibble::tibble(path = path, score = score, subject = subject, disease = disease, age = age, sex = sex)
  }
  candidates <- purrr::map_dfr(files, inspect) |> dplyr::arrange(dplyr::desc(score), path)
  valid <- candidates |> dplyr::filter(is.finite(score), score >= 1900)
  if (nrow(valid) == 0L) stop("No compatible metadata CSV found.")
  valid$path[[1L]]
}

rp_merge_metadata <- function(panel, feature_table, feature_columns = NULL) {
  rp_install_and_load(c("dplyr", "readr", "stringr", "tibble", "purrr"))
  panel <- toupper(panel)
  source_file <- rp_find_metadata_file()
  raw <- readr::read_csv(source_file, show_col_types = FALSE, progress = FALSE)
  subject_col <- rp_first_col(raw, c("subject_id", "dbg_id", "subject"))
  disease_col <- rp_first_col(raw, c("disease_group", "disease_group_model", "disease_group_clinical", "group", "condition", "arm_name"))
  age_col <- rp_first_col(raw, c("age_for_model", "age_raw", "age_clinical", "age_years", "age"))
  sex_col <- rp_first_col(raw, c("sex", "sex_for_clinical_model", "sex_clinical", "gender"))
  metadata <- raw |> dplyr::transmute(
    subject_id = as.character(.data[[subject_col]]),
    disease_group = rp_standardize_disease(.data[[disease_col]]),
    age_raw = suppressWarnings(as.numeric(.data[[age_col]])),
    sex = rp_standardize_sex(.data[[sex_col]])
  ) |> dplyr::filter(!is.na(subject_id), subject_id != "") |>
    dplyr::group_by(subject_id) |> dplyr::summarise(
      disease_group = dplyr::first(stats::na.omit(disease_group), default = NA_character_),
      age_raw = dplyr::first(stats::na.omit(age_raw), default = NA_real_),
      sex = dplyr::first(stats::na.omit(sex), default = NA_character_), .groups = "drop"
    ) |> dplyr::mutate(
      disease_group = factor(disease_group, levels = c("Healthy control", "Cancer patient")),
      age_for_model = ifelse(age_raw >= 18 & age_raw <= 100, age_raw, NA_real_),
      age_group_for_model = factor(dplyr::case_when(
        is.na(age_for_model) ~ NA_character_, age_for_model < 40 ~ "Young_<40",
        age_for_model < 60 ~ "Middle_40_59", TRUE ~ "Older_60plus"
      ), levels = c("Young_<40", "Middle_40_59", "Older_60plus")),
      sex = factor(sex),
      sex_binary = factor(ifelse(as.character(sex) %in% c("Female", "Male"), as.character(sex), NA_character_), levels = c("Female", "Male"))
    )
  merged <- feature_table |> dplyr::mutate(subject_id = as.character(subject_id)) |> dplyr::left_join(metadata, by = "subject_id")
  if (is.null(feature_columns)) feature_columns <- names(merged)[vapply(merged, is.numeric, logical(1)) & grepl("^(pct_|median_|ratio_|score_)", names(merged))]
  list(data = merged, metadata = metadata, feature_columns = feature_columns, source_file = source_file)
}

rp_fit_one <- function(data, outcome, subset_label = "main_valid_age_all_sex_categories", binary_sex = FALSE, event_col = NULL, min_events = 0) {
  d <- data |> dplyr::transmute(
    y = suppressWarnings(as.numeric(.data[[outcome]])),
    disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")),
    age_for_model = suppressWarnings(as.numeric(age_for_model)), sex = factor(as.character(sex)),
    event_count = if (is.null(event_col)) Inf else suppressWarnings(as.numeric(.data[[event_col]]))
  ) |> dplyr::filter(!is.na(y), !is.na(disease_group), !is.na(age_for_model), !is.na(sex), !is.na(event_count), event_count >= min_events)
  if (binary_sex) d <- d |> dplyr::filter(as.character(sex) %in% c("Female", "Male")) |> dplyr::mutate(sex = droplevels(sex))
  if (nrow(d) < 10L || dplyr::n_distinct(d$disease_group) < 2L) return(tibble::tibble(
    subset_label = subset_label, feature = outcome, n_model = nrow(d), n_healthy = sum(d$disease_group == "Healthy control"),
    n_cancer = sum(d$disease_group == "Cancer patient"), healthy_mean = NA_real_, cancer_mean = NA_real_,
    healthy_median = NA_real_, cancer_median = NA_real_, crude_mean_difference_cancer_minus_healthy = NA_real_,
    beta_cancer_vs_healthy = NA_real_, ci_low = NA_real_, ci_high = NA_real_, t_value = NA_real_, p_value = NA_real_,
    model_formula = "y ~ disease_group + age_for_model + sex", direction = NA_character_
  ))
  fit <- stats::lm(y ~ disease_group + age_for_model + sex, data = d)
  ct <- summary(fit)$coefficients
  cn <- grep("^disease_group", rownames(ct), value = TRUE)[1]
  ci <- stats::confint(fit, parm = cn)
  h <- d$y[d$disease_group == "Healthy control"]; c <- d$y[d$disease_group == "Cancer patient"]
  beta <- unname(ct[cn, "Estimate"])
  tibble::tibble(
    subset_label = subset_label, feature = outcome, n_model = stats::nobs(fit), n_healthy = length(h), n_cancer = length(c),
    healthy_mean = mean(h), cancer_mean = mean(c), healthy_median = stats::median(h), cancer_median = stats::median(c),
    crude_mean_difference_cancer_minus_healthy = mean(c) - mean(h), beta_cancer_vs_healthy = beta,
    ci_low = unname(ci[1]), ci_high = unname(ci[2]), t_value = unname(ct[cn, "t value"]),
    p_value = unname(ct[cn, "Pr(>|t|)"]), model_formula = "y ~ disease_group + age_for_model + sex",
    direction = ifelse(beta > 0, "higher_in_cancer", ifelse(beta < 0, "lower_in_cancer", "no_difference"))
  )
}

rp_fit_feature_set <- function(data, module_map, subset_label = "main_valid_age_all_sex_categories", binary_sex = FALSE, event_col = NULL, min_events = 0) {
  result <- purrr::map_dfr(module_map$feature, ~rp_fit_one(data, .x, subset_label, binary_sex, event_col, min_events)) |>
    dplyr::left_join(module_map, by = "feature") |> dplyr::mutate(fdr_all = p.adjust(p_value, "BH")) |>
    dplyr::group_by(module) |> dplyr::mutate(fdr_within_module = p.adjust(p_value, "BH")) |> dplyr::ungroup()
  result |> dplyr::select(subset_label, module, dplyr::everything()) |> dplyr::arrange(p_value)
}

rp_build_scores <- function(data, score_definitions, integrated_name = NULL) {
  scored <- data
  for (name in names(score_definitions)) {
    def <- score_definitions[[name]]
    components <- list()
    for (f in def$positive %||% character()) components[[paste0("plus__", f)]] <- rp_z(scored[[f]])
    for (f in def$negative %||% character()) components[[paste0("minus__", f)]] <- -rp_z(scored[[f]])
    m <- as.data.frame(components, check.names = FALSE)
    scored[[name]] <- rowMeans(m, na.rm = TRUE)
    scored[[name]][rowSums(!is.na(m)) == 0L] <- NA_real_
  }
  if (!is.null(integrated_name)) {
    components <- names(score_definitions)
    m <- as.data.frame(lapply(scored[components], rp_z), check.names = FALSE)
    scored[[integrated_name]] <- rowMeans(m, na.rm = TRUE)
    scored[[integrated_name]][rowSums(!is.na(m)) == 0L] <- NA_real_
  }
  scored
}

rp_fit_scores <- function(scored_data, score_names, event_col = NULL, min_events = 0, binary_sex = FALSE) {
  result <- purrr::map_dfr(score_names, function(name) {
    x <- rp_fit_one(scored_data, name, subset_label = ifelse(binary_sex, "binary_sex_only", "main"), binary_sex = binary_sex, event_col = event_col, min_events = min_events)
    dplyr::rename(x, score = feature)
  }) |> dplyr::mutate(fdr = p.adjust(p_value, "BH")) |> dplyr::arrange(p_value)
  result
}

rp_greedy_match <- function(data, caliper) {
  d <- data |> dplyr::filter(!is.na(age_for_model), as.character(sex) %in% c("Female", "Male"), !is.na(disease_group))
  cancer <- d |> dplyr::filter(disease_group == "Cancer patient") |> dplyr::arrange(sex, age_for_model, subject_id)
  healthy <- d |> dplyr::filter(disease_group == "Healthy control") |> dplyr::arrange(sex, age_for_model, subject_id)
  used <- rep(FALSE, nrow(healthy)); pairs <- vector("list", 0L)
  for (i in seq_len(nrow(cancer))) {
    candidates <- which(!used & as.character(healthy$sex) == as.character(cancer$sex[i]))
    if (length(candidates) == 0L) next
    distances <- abs(healthy$age_for_model[candidates] - cancer$age_for_model[i])
    j_local <- order(distances, healthy$subject_id[candidates])[1]
    j <- candidates[j_local]
    if (distances[j_local] <= caliper) {
      used[j] <- TRUE
      pairs[[length(pairs) + 1L]] <- tibble::tibble(pair_id = length(pairs) + 1L, cancer_row = i, healthy_row = j, abs_age_difference = distances[j_local])
    }
  }
  pair_table <- dplyr::bind_rows(pairs)
  if (nrow(pair_table) == 0L) return(list(data = d[0, ], pairs = pair_table))
  cancer_m <- cancer[pair_table$cancer_row, ]; healthy_m <- healthy[pair_table$healthy_row, ]
  cancer_m$pair_id <- pair_table$pair_id; healthy_m$pair_id <- pair_table$pair_id
  list(data = dplyr::bind_rows(cancer_m, healthy_m), pairs = pair_table)
}

rp_age_sensitivity <- function(data, outcomes, event_col = NULL, min_events = 0, calipers = c(5, 10)) {
  stratified <- purrr::map_dfr(levels(data$age_group_for_model), function(group) {
    subset <- data |> dplyr::filter(age_group_for_model == group)
    purrr::map_dfr(outcomes, function(outcome) rp_fit_one(subset, outcome, paste0("age_group_", group), FALSE, event_col, min_events))
  }) |> dplyr::mutate(fdr = p.adjust(p_value, "BH"))
  matched_results <- list(); matched_data <- list(); match_summaries <- list()
  for (caliper in calipers) {
    m <- rp_greedy_match(data, caliper)
    matched_data[[as.character(caliper)]] <- m$data
    matched_results[[as.character(caliper)]] <- purrr::map_dfr(outcomes, function(outcome) rp_fit_one(
      m$data, outcome, paste0("same_sex_nearest_age_caliper_", caliper, "y"), FALSE, event_col, min_events
    )) |> dplyr::mutate(fdr = p.adjust(p_value, "BH"))
    match_summaries[[as.character(caliper)]] <- tibble::tibble(
      caliper_years = caliper, n_pairs = nrow(m$pairs), n_rows = nrow(m$data),
      mean_abs_age_difference = mean(m$pairs$abs_age_difference, na.rm = TRUE),
      median_abs_age_difference = stats::median(m$pairs$abs_age_difference, na.rm = TRUE),
      max_abs_age_difference = max(m$pairs$abs_age_difference, na.rm = TRUE)
    )
  }
  interactions <- purrr::map_dfr(outcomes, function(outcome) {
    d <- data |> dplyr::transmute(y = as.numeric(.data[[outcome]]), disease_group = factor(as.character(disease_group), levels = c("Healthy control", "Cancer patient")), age_z = as.numeric(scale(age_for_model)), sex = factor(as.character(sex))) |> dplyr::filter(stats::complete.cases(.))
    fit <- stats::lm(y ~ disease_group * age_z + sex, data = d); ct <- summary(fit)$coefficients
    cn <- grep("disease_group.*:age_z", rownames(ct), value = TRUE)[1]
    tibble::tibble(feature = outcome, n_model = stats::nobs(fit), interaction_beta = unname(ct[cn, "Estimate"]), interaction_t = unname(ct[cn, "t value"]), interaction_p = unname(ct[cn, "Pr(>|t|)"]))
  }) |> dplyr::mutate(interaction_fdr = p.adjust(interaction_p, "BH"))
  list(stratified = stratified, matched = dplyr::bind_rows(matched_results), matched_data = matched_data,
       match_summary = dplyr::bind_rows(match_summaries), interactions = interactions)
}

rp_compare_csv <- function(generated, reference, keys, tolerance = 1e-6) {
  if (!file.exists(generated) || !file.exists(reference)) return(tibble::tibble(pass = FALSE, detail = "missing file"))
  g <- readr::read_csv(generated, show_col_types = FALSE); r <- readr::read_csv(reference, show_col_types = FALSE)
  common <- intersect(names(g), names(r)); keys <- intersect(keys, common)
  if (length(keys) > 0L) { g <- g |> dplyr::arrange(dplyr::across(dplyr::all_of(keys))); r <- r |> dplyr::arrange(dplyr::across(dplyr::all_of(keys))) }
  if (nrow(g) != nrow(r)) return(tibble::tibble(pass = FALSE, detail = paste("row mismatch", nrow(g), nrow(r))))
  failures <- character()
  for (col in common) {
    gv <- g[[col]]; rv <- r[[col]]
    if (is.numeric(gv) && is.numeric(rv)) {
      bad <- !(is.na(gv) & is.na(rv)) & (is.na(gv) != is.na(rv) | abs(gv-rv) > tolerance * pmax(1, abs(rv)))
    } else bad <- ifelse(is.na(gv), "<NA>", as.character(gv)) != ifelse(is.na(rv), "<NA>", as.character(rv))
    if (any(bad, na.rm = TRUE)) failures <- c(failures, paste0(col, "[", sum(bad, na.rm = TRUE), "]"))
  }
  tibble::tibble(pass = length(failures) == 0L, detail = ifelse(length(failures) == 0L, "matched", paste(failures, collapse = "; ")))
}
