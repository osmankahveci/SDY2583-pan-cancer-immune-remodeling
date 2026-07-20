# ============================================================
# SDY2583 CP22
# STEP 7B SAFE: FlowJo-style representative B-cell gating figures
#
# Creates:
# - representative healthy and cancer subjects based on integrated
#   humoral B-cell remodeling score
# - B-cell gating and phenotype panels:
#   1) Dump vs CD19
#   2) IgD vs CD27
#   3) IgA vs IgG in switched-memory-like B cells
#   4) CD27 vs CD38
#   5) CD10 vs CD24
#   6) CD24 vs CD39
#
# Output:
# outputs/CP22/05_flow_figures_flowjo_style_TALL_LAYOUT_WHITE_BG
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

pkgs <- c("dplyr", "readr", "stringr", "tibble", "ggplot2", "patchwork")

for (p in pkgs) {
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
  library(ggplot2)
  library(patchwork)
})

# Namespace safety
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
left_join <- dplyr::left_join
group_by <- dplyr::group_by
summarise <- dplyr::summarise
slice_head <- dplyr::slice_head
ungroup <- dplyr::ungroup
case_when <- dplyr::case_when

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP22")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "05_flow_figures_flowjo_style_TALL_LAYOUT_WHITE_BG")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

step4_rdata <- file.path(rdata_dir, "SDY2583_CP22_STEP4_composite_scores.RData")
step1_rdata <- file.path(rdata_dir, "SDY2583_CP22_STEP1_fcs_inventory_marker_QC.RData")

if (!file.exists(step4_rdata)) stop("Step 4 RData not found: ", step4_rdata)
if (!file.exists(step1_rdata)) stop("Step 1 RData not found: ", step1_rdata)

load(step1_rdata)
load(step4_rdata)

if (exists("analysis_df")) {
  cp22_plot_source_df <- analysis_df
} else if (exists("cp22_scores_data")) {
  cp22_plot_source_df <- cp22_scores_data
} else if (exists("cp22_analysis_data_with_scores")) {
  cp22_plot_source_df <- cp22_analysis_data_with_scores
} else {
  cat("\nObjects loaded from Step 4 RData:\n")
  print(ls())
  stop("Step 4 RData içinde analysis_df / cp22_scores_data / cp22_analysis_data_with_scores bulunamadı.")
}

if (!("file_path" %in% names(cp22_plot_source_df))) {
  stop("Step 4 data içinde file_path yok. Representative flow figure için FCS file_path gerekli.")
}

# ------------------------------------------------------------
# 2. Thresholds
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

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------

get_marker_map <- function(ff) {
  pp <- Biobase::pData(flowCore::parameters(ff))
  channels <- as.character(pp$name)
  markers <- as.character(pp$desc)
  idx <- is.na(markers) | markers == ""
  markers[idx] <- channels[idx]
  tibble(parameter_index = seq_along(channels), channel_name = channels, marker_desc = markers)
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
  ch$DUMP <- find_channel_regex_marker(map, "^Viability")
  if (is.na(ch$DUMP)) ch$DUMP <- find_channel_by_channel(map, "BV510-A")
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
  trans <- flowCore::logicleTransform("fixed_logicle", w = 0.5, t = 262144, m = 4.5, a = 0)
  tryCatch({
    tf <- flowCore::transformList(fluor_ch, trans)
    flowCore::transform(ff, tf)
  }, error = function(e) ff)
}

sample_n_safe <- function(df, n = 25000) {
  if (nrow(df) <= n) return(df)
  df[sample(seq_len(nrow(df)), n), , drop = FALSE]
}

make_density_panel <- function(df, x, y, title, subtitle = NULL, xline = NULL, yline = NULL) {
  # FlowJo-style white-background figure:
  # black point cloud + density contours + visible gate lines.
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 0.08, alpha = 0.08) +
    geom_density_2d(linewidth = 0.25) +
    labs(title = title, subtitle = subtitle, x = x, y = y) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8)
    )

  if (!is.null(xline)) p <- p + geom_vline(xintercept = xline, linetype = "dashed")
  if (!is.null(yline)) p <- p + geom_hline(yintercept = yline, linetype = "dashed")
  p
}

extract_plot_data <- function(fp) {
  ff_raw <- flowCore::read.FCS(fp, transformation = FALSE, truncate_max_range = FALSE)
  map <- get_marker_map(ff_raw)
  ch <- build_channel_map(map)

  ff <- ff_raw %>%
    apply_compensation_safely() %>%
    transform_safely()

  ex <- as.data.frame(flowCore::exprs(ff))

  d <- tibble(
    DUMP = as.numeric(ex[[ch$DUMP]]),
    CD19 = as.numeric(ex[[ch$CD19]]),
    CD27 = as.numeric(ex[[ch$CD27]]),
    IgD  = as.numeric(ex[[ch$IgD]]),
    IgM  = as.numeric(ex[[ch$IgM]]),
    IgA  = as.numeric(ex[[ch$IgA]]),
    IgG  = as.numeric(ex[[ch$IgG]]),
    CD38 = as.numeric(ex[[ch$CD38]]),
    CD138 = as.numeric(ex[[ch$CD138]]),
    CD24 = as.numeric(ex[[ch$CD24]]),
    CD10 = as.numeric(ex[[ch$CD10]]),
    CD39 = as.numeric(ex[[ch$CD39]])
  ) %>%
    mutate(
      dump_low = DUMP <= thresholds$DUMP_LOW,
      B_cell = dump_low & CD19 > thresholds$CD19,
      naive_like = B_cell & IgD > thresholds$IgD & CD27 <= thresholds$CD27,
      unswitched_memory_like = B_cell & IgD > thresholds$IgD & CD27 > thresholds$CD27,
      switched_memory_like = B_cell & IgD <= thresholds$IgD & CD27 > thresholds$CD27,
      plasmablast_like = B_cell & CD38 > thresholds$CD38_HIGH & CD27 > thresholds$CD27,
      transitional_like = B_cell & CD10 > thresholds$CD10 & CD24 > thresholds$CD24 & CD38 > thresholds$CD38,
      cd39_regulatory_like = B_cell & CD39 > thresholds$CD39
    )

  d
}

safe_pct <- function(x) round(100 * mean(x, na.rm = TRUE), 1)

# ------------------------------------------------------------
# 4. Select representative subjects
# ------------------------------------------------------------

score_col <- "CP22_integrated_humoral_B_cell_remodeling_score"

rep_candidates <- cp22_plot_source_df %>%
  filter(
    feature_ok == TRUE,
    !is.na(disease_group),
    !is.na(.data[[score_col]]),
    !is.na(file_path),
    n_cd19_b >= 1000
  ) %>%
  mutate(
    disease_group_chr = as.character(disease_group)
  ) %>%
  filter(disease_group_chr %in% c("Healthy control", "Cancer patient"))

group_medians <- rep_candidates %>%
  group_by(disease_group_chr) %>%
  summarise(group_median_score = median(.data[[score_col]], na.rm = TRUE), .groups = "drop")

representative_samples <- rep_candidates %>%
  left_join(group_medians, by = "disease_group_chr") %>%
  mutate(abs_diff = abs(.data[[score_col]] - group_median_score)) %>%
  group_by(disease_group_chr) %>%
  arrange(abs_diff, desc(n_cd19_b)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(integrated_score = .data[[score_col]]) %>%
  select(
    subject_id,
    disease_group = disease_group_chr,
    file_name,
    file_path,
    n_cd19_b,
    integrated_score,
    group_median_score,
    abs_diff
  )

readr::write_csv(
  representative_samples,
  file.path(out_dir, "SDY2583_CP22_representative_samples_STEP7B.csv")
)

# ------------------------------------------------------------
# 5. Build plots
# ------------------------------------------------------------

plots_by_sample <- list()
flow_stats <- list()

for (i in seq_len(nrow(representative_samples))) {

  row <- representative_samples[i, ]
  d <- extract_plot_data(row$file_path)

  label <- paste0(row$disease_group, " | ", row$subject_id)

  all_sample <- sample_n_safe(d, 30000)
  b_sample <- sample_n_safe(d %>% filter(B_cell), 25000)
  switched_sample <- sample_n_safe(d %>% filter(switched_memory_like), 15000)

  stats_i <- tibble(
    subject_id = row$subject_id,
    disease_group = row$disease_group,
    n_total = nrow(d),
    n_B_cell = sum(d$B_cell, na.rm = TRUE),
    pct_B_cell_total = safe_pct(d$B_cell),
    pct_naive_like_within_B = safe_pct(d$naive_like[d$B_cell]),
    pct_unswitched_memory_like_within_B = safe_pct(d$unswitched_memory_like[d$B_cell]),
    pct_switched_memory_like_within_B = safe_pct(d$switched_memory_like[d$B_cell]),
    pct_plasmablast_like_within_B = safe_pct(d$plasmablast_like[d$B_cell]),
    pct_transitional_like_within_B = safe_pct(d$transitional_like[d$B_cell]),
    pct_cd39_regulatory_like_within_B = safe_pct(d$cd39_regulatory_like[d$B_cell])
  )

  flow_stats[[i]] <- stats_i

  p1 <- make_density_panel(
    all_sample,
    "DUMP",
    "CD19",
    title = paste0(label, "\nDump-low / CD19+ B-cell gate"),
    subtitle = paste0("B cells: ", stats_i$pct_B_cell_total, "% of total events"),
    xline = thresholds$DUMP_LOW,
    yline = thresholds$CD19
  )

  p2 <- make_density_panel(
    b_sample,
    "IgD",
    "CD27",
    title = "B-cell memory-state gate",
    subtitle = paste0(
      "Switched-like: ", stats_i$pct_switched_memory_like_within_B,
      "% | Unswitched-like: ", stats_i$pct_unswitched_memory_like_within_B, "%"
    ),
    xline = thresholds$IgD,
    yline = thresholds$CD27
  )

  p3 <- make_density_panel(
    switched_sample,
    "IgA",
    "IgG",
    title = "IgA / IgG structure within switched-memory-like B cells",
    subtitle = "Shown only within IgD−CD27+ switched-memory-like B cells",
    xline = thresholds$IgA,
    yline = thresholds$IgG
  )

  p4 <- make_density_panel(
    b_sample,
    "CD27",
    "CD38",
    title = "Plasmablast-like activation axis",
    subtitle = paste0("Plasmablast-like: ", stats_i$pct_plasmablast_like_within_B, "% of B cells"),
    xline = thresholds$CD27,
    yline = thresholds$CD38_HIGH
  )

  p5 <- make_density_panel(
    b_sample,
    "CD10",
    "CD24",
    title = "Transitional / immature-like B-cell axis",
    subtitle = paste0("CD10+CD24+CD38+ transitional-like: ", stats_i$pct_transitional_like_within_B, "% of B cells"),
    xline = thresholds$CD10,
    yline = thresholds$CD24
  )

  p6 <- make_density_panel(
    b_sample,
    "CD24",
    "CD39",
    title = "CD39+ regulatory-like B-cell axis",
    subtitle = paste0("CD39+ regulatory-like: ", stats_i$pct_cd39_regulatory_like_within_B, "% of B cells"),
    xline = thresholds$CD24,
    yline = thresholds$CD39
  )

  plots_by_sample[[i]] <- (p1 / p2 / p3 / p4 / p5 / p6) +
    plot_annotation(
      title = paste0("CP22 FlowJo-style representative B-cell gating: ", label)
    )

  safe_name <- paste0(
    "FigureS_CP22_FlowJo_style_Bcell_gating_",
    gsub("[^A-Za-z0-9]+", "_", row$disease_group),
    "_",
    row$subject_id
  )

  ggsave(
    file.path(out_dir, paste0(safe_name, ".png")),
    plots_by_sample[[i]],
    width = 7.5,
    height = 16,
    dpi = 300
  )

  ggsave(
    file.path(out_dir, paste0(safe_name, ".pdf")),
    plots_by_sample[[i]],
    width = 7.5,
    height = 16
  )
}

flow_stats_df <- bind_rows(flow_stats)

readr::write_csv(
  flow_stats_df,
  file.path(out_dir, "SDY2583_CP22_representative_flow_stats_STEP7B.csv")
)

# Combined healthy/cancer figure
if (length(plots_by_sample) == 2) {
  combined_plot <- plots_by_sample[[1]] | plots_by_sample[[2]]

  ggsave(
    file.path(out_dir, "FigureS_CP22_FlowJo_style_Bcell_gating_healthy_vs_cancer_combined.png"),
    combined_plot,
    width = 15,
    height = 16,
    dpi = 300
  )

  ggsave(
    file.path(out_dir, "FigureS_CP22_FlowJo_style_Bcell_gating_healthy_vs_cancer_combined.pdf"),
    combined_plot,
    width = 15,
    height = 16
  )
}

save(
  representative_samples,
  flow_stats_df,
  plots_by_sample,
  file = file.path(rdata_dir, "SDY2583_CP22_STEP7B_representative_flow_figures.RData")
)

cat("\n============================================================\n")
cat("SDY2583 CP22 STEP 7B COMPLETE: FLOWJO-STYLE FIGURES\n")
cat("============================================================\n")

cat("\nRepresentative samples:\n")
print(as.data.frame(representative_samples), row.names = FALSE)

cat("\nRepresentative flow stats:\n")
print(as.data.frame(flow_stats_df), row.names = FALSE)

cat("\nFigures saved in:\n")
print(out_dir)

cat("\nFigure files:\n")
print(list.files(out_dir, pattern = "FigureS_CP22_.*\\.(png|pdf)$", full.names = FALSE))

cat("============================================================\n")
