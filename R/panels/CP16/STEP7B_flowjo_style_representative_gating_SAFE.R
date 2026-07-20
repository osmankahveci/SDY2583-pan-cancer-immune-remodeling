# ============================================================
# SDY2583 CP16
# STEP 7B SAFE: FlowJo-style representative gating figures
#
# Purpose:
#   Generate representative real-FCS gating figures for CP16.
#
# Important:
#   - No synthetic events are generated.
#   - All plots are drawn from real CP16 FCS events.
#   - Each gate is saved as a separate large figure.
#   - Multipage PDFs are also created with one gate per page.
#   - The script does NOT cram all gates into one page.
#
# Input:
#   outputs/CP16/11_RData/
#     SDY2583_CP16_STEP4_composite_scores.RData
#
# Output:
#   outputs/CP16/09_flowjo_style_representative_gating
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "ggplot2", "purrr", "tidyr")

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
  library(ggplot2)
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
# 2. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "09_flowjo_style_representative_gating")

individual_dir <- file.path(out_dir, "individual_gate_figures")
multipage_dir <- file.path(out_dir, "multipage_per_sample")
tables_dir <- file.path(out_dir, "source_tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(multipage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

step4_rdata <- file.path(rdata_dir, "SDY2583_CP16_STEP4_composite_scores.RData")

if (!file.exists(step4_rdata)) {
  stop("Step 4 RData bulunamadı: ", step4_rdata)
}

load(step4_rdata)

# Reset paths after load.
analysis_dir <- sd_analysis_dir("CP16")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "09_flowjo_style_representative_gating")
individual_dir <- file.path(out_dir, "individual_gate_figures")
multipage_dir <- file.path(out_dir, "multipage_per_sample")
tables_dir <- file.path(out_dir, "source_tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(multipage_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cp16_scores_data")) stop("cp16_scores_data bulunamadı.")
if (!exists("composite_score_cols")) stop("composite_score_cols bulunamadı.")

integrated_score <- "CP16_integrated_APC_DC_myeloid_remodeling_score"

if (!(integrated_score %in% names(cp16_scores_data))) {
  stop("Integrated CP16 score bulunamadı: ", integrated_score)
}

if (!("file_path" %in% names(cp16_scores_data))) {
  stop("cp16_scores_data içinde file_path yok. Step 2/3/4 objelerini kontrol et.")
}

# ------------------------------------------------------------
# 3. Thresholds used for main representative gate overlays
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

clean_stub <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("_$", "")
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
  ch$FceRI <- find_marker_channel_regex(map, "FceRI|FcER1|FcERI|FCER1")
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

safe_median <- function(x) {
  if (length(x) == 0) return(NA_real_)
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::median(x, na.rm = TRUE))
}

get_xlim_ylim <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) {
    return(list(xlim = c(-1, 5), ylim = c(-1, 5)))
  }

  xq <- stats::quantile(x[ok], probs = c(0.005, 0.995), na.rm = TRUE)
  yq <- stats::quantile(y[ok], probs = c(0.005, 0.995), na.rm = TRUE)

  xpad <- diff(xq) * 0.08
  ypad <- diff(yq) * 0.08

  if (!is.finite(xpad) || xpad == 0) xpad <- 0.5
  if (!is.finite(ypad) || ypad == 0) ypad <- 0.5

  list(
    xlim = c(xq[1] - xpad, xq[2] + xpad),
    ylim = c(yq[1] - ypad, yq[2] + ypad)
  )
}

downsample_plot_df <- function(df, max_points = 60000) {
  df <- df %>%
    filter(is.finite(x), is.finite(y))

  if (nrow(df) <= max_points) return(df)

  set.seed(2583)

  if ("gate_status" %in% names(df)) {
    gated <- df %>% filter(gate_status != "Parent/background")
    bg <- df %>% filter(gate_status == "Parent/background")

    n_gated_keep <- min(nrow(gated), round(max_points * 0.35))
    n_bg_keep <- max_points - n_gated_keep

    out <- bind_rows(
      if (nrow(gated) > 0) gated[sample(seq_len(nrow(gated)), size = n_gated_keep), , drop = FALSE] else gated,
      if (nrow(bg) > 0) bg[sample(seq_len(nrow(bg)), size = min(nrow(bg), n_bg_keep)), , drop = FALSE] else bg
    )

    return(out)
  }

  df[sample(seq_len(nrow(df)), size = max_points), , drop = FALSE]
}

base_flow_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35)
    )
}

make_flow_plot <- function(
    df,
    xlab,
    ylab,
    title,
    subtitle,
    gate_type = c("none", "quadrant", "x_low_y_high", "x_high_y_high", "x_high_y_low"),
    x_thr = NA_real_,
    y_thr = NA_real_,
    gate_label = NULL
) {

  gate_type <- match.arg(gate_type)

  lims <- get_xlim_ylim(df$x, df$y)
  plot_df <- downsample_plot_df(df, max_points = 60000)

  p <- ggplot(plot_df, aes(x = x, y = y)) +
    geom_point(aes(alpha = gate_status), size = 0.35) +
    scale_alpha_manual(
      values = c("Parent/background" = 0.18, "Gate-positive" = 0.65),
      drop = FALSE
    ) +
    guides(alpha = guide_legend(title = NULL, override.aes = list(size = 2))) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
    ) +
    base_flow_theme(12)

  # Density contours are helpful for FlowJo-like appearance.
  # If the subset is very small, silently skip contours.
  if (nrow(plot_df) >= 300) {
    p <- p +
      tryCatch(
        stat_density_2d(
          bins = 7,
          linewidth = 0.22,
          alpha = 0.45,
          show.legend = FALSE
        ),
        error = function(e) NULL
      )
  }

  if (gate_type %in% c("quadrant", "x_low_y_high", "x_high_y_high", "x_high_y_low")) {
    if (is.finite(x_thr)) p <- p + geom_vline(xintercept = x_thr, linetype = "dashed", linewidth = 0.45)
    if (is.finite(y_thr)) p <- p + geom_hline(yintercept = y_thr, linetype = "dashed", linewidth = 0.45)
  }

  if (!is.null(gate_label)) {
    p <- p +
      annotate(
        "label",
        x = lims$xlim[1] + 0.02 * diff(lims$xlim),
        y = lims$ylim[2] - 0.05 * diff(lims$ylim),
        hjust = 0,
        label = gate_label,
        size = 3.2,
        label.size = 0.2
      )
  }

  p
}

save_plot_both <- function(plot_obj, file_stub, width = 6.5, height = 5.5) {
  png_file <- file.path(individual_dir, paste0(file_stub, ".png"))
  pdf_file <- file.path(individual_dir, paste0(file_stub, ".pdf"))

  ggplot2::ggsave(png_file, plot_obj, width = width, height = height, dpi = 320, bg = "white")
  ggplot2::ggsave(pdf_file, plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")

  tibble(
    figure_stub = file_stub,
    png_file = png_file,
    pdf_file = pdf_file
  )
}

# ------------------------------------------------------------
# 5. Select representative samples
# ------------------------------------------------------------

representative_candidates <- cp16_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(.data[[integrated_score]]),
    !is.na(file_path),
    file.exists(file_path),
    disease_group %in% c("Healthy control", "Cancer patient")
  ) %>%
  mutate(
    disease_group_chr = as.character(disease_group)
  )

if (nrow(representative_candidates) == 0) {
  stop("Representative sample seçimi için uygun CP16 dosyası bulunamadı.")
}

group_medians <- representative_candidates %>%
  group_by(disease_group_chr) %>%
  summarise(
    median_integrated_score = median(.data[[integrated_score]], na.rm = TRUE),
    .groups = "drop"
  )

representative_samples <- representative_candidates %>%
  left_join(group_medians, by = "disease_group_chr") %>%
  mutate(
    abs_distance_to_group_median = abs(.data[[integrated_score]] - median_integrated_score),
    representative_role = case_when(
      disease_group_chr == "Healthy control" ~ "Healthy_median_integrated_score",
      disease_group_chr == "Cancer patient" ~ "Cancer_median_integrated_score",
      TRUE ~ disease_group_chr
    )
  ) %>%
  arrange(disease_group_chr, abs_distance_to_group_median, subject_id) %>%
  group_by(disease_group_chr) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(
    representative_role,
    subject_id,
    disease_group,
    age_for_model,
    sex,
    file_name,
    file_path,
    all_of(integrated_score),
    abs_distance_to_group_median
  )

write_csv(representative_samples, file.path(tables_dir, "CP16_representative_samples_STEP7B.csv"))

# ------------------------------------------------------------
# 6. Read, transform, and gate representative samples
# ------------------------------------------------------------

process_representative_sample <- function(row_df) {

  fp <- row_df$file_path[1]

  ff_raw <- flowCore::read.FCS(
    fp,
    transformation = FALSE,
    truncate_max_range = FALSE
  )

  marker_map_file <- get_marker_map(ff_raw)
  ch <- build_channel_map(marker_map_file)

  required <- c("DUMP", "CD45", "HLA_DR", "CD11c", "CD14", "CD16", "CD1c", "CD141", "CLEC9A", "CD123", "FceRI", "CD13")
  missing <- required[vapply(ch[required], function(x) is.na(x), logical(1))]

  if (length(missing) > 0) {
    stop("Missing required markers in ", basename(fp), ": ", paste(missing, collapse = ", "))
  }

  ff <- ff_raw %>%
    apply_compensation_safely() %>%
    transform_safely()

  ex <- flowCore::exprs(ff)

  dat <- tibble(
    DUMP = get_vec(ex, ch$DUMP),
    CD45 = get_vec(ex, ch$CD45),
    HLA_DR = get_vec(ex, ch$HLA_DR),
    CD11c = get_vec(ex, ch$CD11c),
    CD14 = get_vec(ex, ch$CD14),
    CD16 = get_vec(ex, ch$CD16),
    CD1c = get_vec(ex, ch$CD1c),
    CD141 = get_vec(ex, ch$CD141),
    CLEC9A = get_vec(ex, ch$CLEC9A),
    CD123 = get_vec(ex, ch$CD123),
    FceRI = get_vec(ex, ch$FceRI),
    CD13 = get_vec(ex, ch$CD13)
  )

  dat <- dat %>%
    mutate(
      dump_low = DUMP <= thresholds$DUMP_LOW,
      cd45_dump_low = dump_low & CD45 > thresholds$CD45,
      hladr_apc_core = cd45_dump_low & HLA_DR > thresholds$HLA_DR,
      cd11c_hladr_apc_like = hladr_apc_core & CD11c > thresholds$CD11c,

      cd14_mono_like = cd45_dump_low & CD14 > thresholds$CD14,
      cd14_hladr_pos = cd14_mono_like & HLA_DR > thresholds$HLA_DR,
      cd14_hladr_low = cd14_mono_like & HLA_DR <= thresholds$HLA_DR,

      cd1c_apc_like = hladr_apc_core & CD11c > thresholds$CD11c & CD1c > thresholds$CD1c,
      cd1c_cd14low_cdc2_like = cd1c_apc_like & CD14 <= thresholds$CD14,
      cd1c_fceri_apc_like = cd1c_apc_like & FceRI > thresholds$FceRI,

      cd123_hladr_dc_like = hladr_apc_core & CD123 > thresholds$CD123,
      cd123_pdc_like = hladr_apc_core &
        CD123 > thresholds$CD123 &
        CD11c <= thresholds$CD11c &
        CD14 <= thresholds$CD14 &
        CD16 <= thresholds$CD16,
      cd123_cd11c_mixed_apc_like = hladr_apc_core & CD123 > thresholds$CD123 & CD11c > thresholds$CD11c,

      cd141_pos_hladr = hladr_apc_core & CD141 > thresholds$CD141,
      cd141_clec9a_cdc1_like = hladr_apc_core &
        CD11c > thresholds$CD11c &
        CD141 > thresholds$CD141 &
        CLEC9A > thresholds$CLEC9A &
        CD14 <= thresholds$CD14,

      cd13_cd11c_hladr_apc_like = hladr_apc_core & CD13 > thresholds$CD13 & CD11c > thresholds$CD11c
    )

  n_total <- nrow(dat)
  n_cd45_dump_low <- safe_sum(dat$cd45_dump_low)
  n_hladr_core <- safe_sum(dat$hladr_apc_core)
  n_cd14_mono <- safe_sum(dat$cd14_mono_like)
  n_cd1c_apc <- safe_sum(dat$cd1c_apc_like)

  gate_percentages <- tibble(
    representative_role = row_df$representative_role[1],
    subject_id = row_df$subject_id[1],
    disease_group = as.character(row_df$disease_group[1]),
    file_name = row_df$file_name[1],
    total_events = n_total,
    n_cd45_dump_low = n_cd45_dump_low,
    n_hladr_apc_core = n_hladr_core,
    pct_cd45_dump_low_within_total = safe_pct(n_cd45_dump_low, n_total),
    pct_hladr_apc_core_within_cd45_dump_low = safe_pct(n_hladr_core, n_cd45_dump_low),
    pct_cd11c_hladr_apc_like_within_cd45_dump_low = safe_pct(safe_sum(dat$cd11c_hladr_apc_like), n_cd45_dump_low),
    pct_cd11c_hladr_apc_like_within_hladr_core = safe_pct(safe_sum(dat$cd11c_hladr_apc_like), n_hladr_core),
    pct_cd14_mono_like_within_cd45_dump_low = safe_pct(n_cd14_mono, n_cd45_dump_low),
    pct_cd14_hladr_low_within_cd45_dump_low = safe_pct(safe_sum(dat$cd14_hladr_low), n_cd45_dump_low),
    pct_hladr_low_within_cd14_mono_like = safe_pct(safe_sum(dat$cd14_hladr_low), n_cd14_mono),
    pct_cd1c_apc_like_within_cd45_dump_low = safe_pct(n_cd1c_apc, n_cd45_dump_low),
    pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low = safe_pct(safe_sum(dat$cd1c_cd14low_cdc2_like), n_cd45_dump_low),
    pct_cd1c_fceri_apc_like_within_cd45_dump_low = safe_pct(safe_sum(dat$cd1c_fceri_apc_like), n_cd45_dump_low),
    pct_cd123_hladr_dc_like_within_cd45_dump_low = safe_pct(safe_sum(dat$cd123_hladr_dc_like), n_cd45_dump_low),
    pct_cd123_pdc_like_within_cd45_dump_low = safe_pct(safe_sum(dat$cd123_pdc_like), n_cd45_dump_low),
    pct_cd141_clec9a_cdc1_like_within_cd45_dump_low = safe_pct(safe_sum(dat$cd141_clec9a_cdc1_like), n_cd45_dump_low),
    pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low = safe_pct(safe_sum(dat$cd13_cd11c_hladr_apc_like), n_cd45_dump_low)
  )

  list(
    sample_row = row_df,
    events = dat,
    gate_percentages = gate_percentages,
    marker_map = marker_map_file
  )
}

processed_samples <- lapply(seq_len(nrow(representative_samples)), function(i) {
  process_representative_sample(representative_samples[i, ])
})

representative_gate_percentages <- bind_rows(lapply(processed_samples, function(x) x$gate_percentages))
write_csv(representative_gate_percentages, file.path(tables_dir, "CP16_representative_gate_percentages_STEP7B.csv"))

representative_marker_maps <- bind_rows(lapply(processed_samples, function(x) {
  x$marker_map %>%
    mutate(
      subject_id = x$sample_row$subject_id[1],
      representative_role = x$sample_row$representative_role[1],
      file_name = x$sample_row$file_name[1]
    )
}))
write_csv(representative_marker_maps, file.path(tables_dir, "CP16_representative_marker_maps_STEP7B.csv"))

# ------------------------------------------------------------
# 7. Build plot objects for one sample
# ------------------------------------------------------------

build_plots_for_sample <- function(proc) {

  dat <- proc$events
  row_df <- proc$sample_row
  perc <- proc$gate_percentages

  role <- row_df$representative_role[1]
  sid <- row_df$subject_id[1]
  group <- as.character(row_df$disease_group[1])
  age <- row_df$age_for_model[1]
  sex <- as.character(row_df$sex[1])

  title_prefix <- paste0(role, " / ", sid)
  subtitle_prefix <- paste0(group, "; age ", age, "; ", sex, "; real FCS events")

  plots <- list()

  # 01. Dump-low CD45+ primary compartment
  df1 <- dat %>%
    transmute(
      x = DUMP,
      y = CD45,
      gate_status = ifelse(cd45_dump_low, "Gate-positive", "Parent/background")
    )

  plots[["01_Dump_vs_CD45_primary"]] <- make_flow_plot(
    df1,
    xlab = "Dump / Viability channel",
    ylab = "CD45",
    title = paste0(title_prefix, " — CD45+ dump-low primary gate"),
    subtitle = subtitle_prefix,
    gate_type = "x_low_y_high",
    x_thr = thresholds$DUMP_LOW,
    y_thr = thresholds$CD45,
    gate_label = paste0(
      "CD45+ dump-low: ",
      signif(perc$pct_cd45_dump_low_within_total[1], 3),
      "% of total"
    )
  )

  # 02. HLA-DR+ APC core within CD45 dump-low
  parent2 <- dat %>% filter(cd45_dump_low)
  df2 <- parent2 %>%
    transmute(
      x = CD11c,
      y = HLA_DR,
      gate_status = ifelse(hladr_apc_core, "Gate-positive", "Parent/background")
    )

  plots[["02_CD11c_vs_HLA_DR_APC_core"]] <- make_flow_plot(
    df2,
    xlab = "CD11c",
    ylab = "HLA-DR",
    title = paste0(title_prefix, " — HLA-DR+ APC-like core"),
    subtitle = paste0(subtitle_prefix, "; parent: CD45+ dump-low"),
    gate_type = "x_high_y_high",
    x_thr = thresholds$CD11c,
    y_thr = thresholds$HLA_DR,
    gate_label = paste0(
      "HLA-DR+ core: ",
      signif(perc$pct_hladr_apc_core_within_cd45_dump_low[1], 3),
      "% of CD45+ dump-low"
    )
  )

  # 03. CD14 / HLA-DR monocyte-like HLA-DR-low phenotype
  df3 <- parent2 %>%
    transmute(
      x = CD14,
      y = HLA_DR,
      gate_status = ifelse(cd14_hladr_low, "Gate-positive", "Parent/background")
    )

  plots[["03_CD14_vs_HLA_DR_monocyte_HLADR_low"]] <- make_flow_plot(
    df3,
    xlab = "CD14",
    ylab = "HLA-DR",
    title = paste0(title_prefix, " — CD14+HLA-DR-low monocyte-like phenotype"),
    subtitle = paste0(subtitle_prefix, "; parent: CD45+ dump-low"),
    gate_type = "x_high_y_low",
    x_thr = thresholds$CD14,
    y_thr = thresholds$HLA_DR,
    gate_label = paste0(
      "CD14+HLA-DR-low: ",
      signif(perc$pct_cd14_hladr_low_within_cd45_dump_low[1], 3),
      "% of CD45+ dump-low"
    )
  )

  # 04. CD1c / CD14 cDC2-like phenotype within HLA-DR core
  parent4 <- dat %>% filter(hladr_apc_core)
  df4 <- parent4 %>%
    transmute(
      x = CD14,
      y = CD1c,
      gate_status = ifelse(cd1c_cd14low_cdc2_like, "Gate-positive", "Parent/background")
    )

  plots[["04_CD14_vs_CD1c_cDC2_like"]] <- make_flow_plot(
    df4,
    xlab = "CD14",
    ylab = "CD1c",
    title = paste0(title_prefix, " — CD1c+CD14-low cDC2-like phenotype"),
    subtitle = paste0(subtitle_prefix, "; parent: HLA-DR+ APC-like core"),
    gate_type = "x_low_y_high",
    x_thr = thresholds$CD14,
    y_thr = thresholds$CD1c,
    gate_label = paste0(
      "CD1c+CD14-low: ",
      signif(perc$pct_cd1c_cd14low_cdc2_like_within_cd45_dump_low[1], 3),
      "% of CD45+ dump-low"
    )
  )

  # 05. CD1c / FcERI APC-like phenotype
  df5 <- parent4 %>%
    transmute(
      x = CD1c,
      y = FceRI,
      gate_status = ifelse(cd1c_fceri_apc_like, "Gate-positive", "Parent/background")
    )

  plots[["05_CD1c_vs_FcERI_APC_like"]] <- make_flow_plot(
    df5,
    xlab = "CD1c",
    ylab = "FcERI",
    title = paste0(title_prefix, " — CD1c+FcERI+ APC-like phenotype"),
    subtitle = paste0(subtitle_prefix, "; parent: HLA-DR+ APC-like core"),
    gate_type = "x_high_y_high",
    x_thr = thresholds$CD1c,
    y_thr = thresholds$FceRI,
    gate_label = paste0(
      "CD1c+FcERI+: ",
      signif(perc$pct_cd1c_fceri_apc_like_within_cd45_dump_low[1], 3),
      "% of CD45+ dump-low"
    )
  )

  # 06. CD123 / CD11c pDC-like phenotype
  df6 <- parent4 %>%
    transmute(
      x = CD11c,
      y = CD123,
      gate_status = ifelse(cd123_pdc_like, "Gate-positive", "Parent/background")
    )

  plots[["06_CD11c_vs_CD123_pDC_like"]] <- make_flow_plot(
    df6,
    xlab = "CD11c",
    ylab = "CD123",
    title = paste0(title_prefix, " — CD123+CD11c-low pDC-like phenotype"),
    subtitle = paste0(subtitle_prefix, "; parent: HLA-DR+ APC-like core"),
    gate_type = "x_low_y_high",
    x_thr = thresholds$CD11c,
    y_thr = thresholds$CD123,
    gate_label = paste0(
      "CD123+ pDC-like: ",
      signif(perc$pct_cd123_pdc_like_within_cd45_dump_low[1], 3),
      "% of CD45+ dump-low"
    )
  )

  # 07. CD141 / CLEC9A cDC1-like phenotype
  df7 <- parent4 %>%
    transmute(
      x = CD141,
      y = CLEC9A,
      gate_status = ifelse(cd141_clec9a_cdc1_like, "Gate-positive", "Parent/background")
    )

  plots[["07_CD141_vs_CLEC9A_cDC1_like"]] <- make_flow_plot(
    df7,
    xlab = "CD141",
    ylab = "CLEC9A",
    title = paste0(title_prefix, " — CD141+CLEC9A+ cDC1-like phenotype"),
    subtitle = paste0(subtitle_prefix, "; parent: HLA-DR+ APC-like core"),
    gate_type = "x_high_y_high",
    x_thr = thresholds$CD141,
    y_thr = thresholds$CLEC9A,
    gate_label = paste0(
      "CD141+CLEC9A+: ",
      signif(perc$pct_cd141_clec9a_cdc1_like_within_cd45_dump_low[1], 3),
      "% of CD45+ dump-low"
    )
  )

  # 08. CD13 / CD11c myeloid APC-like phenotype
  df8 <- parent4 %>%
    transmute(
      x = CD13,
      y = CD11c,
      gate_status = ifelse(cd13_cd11c_hladr_apc_like, "Gate-positive", "Parent/background")
    )

  plots[["08_CD13_vs_CD11c_APC_myeloid_like"]] <- make_flow_plot(
    df8,
    xlab = "CD13",
    ylab = "CD11c",
    title = paste0(title_prefix, " — CD13+CD11c+HLA-DR+ APC-like myeloid phenotype"),
    subtitle = paste0(subtitle_prefix, "; parent: HLA-DR+ APC-like core"),
    gate_type = "x_high_y_high",
    x_thr = thresholds$CD13,
    y_thr = thresholds$CD11c,
    gate_label = paste0(
      "CD13+CD11c+HLA-DR+: ",
      signif(perc$pct_cd13_cd11c_hladr_apc_like_within_cd45_dump_low[1], 3),
      "% of CD45+ dump-low"
    )
  )

  plots
}

# ------------------------------------------------------------
# 8. Save individual figures and multipage PDFs
# ------------------------------------------------------------

figure_manifest_list <- list()
multipage_manifest_list <- list()

for (idx in seq_along(processed_samples)) {

  proc <- processed_samples[[idx]]
  row_df <- proc$sample_row

  sample_stub <- clean_stub(paste0(row_df$representative_role[1], "_", row_df$subject_id[1]))
  plots <- build_plots_for_sample(proc)

  # Individual files
  for (gate_name in names(plots)) {
    file_stub <- paste0(sample_stub, "__", gate_name)
    figure_manifest_list[[paste0(sample_stub, "_", gate_name)]] <- save_plot_both(
      plots[[gate_name]],
      file_stub,
      width = 6.8,
      height = 5.6
    ) %>%
      mutate(
        representative_role = row_df$representative_role[1],
        subject_id = row_df$subject_id[1],
        disease_group = as.character(row_df$disease_group[1]),
        gate_stub = gate_name
      )
  }

  # Multipage PDF, one gate per page.
  multipage_pdf <- file.path(multipage_dir, paste0(sample_stub, "__multipage_one_gate_per_page.pdf"))
  grDevices::cairo_pdf(multipage_pdf, width = 6.8, height = 5.6)
  for (gate_name in names(plots)) {
    print(plots[[gate_name]])
  }
  dev.off()

  multipage_manifest_list[[sample_stub]] <- tibble(
    representative_role = row_df$representative_role[1],
    subject_id = row_df$subject_id[1],
    disease_group = as.character(row_df$disease_group[1]),
    multipage_pdf = multipage_pdf
  )
}

figure_manifest <- bind_rows(figure_manifest_list)
multipage_manifest <- bind_rows(multipage_manifest_list)

write_csv(figure_manifest, file.path(out_dir, "SDY2583_CP16_STEP7B_flowjo_style_figure_manifest.csv"))
write_csv(multipage_manifest, file.path(out_dir, "SDY2583_CP16_STEP7B_flowjo_style_multipage_manifest.csv"))

# ------------------------------------------------------------
# 9. Save RData
# ------------------------------------------------------------

save(
  representative_samples,
  representative_gate_percentages,
  representative_marker_maps,
  figure_manifest,
  multipage_manifest,
  thresholds,
  file = file.path(rdata_dir, "SDY2583_CP16_STEP7B_flowjo_style_representative_gating.RData")
)

# ------------------------------------------------------------
# 10. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP16 STEP 7B COMPLETE: FLOWJO-STYLE REPRESENTATIVE GATING\n")
cat("============================================================\n")

cat("\nRepresentative samples:\n")
print(as.data.frame(representative_samples), row.names = FALSE)

cat("\nRepresentative gate percentages:\n")
print(as.data.frame(representative_gate_percentages), row.names = FALSE)

cat("\nFigure manifest:\n")
print(as.data.frame(figure_manifest), row.names = FALSE)

cat("\nMultipage PDFs:\n")
print(as.data.frame(multipage_manifest), row.names = FALSE)

cat("\nIndividual figures saved in:\n")
print(individual_dir)

cat("\nMultipage PDFs saved in:\n")
print(multipage_dir)

cat("\nSource tables saved in:\n")
print(tables_dir)

cat("\nNo synthetic events were generated. All panels were plotted from real CP16 FCS events.\n")

cat("============================================================\n")
