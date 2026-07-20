# ============================================================
# SDY2583 CP10
# STEP 7B SAFE: FlowJo-style representative gating figures
#
# IMPORTANT:
# - This script does NOT draw synthetic / schematic flow plots.
# - It reads the actual CP10 FCS files and plots real events.
# - Marker-name tolerant extraction is used because DBG161/DBG177 have
#   shifted channel positions.
# - Each gate is saved as a separate high-resolution figure.
# - No single-page squeezed multi-gate layout is produced.
#
# Outputs:
#   outputs/CP10/09_flowjo_style_representative_gating
#
# Representative samples:
#   1) Healthy control nearest to the healthy median CP10 integrated score
#   2) Cancer patient nearest to the cancer median CP10 integrated score
#
# Phenotype terminology:
#   CD13+CD66b+ granulocyte-like, not definitive neutrophils
#   CCR3+CD66b+ eosinophil-like granulocytic, not definitive eosinophils
#   CD123+HLA-DR+ DC-like, not definitive dendritic cells
#   CD14+HLA-DR-low monocyte-like, not definitive MDSCs
# ============================================================

rm(list = ls())
source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

cran_pkgs <- c("dplyr", "readr", "stringr", "tibble", "ggplot2", "scales")

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
  library(scales)
})

# Namespace safety
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
bind_rows <- dplyr::bind_rows
case_when <- dplyr::case_when
n <- dplyr::n

set.seed(2583)

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

analysis_dir <- sd_analysis_dir("CP10")
rdata_dir <- file.path(analysis_dir, "11_RData")
out_dir <- file.path(analysis_dir, "09_flowjo_style_representative_gating")
individual_dir <- file.path(out_dir, "individual_gate_figures")
multipage_dir <- file.path(out_dir, "multipage_per_sample")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(multipage_dir, recursive = TRUE, showWarnings = FALSE)

step1_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP1_fcs_inventory_marker_QC.RData")
step4_rdata <- file.path(rdata_dir, "SDY2583_CP10_STEP4_composite_scores.RData")

if (!file.exists(step1_rdata)) stop("Step 1 RData bulunamadı: ", step1_rdata)
if (!file.exists(step4_rdata)) stop("Step 4 RData bulunamadı: ", step4_rdata)

load(step1_rdata)
load(step4_rdata)

if (!exists("cp10_scores_data")) stop("cp10_scores_data bulunamadı.")
if (!exists("fcs_files")) stop("fcs_files bulunamadı.")

# ------------------------------------------------------------
# 3. Main thresholds, same as CP10 Step 2
# ------------------------------------------------------------

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

threshold_table <- tibble(
  marker_or_gate = names(thresholds),
  threshold = unlist(thresholds)
)

write_csv(
  threshold_table,
  file.path(out_dir, "SDY2583_CP10_STEP7B_gate_thresholds_used.csv")
)

# ------------------------------------------------------------
# 4. Helper functions: marker-name tolerant FCS reading
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
  ch$VIABILITY <- find_marker_channel_regex(map, "^Viability$")

  if (is.na(ch$VIABILITY)) {
    ch$VIABILITY <- find_marker_channel_regex(map, "Viability|Live|Dead")
  }

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

safe_pct <- function(mask, denom) {
  den <- sum(denom, na.rm = TRUE)
  if (is.na(den) || den <= 0) return(NA_real_)
  100 * sum(mask & denom, na.rm = TRUE) / den
}

fmt_pct <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.1f%%", x))
}

# ------------------------------------------------------------
# 5. Select representative samples
# ------------------------------------------------------------

score_for_rep <- "CP10_integrated_myeloid_granulocytic_remodeling_score"

if (!(score_for_rep %in% names(cp10_scores_data))) {
  stop("Representative score not found: ", score_for_rep)
}

candidate_reps <- cp10_scores_data %>%
  filter(
    model_ready_age_sex == TRUE,
    !is.na(.data[[score_for_rep]]),
    !is.na(file_path),
    file.exists(file_path),
    disease_group %in% c("Healthy control", "Cancer patient")
  ) %>%
  mutate(
    disease_group_chr = as.character(disease_group),
    rep_score = as.numeric(.data[[score_for_rep]])
  )

if (nrow(candidate_reps) == 0) {
  stop("Temsilci örnek seçilecek uygun kayıt bulunamadı.")
}

healthy_median <- candidate_reps %>%
  filter(disease_group_chr == "Healthy control") %>%
  summarise(x = median(rep_score, na.rm = TRUE)) %>%
  pull(x)

cancer_median <- candidate_reps %>%
  filter(disease_group_chr == "Cancer patient") %>%
  summarise(x = median(rep_score, na.rm = TRUE)) %>%
  pull(x)

healthy_rep <- candidate_reps %>%
  filter(disease_group_chr == "Healthy control") %>%
  mutate(abs_dist = abs(rep_score - healthy_median)) %>%
  arrange(abs_dist, subject_id) %>%
  slice(1) %>%
  mutate(representative_role = "Healthy_median_integrated_score")

cancer_rep <- candidate_reps %>%
  filter(disease_group_chr == "Cancer patient") %>%
  mutate(abs_dist = abs(rep_score - cancer_median)) %>%
  arrange(abs_dist, subject_id) %>%
  slice(1) %>%
  mutate(representative_role = "Cancer_median_integrated_score")

representative_samples <- bind_rows(healthy_rep, cancer_rep) %>%
  select(
    representative_role,
    subject_id,
    disease_group,
    age_for_model,
    sex,
    file_name,
    file_path,
    rep_score
  )

write_csv(
  representative_samples,
  file.path(out_dir, "SDY2583_CP10_STEP7B_representative_samples.csv")
)

# ------------------------------------------------------------
# 6. Read and prepare one representative sample
# ------------------------------------------------------------

read_prepare_sample <- function(fp) {

  ff_raw <- flowCore::read.FCS(
    fp,
    transformation = FALSE,
    truncate_max_range = FALSE
  )

  marker_map <- get_marker_map(ff_raw)
  ch <- build_channel_map(marker_map)

  required <- c(
    "VIABILITY", "CD45", "CD3", "CD19", "CD56",
    "CD14", "HLA_DR", "CD11c", "CD13", "CD66b", "CCR3", "CD123"
  )

  missing <- required[vapply(ch[required], function(x) is.na(x), logical(1))]
  if (length(missing) > 0) {
    stop("Missing channels in ", basename(fp), ": ", paste(missing, collapse = ", "))
  }

  ff <- ff_raw %>%
    apply_compensation_safely() %>%
    transform_safely()

  ex <- flowCore::exprs(ff)

  df <- tibble(
    EventID = seq_len(nrow(ex)),
    Viability = get_vec(ex, ch$VIABILITY),
    CD45 = get_vec(ex, ch$CD45),
    CD3 = get_vec(ex, ch$CD3),
    CD19 = get_vec(ex, ch$CD19),
    CD56 = get_vec(ex, ch$CD56),
    CD14 = get_vec(ex, ch$CD14),
    HLA_DR = get_vec(ex, ch$HLA_DR),
    CD11c = get_vec(ex, ch$CD11c),
    CD13 = get_vec(ex, ch$CD13),
    CD66b = get_vec(ex, ch$CD66b),
    CCR3 = get_vec(ex, ch$CCR3),
    CD123 = get_vec(ex, ch$CD123)
  ) %>%
    mutate(
      viable = Viability <= thresholds$VIABILITY_LOW,
      cd45_pos = viable & CD45 > thresholds$CD45,
      t_like = cd45_pos & CD3 > thresholds$CD3,
      b_like = cd45_pos & CD19 > thresholds$CD19,
      nk_like = cd45_pos & CD3 <= thresholds$CD3 & CD56 > thresholds$CD56,
      lymphoid_like_any = cd45_pos & (CD3 > thresholds$CD3 | CD19 > thresholds$CD19 | CD56 > thresholds$CD56),
      cd14_mono_like = cd45_pos & CD14 > thresholds$CD14,
      cd14_hladr_low = cd14_mono_like & HLA_DR <= thresholds$HLA_DR,
      cd14_hladr_pos = cd14_mono_like & HLA_DR > thresholds$HLA_DR,
      cd14_cd11c_hladr_pos = cd14_mono_like & CD11c > thresholds$CD11c & HLA_DR > thresholds$HLA_DR,
      cd11c_hladr_apc_like = cd45_pos & CD11c > thresholds$CD11c & HLA_DR > thresholds$HLA_DR,
      cd66b_gran_like = cd45_pos & CD66b > thresholds$CD66b,
      cd13_cd66b_gran_like = cd45_pos & CD13 > thresholds$CD13 & CD66b > thresholds$CD66b,
      ccr3_pos = cd45_pos & CCR3 > thresholds$CCR3,
      ccr3_cd66b_eosinophil_like = cd45_pos & CCR3 > thresholds$CCR3 & CD66b > thresholds$CD66b,
      ccr3_cd13_cd66b_gran_like = cd45_pos & CCR3 > thresholds$CCR3 & CD13 > thresholds$CD13 & CD66b > thresholds$CD66b,
      cd123_pos = cd45_pos & CD123 > thresholds$CD123,
      cd123_hladr_dc_like = cd45_pos & CD123 > thresholds$CD123 & HLA_DR > thresholds$HLA_DR,
      cd123_hladr_cd11c_neg_pdc_like =
        cd45_pos &
        CD123 > thresholds$CD123 &
        HLA_DR > thresholds$HLA_DR &
        CD11c <= thresholds$CD11c &
        CD14 <= thresholds$CD14 &
        CD3 <= thresholds$CD3 &
        CD19 <= thresholds$CD19
    )

  list(
    data = df,
    marker_map = marker_map,
    channels = ch
  )
}

# ------------------------------------------------------------
# 7. FlowJo-style plotting functions
# ------------------------------------------------------------

thin_for_plot <- function(df, max_events = 120000) {
  # Do not use .data[[1]] here; .data pronoun requires column names.
  # Filtering for finite x/y already happens before this function is called.
  if (nrow(df) <= max_events) return(df)
  df[sample(seq_len(nrow(df)), max_events), , drop = FALSE]
}

get_limits <- function(v, threshold_value = NULL) {
  v <- v[is.finite(v)]
  if (length(v) < 10) return(c(0, 5))

  q <- as.numeric(stats::quantile(v, probs = c(0.001, 0.999), na.rm = TRUE))

  if (!is.null(threshold_value) && is.finite(threshold_value)) {
    q[1] <- min(q[1], threshold_value - 0.35)
    q[2] <- max(q[2], threshold_value + 0.35)
  }

  if (!is.finite(q[1]) || !is.finite(q[2]) || q[1] == q[2]) {
    q <- range(v, na.rm = TRUE)
  }

  pad <- 0.03 * diff(q)
  c(q[1] - pad, q[2] + pad)
}

theme_flowjo <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1),
      legend.position = "none"
    )
}

make_flow_plot <- function(
    df,
    xvar,
    yvar,
    xlab,
    ylab,
    xthr = NULL,
    ythr = NULL,
    title,
    subtitle,
    label_text,
    max_events = 120000
) {

  plot_df <- df %>%
    filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))

  n_total_for_panel <- nrow(plot_df)

  if (n_total_for_panel == 0) {
    stop("No finite events for plot: ", title)
  }

  plot_df <- thin_for_plot(plot_df, max_events = max_events)

  xlim <- get_limits(plot_df[[xvar]], xthr)
  ylim <- get_limits(plot_df[[yvar]], ythr)

  label_x <- xlim[1] + 0.04 * diff(xlim)
  label_y <- ylim[2] - 0.05 * diff(ylim)

  p <- ggplot(plot_df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_bin2d(bins = 190) +
    scale_fill_gradientn(
      colors = c("white", "grey85", "green3", "yellow", "orange", "red", "darkred"),
      trans = "sqrt"
    ) +
    stat_density_2d(color = "black", linewidth = 0.22, alpha = 0.45, bins = 7) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(
      title = title,
      subtitle = paste0(subtitle, " | real FCS events shown: ", format(nrow(plot_df), big.mark = ","), " of ", format(n_total_for_panel, big.mark = ",")),
      x = xlab,
      y = ylab
    ) +
    annotate(
      "label",
      x = label_x,
      y = label_y,
      hjust = 0,
      vjust = 1,
      label = label_text,
      size = 3.5,
      label.size = 0.25,
      fill = "white",
      alpha = 0.92
    ) +
    theme_flowjo(base_size = 12)

  if (!is.null(xthr) && is.finite(xthr)) {
    p <- p + geom_vline(xintercept = xthr, linewidth = 0.45, linetype = "dashed")
  }

  if (!is.null(ythr) && is.finite(ythr)) {
    p <- p + geom_hline(yintercept = ythr, linewidth = 0.45, linetype = "dashed")
  }

  p
}

save_gate_plot <- function(plot_obj, sample_stub, gate_stub, width = 7.5, height = 7.5) {

  pdf_file <- file.path(individual_dir, paste0(sample_stub, "__", gate_stub, ".pdf"))
  png_file <- file.path(individual_dir, paste0(sample_stub, "__", gate_stub, ".png"))

  ggsave(pdf_file, plot_obj, width = width, height = height)
  ggsave(png_file, plot_obj, width = width, height = height, dpi = 300)

  tibble(
    sample_stub = sample_stub,
    gate_stub = gate_stub,
    pdf_file = pdf_file,
    png_file = png_file
  )
}

clean_stub <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# ------------------------------------------------------------
# 8. Build plots for one sample
# ------------------------------------------------------------

build_sample_plots <- function(sample_row) {

  sample_row <- as_tibble(sample_row)
  sid <- sample_row$subject_id[1]
  role <- sample_row$representative_role[1]
  group <- as.character(sample_row$disease_group[1])
  fp <- sample_row$file_path[1]

  sample_stub <- clean_stub(paste(role, sid, sep = "__"))

  prep <- read_prepare_sample(fp)
  df <- prep$data

  denom_all <- rep(TRUE, nrow(df))
  denom_cd45 <- df$cd45_pos

  n_all <- sum(denom_all, na.rm = TRUE)
  n_cd45 <- sum(denom_cd45, na.rm = TRUE)

  subtitle_base <- paste0(
    sid, " | ", group,
    " | integrated score=", sprintf("%.3f", sample_row$rep_score[1])
  )

  gate_summary <- tibble(
    representative_role = role,
    subject_id = sid,
    disease_group = group,
    file_name = basename(fp),
    total_events = n_all,
    viable_pct_total = safe_pct(df$viable, denom_all),
    viable_cd45_pct_total = safe_pct(df$cd45_pos, denom_all),
    n_viable_cd45 = n_cd45,
    pct_t_like_within_cd45 = safe_pct(df$t_like, denom_cd45),
    pct_b_like_within_cd45 = safe_pct(df$b_like, denom_cd45),
    pct_nk_like_within_cd45 = safe_pct(df$nk_like, denom_cd45),
    pct_lymphoid_like_any_within_cd45 = safe_pct(df$lymphoid_like_any, denom_cd45),
    pct_cd14_mono_like_within_cd45 = safe_pct(df$cd14_mono_like, denom_cd45),
    pct_cd14_hladr_low_within_cd45 = safe_pct(df$cd14_hladr_low, denom_cd45),
    pct_cd14_hladr_pos_within_cd45 = safe_pct(df$cd14_hladr_pos, denom_cd45),
    pct_cd11c_hladr_apc_like_within_cd45 = safe_pct(df$cd11c_hladr_apc_like, denom_cd45),
    pct_cd13_cd66b_gran_like_within_cd45 = safe_pct(df$cd13_cd66b_gran_like, denom_cd45),
    pct_cd66b_gran_like_within_cd45 = safe_pct(df$cd66b_gran_like, denom_cd45),
    pct_ccr3_cd66b_eosinophil_like_within_cd45 = safe_pct(df$ccr3_cd66b_eosinophil_like, denom_cd45),
    pct_cd123_hladr_dc_like_within_cd45 = safe_pct(df$cd123_hladr_dc_like, denom_cd45),
    pct_cd123_hladr_cd11c_neg_pdc_like_within_cd45 = safe_pct(df$cd123_hladr_cd11c_neg_pdc_like, denom_cd45)
  )

  plot_list <- list()

  # 01 Viability vs CD45
  label_01 <- paste0(
    "Gate: Viability-low CD45+\n",
    "Viability-low: ", fmt_pct(safe_pct(df$viable, denom_all)), "\n",
    "Viability-low CD45+: ", fmt_pct(safe_pct(df$cd45_pos, denom_all))
  )

  plot_list[["01_Viability_vs_CD45"]] <- make_flow_plot(
    df = df,
    xvar = "Viability",
    yvar = "CD45",
    xlab = "Viability",
    ylab = "CD45",
    xthr = thresholds$VIABILITY_LOW,
    ythr = thresholds$CD45,
    title = "Primary leukocyte gate",
    subtitle = subtitle_base,
    label_text = label_01,
    max_events = 140000
  )

  # 02 CD3 vs CD19
  label_02 <- paste0(
    "Denominator: Viability-low CD45+\n",
    "CD3+ T-cell-like: ", fmt_pct(safe_pct(df$t_like, denom_cd45)), "\n",
    "CD19+ B-cell-like: ", fmt_pct(safe_pct(df$b_like, denom_cd45))
  )

  plot_list[["02_CD3_vs_CD19_within_CD45"]] <- make_flow_plot(
    df = df %>% filter(cd45_pos),
    xvar = "CD3",
    yvar = "CD19",
    xlab = "CD3",
    ylab = "CD19",
    xthr = thresholds$CD3,
    ythr = thresholds$CD19,
    title = "T-cell-like and B-cell-like compartments",
    subtitle = subtitle_base,
    label_text = label_02,
    max_events = 120000
  )

  # 03 CD3 vs CD56
  label_03 <- paste0(
    "Denominator: Viability-low CD45+\n",
    "CD3- CD56+ NK-like: ", fmt_pct(safe_pct(df$nk_like, denom_cd45)), "\n",
    "CD3+ T-cell-like: ", fmt_pct(safe_pct(df$t_like, denom_cd45))
  )

  plot_list[["03_CD3_vs_CD56_within_CD45"]] <- make_flow_plot(
    df = df %>% filter(cd45_pos),
    xvar = "CD3",
    yvar = "CD56",
    xlab = "CD3",
    ylab = "CD56",
    xthr = thresholds$CD3,
    ythr = thresholds$CD56,
    title = "NK-like and T-cell-like compartments",
    subtitle = subtitle_base,
    label_text = label_03,
    max_events = 120000
  )

  # 04 CD13 vs CD66b
  label_04 <- paste0(
    "Denominator: Viability-low CD45+\n",
    "CD66b+ granulocyte-like: ", fmt_pct(safe_pct(df$cd66b_gran_like, denom_cd45)), "\n",
    "CD13+CD66b+ granulocyte-like: ", fmt_pct(safe_pct(df$cd13_cd66b_gran_like, denom_cd45))
  )

  plot_list[["04_CD13_vs_CD66b_granulocyte_like"]] <- make_flow_plot(
    df = df %>% filter(cd45_pos),
    xvar = "CD13",
    yvar = "CD66b",
    xlab = "CD13",
    ylab = "CD66b",
    xthr = thresholds$CD13,
    ythr = thresholds$CD66b,
    title = "Granulocyte-like axis",
    subtitle = subtitle_base,
    label_text = label_04,
    max_events = 120000
  )

  # 05 CD14 vs HLA-DR
  label_05 <- paste0(
    "Denominator: Viability-low CD45+\n",
    "CD14+ monocyte-like: ", fmt_pct(safe_pct(df$cd14_mono_like, denom_cd45)), "\n",
    "CD14+HLA-DR-low: ", fmt_pct(safe_pct(df$cd14_hladr_low, denom_cd45)), "\n",
    "CD14+HLA-DR+: ", fmt_pct(safe_pct(df$cd14_hladr_pos, denom_cd45))
  )

  plot_list[["05_CD14_vs_HLA_DR_monocyte_like"]] <- make_flow_plot(
    df = df %>% filter(cd45_pos),
    xvar = "CD14",
    yvar = "HLA_DR",
    xlab = "CD14",
    ylab = "HLA-DR",
    xthr = thresholds$CD14,
    ythr = thresholds$HLA_DR,
    title = "CD14/HLA-DR monocyte-like phenotypes",
    subtitle = subtitle_base,
    label_text = label_05,
    max_events = 120000
  )

  # 06 CD11c vs HLA-DR
  label_06 <- paste0(
    "Denominator: Viability-low CD45+\n",
    "CD11c+HLA-DR+ APC-like myeloid: ", fmt_pct(safe_pct(df$cd11c_hladr_apc_like, denom_cd45)), "\n",
    "CD14+CD11c+HLA-DR+: ", fmt_pct(safe_pct(df$cd14_cd11c_hladr_pos, denom_cd45))
  )

  plot_list[["06_CD11c_vs_HLA_DR_APC_like"]] <- make_flow_plot(
    df = df %>% filter(cd45_pos),
    xvar = "CD11c",
    yvar = "HLA_DR",
    xlab = "CD11c",
    ylab = "HLA-DR",
    xthr = thresholds$CD11c,
    ythr = thresholds$HLA_DR,
    title = "APC-like myeloid axis",
    subtitle = subtitle_base,
    label_text = label_06,
    max_events = 120000
  )

  # 07 CCR3 vs CD66b
  label_07 <- paste0(
    "Denominator: Viability-low CD45+\n",
    "CCR3+CD66b+ eosinophil-like: ", fmt_pct(safe_pct(df$ccr3_cd66b_eosinophil_like, denom_cd45)), "\n",
    "CCR3+CD13+CD66b+: ", fmt_pct(safe_pct(df$ccr3_cd13_cd66b_gran_like, denom_cd45))
  )

  plot_list[["07_CCR3_vs_CD66b_eosinophil_like"]] <- make_flow_plot(
    df = df %>% filter(cd45_pos),
    xvar = "CCR3",
    yvar = "CD66b",
    xlab = "CD193 / CCR3",
    ylab = "CD66b",
    xthr = thresholds$CCR3,
    ythr = thresholds$CD66b,
    title = "CCR3/CD193-associated eosinophil-like granulocytic axis",
    subtitle = subtitle_base,
    label_text = label_07,
    max_events = 120000
  )

  # 08 CD123 vs HLA-DR
  label_08 <- paste0(
    "Denominator: Viability-low CD45+\n",
    "CD123+HLA-DR+ DC-like: ", fmt_pct(safe_pct(df$cd123_hladr_dc_like, denom_cd45)), "\n",
    "CD123+HLA-DR+CD11c-CD14- pDC-like: ", fmt_pct(safe_pct(df$cd123_hladr_cd11c_neg_pdc_like, denom_cd45))
  )

  plot_list[["08_CD123_vs_HLA_DR_DC_like"]] <- make_flow_plot(
    df = df %>% filter(cd45_pos),
    xvar = "CD123",
    yvar = "HLA_DR",
    xlab = "CD123",
    ylab = "HLA-DR",
    xthr = thresholds$CD123,
    ythr = thresholds$HLA_DR,
    title = "CD123/HLA-DR DC-like axis",
    subtitle = subtitle_base,
    label_text = label_08,
    max_events = 120000
  )

  # Save individual figures
  manifest_rows <- bind_rows(lapply(names(plot_list), function(gate_stub) {
    save_gate_plot(plot_list[[gate_stub]], sample_stub, gate_stub)
  }))

  # Save multipage PDF, one gate per page
  multipage_pdf <- file.path(multipage_dir, paste0(sample_stub, "__multipage_one_gate_per_page.pdf"))

  grDevices::pdf(multipage_pdf, width = 7.5, height = 7.5, onefile = TRUE)
  for (gate_stub in names(plot_list)) {
    print(plot_list[[gate_stub]])
  }
  grDevices::dev.off()

  manifest_rows <- manifest_rows %>%
    mutate(
      representative_role = role,
      subject_id = sid,
      disease_group = group,
      file_name = basename(fp),
      multipage_pdf = multipage_pdf
    )

  list(
    manifest = manifest_rows,
    gate_summary = gate_summary,
    marker_map = prep$marker_map %>%
      mutate(
        representative_role = role,
        subject_id = sid,
        file_name = basename(fp)
      )
  )
}

# ------------------------------------------------------------
# 9. Generate all representative plots
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 7B STARTED: FLOWJO-STYLE REAL FCS FIGURES\n")
cat("============================================================\n")

cat("\nRepresentative samples:\n")
print(as.data.frame(representative_samples), row.names = FALSE)

results <- vector("list", nrow(representative_samples))

for (i in seq_len(nrow(representative_samples))) {
  cat("\nProcessing representative sample ", i, " of ", nrow(representative_samples), ": ",
      representative_samples$subject_id[i], "\n", sep = "")
  results[[i]] <- build_sample_plots(representative_samples[i, ])
}

figure_manifest <- bind_rows(lapply(results, function(x) x$manifest))
gate_summary <- bind_rows(lapply(results, function(x) x$gate_summary))
marker_map_representatives <- bind_rows(lapply(results, function(x) x$marker_map))

# ------------------------------------------------------------
# 10. Save summary outputs
# ------------------------------------------------------------

write_csv(
  figure_manifest,
  file.path(out_dir, "SDY2583_CP10_STEP7B_flowjo_style_figure_manifest.csv")
)

write_csv(
  gate_summary,
  file.path(out_dir, "SDY2583_CP10_STEP7B_representative_gate_percentages.csv")
)

write_csv(
  marker_map_representatives,
  file.path(out_dir, "SDY2583_CP10_STEP7B_representative_marker_maps.csv")
)

save(
  representative_samples,
  figure_manifest,
  gate_summary,
  marker_map_representatives,
  thresholds,
  file = file.path(rdata_dir, "SDY2583_CP10_STEP7B_flowjo_style_representative_gating.RData")
)

# ------------------------------------------------------------
# 11. Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SDY2583 CP10 STEP 7B COMPLETE: FLOWJO-STYLE REAL FCS FIGURES\n")
cat("============================================================\n")

cat("\nRepresentative samples:\n")
print(as.data.frame(representative_samples), row.names = FALSE)

cat("\nRepresentative gate percentages:\n")
print(as.data.frame(gate_summary), row.names = FALSE)

cat("\nFigure manifest:\n")
print(as.data.frame(figure_manifest), row.names = FALSE)

cat("\nIndividual figures saved in:\n")
print(individual_dir)

cat("\nMultipage PDFs saved in:\n")
print(multipage_dir)

cat("\nNo synthetic events were generated. All panels were plotted from real CP10 FCS events.\n")

cat("============================================================\n")
