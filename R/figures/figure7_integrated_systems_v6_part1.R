############################################################
# SDY2583 final integrated systems figure (v6; locked manuscript version)
# Clean academic layout with A-B-C left-to-right order
# v6 redesigns the cross-panel network panel for readability
############################################################

rm(list = ls())

############################
# 0) Packages
############################
packages_needed <- c(
  "tidyverse",
  "emmeans",
  "sandwich",
  "patchwork",
  "ggrepel",
  "scales",
  "grid"
)

for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(emmeans)
  library(sandwich)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
})

############################
# 1) User-facing options
############################
INPUT_FILE <- path.expand(Sys.getenv(
  "SDY2583_ALL10_MATRIX",
  unset = "SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv"
))
OUTPUT_DIR <- path.expand(Sys.getenv(
  "SDY2583_FIGURE7_OUT_DIR",
  unset = file.path("outputs", "Figure7_integrated_systems")
))

SHOW_FIGURE_TITLE <- TRUE
PNG_DPI <- 600

FIG_WIDTH  <- 22.0
FIG_HEIGHT <- 8.8

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

############################
# 2) Helper functions
############################
pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

to_clean_chr <- function(x) {
  y <- as.character(x)
  y <- trimws(y)
  y[y %in% c("", "NA", "NaN", "NULL", "None")] <- NA
  y
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

z_within <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

residualize <- function(data, ycol, covars) {
  use_cols <- c(ycol, covars)
  work <- data %>%
    mutate(.row_id_internal = dplyr::row_number()) %>%
    select(.row_id_internal, all_of(use_cols))

  dat <- work[complete.cases(work[, use_cols, drop = FALSE]), , drop = FALSE]

  if (nrow(dat) < 20) return(rep(NA_real_, nrow(data)))
  if (length(unique(dat[[ycol]])) < 2) return(rep(NA_real_, nrow(data)))

  form <- as.formula(
    paste(ycol, "~", paste(covars, collapse = " + "))
  )
  fit <- lm(form, data = dat)

  out <- rep(NA_real_, nrow(data))
  out[dat$.row_id_internal] <- resid(fit)
  out
}

ci95_mean <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) {
    return(c(mean = mean(x), lower = NA_real_, upper = NA_real_, n = n))
  }
  m  <- mean(x)
  se <- sd(x) / sqrt(n)
  tcrit <- qt(0.975, df = n - 1)
  c(mean = m, lower = m - tcrit * se, upper = m + tcrit * se, n = n)
}

choose_score <- function(df, panel, preferred_patterns = NULL) {
  panel_cols <- names(df)[grepl(paste0("^", panel, "_"), names(df))]
  if (length(panel_cols) == 0) return(NA_character_)

  integrated_cols <- panel_cols[grepl("integrated", panel_cols, ignore.case = TRUE)]
  if (length(integrated_cols) == 0) integrated_cols <- panel_cols

  if (!is.null(preferred_patterns) && length(preferred_patterns) > 0) {
    cand <- integrated_cols
    for (pat in preferred_patterns) {
      keep <- cand[grepl(pat, cand, ignore.case = TRUE)]
      if (length(keep) > 0) cand <- keep
    }
    if (length(cand) > 0) {
      nm <- sapply(cand, function(x) sum(!is.na(df[[x]])))
      return(cand[which.max(nm)])
    }
  }

  nm <- sapply(integrated_cols, function(x) sum(!is.na(df[[x]])))
  integrated_cols[which.max(nm)]
}

safe_write_csv <- function(x, file) {
  readr::write_csv(x, file)
}

############################
# 3) Read data and detect key columns
############################
if (!file.exists(INPUT_FILE)) {
  stop("Input CSV not found in working directory: ", INPUT_FILE)
}

raw <- readr::read_csv(INPUT_FILE, show_col_types = FALSE, progress = FALSE)

subject_col <- pick_first_existing(raw, c("subject_id", "subject_accession", "participant_id"))
age_col     <- pick_first_existing(raw, c("age_clinical", "age_for_model", "age", "age_raw"))
sex_col     <- pick_first_existing(raw, c("sex_clinical", "sex_for_model", "sex"))
disease_col <- pick_first_existing(raw, c("disease_group_model", "disease_group", "group", "disease"))
immuno_col  <- pick_first_existing(raw, c("original_cluster", "immunotype", "cluster", "immunotype_group"))

if (is.na(subject_col)) stop("No subject ID column found.")
if (is.na(age_col))     stop("No age column found.")
if (is.na(sex_col))     stop("No sex column found.")
if (is.na(disease_col)) stop("No disease group column found.")
if (is.na(immuno_col))  stop("No immunotype / cluster column found.")

############################
# 4) Define principal score map
############################
principal_map <- tibble::tribble(
  ~panel, ~domain, ~preferred_patterns,                           ~display_label,              ~network_label,           ~domain_color,
  "CP24", "CD8 / T-cell",      "cd8",                            "CP24\nCD8",                "CP24\nCD8",             "#1565C0",
  "CP7",  "CD8 / T-cell",      "checkpoint",                     "CP7\nCheckpoint",          "CP7\nCheckpoint",       "#1565C0",
  "CP28", "T/NK interface",    "t.*nk|tnk|interface",            "CP28\nT/NK interface",     "CP28\nT/NK",            "#FF6F2C",
  "CP26", "T/NK interface",    "nk",                             "CP26\nNK",                 "CP26\nNK",              "#FF6F2C",
  "CP8",  "CD4 / regulatory",  "cd4|helper|regulatory",          "CP8\nCD4 helper/reg",      "CP8\nCD4 helper",       "#4E9A06",
  "CP25", "CD4 / regulatory",  "regulatory|checkpoint|cd4",      "CP25\nCD4 reg-check",      "CP25\nCD4 reg-check",   "#4E9A06",
  "CP22", "Humoral B-cell",    "humoral|b.cell|bcell|ig",        "CP22\nHumoral B-cell",     "CP22\nB-cell",          "#8E7CC3",
  "CP10", "Myeloid / APC / monocyte", "myeloid|granul",          "CP10\nMyeloid/gran",       "CP10\nMyeloid/gran",    "#8E24AA",
  "CP16", "Myeloid / APC / monocyte", "apc|dc",                  "CP16\nAPC/DC-like",        "CP16\nAPC/DC-like",     "#8E24AA",
  "CP23", "Myeloid / APC / monocyte", "monocyte|macrophage|mono","CP23\nMono/mac-like",      "CP23\nMono/mac",        "#8E24AA"
)

principal_map$score_col <- NA_character_

for (i in seq_len(nrow(principal_map))) {
  pats <- unlist(strsplit(principal_map$preferred_patterns[i], "\\|"))
  principal_map$score_col[i] <- choose_score(raw, principal_map$panel[i], preferred_patterns = pats)
}

if (any(is.na(principal_map$score_col))) {
  missing_panels <- principal_map$panel[is.na(principal_map$score_col)]
  stop("Could not identify principal score(s) for panel(s): ",
       paste(missing_panels, collapse = ", "))
}

############################
# 5) Analysis dataset
############################
dat <- raw %>%
  transmute(
    subject_id = to_clean_chr(.data[[subject_col]]),
    age        = to_num(.data[[age_col]]),
    sex        = factor(to_clean_chr(.data[[sex_col]])),
    disease    = factor(to_clean_chr(.data[[disease_col]])),
    immunotype = factor(
      to_clean_chr(.data[[immuno_col]]),
      levels = c("G1_Naive", "G2_Primed", "G3_Progressive", "G4_Chronic", "G5_Suppressive")
    )
  )

for (i in seq_len(nrow(principal_map))) {
  dat[[principal_map$panel[i]]] <- to_num(raw[[principal_map$score_col[i]]])
}

score_panels <- principal_map$panel
complete_vars <- c("subject_id", "age", "sex", "disease", "immunotype", score_panels)

analysis_df <- dat %>%
  filter(complete.cases(across(all_of(complete_vars))))

cat("Complete-case analysis population: n =", nrow(analysis_df), "\n\n")

############################
# 6) PANEL A - Immunotype adjusted profiles
############################
emm_list <- list()

for (i in seq_len(nrow(principal_map))) {
  pn <- principal_map$panel[i]

  fit_df <- analysis_df %>%
    select(all_of(c("immunotype", "age", "sex", "disease", pn))) %>%
    rename(score = all_of(pn))

  fit <- lm(score ~ immunotype + age + sex + disease, data = fit_df)
  vc <- sandwich::vcovHC(fit, type = "HC3")

  emm <- emmeans::emmeans(fit, ~ immunotype, vcov. = vc) %>%
    as.data.frame() %>%
    mutate(
      panel = pn,
      score_col = principal_map$score_col[principal_map$panel == pn],
      domain = principal_map$domain[principal_map$panel == pn],
      display_label = principal_map$display_label[principal_map$panel == pn]
    )

  emm_list[[pn]] <- emm
}

emm_df <- bind_rows(emm_list) %>%
  group_by(panel) %>%
  mutate(profile_z = z_within(emmean)) %>%
  ungroup()

immu_pretty <- c(
  G1_Naive = "G1 Naive",
  G2_Primed = "G2 Primed",
  G3_Progressive = "G3 Progressive",
  G4_Chronic = "G4 Chronic",
  G5_Suppressive = "G5 Suppressive"
)

emm_df <- emm_df %>%
  mutate(
    immunotype_pretty = recode(as.character(immunotype), !!!immu_pretty),
    display_label = factor(display_label, levels = principal_map$display_label),
    immunotype_pretty = factor(
      immunotype_pretty,
      levels = rev(c("G1 Naive", "G2 Primed", "G3 Progressive", "G4 Chronic", "G5 Suppressive"))
    )
  )

safe_write_csv(
  emm_df %>%
    select(immunotype, panel, score_col, domain, emmean, SE, lower.CL, upper.CL, profile_z),
  file.path(OUTPUT_DIR, "immunotype_adjusted_marginal_means_principal_scores.csv")
)

pA <- ggplot(emm_df, aes(x = display_label, y = immunotype_pretty, fill = profile_z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", profile_z)), size = 3.0, color = "grey20") +
  scale_fill_gradient2(
    low = "#6A5ACD",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    oob = squish,
    name = "Adjusted\nprofile z"
  ) +
  labs(
    title = "Biological resolution of the five immunotypes",
    subtitle = "Adjusted marginal means",
    x = NULL,
    y = NULL,
    tag = "A"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, colour = "grey35"),
    axis.text.x = element_text(size = 8, angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 10, 8, 10)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    barwidth = unit(4.2, "cm"),
    barheight = unit(0.45, "cm")
  ))

