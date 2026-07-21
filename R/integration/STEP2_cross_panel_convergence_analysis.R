# SDY2583 cross-panel convergence analysis
# Participant-level matrices and outputs remain local and are not committed.

rm(list = ls())
root <- Sys.getenv("SDY2583_REPO_ROOT", unset = ".")
source(file.path(root, "R", "shared", "bootstrap.R"))

pkgs <- c("dplyr", "readr", "tidyr", "purrr", "stringr", "tibble",
          "ggplot2", "igraph", "ggraph", "ggrepel", "scales")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr)
  library(stringr); library(tibble); library(ggplot2)
  library(igraph); library(ggraph); library(ggrepel); library(scales)
})

pick <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) hit[1] else NA_character_
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
chr <- function(x) {
  y <- trimws(as.character(x)); y[y %in% c("", "NA", "NULL", "None")] <- NA_character_; y
}
disease_group <- function(x) {
  y <- tolower(chr(x))
  case_when(str_detect(y, "cancer|patient|tumou?r|malignan") ~ "Cancer",
            str_detect(y, "healthy|control") ~ "Healthy", TRUE ~ NA_character_)
}

residualize <- function(dat, score, covars, min_n) {
  z <- data.frame(y = num(dat[[score]]))
  for (v in covars[covars %in% names(dat)]) z[[v]] <- dat[[v]]
  keep <- is.finite(z$y) & complete.cases(z)
  out <- rep(NA_real_, nrow(dat))
  if (sum(keep) < min_n || sd(z$y[keep]) == 0) return(out)
  usable <- setdiff(names(z), "y")
  usable <- usable[vapply(usable, function(v) length(unique(z[[v]][keep])) > 1, logical(1))]
  f <- as.formula(if (length(usable)) paste("y ~", paste(usable, collapse = "+")) else "y ~ 1")
  fit <- tryCatch(lm(f, data = z[keep, , drop = FALSE]), error = function(e) NULL)
  if (!is.null(fit)) out[keep] <- residuals(fit)
  out
}

cor_pairs <- function(dat, registry, cohort, min_n) {
  pairs <- combn(registry$score, 2, simplify = FALSE)
  out <- map_dfr(pairs, function(p) {
    x <- num(dat[[p[1]]]); y <- num(dat[[p[2]]]); ok <- is.finite(x) & is.finite(y)
    n <- sum(ok)
    if (n < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
      return(tibble(cohort, score_1 = p[1], score_2 = p[2], n_pair = n,
                    rho = NA_real_, p_value = NA_real_))
    }
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    tibble(cohort, score_1 = p[1], score_2 = p[2], n_pair = n,
           rho = unname(ct$estimate), p_value = ct$p.value)
  })
  key <- registry %>% select(score, panel, label)
  out %>%
    left_join(key, by = c("score_1" = "score")) %>% rename(panel_1 = panel, label_1 = label) %>%
    left_join(key, by = c("score_2" = "score")) %>% rename(panel_2 = panel, label_2 = label) %>%
    mutate(cross_panel = panel_1 != panel_2,
           panel_low = pmin(panel_1, panel_2), panel_high = pmax(panel_1, panel_2),
           panel_pair = paste(panel_low, panel_high, sep = "–"),
           pair_id = paste(pmin(score_1, score_2), pmax(score_1, score_2), sep = "||"),
           q_all = p.adjust(p_value, method = "BH")) %>%
    group_by(cross_panel) %>% mutate(q_scope = p.adjust(p_value, method = "BH")) %>% ungroup()
}

# Input: explicit file first; otherwise latest integrated matrix under SDY2583_INTEGRATED_DIR.
explicit <- path.expand(Sys.getenv("SDY2583_CROSS_PANEL_MATRIX_FILE", unset = ""))
dir_in <- sd_integrated_dir()
files <- if (nzchar(explicit) && file.exists(explicit)) explicit else list.files(
  dir_in, pattern = "SDY2583_integrated_clinical_immune_score_matrix_.*\\.csv$",
  recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("Set SDY2583_CROSS_PANEL_MATRIX_FILE to the ALL10 subject-level CSV.")
priority <- 100 * grepl("ALL10", basename(files), ignore.case = TRUE) +
  80 * grepl("with_CP23", basename(files), ignore.case = TRUE) +
  60 * grepl("ALL_PANELS", basename(files), ignore.case = TRUE) -
  100 * grepl("summary|manifest|dictionary", basename(files), ignore.case = TRUE)
input_file <- files[order(priority, file.info(files)$mtime, decreasing = TRUE)][1]
out_dir <- path.expand(Sys.getenv("SDY2583_CROSS_PANEL_OUT_DIR",
                                  unset = file.path(dir_in, "09_cross_panel_convergence")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
min_n <- as.integer(Sys.getenv("SDY2583_CROSS_PANEL_MIN_N", unset = "50"))
if (is.na(min_n) || min_n < 20) min_n <- 50L

raw <- read_csv(input_file, show_col_types = FALSE, progress = FALSE)
subject_col <- pick(raw, c("subject_id", "subject_accession", "participant_id"))
age_col <- pick(raw, c("age_for_model", "age_clinical", "age", "age_raw"))
sex_col <- pick(raw, c("sex_for_model", "sex_for_clinical_model", "sex"))
dis_col <- pick(raw, c("disease_group_model", "disease_group_for_clinical",
                       "disease_group", "disease", "group"))
if (is.na(subject_col)) stop("Subject identifier not found.")

score_cols <- grep("^CP[0-9]+_.*_score$|^CP[0-9]+_integrated_", names(raw), value = TRUE)
score_cols <- score_cols[vapply(score_cols, function(s) {
  x <- num(raw[[s]]); sum(is.finite(x)) >= min_n && length(unique(x[is.finite(x)])) > 1
}, logical(1))]
if (length(score_cols) < 2) stop("Fewer than two analyzable score columns found.")

dat <- raw %>% transmute(
  subject_id = as.character(.data[[subject_col]]),
  age_model = if (!is.na(age_col)) num(.data[[age_col]]) else NA_real_,
  sex_model = if (!is.na(sex_col)) factor(chr(.data[[sex_col]])) else factor(NA_character_),
  disease_model = if (!is.na(dis_col)) factor(disease_group(.data[[dis_col]])) else factor(NA_character_),
  across(all_of(score_cols), num)
) %>% filter(!is.na(subject_id), subject_id != "") %>% distinct(subject_id, .keep_all = TRUE)

registry <- tibble(score = score_cols) %>% mutate(
  panel = str_extract(score, "^CP[0-9]+"),
  short_label = score %>% str_remove("^CP[0-9]+_") %>% str_remove("_score$") %>%
    str_replace_all("_", " ") %>% str_squish(),
  label = paste0(panel, ": ", short_label),
  priority = case_when(str_detect(str_to_lower(score), "integrated") ~ 1L,
                       str_detect(str_to_lower(score), "core") ~ 2L,
                       str_detect(str_to_lower(score), "composition") ~ 3L, TRUE ~ 4L),
  n_available = vapply(score, function(s) sum(is.finite(dat[[s]])), integer(1)),
  missing_fraction = 1 - n_available / nrow(dat)
) %>% arrange(as.integer(str_remove(panel, "CP")), priority, score)
principal <- registry %>% group_by(panel) %>% arrange(priority, desc(n_available), .by_group = TRUE) %>%
  slice(1) %>% ungroup()

full <- dat
for (s in score_cols) full[[s]] <- residualize(dat, s, c("age_model", "sex_model", "disease_model"), min_n)
cancer <- dat %>% filter(disease_model == "Cancer")
for (s in score_cols) cancer[[s]] <- residualize(cancer, s, c("age_model", "sex_model"), min_n)
healthy <- dat %>% filter(disease_model == "Healthy")
for (s in score_cols) healthy[[s]] <- residualize(healthy, s, c("age_model", "sex_model"), min_n)

r_full <- cor_pairs(full, registry, "All: age/sex/disease adjusted", min_n)
r_cancer <- cor_pairs(cancer, registry, "Cancer: age/sex adjusted", min_n)
r_healthy <- cor_pairs(healthy, registry, "Healthy: age/sex adjusted", min_n)
all_results <- bind_rows(r_full, r_cancer, r_healthy)
cross <- all_results %>% filter(cross_panel, is.finite(rho), is.finite(p_value))

panel_summary <- cross %>% group_by(cohort, panel_low, panel_high, panel_pair) %>% summarise(
  n_score_pairs = n(), median_rho = median(rho), median_abs_rho = median(abs(rho)),
  q25_rho = quantile(rho, .25), q75_rho = quantile(rho, .75),
  proportion_positive = mean(rho > 0), n_fdr = sum(q_scope < .05),
  n_positive_fdr = sum(q_scope < .05 & rho > 0), proportion_fdr = mean(q_scope < .05),
  .groups = "drop")

axes <- tribble(
  ~axis, ~panel_a, ~panel_b,
  "CD8 differentiation/checkpoint", "CP7", "CP24",
  "CD8/T-cell interface", "CP7", "CP28",
  "CD8/T-cell interface", "CP24", "CP28",
  "CD4 helper/regulatory", "CP8", "CP25",
  "Myeloid/APC", "CP10", "CP16",
  "Myeloid/monocyte", "CP10", "CP23",
  "APC/monocyte", "CP16", "CP23",
  "NK/T–NK interface", "CP26", "CP28") %>%
  mutate(panel_low = pmin(panel_a, panel_b), panel_high = pmax(panel_a, panel_b))
prespecified <- panel_summary %>% inner_join(axes %>% select(axis, panel_low, panel_high),
                                             by = c("panel_low", "panel_high")) %>%
  select(axis, cohort, everything())

robust <- r_full %>% filter(cross_panel) %>%
  select(pair_id, score_1, score_2, panel_1, panel_2, n_full = n_pair,
         rho_full = rho, q_full = q_scope) %>%
  left_join(r_cancer %>% filter(cross_panel) %>%
              select(pair_id, n_cancer = n_pair, rho_cancer = rho, q_cancer = q_scope), by = "pair_id") %>%
  left_join(r_healthy %>% filter(cross_panel) %>%
              select(pair_id, n_healthy = n_pair, rho_healthy = rho, q_healthy = q_scope), by = "pair_id") %>%
  mutate(same_direction_cancer = sign(rho_full) == sign(rho_cancer),
         same_direction_healthy = sign(rho_full) == sign(rho_healthy),
         robust_positive = q_full < .05 & rho_full > 0 & same_direction_cancer &
           (is.na(same_direction_healthy) | same_direction_healthy)) %>%
  arrange(desc(robust_positive), q_full, desc(abs(rho_full)))

# Panel-level heatmap.
panel_levels <- registry %>% distinct(panel) %>%
  mutate(n = as.integer(str_remove(panel, "CP"))) %>% arrange(n) %>% pull(panel)
plot_panel <- panel_summary %>% filter(cohort == "All: age/sex/disease adjusted") %>%
  mutate(panel_low = factor(panel_low, panel_levels), panel_high = factor(panel_high, panel_levels))
p1 <- ggplot(plot_panel, aes(panel_low, panel_high, fill = median_rho)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", median_rho)), size = 3) +
  scale_fill_gradient2(midpoint = 0, limits = c(-1, 1), oob = squish) + coord_fixed() +
  labs(title = "Cross-panel immune-remodeling convergence",
       subtitle = "Median adjusted Spearman correlation across score pairs",
       x = NULL, y = NULL, fill = "Median rho") + theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "SDY2583_cross_panel_panel_level_heatmap.png"), p1,
       width = 8.5, height = 7.5, dpi = 500)
ggsave(file.path(out_dir, "SDY2583_cross_panel_panel_level_heatmap.pdf"), p1,
       width = 8.5, height = 7.5)

# One principal integrated/core score per panel.
ps <- principal$score
sym <- r_full %>% filter(score_1 %in% ps, score_2 %in% ps) %>%
  select(score_1, score_2, rho, q_scope)
sym <- bind_rows(sym, sym %>% transmute(score_1 = score_2, score_2 = score_1, rho, q_scope),
                 tibble(score_1 = ps, score_2 = ps, rho = 1, q_scope = 0)) %>% distinct(score_1, score_2)
labels <- principal %>% select(score, label)
sym <- sym %>% left_join(labels, by = c("score_1" = "score")) %>% rename(label_1 = label) %>%
  left_join(labels, by = c("score_2" = "score")) %>% rename(label_2 = label) %>%
  mutate(label_1 = factor(label_1, principal$label), label_2 = factor(label_2, rev(principal$label)))
p2 <- ggplot(sym, aes(label_1, label_2, fill = rho)) + geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.7) +
  scale_fill_gradient2(midpoint = 0, limits = c(-1, 1), oob = squish) + coord_fixed() +
  labs(title = "Principal cross-panel remodeling axes",
       subtitle = "Adjusted residual correlations; one automatically selected score per panel",
       x = NULL, y = NULL, fill = "Spearman rho") + theme_minimal(base_size = 9.5) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "SDY2583_principal_score_convergence_heatmap.png"), p2,
       width = 10, height = 9, dpi = 500)
ggsave(file.path(out_dir, "SDY2583_principal_score_convergence_heatmap.pdf"), p2,
       width = 10, height = 9)

edges <- r_full %>% filter(score_1 %in% ps, score_2 %in% ps, q_scope < .05, abs(rho) >= .20) %>%
  transmute(from = score_1, to = score_2, rho, abs_rho = abs(rho),
            direction = ifelse(rho >= 0, "Positive", "Negative"))
nodes <- principal %>% transmute(name = score, panel, label)
if (nrow(edges)) {
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  set.seed(2583)
  p3 <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = abs_rho, linetype = direction), alpha = .65) +
    geom_node_point(aes(shape = panel), size = 4) +
    geom_node_text(aes(label = label), repel = TRUE, size = 3) +
    scale_edge_width(range = c(.4, 2.5)) + theme_void() +
    labs(title = "Network of principal remodeling axes",
         subtitle = "Edges: adjusted FDR < 0.05 and |rho| >= 0.20",
         edge_width = "|rho|", edge_linetype = "Direction", shape = "Panel")
  ggsave(file.path(out_dir, "SDY2583_principal_score_convergence_network.png"), p3,
         width = 10, height = 8, dpi = 500)
  ggsave(file.path(out_dir, "SDY2583_principal_score_convergence_network.pdf"), p3,
         width = 10, height = 8)
}

write_csv(registry, file.path(out_dir, "SDY2583_cross_panel_score_registry.csv"))
write_csv(principal, file.path(out_dir, "SDY2583_principal_score_selection.csv"))
write_csv(all_results, file.path(out_dir, "SDY2583_all_pairwise_score_correlations.csv"))
write_csv(cross, file.path(out_dir, "SDY2583_cross_panel_score_correlations.csv"))
write_csv(panel_summary, file.path(out_dir, "SDY2583_panel_pair_convergence_summary.csv"))
write_csv(prespecified, file.path(out_dir, "SDY2583_prespecified_axis_convergence_summary.csv"))
write_csv(robust, file.path(out_dir, "SDY2583_cross_panel_convergence_robustness.csv"))
write_csv(tibble(input_file = normalizePath(input_file), n_subjects = nrow(dat),
                 n_scores = length(score_cols), n_panels = n_distinct(registry$panel),
                 n_cancer = sum(dat$disease_model == "Cancer", na.rm = TRUE),
                 n_healthy = sum(dat$disease_model == "Healthy", na.rm = TRUE), min_pairwise_n = min_n),
          file.path(out_dir, "SDY2583_cross_panel_convergence_manifest.csv"))
writeLines(c(
  "SDY2583 CROSS-PANEL CONVERGENCE ANALYSIS",
  paste("Input:", input_file), paste("Subjects:", nrow(dat)), paste("Scores:", length(score_cols)),
  paste("Panels:", n_distinct(registry$panel)),
  paste("Positive FDR-significant cross-panel pairs:",
        sum(r_full$cross_panel & r_full$q_scope < .05 & r_full$rho > 0, na.rm = TRUE)),
  paste("Robust positive pairs:", sum(robust$robust_positive, na.rm = TRUE))),
  file.path(out_dir, "SDY2583_cross_panel_convergence_summary.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "SDY2583_cross_panel_convergence_sessionInfo.txt"))
cat("Cross-panel convergence analysis complete. Outputs: ", out_dir, "\n", sep = "")
