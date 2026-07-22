############################
# 7) PANEL B - Cross-panel convergence network
############################
# Use pairwise-available data to match the main principal-score analysis.
network_base <- dat %>%
  filter(complete.cases(across(all_of(c("age", "sex", "disease")))))

for (pn in score_panels) {
  network_base[[pn]] <- residualize(network_base, pn, c("age", "sex", "disease"))
}

pairs <- combn(score_panels, 2, simplify = FALSE)

cor_tbl <- purrr::map_dfr(pairs, function(pr) {
  x <- network_base[[pr[1]]]
  y <- network_base[[pr[2]]]
  ok <- is.finite(x) & is.finite(y)
  n_pair <- sum(ok)
  if (n_pair < 20) {
    return(tibble(panel_1 = pr[1], panel_2 = pr[2], n = n_pair, rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(
    panel_1 = pr[1],
    panel_2 = pr[2],
    n = n_pair,
    rho = unname(ct$estimate),
    p = ct$p.value
  )
}) %>%
  mutate(q = p.adjust(p, method = "BH"), significant = q < 0.05)

safe_write_csv(cor_tbl, file.path(OUTPUT_DIR, "principal_score_adjusted_spearman_correlations.csv"))

n_sig_pairs <- sum(cor_tbl$significant, na.rm = TRUE)

# Manual node layout for a cleaner academic network panel
node_tbl <- tibble::tribble(
  ~panel, ~x,  ~y,
  "CP24", 1.15, 3.20,
  "CP7",  1.15, 2.18,
  "CP28", 3.10, 3.20,
  "CP26", 3.10, 2.18,
  "CP8",  5.55, 3.20,
  "CP25", 5.55, 2.18,
  "CP22", 7.85, 2.75,
  "CP10", 2.25, 0.72,
  "CP16", 4.15, 0.30,
  "CP23", 6.15, 0.72
) %>%
  left_join(principal_map %>% select(panel, network_label, domain_color), by = "panel")

# Manually controlled edges so lines do not run through labels
edge_manual <- tibble::tribble(
  ~panel_1, ~panel_2, ~x,   ~y,   ~xend, ~yend, ~label_x, ~label_y,
  "CP24", "CP28",  1.55, 3.20, 2.72,  3.20,  2.10,     3.42,
  "CP7",  "CP24",  1.15, 2.50, 1.15,  2.88,  0.58,     2.72,
  "CP7",  "CP28",  1.55, 2.35, 2.72,  2.88,  2.08,     2.73,
  "CP28", "CP26",  3.10, 2.88, 3.10,  2.50,  3.48,     2.70,
  "CP8",  "CP25",  5.55, 2.88, 5.55,  2.50,  5.93,     2.70,
  "CP10", "CP23",  2.68, 0.72, 5.72,  0.72,  4.18,     0.96,
  "CP16", "CP23",  4.55, 0.40, 5.72,  0.62,  5.18,     0.48
)

# Attach actual rho values from the analysis table
edge_tbl <- edge_manual %>%
  left_join(cor_tbl %>% select(panel_1, panel_2, rho, q, significant), by = c("panel_1", "panel_2"))

# Handle potentially reversed pair ordering
missing_idx <- which(is.na(edge_tbl$rho))
if (length(missing_idx) > 0) {
  rev_lookup <- cor_tbl %>%
    transmute(panel_1 = panel_2, panel_2 = panel_1, rho, q, significant)
  edge_tbl2 <- edge_tbl[missing_idx, ] %>%
    select(panel_1, panel_2) %>%
    left_join(rev_lookup, by = c("panel_1", "panel_2"))
  edge_tbl$rho[missing_idx] <- edge_tbl2$rho
  edge_tbl$q[missing_idx] <- edge_tbl2$q
  edge_tbl$significant[missing_idx] <- edge_tbl2$significant
}

edge_tbl <- edge_tbl %>%
  filter(significant) %>%
  mutate(
    edge_lab = sprintf("ρ = %.3f", rho),
    line_w = scales::rescale(abs(rho), to = c(1.0, 2.2))
  )

# Domain headers
network_headers <- tibble::tribble(
  ~label,                      ~x,   ~y,   ~col,
  "CD8 / T-cell",             1.15, 3.72, "#1565C0",
  "T/NK interface",           3.10, 3.72, "#FF6F2C",
  "CD4 / regulatory",         5.55, 3.72, "#4E9A06",
  "Humoral B-cell",           7.85, 3.25, "#8E7CC3",
  "Myeloid / APC / monocyte", 4.15, 1.22, "#8E24AA"
)

pB <- ggplot() +
  geom_segment(
    data = edge_tbl,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = line_w),
    color = alpha("#5A5A5A", 0.90),
    lineend = "round"
  ) +
  geom_label(
    data = edge_tbl,
    aes(x = label_x, y = label_y, label = edge_lab),
    size = 3.3,
    fill = alpha("white", 0.95),
    color = "#333333",
    label.size = 0,
    label.padding = unit(0.10, "lines")
  ) +
  geom_label(
    data = node_tbl,
    aes(x = x, y = y, label = network_label),
    fill = alpha("white", 0.96),
    color = node_tbl$domain_color,
    label.size = 0.8,
    label.r = unit(0.12, "lines"),
    label.padding = unit(0.18, "lines"),
    size = 3.9,
    fontface = "bold",
    lineheight = 0.93
  ) +
  geom_text(
    data = network_headers,
    aes(x = x, y = y, label = label),
    color = network_headers$col,
    size = 4.7,
    fontface = "bold"
  ) +
  annotate(
    "label",
    x = 4.05, y = -0.18,
    label = paste0(
      n_sig_pairs,
      "/45 principal-score pairs FDR-significant; directions preserved in cancer-only and healthy-only analyses"
    ),
    size = 3.8,
    fill = alpha("white", 0.96),
    color = "#333333",
    label.size = 0.6,
    label.r = unit(0.16, "lines"),
    label.padding = unit(0.20, "lines")
  ) +
  scale_linewidth_identity() +
  coord_cartesian(xlim = c(0.15, 8.55), ylim = c(-0.35, 4.05), clip = "off") +
  labs(
    title = "Cross-panel convergence network",
    subtitle = "Selected FDR-supported relationships among the 10 principal integrated scores",
    x = NULL, y = NULL, tag = "B"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, colour = "grey35"),
    plot.margin = margin(8, 12, 8, 12)
  )

############################
# 8) PANEL C - PCA centroids
############################
pca_df <- analysis_df %>%
  select(immunotype, age, sex, disease, all_of(score_panels))

pca_resid <- pca_df
for (pn in score_panels) {
  pca_resid[[pn]] <- residualize(pca_df, pn, c("age", "sex", "disease"))
}

pca_mat <- pca_resid %>% select(all_of(score_panels)) %>% as.matrix()
pca_mat <- scale(pca_mat)

pc_fit <- prcomp(pca_mat, center = FALSE, scale. = FALSE)
pc_var <- 100 * (pc_fit$sdev^2 / sum(pc_fit$sdev^2))

pc_scores <- as.data.frame(pc_fit$x[, 1:2])
names(pc_scores) <- c("PC1", "PC2")
pc_scores$immunotype <- pca_df$immunotype

pc1_load <- pc_fit$rotation[, 1]
if (mean(pc1_load, na.rm = TRUE) < 0) {
  pc_scores$PC1 <- -pc_scores$PC1
  pc_fit$rotation[, 1] <- -pc_fit$rotation[, 1]
}

myeloid_panels <- c("CP10", "CP16", "CP23", "CP26")
cd8_tnk_panels <- c("CP24", "CP28", "CP7", "CP25")
pc2_load <- pc_fit$rotation[, 2]
myeloid_mean <- mean(pc2_load[rownames(pc_fit$rotation) %in% myeloid_panels], na.rm = TRUE)
cd8_mean     <- mean(pc2_load[rownames(pc_fit$rotation) %in% cd8_tnk_panels], na.rm = TRUE)

if ((myeloid_mean - cd8_mean) < 0) {
  pc_scores$PC2 <- -pc_scores$PC2
  pc_fit$rotation[, 2] <- -pc_fit$rotation[, 2]
}

centroids <- pc_scores %>%
  group_by(immunotype) %>%
  summarise(
    PC1_mean = ci95_mean(PC1)[["mean"]],
    PC1_low  = ci95_mean(PC1)[["lower"]],
    PC1_high = ci95_mean(PC1)[["upper"]],
    PC2_mean = ci95_mean(PC2)[["mean"]],
    PC2_low  = ci95_mean(PC2)[["lower"]],
    PC2_high = ci95_mean(PC2)[["upper"]],
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = recode(as.character(immunotype), !!!immu_pretty))

safe_write_csv(centroids, file.path(OUTPUT_DIR, "disease_adjusted_pca_immunotype_centroids.csv"))

loadings_df <- as.data.frame(pc_fit$rotation[, 1:2])
loadings_df$panel <- rownames(loadings_df)
loadings_df <- loadings_df %>%
  left_join(principal_map %>% select(panel, score_col, domain, display_label), by = "panel")

safe_write_csv(
  tibble(
    panel = loadings_df$panel,
    score_col = loadings_df$score_col,
    domain = loadings_df$domain,
    display_label = loadings_df$display_label,
    PC1_loading = loadings_df$PC1,
    PC2_loading = loadings_df$PC2
  ),
  file.path(OUTPUT_DIR, "disease_adjusted_pca_loadings.csv")
)

pC <- ggplot(centroids, aes(x = PC1_mean, y = PC2_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey75") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey75") +
  geom_segment(aes(x = PC1_low, xend = PC1_high, y = PC2_mean, yend = PC2_mean),
               linewidth = 0.8, color = "grey50") +
  geom_segment(aes(x = PC1_mean, xend = PC1_mean, y = PC2_low, yend = PC2_high),
               linewidth = 0.8, color = "grey50") +
  geom_point(aes(color = immunotype), size = 4.2, show.legend = FALSE) +
  ggrepel::geom_text_repel(
    aes(label = label, color = immunotype),
    size = 4.0,
    fontface = "bold",
    show.legend = FALSE,
    max.overlaps = Inf,
    seed = 2583,
    segment.color = NA
  ) +
  annotate(
    "text",
    x = max(centroids$PC1_high, na.rm = TRUE) * 0.74,
    y = max(centroids$PC2_high, na.rm = TRUE) * 1.12,
    label = "Myeloid / monocyte / NK pole",
    color = "#8E24AA",
    fontface = "bold",
    size = 4.1
  ) +
  annotate(
    "text",
    x = max(centroids$PC1_high, na.rm = TRUE) * 0.30,
    y = min(centroids$PC2_low, na.rm = TRUE) * 1.16,
    label = "CD8 / T-NK pole",
    color = "#1565C0",
    fontface = "bold",
    size = 4.1
  ) +
  scale_color_manual(
    values = c(
      "G1_Naive" = "#80B1D3",
      "G2_Primed" = "#7FC97F",
      "G3_Progressive" = "#BC80BD",
      "G4_Chronic" = "#1565C0",
      "G5_Suppressive" = "#D73027"
    )
  ) +
  labs(
    title = "Global remodeling dimensions",
    subtitle = "Immunotype centroids with 95% confidence intervals",
    x = paste0("PC1: global remodeling burden (", sprintf("%.1f", pc_var[1]), "%)"),
    y = paste0("PC2: remodeling configuration (", sprintf("%.1f", pc_var[2]), "%)"),
    tag = "C"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, colour = "grey35"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(8, 12, 8, 12)
  )

