############################
# 9) Combine panels in A-B-C order
############################
main_title <- "Coordinated cross-panel immune remodeling resolves the five immunotypes along shared and divergent systemic axes"

final_plot <- pA + pB + pC +
  patchwork::plot_layout(widths = c(1.08, 1.52, 1.10))

if (SHOW_FIGURE_TITLE) {
  final_plot <- final_plot +
    patchwork::plot_annotation(
      title = main_title,
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
      )
    )
}

############################
# 10) Save figure
############################
png_file <- file.path(OUTPUT_DIR, "Figure_integrated_cross_panel_immunotype_PCA_v6_600dpi.png")
pdf_file <- file.path(OUTPUT_DIR, "Figure_integrated_cross_panel_immunotype_PCA_v6_vector.pdf")

ggsave(
  filename = png_file,
  plot = final_plot,
  width = FIG_WIDTH,
  height = FIG_HEIGHT,
  dpi = PNG_DPI,
  limitsize = FALSE,
  bg = "white"
)

if (capabilities("cairo")) {
  ggsave(
    filename = pdf_file,
    plot = final_plot,
    width = FIG_WIDTH,
    height = FIG_HEIGHT,
    device = cairo_pdf,
    limitsize = FALSE,
    bg = "white"
  )
} else {
  ggsave(
    filename = pdf_file,
    plot = final_plot,
    width = FIG_WIDTH,
    height = FIG_HEIGHT,
    device = "pdf",
    limitsize = FALSE,
    bg = "white"
  )
}

############################
# 11) Write figure metrics / caption helper
############################
caption_lines <- c(
  "Figure X. Coordinated cross-panel remodeling resolves the five peripheral immunotypes along shared and divergent systemic immune axes.",
  "",
  "(A) Adjusted immunotype profiles across the ten principal integrated remodeling scores.",
  "Estimated marginal means were derived from models adjusted for age, sex, and disease group and standardized within score for visualization.",
  "",
  "(B) Network of the ten principal panel-level scores.",
  "Edges display selected strong, FDR-supported, biologically interpretable adjusted Spearman correlations.",
  paste0("Overall, ", n_sig_pairs, " of 45 principal-score pairs were FDR-significant."),
  "",
  "(C) Disease-adjusted PCA centroid map for the five immunotypes.",
  paste0("PC1 explains ", sprintf('%.1f', pc_var[1]), "% and PC2 explains ", sprintf('%.1f', pc_var[2]), "% of variance."),
  "PC1 represents global multi-compartment remodeling burden.",
  "PC2 separates a CD8/T-NK-oriented configuration from a myeloid/monocyte/NK-oriented configuration."
)

writeLines(caption_lines, con = file.path(OUTPUT_DIR, "Figure_integrated_caption_and_metrics_v6.txt"))

############################
# 12) Final console output
############################
cat("Completed.\n\n")
cat("PNG: ", normalizePath(png_file), "\n", sep = "")
cat("PDF: ", normalizePath(pdf_file), "\n", sep = "")
cat("Supporting tables: ", normalizePath(OUTPUT_DIR), "\n", sep = "")

invisible(final_plot)
