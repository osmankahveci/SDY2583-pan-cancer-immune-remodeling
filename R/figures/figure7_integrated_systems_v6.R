# Locked v6 integrated manuscript figure.
# Run from the repository root after setting SDY2583_ALL10_MATRIX.

input_path <- path.expand(Sys.getenv(
  "SDY2583_ALL10_MATRIX",
  unset = "SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv"
))
if (!file.exists(input_path)) stop("Input CSV not found: ", input_path)

# The locked v6 source predates the final clinical-column aliases. Build a
# temporary local copy with compatible aliases; no participant-level data are
# written inside the repository.
figure_input <- utils::read.csv(
  input_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
if ("age_for_clinical_model" %in% names(figure_input) &&
    !"age_for_model" %in% names(figure_input)) {
  figure_input$age_for_model <- figure_input$age_for_clinical_model
}
if ("sex_for_clinical_model" %in% names(figure_input) &&
    !"sex_for_model" %in% names(figure_input)) {
  figure_input$sex_for_model <- figure_input$sex_for_clinical_model
}

temporary_input <- tempfile(pattern = "SDY2583_figure7_", fileext = ".csv")
utils::write.csv(figure_input, temporary_input, row.names = FALSE, na = "")
Sys.setenv(
  SDY2583_FIGURE7_ORIGINAL_INPUT = input_path,
  SDY2583_FIGURE7_TEMP_INPUT = temporary_input,
  SDY2583_ALL10_MATRIX = temporary_input
)

source(file.path("R", "figures", "figure7_integrated_systems_v6_part1.R"))
source(file.path("R", "figures", "figure7_integrated_systems_v6_part2.R"))
source(file.path("R", "figures", "figure7_integrated_systems_v6_part3.R"))

unlink(Sys.getenv("SDY2583_FIGURE7_TEMP_INPUT"))
Sys.setenv(SDY2583_ALL10_MATRIX = Sys.getenv("SDY2583_FIGURE7_ORIGINAL_INPUT"))
Sys.unsetenv(c("SDY2583_FIGURE7_ORIGINAL_INPUT", "SDY2583_FIGURE7_TEMP_INPUT"))
