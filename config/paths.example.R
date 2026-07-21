# Copy this file to config/paths.R, edit local values, and source it before
# running scripts. config/paths.R is ignored by Git.
Sys.setenv(
  SDY2583_REPO_ROOT = normalizePath(".", mustWork = FALSE),
  SDY2583_IMMPORT_DOWNLOAD_DIR = "/path/to/unpacked/ImmPort/SDY2583",
  SDY2583_METADATA_DIR = "/path/to/derived/metadata",
  SDY2583_INTEGRATED_DIR = "/path/to/derived/clinical_integration",
  SDY2583_INTEGRATION_RDATA_DIR = "/path/to/derived/integration_rdata",
  # Optional: use an existing compatible subject-level matrix instead of
  # rebuilding it from the ImmPort tabular package.
  SDY2583_METADATA_MATRIX_FILE = "/path/to/SDY2583_subject_FCS_metadata_matrix.csv",
  # Optional local cancer-subgroup/therapy annotation matrix.
  SDY2583_CLINICAL_ANNOTATION_FILE = "/path/to/SDY2583_clinical_annotation.csv",
  # Cross-panel convergence input/output. The input should be the subject-level
  # ALL10 matrix, for example ALL10_with_CP23.csv.
  SDY2583_CROSS_PANEL_MATRIX_FILE = "/path/to/SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv",
  SDY2583_CROSS_PANEL_OUT_DIR = "/path/to/derived/clinical_integration/09_cross_panel_convergence",
  SDY2583_CROSS_PANEL_MIN_N = "50",

  SDY2583_CP7_FCS_DIR = "/path/to/CP7/fcs",
  SDY2583_CP8_FCS_DIR = "/path/to/CP8/fcs",
  SDY2583_CP10_FCS_DIR = "/path/to/CP10/fcs",
  SDY2583_CP16_FCS_DIR = "/path/to/CP16/fcs",
  SDY2583_CP22_FCS_DIR = "/path/to/CP22/fcs",
  SDY2583_CP23_FCS_DIR = "/path/to/CP23/fcs",
  SDY2583_CP24_FCS_DIR = "/path/to/CP24/fcs",
  SDY2583_CP25_FCS_DIR = "/path/to/CP25/fcs",
  SDY2583_CP26_FCS_DIR = "/path/to/CP26/fcs",
  SDY2583_CP28_FCS_DIR = "/path/to/CP28/fcs",

  SDY2583_CP7_ANALYSIS_DIR = "/path/to/outputs/CP7",
  SDY2583_CP8_ANALYSIS_DIR = "/path/to/outputs/CP8",
  SDY2583_CP10_ANALYSIS_DIR = "/path/to/outputs/CP10",
  SDY2583_CP16_ANALYSIS_DIR = "/path/to/outputs/CP16",
  SDY2583_CP22_ANALYSIS_DIR = "/path/to/outputs/CP22",
  SDY2583_CP23_ANALYSIS_DIR = "/path/to/outputs/CP23",
  SDY2583_CP24_ANALYSIS_DIR = "/path/to/outputs/CP24",
  SDY2583_CP25_ANALYSIS_DIR = "/path/to/outputs/CP25",
  SDY2583_CP26_ANALYSIS_DIR = "/path/to/outputs/CP26",
  SDY2583_CP28_ANALYSIS_DIR = "/path/to/outputs/CP28",

  # Optional archived-output roots for deep table-level validation.
  SDY2583_CP7_REFERENCE_DIR = "/path/to/archive/CP7",
  SDY2583_CP8_REFERENCE_DIR = "/path/to/archive/CP8",
  SDY2583_CP10_REFERENCE_DIR = "/path/to/archive/CP10",
  SDY2583_CP16_REFERENCE_DIR = "/path/to/archive/CP16",
  SDY2583_CP22_REFERENCE_DIR = "/path/to/archive/CP22",
  SDY2583_CP23_REFERENCE_DIR = "/path/to/archive/CP23",
  SDY2583_CP24_REFERENCE_DIR = "/path/to/archive/CP24",
  SDY2583_CP25_REFERENCE_DIR = "/path/to/archive/CP25",
  SDY2583_CP26_REFERENCE_DIR = "/path/to/archive/CP26",
  SDY2583_CP28_REFERENCE_DIR = "/path/to/archive/CP28"
)
