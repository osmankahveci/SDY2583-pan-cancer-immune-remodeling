# Copy this file to config/paths.R, edit the values, and source it before
# running a panel script. config/paths.R is ignored by Git.

Sys.setenv(
  SDY2583_REPO_ROOT = normalizePath(".", mustWork = FALSE),
  SDY2583_CP10_FCS_DIR = "/path/to/CP10/fcs",
  SDY2583_CP16_FCS_DIR = "/path/to/CP16/fcs",
  SDY2583_CP22_FCS_DIR = "/path/to/CP22/fcs",
  SDY2583_CP23_FCS_DIR = "/path/to/CP23/fcs",
  SDY2583_CP24_FCS_DIR = "/path/to/CP24/fcs",
  SDY2583_CP10_ANALYSIS_DIR = "/path/to/outputs/CP10",
  SDY2583_CP16_ANALYSIS_DIR = "/path/to/outputs/CP16",
  SDY2583_CP22_ANALYSIS_DIR = "/path/to/outputs/CP22",
  SDY2583_CP23_ANALYSIS_DIR = "/path/to/outputs/CP23",
  SDY2583_IMMPORT_DOWNLOAD_DIR = "/path/to/unpacked/ImmPort/SDY2583",
  SDY2583_METADATA_DIR = "/path/to/derived/metadata",
  SDY2583_INTEGRATED_DIR = "/path/to/derived/clinical_integration",
  SDY2583_INTEGRATION_RDATA_DIR = "/path/to/derived/integration_rdata"
)
