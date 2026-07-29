#!/usr/bin/env Rscript

# Run selected representative supplementary flow-cytometry atlases.
# Raw FCS files remain local and are supplied through environment variables.

panel_scripts <- c(
  CP7  = "CP7_supplementary_flow_atlas.R",
  CP8  = "CP8_supplementary_flow_atlas.R",
  CP10 = "CP10_supplementary_flow_atlas.R",
  CP16 = "CP16_supplementary_flow_atlas.R",
  CP22 = "CP22_supplementary_flow_atlas.R",
  CP23 = "CP23_supplementary_flow_atlas.R",
  CP24 = "CP24_supplementary_flow_atlas.R",
  CP25 = "CP25_supplementary_flow_atlas.R",
  CP26 = "CP26_supplementary_flow_atlas.R",
  CP28 = "CP28_supplementary_flow_atlas.R"
)

requested <- Sys.getenv(
  "SDY2583_FLOW_PANELS",
  unset = paste(names(panel_scripts), collapse = ",")
)

selected <- unique(trimws(strsplit(requested, ",", fixed = TRUE)[[1]]))
selected <- selected[nzchar(selected)]

unknown <- setdiff(selected, names(panel_scripts))
if (length(unknown) > 0L) {
  stop(
    "Unknown panel(s) in SDY2583_FLOW_PANELS: ",
    paste(unknown, collapse = ", ")
  )
}

repo_root <- normalizePath(getwd(), mustWork = TRUE)
figure_dir <- file.path(
  repo_root,
  "R",
  "figures",
  "supplementary_flow",
  "RECONSTRUCTED"
)

if (!dir.exists(figure_dir)) {
  stop(
    "Supplementary flow script directory was not found. ",
    "Run this command from the repository root:\n",
    figure_dir
  )
}

status <- data.frame(
  panel = selected,
  script = unname(panel_scripts[selected]),
  status = NA_character_,
  message = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(selected)) {
  panel <- selected[i]
  script_path <- file.path(figure_dir, panel_scripts[[panel]])

  message("Running ", panel, ": ", script_path)

  result <- tryCatch(
    {
      sys.source(
        script_path,
        envir = new.env(parent = globalenv()),
        chdir = FALSE
      )
      list(status = "completed", message = "")
    },
    error = function(e) {
      list(
        status = "failed",
        message = conditionMessage(e)
      )
    }
  )

  status$status[i] <- result$status
  status$message[i] <- result$message

  if (identical(result$status, "failed")) {
    print(status, row.names = FALSE)
    stop(
      "Supplementary flow workflow failed for ",
      panel,
      ": ",
      result$message
    )
  }
}

print(status, row.names = FALSE)
message("Selected supplementary flow workflows completed.")
