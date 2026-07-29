#!/usr/bin/env Rscript

# Run selected representative supplementary flow-cytometry atlases.
# Raw FCS files remain local and are supplied through environment variables.

available_panels <- c(
  "CP7", "CP8", "CP10", "CP16", "CP22",
  "CP23", "CP24", "CP25", "CP26", "CP28"
)

requested <- Sys.getenv(
  "SDY2583_FLOW_PANELS",
  unset = paste(available_panels, collapse = ",")
)
selected <- unique(trimws(strsplit(requested, ",", fixed = TRUE)[[1]]))
selected <- selected[nzchar(selected)]

unknown <- setdiff(selected, available_panels)
if (length(unknown) > 0L) {
  stop(
    "Unknown panel(s) in SDY2583_FLOW_PANELS: ",
    paste(unknown, collapse = ", ")
  )
}

master_script <- file.path(
  getwd(),
  "R",
  "figures",
  "supplementary_flow",
  "RECONSTRUCTED",
  "run_supplementary_flow_atlas.R"
)
if (!file.exists(master_script)) {
  stop(
    "Supplementary flow master script was not found. ",
    "Run this command from the repository root:\n",
    master_script
  )
}

status <- data.frame(
  panel = selected,
  status = NA_character_,
  message = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(selected)) {
  panel <- selected[i]
  message("Running supplementary flow atlas for ", panel)
  Sys.setenv(SDY2583_FLOW_PANEL = panel)

  result <- tryCatch(
    {
      sys.source(
        master_script,
        envir = new.env(parent = globalenv()),
        chdir = FALSE
      )
      list(status = "completed", message = "")
    },
    error = function(e) {
      list(status = "failed", message = conditionMessage(e))
    }
  )

  status$status[i] <- result$status
  status$message[i] <- result$message

  if (identical(result$status, "failed")) {
    print(status, row.names = FALSE)
    Sys.unsetenv("SDY2583_FLOW_PANEL")
    stop(
      "Supplementary flow workflow failed for ",
      panel,
      ": ",
      result$message
    )
  }
}

Sys.unsetenv("SDY2583_FLOW_PANEL")
print(status, row.names = FALSE)
message("Selected supplementary flow workflows completed.")
