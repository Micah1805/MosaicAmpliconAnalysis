# Master script to run all analysis steps in sequence with logging and error handling

log_file <- file("./output/run_pipeline_2.log", open = "wt")
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg), file = log_file, append = TRUE)
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

scripts <- list(
 # "./scripts/01_generate_asv_tables_prokaryotes_JH.R",
# "./scripts/02_generate_asv_tables_eukaryotes_JH.R",
  "./scripts/03_normalize_abundance_prokaryotes_JH.R",
  "./scripts/04_normalize_abundance_eukaryotes_JH.R",
  "./scripts/05_prepare_joint_tables_JH.R",
  "./scripts/06_create_taxa_environment_tables_JH.R",
  "./scripts/07_prepare_rda_environmental_data_JH.R"
)

log_message("Pipeline started.")

for (script in scripts) {
  log_message(paste("Starting script:", script))
  tryCatch(
    {
      source(script)
      log_message(paste("Finished script:", script))
    },
    error = function(e) {
      log_message(paste("ERROR in script:", script, ":", e$message))
    }
  )
}

log_message("Pipeline finished.")
close(log_file)
