# Master script to run all analysis steps in sequence with logging and error handling

log_file <- file("./output/run_pipeline.log", open = "wt")
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg), file = log_file, append = TRUE)
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

scripts <- list(
  "./scripts/01_generate_asv_tables_prokaryotes_JH.R",
  "./scripts/02_generate_asv_tables_eukaryotes_JH.R",
  "./scripts/03_normalize_abundance_prokaryotes_JH.R",
  "./scripts/04_normalize_abundance_eukaryotes_JH.R",
  "./scripts/05_prepare_joint_tables_JH.R",
  "./scripts/06_create_taxa_environment_tables_JH.R",
  "./scripts/07_prepare_rda_environmental_data_JH.R",
  "./scripts/08_Create_microeco_objects.R",
  "./scripts/09_Calculate_and_plot_diversity.R",
  "./scripts/10_Create_a_heatmap_showing_the_beta-diversity.R",
  "./scripts/11_Change_in_species_composition_over_time.R",
  "./scripts/12_Change_in_species_composition_removed_taxa.R",
  "./scripts/13_Redundancy_analysis_(RDA)_with_interpolated_environmental_data.R"
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
