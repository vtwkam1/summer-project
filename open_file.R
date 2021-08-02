# Easy way to open specific scenario files

library(tidyverse)

### Write function
open_file <- function(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, file_ending, output_folder) {
    
  file_specific <- file_ending %>% str_replace("\\.", "\\.")
  
  scenario_pattern <- sprintf("^r%s_transmiss%s_crossimm%s_seed%s_vacc%s_start_%s", r0_strain1, transmiss, crossimm, seed, vacc, vacc_start) %>% str_replace_all("\\.", "pt")
  
  output_files <- list.files(file.path(".", output_folder))
  
  file_required <- str_subset(output_files, scenario_pattern) %>% str_subset(file_specific)
  
  current_dir <- getwd()
  
  shell(file.path(current_dir, output_folder, file_required), wait=F)
}

### Run function
output_folder <- "seir_model_output"

file_ending <- "plot_trim.png"

r0_strain1 <- 1.5
transmiss <- 1.25
crossimm <- 0.5
seed <- 35
vacc <- 0
vacc_start <- 0

open_file(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, file_ending, output_folder)
