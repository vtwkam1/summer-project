# Easy way to open specific scenario files

library(tidyverse)

### Write function
open_file <- function(r0_strain1, transmiss, crossimm, preinf, seed, vacc, vacc_start, file_ending, output_folder) {
    
  file_specific <- file_ending %>% str_replace("\\.", "\\\\\\.")
  
  scenario_pattern <- sprintf("^r%s_transmiss%s_crossimm%s_seed%s_preinf%s_vacc%s_start_%s", r0_strain1, transmiss, crossimm, seed, preinf, vacc, vacc_start) %>% str_replace_all("\\.", "pt")
  
  output_files <- list.files(file.path(".", output_folder))
  
  file_required <- str_subset(output_files, scenario_pattern) %>% str_subset(file_specific)
  
  current_dir <- getwd()
  
  shell(file.path(current_dir, output_folder, file_required), wait=F)
}

### Run function
output_folder <- "seir_model_output"

file_ending <- "plot_full.png"

r0_strain1 <- 1.5
transmiss <- 1.06
crossimm <- 1
preinf <- 2.5
seed <- 301.16
vacc <- 1
vacc_start <- 0.1111

open_file(r0_strain1, transmiss, crossimm, preinf, seed, vacc, vacc_start, file_ending, output_folder)

# Move selected runs to Archive
selected_scenarios <- str_c("^", run_longer$scenario)

output_files <- list.files(file.path(".", output_folder))

files_required <- map(selected_scenarios, str_subset, string = output_files) %>%  unlist(use.names=F)

move_file <- function(x) {
  file.rename(from = file.path(".", output_folder, x),
              to = file.path(".", "Archive", "seir_model_output_cutshort", x)) 
  }

lapply(files_required, move_file)