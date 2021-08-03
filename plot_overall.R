################ Change vaccstart str extract

library(tidyverse)
# library(RColorBrewer)

theme_set(
    theme_gray(base_size = 12)
)

output_folder <- "seir_model_output_transmiss"

end <- read.csv(file.path(".", output_folder, "end.csv"))

end %>% 
    count(match_set)

# Extract parameters
end <- end %>%
    mutate(scenario_pt = str_replace_all(scenario, "pt", "."),
           file = str_c(scenario, sys_time, sep="_")) %>% 
    mutate(
        r0_strain1 = str_extract(scenario_pt, "(?<=r)\\d*\\.?\\d*"),
        crossimm = str_extract(scenario_pt, "(?<=crossimm)\\d*\\.?\\d*"),
        seed = str_extract(scenario_pt, "(?<=seed)\\d*"),
        vacc = str_extract(scenario_pt, "(?<=vacc)\\d*\\.?\\d*"),
        vacc_start = str_extract(scenario_pt, "(?<=start_)\\d*\\.?\\d*"),
        preinf = str_extract(scenario_pt, "(?<=preinf)\\d*\\.?\\d*")) %>% 
    mutate(transmiss = case_when(preinf!="2.5" ~ str_extract(scenario_pt, "(?<=transmiss)\\d*\\.?\\d*"),
                              preinf=="2.5" ~ as.character(higher_r0/as.numeric(r0_strain1))))


# Combos
combos <- end %>% 
    select(r0_strain1, transmiss=higher_r0_rel, crossimm, seed, preinf=higher_r0_preinf, match_set=scenario) %>%
    mutate(vacc=1,
           vacc_start30 = 0.3,
           vacc_start60 = 0.6) %>% 
    pivot_longer(starts_with("vacc_start"), names_to = "name", values_to = "vacc_start") %>% 
    select(!name) %>% 
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric))
    
write.csv(combos, file.path(".", output_folder, "vacc_combos.csv"), row.names=F)
# Combos
# combos <- end %>% 
#     select(r0_strain1, transmiss=higher_r0_rel, crossimm, seed, vacc, vacc_start, preinf=higher_r0_preinf, match_set=scenario) %>% 
#     rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
#     mutate(across(!match_set, as.numeric))
# 
# write.csv(combos, file.path(".", output_folder, "more_transmiss_combos.csv"), row.names=F)

## Check missing
# r0_strain1_range <- c(1.5, 3)
# transmiss_range <- seq(1, 2, 0.25)
# crossimm_range <- c(0.5, 0.75, 0.95, 1)
# seed_range <- c(0, 35, 70)
# vacc_range <- c(0, 1)
# vacc_start_range <- c(0, 35)
# 
# combos <- crossing(r0_strain1_range, transmiss_range, crossimm_range, seed_range, vacc_range, vacc_start_range) %>% 
#     filter(!(vacc_range==0 & vacc_start_range==35)) %>% 
#     mutate(scenario = sprintf("r%s_transmiss%s_crossimm%s_seed%s_vacc%s_start_%s", r0_strain1_range, transmiss_range, crossimm_range, seed_range, vacc_range, vacc_start_range))

# missing <- combos$scenario[!(combos$scenario %in% end$scenario_pt)]
# nrows(missing)

### Remove incorrect files
# missing_pt <- str_replace_all(missing, "\\.", "pt")
# 
# folder_files <- list.files(file.path(".", output_folder))
# 
# files_remove <- folder_files[str_detect(folder_files, paste(missing_pt, collapse = '|'))]
# 
# files_remove <- str_c("./seir_model_output/", files_remove)
# 
# file.remove(files_remove)

# Plot epidemic size
# end %>% 
#     select(scenario, strain1_only, bothstrains, strain2_only) %>% 
#     pivot_longer(c(!scenario), names_to = "strain", values_to = "count") %>%
#     ggplot(aes(x=scenario, y=count, fill=strain)) + 
#     geom_bar(position="stack", stat="identity") +
#     # theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
#     scale_y_continuous(labels = scales::label_comma()) +
#     coord_flip()

for (i in unique(end$vacc)){
    
    graph_vacc <- 1
    graph_vacc_start <- 35
    graph_r0 <- 3
    
    end %>% 
        filter(vacc==graph_vacc, vacc_start==graph_vacc_start, r0_strain1==graph_r0) %>% 
        select(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, strain1_only, bothstrains, strain2_only) %>% 
        pivot_longer(c(strain1_only, bothstrains, strain2_only), names_to = "strain", values_to = "count") %>% 
        ggplot(aes(x=seed, y=count, fill=strain)) +
        geom_bar(position="stack", stat="identity") +
        facet_grid(cols = vars(transmiss),
                   rows = vars(crossimm),
                   labeller=label_both) +
        scale_y_continuous(labels = scales::label_comma()) +
        labs(caption=sprintf("Vaccination = %s, start day = %s\nStrain 1 R0 = %s", graph_vacc, graph_vacc_start, graph_r0)) +
        theme(plot.caption = element_text(hjust = 1))
    
    graph_name <- file.path(".", output_folder, sprintf("end_vacc%s_start%s_r%s.png", graph_vacc, graph_vacc_start, graph_r0))
    
    ggsave(graph_name, width = 25, height = 15, units = "cm", dpi = 300)
    
}

# Peak timing and peak new infectious cases
end <- end %>%
    mutate(file = str_c(scenario, "_", sys_time))

for (i in 310:length(end$file)) {
    filename <- sprintf("%s_matched.csv", end$file[i])
    match_growth <- read.csv(file.path(".", output_folder, filename))
    
    end$higher_r0[i] <- match_growth$higher_r0[match_growth$name=="r0"]
    
    end$shorter_preinf[i] <- match_growth$shorter_preinf[match_growth$name=="preinf_period"]
    
    end$shorter_preinf_gentime[i] <- match_growth$shorter_preinf[match_growth$name=="gen_time"]
    
    end$reduc_crossimm[i] <- match_growth$reduc_crossimm[match_growth$name=="crossimm"]
    
    end$reduc_vacceff[i] <- match_growth$reduc_vacceff[match_growth$name=="vacceff_i"]
    
    end$sus2[i] <- match_growth$higher_r0[match_growth$name=="sus2"]
    
    end$growth_rate2_intro[i] <- match_growth$higher_r0[match_growth$name=="growth_rate"]
    
}

write.csv(end, file.path(".", output_folder, "end_calc.csv"), row.names=F)

# Plot peak cases
graph_r0 <- 3

end %>% 
    filter(r0_strain1 == graph_r0) %>% 
    select(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, peak_new_I) %>%
    ggplot(aes(x=seed, y=peak_new_I, color=vacc, shape=vacc_start)) +
        geom_point(position=position_dodge(width=0.6)) +
    facet_grid(cols = vars(transmiss),
               rows = vars(crossimm),
               labeller=label_both) +
    scale_y_continuous(breaks=seq(0,5000000,500000),
                       limits=c(0,5000000),
                       labels = scales::label_comma()) +
    labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
    theme(plot.caption = element_text(hjust = 1))

graph_name <- file.path(".", output_folder, sprintf("end_peakinf_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 15, units = "cm", dpi = 300)

# Peak day
end %>% 
    filter(r0_strain1 == graph_r0) %>% 
    select(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, peak_time) %>%
    ggplot(aes(x=seed, y=peak_time, color=vacc, shape=vacc_start)) +
    geom_point(position=position_dodge(width=0.6)) +
    facet_grid(cols = vars(transmiss),
               rows = vars(crossimm),
               labeller=label_both) +
    scale_y_continuous(breaks=seq(0,500,100),
                       limits=c(0,450),
                       labels = scales::label_comma()) +
    labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
    theme(plot.caption = element_text(hjust = 1))

graph_name <- file.path(".", output_folder, sprintf("end_peaktime_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 15, units = "cm", dpi = 300)

# Epidemic size against susceptibles at intro
end %>% ggplot(aes(x=sus2, y=allinf, color=r0_strain1, shape=vacc)) +
    geom_point() +
    facet_grid(cols = vars(transmiss),
               rows = vars(crossimm),
               labeller=label_both) 

graph_r0 <- 1.5

end %>% filter(r0_strain1 == graph_r0) %>% 
    ggplot(aes(x=seed, y=shorter_preinf, color=vacc_start, shape=vacc)) +
    geom_point() +
    facet_grid(cols = vars(transmiss),
               rows = vars(crossimm),
               labeller=label_both) +
    geom_hline(yintercept=0)

