################ Specify factor with scale manual

library(tidyverse)
library(ggrepel)
library(cowplot)
# library(RColorBrewer)

theme_set(
    theme_gray(base_size = 12)
)

output_folder <- "seir_model_output"

end <- read.csv(file.path(".", output_folder, "end.csv"))

# Remove if running newest run_model_loop_vaccrate_preinf
population <- c(total_pop = 1000000)

# Extract parameters
end <- end %>%
    mutate(scenario_pt = str_replace_all(scenario, "pt", "."),
           file = str_c(scenario, sys_time, sep="_")) %>% 
    mutate(
        r0_strain1 = str_extract(scenario_pt, "(?<=r)\\d*\\.?\\d*"),
        transmiss = str_extract(scenario_pt, "(?<=transmiss)\\d*\\.?\\d*"),
        crossimm = str_extract(scenario_pt, "(?<=crossimm)\\d*\\.?\\d*"),
        seed = str_extract(scenario_pt, "(?<=seed)\\d*\\.?\\d*"),
        vacc = str_extract(scenario_pt, "(?<=vacc)\\d*\\.?\\d*"),
        vacc_start = str_extract(scenario_pt, "(?<=start_)\\d*\\.?\\d*"),
        preinf = str_extract(scenario_pt, "(?<=preinf)\\d*\\.?\\d*")) %>% 
    mutate(across(c(r0_strain1, crossimm, seed, vacc, vacc_start, preinf, transmiss), as.numeric)) %>% 
    mutate(match_set_pt = str_replace_all(match_set, "pt", ".")) %>% 
    mutate(match_preinf = str_extract(match_set_pt, "(?<=preinf)\\d*\\.?\\d*")) %>% 
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "More transmissable (Pre-infectious 0.5 days)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "More transmissable (Pre-infectious 1.5 days)",
                             preinf==0.5 ~ "Pre-infectious 0.5 days",
                             preinf==1.5 ~ "Pre-infectious 1.5 days",
                             preinf==2.5 & transmiss==1 ~ "Resident only")) %>% 
    mutate(label = fct_relevel(label,
                               c("Resident only",
                                 "Pre-infectious 0.5 days",
                                 "More transmissable (Pre-infectious 0.5 days)",
                                 "Pre-infectious 1.5 days",
                                 "More transmissable (Pre-infectious 1.5 days)"
                               ))) %>% 
    mutate(r0_strain2 = case_when(transmiss==1 ~ r0_strain1,
                                  transmiss > 1 ~ higher_r0),
           proppop_peak_new_I = peak_new_I/population[["total_pop"]])

legend_colour <- c("Resident only" = "#6A3D9A",
                   "Pre-infectious 0.5 days" = "#FB9A99",
                   "More transmissable (Pre-infectious 0.5 days)" = "#E31A1C",
                   "Pre-infectious 1.5 days" = "#A6CEE3",
                   "More transmissable (Pre-infectious 1.5 days)" = "#1F78B4")

# write.csv(end, file.path(".", output_folder, "end_mod.csv", row.names=F))

# Count match sets
end %>% count(match_set) # 99 unique match_sets
# [(3 R0 strain1 x 4 crossimm x 2 preinf x 4 seed x 3 vacc_start) x 2 match transmiss] + (3 R0 strain 1 x 3 vacc_start)
# (288x2) + 9 = 585 permutations

# Extract r0s -------
## Check if matching for more transmiss done correctly
end %>%
    filter(transmiss > 1) %>% 
    group_by(match_set_pt) %>% 
    summarise(r0_strain1 = n_distinct(r0_strain1),
              more_transmiss = n_distinct(higher_r0/r0_strain1),
              more_transmiss_label = n_distinct(transmiss),
              r0_strain2 = n_distinct(r0_strain2),
              match_preinf = n_distinct(as.numeric(match_preinf)),
              seed = n_distinct(seed),
              crossimm = n_distinct(crossimm)) %>% 
    summarise(across(!match_set_pt, ~ sum(.x!=1, na.rm = TRUE)))

matched_r0 <- end %>%
    filter(transmiss > 1) %>% 
    group_by(match_set_pt) %>% 
    summarise(r0_strain1 = mean(r0_strain1),
              more_transmiss = mean(higher_r0/r0_strain1),
              more_transmiss_label = mean(transmiss),
              r0_strain2 = mean(r0_strain2),
              match_preinf = mean(as.numeric(match_preinf)),
              seed = mean(seed),
              crossimm = mean(crossimm)) %>% 
    mutate(crossimm = factor(sprintf("%s%%", crossimm*100)),
           match_preinf = factor(sprintf("%s days", match_preinf)),
           seed = factor(round(seed))) %>% 
    mutate(crossimm = fct_shift(crossimm, 1))
    

# (matched_r0, file.path(".", output_folder, "matched_r0.csv"), row.names = F)

## Plot matched R0s
more_transmiss_r1pt5 <- matched_r0 %>% 
    filter(r0_strain1==1.5) %>%
    ggplot(aes(x=match_preinf, y=r0_strain2, group=crossimm, colour=seed, shape=crossimm)) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_hline(yintercept=1.5) +
    labs(title = bquote("Resident strain" ~ R[0] == 1.5),
         x = "",
         y = "",
         shape = "Cross-immunity",
         colour = "Variant seed day") +
    scale_colour_brewer(palette="YlOrRd") +
    theme_dark() +
    scale_shape(guide = 'none')

more_transmiss_r2pt5 <- matched_r0 %>% 
    filter(r0_strain1==2.5) %>%
    ggplot(aes(x=match_preinf, y=r0_strain2, group=crossimm, colour=seed, shape=crossimm)) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_hline(yintercept=2.5) +
    labs(title = bquote("Resident strain" ~ R[0] == 2.5),
         x = "",
         y = bquote("More transmissible variant" ~~ R[0]),
         shape = "Cross-immunity",
         colour = "Variant seed day") +
    scale_colour_brewer(palette="YlOrRd") +
    theme_dark() +
    scale_shape(guide = 'none')

more_transmiss_r4 <- matched_r0 %>% 
    filter(r0_strain1==4) %>%
    ggplot(aes(x=match_preinf, y=r0_strain2, group=crossimm, colour=seed, shape=crossimm)) +
    geom_point(position=position_dodge(width=0.5)) +
    geom_hline(yintercept=4) +
    labs(title = bquote("Resident strain" ~ R[0] == 4),
         x = "Pre-infectious period matched to",
         y = "",
         shape = "Cross-immunity",
         colour = "Variant seed day") +
    scale_colour_brewer(palette="YlOrRd") +
    theme_dark() 

more_transmiss_graph <- plot_grid(
    more_transmiss_r1pt5, more_transmiss_r2pt5, more_transmiss_r4,
    ncol=1,
    align = "v"
    )

ggsave(file.path(".", output_folder, "more_transmiss_r0s.png"), width = 20, height = 20, units = "cm", dpi = 300)

# Check which ended prematurely ----
run_longer <- end %>% filter(I > 5) ## Removed and reran in other script

file.copy(file.path(".", output_folder, "end.csv"), file.path(".", output_folder, "end_copy.csv")) # Save copy of end.csv

end <- end %>% filter(I <= 5) %>% select(scenario:match_set) # Remove scenarios which need rerunning from end and drop mutated columns

write.csv(end, file.path(".", output_folder, "end.csv", row.names=F))

## Create combos for those which need rerunning
combos <- run_longer %>% 
    select(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, preinf, match_set, higher_r0) %>%
    mutate(transmiss = case_when(transmiss==1 ~ transmiss,
                                 transmiss > 1 ~ higher_r0/r0_strain1)) %>% 
    select(!higher_r0) %>% 
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric))

write.csv(combos, file.path(".", output_folder, "run_longer.csv"), row.names = F)


# Create combos with matched growth rate, higher R0 ----
combos <- end %>%
    filter(seed != 1200) %>% 
    select(r0_strain1, transmiss=higher_r0_rel, crossimm, seed, vacc, vacc_start, preinf=higher_r0_preinf, match_set=scenario) %>%
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric))

write.csv(combos, file.path(".", output_folder, "more_transmiss_combos.csv"), row.names=F)

# Combos for vaccination ----
combos_preinf_vacc <- end %>%
    filter(seed != 1200, transmiss==1) %>%
    select(r0_strain1, transmiss, crossimm, seed, preinf, match_set) %>%
    mutate(vacc=1,
           vacc_start1 = case_when(r0_strain1==1.5 ~ 0.1111,
                                   r0_strain1==2.5 ~ 0.2,
                                   r0_strain1==4 ~ 0.25),
           vacc_start2 = case_when(r0_strain1==1.5 ~ 0.3333,
                                    r0_strain1==2.5 ~ 0.45,
                                    r0_strain1==4 ~ 0.5625)) %>%
    pivot_longer(starts_with("vacc_start"), names_to = "name", values_to = "vacc_start") %>%
    select(!name) %>%
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric)) %>% 
    relocate(match_set, .after = last_col())

combos_transmiss_vacc <- end %>%
    filter(seed != 1200, transmiss>1) %>%
    select(r0_strain1, transmiss, higher_r0, crossimm, seed, preinf, match_set) %>%
    mutate(transmiss = higher_r0/r0_strain1) %>% 
    mutate(vacc=1,
           vacc_start1 = case_when(r0_strain1==1.5 ~ 0.1111,
                                   r0_strain1==2.5 ~ 0.2,
                                   r0_strain1==4 ~ 0.25),
           vacc_start2 = case_when(r0_strain1==1.5 ~ 0.3333,
                                   r0_strain1==2.5 ~ 0.45,
                                   r0_strain1==4 ~ 0.5625)) %>%
    pivot_longer(starts_with("vacc_start"), names_to = "name", values_to = "vacc_start") %>%
    select(!c(name, higher_r0)) %>%
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric)) %>% 
    relocate(match_set, .after = last_col())

combos <- rbind(combos_preinf_vacc, combos_transmiss_vacc)

write.csv(combos, file.path(".", output_folder, "vacc_combos.csv"), row.names=F)



## Check missing ----
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

end %>%
    filter(r0_strain1 == 4, vacc_start==0, seed!=1200) %>% 
    ggplot(aes(x=scenario, y=allinf, fill=preinf)) +
        geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### Total size ------
graph_r0 <- 4
graph_vaccstart_range <- sort(unique(end$vacc_start[end$r0_strain1==graph_r0]))

allinf_strain1only_vacc0 <- end$allinf[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[1]]
allinf_strain1only_vacc1 <- end$allinf[end$r0_strain1== graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[2]]
allinf_strain1only_vacc2 <- end$allinf[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[3]]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_allinf = case_when(vacc_start==graph_vaccstart_range[1] ~ allinf/allinf_strain1only_vacc0,
                                      vacc_start==graph_vaccstart_range[2] ~ allinf/allinf_strain1only_vacc1,
                                      vacc_start==graph_vaccstart_range[3] ~ allinf/allinf_strain1only_vacc2)) %>% 
    mutate(vacc_start = sprintf("Vaccine coverage\n%s%%", vacc_start*100),
           crossimm = sprintf("Cross-immunity\n%s%%", crossimm*100),
           seed = round(seed)) %>% 
    mutate(seed = factor(replace(seed, seed==1200, "NA"))) %>%
    mutate(crossimm = fct_relevel(crossimm,
                                  c("Cross-immunity\n25%",
                                    "Cross-immunity\n50%",
                                    "Cross-immunity\n75%",
                                    "Cross-immunity\n100%")),
           seed = fct_relevel(seed, "NA"))

## Relative, final number infected
graph_table %>% 
    filter(seed!="NA") %>% 
    ggplot(aes(x=seed, y=rel_allinf, colour=label, group=match_set)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point(position=position_dodge(width=0.5)) +
    facet_grid(rows=vars(crossimm),
               cols=vars(vacc_start)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Novel variant seed time (days)",
         y = "Final number infected, relative to resident only",
         colour = "Variant characteristic") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour)

graph_name <- file.path(".", output_folder, sprintf("rel_allinf_r%s.png", graph_r0))
# graph_name <- file.path(".", output_folder, sprintf("rel_allinf_r%s_log.png", graph_r0))

ggsave(graph_name, width = 30, height = 16, units = "cm", dpi = 300)

## Final number infected as percentage of population
graph_table %>%
    ggplot(aes(x=seed, y=proppop_allinf, colour=label, group=match_set)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
        x = "Novel variant seed time (days)",
        y = "Final number infected (% population)",
        colour = "Variant characteristic") +
        scale_colour_manual(values=legend_colour) +
    scale_y_continuous(labels = scales::label_percent(accuracy=1),
                       limits = c(0, NA))
    
graph_name <- file.path(".", output_folder, sprintf("allinf_r%s.png", graph_r0))

ggsave(graph_name, width = 30, height = 16, units = "cm", dpi = 300)



# Peak new cases -----------------------------
graph_r0 <- 4
graph_vaccstart_range <- sort(unique(end$vacc_start[end$r0_strain1==graph_r0]))

peakinf_strain1only_vacc0 <- end$peak_new_I[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[1]]
peakinf_strain1only_vacc1 <- end$peak_new_I[end$r0_strain1== graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[2]]
peakinf_strain1only_vacc2 <- end$peak_new_I[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[3]]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_peak_new_I = case_when(vacc_start==graph_vaccstart_range[1] ~ peak_new_I/peakinf_strain1only_vacc0,
                                  vacc_start==graph_vaccstart_range[2] ~ peak_new_I/peakinf_strain1only_vacc1,
                                  vacc_start==graph_vaccstart_range[3] ~ peak_new_I/peakinf_strain1only_vacc2)) %>% 
    mutate(vacc_start = sprintf("Vaccine coverage\n%s%%", vacc_start*100),
           crossimm = sprintf("Cross-immunity\n%s%%", crossimm*100),
           seed = round(seed)) %>% 
    mutate(seed = factor(replace(seed, seed==1200, "NA"))) %>%
    mutate(crossimm = fct_relevel(crossimm,
                                  c("Cross-immunity\n25%",
                                    "Cross-immunity\n50%",
                                    "Cross-immunity\n75%",
                                    "Cross-immunity\n100%")),
           seed = fct_relevel(seed, "NA"))

# Relative, peak cases
graph_table %>% 
    filter(seed!="NA") %>% 
    ggplot(aes(x=seed, y=rel_peak_new_I, colour=label, group=match_set)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point(position=position_dodge(width=0.5)) +
    facet_grid(rows=vars(crossimm),
               cols=vars(vacc_start)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Novel variant seed time (days)",
         y = "Peak cases, relative to resident only",
         colour = "Variant characteristic") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour)

graph_name <- file.path(".", output_folder, sprintf("rel_peak_new_I_r%s.png", graph_r0))

ggsave(graph_name, width = 30, height = 16, units = "cm", dpi = 300)

# Peak cases
graph_table %>% 
    ggplot(aes(x=seed, y=proppop_peak_new_I, colour=label, group=match_set)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Novel variant seed time (days)",
         y = "Peak daily cases (% population)",
         colour = "Variant characteristic") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour) +
    scale_y_continuous(labels = scales::label_percent(accuracy=0.1),
                       limits=c(0,NA))
    
graph_name <- file.path(".", output_folder, sprintf("peak_new_I_r%s.png", graph_r0))

ggsave(graph_name, width = 30, height = 16, units = "cm", dpi = 300)

# Peak time -----
graph_r0 <- 4
graph_vaccstart_range <- sort(unique(end$vacc_start[end$r0_strain1==graph_r0]))

peaktime_strain1only_vacc0 <- end$peak_time[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[1]]
peaktime_strain1only_vacc1 <- end$peak_time[end$r0_strain1== graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[2]]
peaktime_strain1only_vacc2 <- end$peak_time[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == graph_vaccstart_range[3]]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_peaktime = case_when(vacc_start==graph_vaccstart_range[1] ~ peak_time-peaktime_strain1only_vacc0,
                                   vacc_start==graph_vaccstart_range[2] ~ peak_time-peaktime_strain1only_vacc1,
                                   vacc_start==graph_vaccstart_range[3] ~ peak_time-peaktime_strain1only_vacc2)) %>% 
    mutate(vacc_start = sprintf("Vaccine coverage\n%s%%", vacc_start*100),
           crossimm = sprintf("Cross-immunity\n%s%%", crossimm*100),
           seed = round(seed)) %>% 
    mutate(seed = factor(replace(seed, seed==1200, "NA"))) %>%
    mutate(crossimm = fct_relevel(crossimm,
                                  c("Cross-immunity\n25%",
                                    "Cross-immunity\n50%",
                                    "Cross-immunity\n75%",
                                    "Cross-immunity\n100%")),
           seed = fct_relevel(seed, "NA"))

# Relative peak time
graph_table %>% 
    filter(seed!="NA") %>% 
    ggplot(aes(x=seed, y=rel_peaktime, colour=label, group=match_set)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point(position=position_dodge(width=0.5)) +
    facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Novel variant seed time (days)",
         y = "Peak day, relative to resident only",
         colour = "Variant characteristic") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour) 

graph_name <- file.path(".", output_folder, sprintf("rel_peaktime_r%s.png", graph_r0))

ggsave(graph_name, width = 30, height = 16, units = "cm", dpi = 300)

# Peak time
graph_table %>% 
    ggplot(aes(x=seed, y=peak_time, colour=label, group=match_set)) +
    geom_point(position=position_dodge(width=0.5)) +
    facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Novel variant seed time (days)",
         y = "Peak day",
         colour = "Variant characteristic") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour) +
    scale_y_continuous(labels = scales::label_comma())

graph_name <- file.path(".", output_folder, sprintf("peaktime_r%s.png", graph_r0))

ggsave(graph_name, width = 30, height = 16, units = "cm", dpi = 300)


###### Cumulative --------
graph_r0 <- 4

graph_table <- end %>% 
    filter(r0_strain1==graph_r0) %>% 
    select(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, preinf, match_preinf, label, proppop_strain1_only, proppop_bothstrains, proppop_strain2_only) %>%
    mutate(crossimm_vacc = factor(sprintf("Cross-immunity %s%%\nVaccine coverage %s%%", crossimm*100, vacc_start*100)),
           seed_label = factor(case_when(seed==1200 ~ "NA", 
                               seed!=1200 ~ sprintf("Seed day %s", round(seed))))) %>%
    mutate(crossimm_vacc = fct_shift(crossimm_vacc, 3)) %>% 
    pivot_longer(c(proppop_strain1_only, proppop_bothstrains, proppop_strain2_only), names_to = "strain", values_to = "count") %>% 
    mutate(strain = fct_recode(strain,
                               "Both strains" = "proppop_bothstrains",
                               "Resident only" = "proppop_strain1_only",
                               "Variant only" = "proppop_strain2_only"),
           label = fct_recode(label,
                              "More transmissable\n(Pre-infectious 0.5 days)" = "More transmissable (Pre-infectious 0.5 days)",
                              "More transmissable\n(Pre-infectious 1.5 days)" = "More transmissable (Pre-infectious 1.5 days)"))

ggplot(graph_table, aes(x=label, y=count, fill=strain)) +
    geom_bar(position="stack", stat="identity") +
    facet_grid(cols = vars(seed_label),
               rows = vars(crossimm_vacc),
               scale = "free_x", 
               space = "free_x") +
    scale_y_continuous(labels = scales::label_percent(accuracy=1),
                       breaks = scales::breaks_extended()) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Variant characteristic",
         y = "% population",
         fill = "Infecting strain") +    
    theme(plot.caption = element_text(hjust = 1),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          strip.text.y.right = element_text(angle = 0))

graph_name <- file.path(".", output_folder, sprintf("cumulative_r%s.png", graph_r0))

ggsave(graph_name, width = 24, height = 24, units = "cm", dpi = 300)
    
#### Plotting new infections ----
# graph_r0 <- 4
# graph_seed <- 0
# 
# plot_scenarios <- end %>% filter(r0_strain1==graph_r0, seed==graph_seed)
# 
# plot_scenarios <- rbind(plot_scenarios, end[end$seed==1200 & end$r0_strain1==graph_r0,])
# 
# daily_inf_table <- data.frame(scenario=numeric(), day=numeric(), daily_new_I1=numeric(), daily_new_I2=numeric(), daily_new_I=numeric(), prop_new_I1=numeric(), prop_new_I2=numeric(), r0_strain1=numeric(), crossimm=numeric(), seed=numeric(), vacc=numeric(), vacc_start=numeric(), preinf=numeric(), transmiss=numeric(), match_preinf=numeric())
# 
# for (i in 1:length(plot_scenarios$file)) {
#     filename <- sprintf("%s_dailyinf.csv", plot_scenarios$file[i])
# 
#     table <- read.csv(file.path(".", output_folder, filename)) %>% 
#         mutate(scenario=plot_scenarios$scenario[i],
#                .before=day) %>% 
#         mutate(r0_strain1=plot_scenarios$r0_strain1[i], 
#                crossimm=plot_scenarios$crossimm[i], 
#                seed=plot_scenarios$seed[i], 
#                vacc=plot_scenarios$vacc[i], 
#                vacc_start=plot_scenarios$vacc_start[i], 
#                preinf=plot_scenarios$preinf[i], 
#                transmiss=plot_scenarios$transmiss[i], 
#                match_preinf=plot_scenarios$match_preinf[i])
#     
#     daily_inf_table <- rbind(daily_inf_table, table)
# }
# 
# daily_inf_table %>% 
#     ggplot(aes(x=day, y=daily_new_I, colour=label)) +
#     geom_line(size=0.8) +
#     facet_grid(cols=vars(crossimm),
#                rows=vars(vacc_start)) +
#     scale_colour_manual(values=c(
#         "#1F78B4",
#         "#E31A1C",
#         "#A6CEE3",
#         "#FB9A99",
#         "#6A3D9A"
#     )) +
#     scale_x_continuous(limits=c(0,300)) +
#     geom_vline(xintercept = graph_seed, linetype="dashed", size=0.5) +
#     scale_y_continuous(labels = scales::label_comma())

#### Plotting new infections, all seeds----
graph_r0 <- 4
graph_limit <- 250

plot_scenarios <- end %>% filter(r0_strain1==graph_r0)

# plot_scenarios <- rbind(plot_scenarios, end[end$seed==1200 & end$r0_strain1==graph_r0,])

daily_inf_table <- data.frame(scenario=numeric(), day=numeric(), daily_new_I1=numeric(), daily_new_I2=numeric(), daily_new_I=numeric(), prop_new_I1=numeric(), prop_new_I2=numeric(), r0_strain1=numeric(), crossimm=numeric(), seed=numeric(), vacc=numeric(), vacc_start=numeric(), preinf=numeric(), transmiss=numeric(), match_preinf=numeric())

for (i in 1:length(plot_scenarios$file)) {
    filename <- sprintf("%s_dailyinf.csv", plot_scenarios$file[i])
    
    table <- read.csv(file.path(".", output_folder, filename)) %>% 
        mutate(scenario=plot_scenarios$scenario[i],
               .before=day) %>% 
        mutate(r0_strain1=plot_scenarios$r0_strain1[i], 
               crossimm=plot_scenarios$crossimm[i], 
               seed=plot_scenarios$seed[i], 
               vacc=plot_scenarios$vacc[i], 
               vacc_start=plot_scenarios$vacc_start[i], 
               preinf=plot_scenarios$preinf[i], 
               transmiss=plot_scenarios$transmiss[i], 
               match_preinf=plot_scenarios$match_preinf[i])
    
    daily_inf_table <- rbind(daily_inf_table, table)
}

graph_table <- daily_inf_table %>% 
    mutate(seed = if_else(seed==1200, 0, round(seed))) %>% 
    mutate(crossimm_seed = sprintf("Seed day %s\nCross-immunity %s%%", seed, crossimm*100)) %>% 
    mutate(crossimm_seed = fct_relevel(crossimm_seed, unique(str_sort(crossimm_seed, numeric = T)))) %>% 
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "More transmissable (Pre-infectious 0.5 days)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "More transmissable (Pre-infectious 1.5 days)",
                             preinf==0.5 ~ "Pre-infectious 0.5 days",
                             preinf==1.5 ~ "Pre-infectious 1.5 days",
                             preinf==2.5 & transmiss==1 ~ "Resident only")) %>% 
    mutate(label = fct_relevel(label,
                               c("Resident only",
                                 "Pre-infectious 0.5 days",
                                 "More transmissable (Pre-infectious 0.5 days)",
                                 "Pre-infectious 1.5 days",
                                 "More transmissable (Pre-infectious 1.5 days)"
                               ))) %>% 
    mutate(vacc_label = factor(sprintf("Vaccine coverage %s%%", vacc_start*100)))


seed_line <- graph_table %>% 
    group_by(crossimm_seed) %>% 
    summarise(seed_line = mean(seed))
    
ggplot(graph_table, aes(x=day, y=proppop_new_I, colour=label)) +
    geom_line() +
    facet_grid(rows=vars(crossimm_seed),
               cols=vars(vacc_label)) +
    scale_colour_manual(values=legend_colour) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_continuous(labels = scales::label_percent(accuracy=0.1)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day",
         y = "Daily new cases (% population)",
         colour = "Variant characteristic") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0))
    

graph_name <- file.path(".", output_folder, sprintf("newinf_r%s.png", graph_r0))

ggsave(graph_name, width = 28, height = 25, units = "cm", dpi = 300)

#### Plotting susceptibles and R1, all seeds----
graph_r0 <- 2.5
graph_limit <- 500

plot_scenarios <- end %>% filter(r0_strain1==graph_r0)

# plot_scenarios <- rbind(plot_scenarios, end[end$seed==1200 & end$r0_strain1==graph_r0,])

S_R1 <- data.frame(scenario=numeric(), time=numeric(), proppop_S=numeric(), proppop_R1=numeric(), proppop_V=numeric(), r0_strain1=numeric(), crossimm=numeric(), seed=numeric(), vacc=numeric(), vacc_start=numeric(), preinf=numeric(), transmiss=numeric(), match_preinf=numeric())

for (i in 1:length(plot_scenarios$file)) {
    filename <- sprintf("%s_calc.csv", plot_scenarios$file[i])
    
    table <- read.csv(file.path(".", output_folder, filename)) %>%
        select(time, proppop_S, proppop_R1, proppop_V) %>% 
        filter(time %in% seq(0, max(time), 1)) %>% 
        mutate(scenario=plot_scenarios$scenario[i],
               .before=time) %>% 
        mutate(r0_strain1=plot_scenarios$r0_strain1[i], 
               crossimm=plot_scenarios$crossimm[i], 
               seed=plot_scenarios$seed[i], 
               vacc=plot_scenarios$vacc[i], 
               vacc_start=plot_scenarios$vacc_start[i], 
               preinf=plot_scenarios$preinf[i], 
               transmiss=plot_scenarios$transmiss[i], 
               match_preinf=plot_scenarios$match_preinf[i])
    
    S_R1 <- rbind(S_R1, table)
}

write.csv(S_R1, file.path(".", output_folder, sprintf("S_R1_r%s.csv", graph_r0)), row.names=F)

graph_table <- S_R1 %>% 
    mutate(seed = if_else(seed==1200, 0, round(seed))) %>% 
    mutate(crossimm_seed = sprintf("Seed day %s\nCross-immunity %s%%", seed, crossimm*100)) %>% 
    mutate(crossimm_seed = fct_relevel(crossimm_seed, unique(str_sort(crossimm_seed, numeric = T)))) %>% 
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "More transmissable (Pre-infectious 0.5 days)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "More transmissable (Pre-infectious 1.5 days)",
                             preinf==0.5 ~ "Pre-infectious 0.5 days",
                             preinf==1.5 ~ "Pre-infectious 1.5 days",
                             preinf==2.5 & transmiss==1 ~ "Resident only")) %>% 
    mutate(label = fct_relevel(label,
                               c("Resident only",
                                 "Pre-infectious 0.5 days",
                                 "More transmissable (Pre-infectious 0.5 days)",
                                 "Pre-infectious 1.5 days",
                                 "More transmissable (Pre-infectious 1.5 days)"
                               ))) %>% 
    mutate(vacc_label = factor(sprintf("Vaccine coverage %s%%", vacc_start*100))) %>%
    pivot_longer(starts_with("proppop"), names_to="compartment", values_to="proportion") %>% 
    mutate(compartment = fct_recode(compartment,
                                    "Recovered from resident strain" = "proppop_R1",
                                    "Susceptible" = "proppop_S",
                                    "Vaccinated and uninfected" = "proppop_V"))


seed_line <- graph_table %>% 
    group_by(crossimm_seed) %>% 
    summarise(seed_line = mean(seed))

graph_table %>% 
    filter(compartment!="Vaccinated and uninfected") %>% 
    ggplot(aes(x=time, y=proportion, colour=label, linetype=compartment)) +
    geom_line() +
    facet_grid(rows=vars(crossimm_seed),
               cols=vars(vacc_label)) +
    scale_colour_manual(values=legend_colour) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_continuous(labels = scales::label_percent(accuracy=0.1)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day",
         y = "% population",
         colour = "Variant characteristic",
         linetype = "Compartment") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0)) 


graph_name <- file.path(".", output_folder, sprintf("S_R1_r%s.png", graph_r0))

ggsave(graph_name, width = 28, height = 28, units = "cm", dpi = 300)

###  -----
graph_name <- file.path(".", output_folder, sprintf("end_vacc%s_start%s_r%s.png", graph_vacc, graph_vacc_start, graph_r0))

ggsave(graph_name, width = 25, height = 15, units = "cm", dpi = 300)

for (i in unique(end$vacc)){
    
    graph_vacc <- 0
    graph_vacc_start <- 0
    graph_r0 <- 1.5
    
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

