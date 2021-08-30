################ Specify factor with scale manual

library(tidyverse)
library(ggrepel)
library(cowplot)
library(scales)
library(ggtext)
library(Rcpp)
library(patchwork)
# library(RColorBrewer)

theme_set(
    theme_gray(base_size = 11)
)

output_folder <- "seir_model_output"

end <- read.csv(file.path(".", output_folder, "end.csv"))

resident_seedproxy <- 1500

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
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "More transmissable (r[init] matched to d[E] 0.5)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "More transmissable (r[init] matched to d[E] 1.5)",
                             preinf==0.5 ~ "Pre-infectious (d[E]) 0.5 days",
                             preinf==1.5 ~ "Pre-infectious (d[E]) 1.5 days",
                             preinf==2.5 & transmiss==1 ~ "Resident strain only")) %>% 
    mutate(label = fct_relevel(label,
                               c("Resident strain only",
                                 "Pre-infectious (d[E]) 0.5 days",
                                 "More transmissable (r[init] matched to d[E] 0.5)",
                                 "Pre-infectious (d[E]) 1.5 days",
                                 "More transmissable (r[init] matched to d[E] 1.5)"
                               ))) %>% 
    mutate(r0_strain2 = case_when(transmiss==1 ~ r0_strain1,
                                  transmiss > 1 ~ higher_r0)) %>% 
    mutate(vacc_cov = factor(sprintf("Vaccine coverage\n%s%%", vacc_start*100)),
           crossimm_label = factor(sprintf("Cross-immunity\n%s%%", crossimm*100)),
           crossimm_perc = factor(sprintf("%s%%", crossimm*100)),
           seed_round = round(seed)) %>% 
    mutate(seed_round = factor(replace(seed_round, seed_round==resident_seedproxy, "NA"))) %>%
    mutate(seed_round = fct_relevel(seed_round, unique(str_sort(seed_round, numeric = T)))) %>% 
    mutate(crossimm_label = fct_relevel(crossimm_label, unique(str_sort(crossimm_label, numeric = T))),
           crossimm_perc = fct_shift(crossimm_perc, 1),
           seed_round = fct_relevel(seed_round, "NA"))

legend_colour <- c("Resident strain only" = "#008000",
                   "Pre-infectious (d[E]) 0.5 days" = "#ff6262",
                   "More transmissable (r[init] matched to d[E] 0.5)" = "#cc0000",
                   "Pre-infectious (d[E]) 1.5 days" = "#7676ff",
                   "More transmissable (r[init] matched to d[E] 1.5)" = "#0000b3")

variant_label <- c("Resident strain only" = "Resident strain only",
                   "Pre-infectious (d[E]) 0.5 days" = bquote("Pre-infectious" ~ (d[E]) ~ "0.5 days"),
                   "More transmissable (r[init] matched to d[E] 0.5)" = bquote("More transmissable" ~ (r[init] ~ matched ~ to ~ d[E] ~ 0.5)),
                   "Pre-infectious (d[E]) 1.5 days" = bquote("Pre-infectious" ~ (d[E]) ~ "0.5 days"),
                   "More transmissable (r[init] matched to d[E] 1.5)" = bquote("More transmissable" ~ (r[init] ~ matched ~ to ~ d[E] ~ 0.5)))

legend_shape <- c("25%" = 16,
                  "75%" = 17,
                  "100%" = 15)

legend_linetype <- c("25%" = "solid",
                     "75%" = "longdash",
                     "100%" = "dotted")


# write.csv(end, file.path(".", output_folder, "end_mod.csv", row.names=F))

# Count match sets
end %>% count(match_set) %>% count(n) # 75 unique match_sets, [(3 R0 strain1 x 3 crossimm x 2 preinf x 4 seed)] + 3 R0 strain 1
# [(3 R0 strain1 x 3 crossimm x 2 preinf x 4 seed x 3 vacc_start) x 2 match transmiss] + (3 R0 strain 1 x 3 vacc_start)
# (216x2) + 9 = 441 permutations


# Create combos with matched growth rate, higher R0 ----
combos <- end %>%
    filter(seed != resident_seedproxy) %>% 
    select(r0_strain1, transmiss=higher_r0_rel, crossimm, seed, vacc, vacc_start, preinf=higher_r0_preinf, match_set=scenario) %>%
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric))

write.csv(combos, file.path(".", output_folder, "more_transmiss_combos.csv"), row.names=F)

# Combos for vaccination ----
combos_preinf_vacc <- end %>%
    filter(seed != resident_seedproxy, transmiss==1) %>%
    select(r0_strain1, transmiss, crossimm, seed, preinf, match_set) %>%
    mutate(vacc=1,
           vacc_start1 = case_when(r0_strain1==1.5 ~ 0.1667,
                                   r0_strain1==2.5 ~ 0.3,
                                   r0_strain1==4 ~ 0.375),
           vacc_start2 = case_when(r0_strain1==1.5 ~ 0.3334,
                                    r0_strain1==2.5 ~ 0.6,
                                    r0_strain1==4 ~ 0.75)) %>%
    pivot_longer(starts_with("vacc_start"), names_to = "name", values_to = "vacc_start") %>%
    select(!name) %>%
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric)) %>% 
    relocate(match_set, .after = last_col())

combos_transmiss_vacc <- end %>%
    filter(seed != resident_seedproxy, transmiss>1) %>%
    select(r0_strain1, transmiss, higher_r0, crossimm, seed, preinf, match_set) %>%
    mutate(transmiss = higher_r0/r0_strain1) %>% 
    mutate(vacc=1,
           vacc_start1 = case_when(r0_strain1==1.5 ~ 0.1667,
                                   r0_strain1==2.5 ~ 0.3,
                                   r0_strain1==4 ~ 0.375),
           vacc_start2 = case_when(r0_strain1==1.5 ~ 0.3334,
                                   r0_strain1==2.5 ~ 0.6,
                                   r0_strain1==4 ~ 0.75)) %>%
    pivot_longer(starts_with("vacc_start"), names_to = "name", values_to = "vacc_start") %>%
    select(!c(name, higher_r0)) %>%
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric)) %>% 
    relocate(match_set, .after = last_col())

combos <- rbind(combos_preinf_vacc, combos_transmiss_vacc)

combos %>% count(match_set) %>% count(n)

write.csv(combos, file.path(".", output_folder, "vacc_combos.csv"), row.names=F)

# Check which ended prematurely ----
run_longer <- end %>% filter(I > 2) ## Removed and reran in other script

file.copy(file.path(".", output_folder, "end.csv"), file.path(".", output_folder, "end_copy.csv")) # Save copy of end.csv

end <- end %>% filter(I <= 2) %>% select(scenario:match_set) # Remove scenarios which need rerunning from end and drop mutated columns

write.csv(end, file.path(".", output_folder, "end.csv"), row.names=F)

## Create combos for those which need rerunning
combos <- run_longer %>% 
    select(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, preinf, match_set, higher_r0) %>%
    mutate(transmiss = case_when(transmiss==1 ~ transmiss,
                                 transmiss > 1 ~ higher_r0/r0_strain1)) %>% 
    select(!higher_r0) %>% 
    rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
    mutate(across(!match_set, as.numeric))

write.csv(combos, file.path(".", output_folder, "run_longer.csv"), row.names = F)


# Extract more transmissible variant r0s -------
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


write.csv(matched_r0, file.path(".", output_folder, "matched_r0.csv"), row.names = F)

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
    scale_shape(guide = 'none') +
    scale_y_continuous(sec.axis = sec_axis(~ . / 1.5))

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
    scale_shape(guide = 'none') +
    scale_y_continuous(sec.axis = sec_axis(~ . / 2.5,
                                           name="Relative transmissibility"))


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
    theme_dark() +
    scale_y_continuous(sec.axis = sec_axis(~ . / 4))

plot_grid(
    more_transmiss_r1pt5, more_transmiss_r2pt5, more_transmiss_r4,
    ncol=1,
    align = "v"
)

ggsave(file.path(".", output_folder, "more_transmiss_r0s.png"), width = 20, height = 20, units = "cm", dpi = 300)


#### Total size ------
graph_r0 <- 2.5
graph_vaccstart_range <- sort(unique(end$vacc_start[end$r0_strain1==graph_r0]))

allinf_strain1only_vacc0 <- end$allinf[end$r0_strain1==graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[1]]
allinf_strain1only_vacc1 <- end$allinf[end$r0_strain1== graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[2]]
allinf_strain1only_vacc2 <- end$allinf[end$r0_strain1==graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[3]]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_allinf = case_when(vacc_start==graph_vaccstart_range[1] ~ allinf/allinf_strain1only_vacc0,
                                      vacc_start==graph_vaccstart_range[2] ~ allinf/allinf_strain1only_vacc1,
                                      vacc_start==graph_vaccstart_range[3] ~ allinf/allinf_strain1only_vacc2)) %>% 
    mutate(seed_round = fct_drop(seed_round)) %>% 
    mutate(seed_dodge = factor(seed_round, labels=seq(1,length(levels(seed_round))))) %>% 
    mutate(seed_dodge = case_when(crossimm_perc==levels(crossimm_perc)[1] ~ as.numeric(seed_dodge)-0.2,
                                  crossimm_perc==levels(crossimm_perc)[2] | seed==resident_seedproxy ~ as.numeric(seed_dodge),
                                  crossimm_perc==levels(crossimm_perc)[3] ~ as.numeric(seed_dodge)+0.2),
           line_group = str_c(crossimm_perc, label))


# Relative, final number infected (collapsed)
graph_table %>% 
    filter(seed!=resident_seedproxy) %>% 
    ggplot(aes(x=seed_round, y=rel_allinf, colour=label, group=crossimm_perc, shape = crossimm_perc)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final number population infected with at least one strain,\nrelative to resident strain only",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape)

graph_table %>% 
    filter(seed!=resident_seedproxy) %>%
    ggplot(aes(x=seed_dodge, y=rel_allinf, colour=label, group=crossimm_perc, shape = crossimm_perc)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +    
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final number infected with at least one strain,\nrelative to resident strain only",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("rel_allinf_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

### Final number infected as percentage of population, collapsed

graph_table %>%
    ggplot(aes(x=seed_round, y=proppop_allinf, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final % population infected with at least one strain",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent(),
                       limits = c(0, NA)) +
    scale_shape_manual(values=legend_shape)

graph_table %>%
    ggplot(aes(x=seed_dodge, y=proppop_allinf, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final % population infected with at least one strain",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent(),
                       limits = c(0, NA)) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)


graph_name <- file.path(".", output_folder, sprintf("allinf_ collapsed_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Cumulative new infections, relative to population size, collapsed
graph_table %>%
    ggplot(aes(x=seed_round, y=proppop_Cum_Inf, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Overall infections, as % population size",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent(),
                       limits = c(0,NA)
                       # limits = c(0, round(max(graph_table$proppop_Cum_Inf) + 0.1, 1))
                       ) +
    scale_shape_manual(values=legend_shape)

graph_table %>%
    ggplot(aes(x=seed_dodge, y=proppop_Cum_Inf, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Overall infections, as % population size",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent(),
                       limits = c(0,NA)
                       # limits = c(0, round(max(graph_table$proppop_Cum_Inf) + 0.1, 1))
                       ) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("cuminf_ collapsed_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

### Reinfections, collapsed
graph_table %>%
    ggplot(aes(x=seed_round, y=proppop_bothstrains, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final % population re-infected",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent()) +
    scale_shape_manual(values=legend_shape)

graph_table %>%
    ggplot(aes(x=seed_dodge, y=proppop_bothstrains, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final % population re-infected",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent()) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("strain2only_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Second strain only
graph_table %>%
    ggplot(aes(x=seed_dodge, y=proppop_strain2_only, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final % population infected only by novel variant",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent()) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("strain2_only_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Second strain - Either reinfection or second strain only
graph_table %>%
    mutate(proppop_strain2 = proppop_strain2_only + proppop_bothstrains) %>% 
    ggplot(aes(x=seed_dodge, y=proppop_strain2, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final % population infected at least once by novel variant\n(primary infections and reinfections)",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent()) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("proppopstrain2_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Resident strain only
graph_table %>%
    ggplot(aes(x=seed_dodge, y=proppop_strain1_only, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Final % population infected only by resident strain",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent()) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("strain1_only_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Proportion infections which are vaccine breakthrough
graph_table %>%
    filter(vacc!=0) %>% 
    ggplot(aes(x=seed_dodge, y=prop_Cum_Infv, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "% infections among vaccinated",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent()) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("propinf_ vacc_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

# Peak new cases -----------------------------
graph_r0 <- 2.5
graph_vaccstart_range <- sort(unique(end$vacc_start[end$r0_strain1==graph_r0]))

peakinf_strain1only_vacc0 <- end$peak_new_I[end$r0_strain1==graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[1]]
peakinf_strain1only_vacc1 <- end$peak_new_I[end$r0_strain1== graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[2]]
peakinf_strain1only_vacc2 <- end$peak_new_I[end$r0_strain1==graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[3]]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_peak_new_I = case_when(vacc_start==graph_vaccstart_range[1] ~ peak_new_I/peakinf_strain1only_vacc0,
                                  vacc_start==graph_vaccstart_range[2] ~ peak_new_I/peakinf_strain1only_vacc1,
                                  vacc_start==graph_vaccstart_range[3] ~ peak_new_I/peakinf_strain1only_vacc2)) %>% 
    mutate(seed_round = fct_drop(seed_round)) %>% 
    mutate(seed_dodge = factor(seed_round, labels=seq(1,length(levels(seed_round))))) %>% 
    mutate(seed_dodge = case_when(crossimm_perc==levels(crossimm_perc)[1] ~ as.numeric(seed_dodge)-0.2,
                                  crossimm_perc==levels(crossimm_perc)[2] | seed==resident_seedproxy ~ as.numeric(seed_dodge),
                                  crossimm_perc==levels(crossimm_perc)[3] ~ as.numeric(seed_dodge)+0.2),
           line_group = str_c(crossimm_perc, label))

## Relative, peak cases, collapsed
graph_table %>% 
    filter(seed!=resident_seedproxy) %>% 
    ggplot(aes(x=seed_round, y=rel_peak_new_I, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak daily cases, relative to resident strain only",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape)

graph_table %>% 
    filter(seed!=resident_seedproxy) %>% 
    ggplot(aes(x=seed_dodge, y=rel_peak_new_I, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak daily cases, relative to resident strain only",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("rel_peak_new_I_collapse_r%s_log.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

# Peak cases, collapsed
graph_table %>% 
    ggplot(aes(x=seed_round, y=proppop_peak_new_I, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak daily cases, as % population",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent(accuracy=0.1),
                       limits=c(0,NA)) +
    scale_shape_manual(values=legend_shape)

graph_table %>% 
    ggplot(aes(x=seed_dodge, y=proppop_peak_new_I, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak daily cases, as % population",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent(accuracy=0.1),
                       limits=c(0,NA)) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)


graph_name <- file.path(".", output_folder, sprintf("peak_new_I_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Peak resident and peak novel variant
graph_table %>% 
    pivot_longer(c(proppop_peak_new_I1, proppop_peak_new_I2), names_to = "strain", values_to = "count") %>% 
    mutate(strain = fct_recode(strain,
                               "Resident strain" = "proppop_peak_new_I1",
                               "Novel variant" = "proppop_peak_new_I2")) %>% 
    ggplot(aes(x=seed_dodge, y=count, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(rows=vars(vacc_cov),
               cols=vars(strain)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak daily cases, as % population",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_percent(accuracy=0.1),
                       limits=c(0,NA)) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("peak_new_I1_I2rows_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 18, units = "cm", dpi = 300)

# Peak time -----
graph_r0 <- 2.5
graph_vaccstart_range <- sort(unique(end$vacc_start[end$r0_strain1==graph_r0]))

peaktime_strain1only_vacc0 <- end$peak_time[end$r0_strain1==graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[1]]
peaktime_strain1only_vacc1 <- end$peak_time[end$r0_strain1== graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[2]]
peaktime_strain1only_vacc2 <- end$peak_time[end$r0_strain1==graph_r0 & end$seed==resident_seedproxy & end$vacc_start == graph_vaccstart_range[3]]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_peaktime = case_when(vacc_start==graph_vaccstart_range[1] ~ peak_time-peaktime_strain1only_vacc0,
                                   vacc_start==graph_vaccstart_range[2] ~ peak_time-peaktime_strain1only_vacc1,
                                   vacc_start==graph_vaccstart_range[3] ~ peak_time-peaktime_strain1only_vacc2)) %>% 
    mutate(seed_round = fct_drop(seed_round)) %>% 
    mutate(seed_dodge = factor(seed_round, labels=seq(1,length(levels(seed_round))))) %>% 
    mutate(seed_dodge = case_when(crossimm_perc==levels(crossimm_perc)[1] ~ as.numeric(seed_dodge)-0.2,
                                  crossimm_perc==levels(crossimm_perc)[2] | seed==resident_seedproxy ~ as.numeric(seed_dodge),
                                  crossimm_perc==levels(crossimm_perc)[3] ~ as.numeric(seed_dodge)+0.2),
           line_group = str_c(crossimm_perc, label))

## Relative peak time, collapse
graph_table %>% 
    filter(seed!="NA") %>% 
    ggplot(aes(x=seed_round, y=rel_peaktime, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak day, in days from resident strain only peak",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape)
    

graph_table %>% 
    filter(seed!="NA") %>% 
    ggplot(aes(x=seed_dodge, y=rel_peaktime, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_hline(yintercept=1, size=0.4) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak day, in days from resident strain only peak",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)


graph_name <- file.path(".", output_folder, sprintf("rel_peaktime_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Peak time, collapse
graph_table %>% 
    ggplot(aes(x=seed_round, y=peak_time, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point(position=position_dodge(width=0.8)) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak day",
         colour = "Variant characteristic",
         shape = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_y_continuous(labels = label_comma()) +
    scale_shape_manual(values = legend_shape)

resident_peak <- graph_table %>% 
    filter(seed==resident_seedproxy) %>% 
    select(vacc_cov, peak_time)

graph_table %>% 
    ggplot(aes(x=seed_dodge, y=peak_time, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_hline(data=resident_peak,
               aes(yintercept=peak_time),
               size = 1,
               colour=legend_colour[["Resident strain only"]],
               alpha = 0.3) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(cols=vars(vacc_cov)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak day",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("peaktime_collapse2_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 14, units = "cm", dpi = 300)

## Peak time for resident and novel variant
graph_table %>% 
    pivot_longer(c(peak_time_new_I1, peak_time_new_I2), names_to = "strain", values_to = "count") %>% 
    mutate(strain = fct_recode(strain,
                               "Resident strain" = "peak_time_new_I1",
                               "Novel variant" = "peak_time_new_I2")) %>% 
    ggplot(aes(x=seed_dodge, y=count, colour=label, group=crossimm_perc, shape=crossimm_perc)) +
    geom_point() +
    geom_line(aes(group=line_group, linetype=crossimm_perc), 
              size = 0.5,
              alpha = 0.4) +
    facet_grid(rows=vars(vacc_cov),
               cols=vars(strain)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day variant seeded",
         y = "Peak day",
         colour = "Variant characteristic",
         shape = "Cross-immunity",
         linetype = "Cross-immunity") +
    # scale_y_log10() +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_shape_manual(values=legend_shape) +
    scale_x_continuous(labels = function(x) graph_table %>% pull(seed_round) %>% fct_drop() %>% levels() %>% .[x]) +
    scale_linetype_manual(values=legend_linetype) +
    theme(strip.text.y.right = element_text(angle = 0))

graph_name <- file.path(".", output_folder, sprintf("peak_time_I1_I2_collapse_r%s.png", graph_r0))

ggsave(graph_name, width = 25, height = 18, units = "cm", dpi = 300)

###### Cumulative --------
graph_r0 <- 2.5

graph_vacccov <- levels(fct_drop(end$vacc_cov[end$r0_strain1==graph_r0]))

graph_table <- end %>% 
    filter(r0_strain1==graph_r0) %>% 
    select(r0_strain1, transmiss, crossimm, crossimm_label, seed, vacc, vacc_start, vacc_cov, preinf, match_preinf, label, proppop_strain1_only, proppop_bothstrains, proppop_strain2_only) %>%
    mutate(crossimm_vacc = factor(sprintf("Cross-immunity %s%%\nVaccine coverage %s%%", crossimm*100, vacc_start*100)),
           seed_label = factor(case_when(seed==resident_seedproxy ~ "NA", 
                               seed!=resident_seedproxy ~ sprintf("Seed day %s", round(seed))))) %>%
    mutate(crossimm_vacc = fct_shift(crossimm_vacc, 3),
           seed_label = fct_relevel(seed_label, unique(str_sort(seed_label, numeric = T)))) %>% 
    pivot_longer(c(proppop_strain1_only, proppop_bothstrains, proppop_strain2_only), names_to = "strain", values_to = "count") %>% 
    mutate(strain = fct_recode(strain,
                               "Both strains" = "proppop_bothstrains",
                               "Resident strain only" = "proppop_strain1_only",
                               "Variant only" = "proppop_strain2_only"),
           label_new = label) %>% 
    mutate(label_new = fct_recode(label_new,
                              "More transmissable<br>(r[init] matched to d[E] 0.5)" = "More transmissable (r[init] matched to d[E] 0.5)",
                              "More transmissable<br>(r[init] matched to d[E] 1.5)" = "More transmissable (r[init] matched to d[E] 1.5)")) %>% 
    mutate(colour_label = paste("<span style = 'color: ",
                                legend_colour[label],
                                ";'>",
                                gsub("\\]", "</sub>", gsub("\\[", "<sub>", label_new)),
                                "</span>", sep = ""),
           colour_label = fct_reorder(colour_label, as.numeric(label_new)))

# ggplot(graph_table, aes(x=colour_label, y=count, fill=strain)) +
#     geom_bar(position="stack", stat="identity") +
#     facet_grid(cols = vars(seed_label),
#                rows = vars(crossimm_vacc),
#                scale = "free_x", 
#                space = "free_x") +
#     scale_y_continuous(labels = label_percent(accuracy=1),
#                        breaks = breaks_extended()) +
#     labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
#          x = "Variant characteristic",
#          y = "Final % population infected",
#          fill = "Infecting strain") +    
#     theme(plot.caption = element_text(hjust = 1),
#           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#           strip.text.y.right = element_text(angle = 0))
# 
# cum_vacc1 <- graph_table %>% 
#     filter(vacc_cov==graph_vacccov[1]) %>% 
#     ggplot(aes(x=colour_label, y=count, fill=strain)) +
#     geom_bar(position="stack", stat="identity") +
#     facet_grid(cols = vars(seed_label),
#                rows = vars(crossimm_label),
#                scale = "free_x", 
#                space = "free_x") +
#     scale_y_continuous(labels = label_percent(accuracy=1),
#                        limits = c(0,1)) +
#     labs(x = "",
#          y = "",
#          fill = "Infecting strain",
#          title = str_replace(graph_vacccov[1], "\n", " ")) +    
#     theme(plot.caption = element_text(hjust = 1),
#           axis.text.x = element_blank(),
#           axis.ticks.x=element_blank(),
#           strip.text.y.right = element_text(angle = 0))
# 
# cum_vacc2 <- graph_table %>% 
#     filter(vacc_cov==graph_vacccov[2]) %>% 
#     ggplot(aes(x=colour_label, y=count, fill=strain)) +
#     geom_bar(position="stack", stat="identity") +
#     facet_grid(cols = vars(seed_label),
#                rows = vars(crossimm_label),
#                scale = "free_x", 
#                space = "free_x") +
#     scale_y_continuous(labels = label_percent(accuracy=1),
#                        limits = c(0,1)) +
#     labs(x = "",
#          y = "Final % population infected",
#          fill = "Infecting strain",
#          title = str_replace(graph_vacccov[2], "\n", " ")) +    
#     theme(plot.caption = element_text(hjust = 1),
#           axis.text.x = element_blank(),
#           axis.ticks.x=element_blank(),
#           strip.text.x = element_blank(),
#           strip.text.y.right = element_text(angle = 0))
# 
# 
# cum_vacc3 <- graph_table %>% 
#     filter(vacc_cov==graph_vacccov[3]) %>% 
#     ggplot(aes(x=colour_label, y=count, fill=strain)) +
#     geom_bar(position="stack", stat="identity") +
#     facet_grid(cols = vars(seed_label),
#                rows = vars(crossimm_label),
#                scale = "free_x", 
#                space = "free_x") +
#     scale_y_continuous(labels = label_percent(accuracy=1),
#                        limits = c(0,1)) +
#     labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
#          x = "Variant characteristic",
#          y = "",
#          fill = "Infecting strain",
#          title = str_replace(graph_vacccov[3], "\n", " ")) +    
#     theme(plot.caption = element_text(hjust = 1),
#           axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
#           strip.text.x = element_blank(),
#           strip.text.y.right = element_text(angle = 0))
# 
# cum_vacc1 / cum_vacc2 / cum_vacc3 + plot_layout(guides = 'collect')
# 
# graph_name <- file.path(".", output_folder, sprintf("cumulative_r%s.png", graph_r0))
# 
# ggsave(graph_name, width = 24, height = 40, units = "cm", dpi = 300)

## Cumulative dodged
cum_vacc1 <- graph_table %>% 
    filter(vacc_cov==graph_vacccov[1]) %>% 
    ggplot(aes(x=colour_label, y=count, fill=strain)) +
    geom_bar(position="dodge", stat="identity") +
    facet_grid(cols = vars(seed_label),
               rows = vars(crossimm_label),
               scale = "free_x", 
               space = "free_x") +
    scale_y_continuous(labels = label_percent(accuracy=1),
                       limits = c(0,1)) +
    labs(x = "",
         y = "",
         fill = "Infecting strain",
         title = str_replace(graph_vacccov[1], "\n", " ")) +    
    theme(plot.caption = element_text(hjust = 1),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.y.right = element_text(angle = 0))

cum_vacc2 <- graph_table %>% 
    filter(vacc_cov==graph_vacccov[2]) %>% 
    ggplot(aes(x=colour_label, y=count, fill=strain)) +
    geom_bar(position="dodge", stat="identity") +
    facet_grid(cols = vars(seed_label),
               rows = vars(crossimm_label),
               scale = "free_x", 
               space = "free_x") +
    scale_y_continuous(labels = label_percent(accuracy=1),
                       limits = c(0,1)) +
    labs(x = "",
         y = "Final % population infected",
         fill = "Infecting strain",
         title = str_replace(graph_vacccov[2], "\n", " ")) +    
    theme(plot.caption = element_text(hjust = 1),
          axis.text.x = element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_blank(),
          strip.text.y.right = element_text(angle = 0))


cum_vacc3 <- graph_table %>% 
    filter(vacc_cov==graph_vacccov[3]) %>% 
    ggplot(aes(x=colour_label, y=count, fill=strain)) +
    geom_bar(position="dodge", stat="identity") +
    facet_grid(cols = vars(seed_label),
               rows = vars(crossimm_label),
               scale = "free_x", 
               space = "free_x") +
    scale_y_continuous(labels = label_percent(accuracy=1),
                       limits = c(0,1)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Variant characteristic",
         y = "",
         fill = "Infecting strain",
         title = str_replace(graph_vacccov[3], "\n", " ")) +    
    theme(plot.caption = element_text(hjust = 1),
          axis.text.x = element_markdown(angle = 90, vjust = 0.5, hjust=1),
          strip.text.x = element_blank(),
          strip.text.y.right = element_text(angle = 0))

cum_vacc1 / cum_vacc2 / cum_vacc3 + plot_layout(guides = 'collect')

graph_name <- file.path(".", output_folder, sprintf("cumulative_dodge_r%s.png", graph_r0))

ggsave(graph_name, width = 30, height = 40, units = "cm", dpi = 300)

#### Plotting new infections, all seeds----

# legend_colour[["Pre-infectious 0.5 days"]] <- "#ff8989"
# legend_colour[["Pre-infectious 1.5 days"]] <- "#8989ff"
# names(legend_linetype) <- str_c("Cross-immunity", names(legend_linetype), sep=" ")


graph_r0 <- 2.5
graph_limit <- 400

plot_scenarios <- end %>% filter(r0_strain1==graph_r0)

# plot_scenarios <- rbind(plot_scenarios, end[end$seed==resident_seedproxy & end$r0_strain1==graph_r0,])

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
    mutate(seed = if_else(seed==resident_seedproxy, 0, round(seed))) %>% 
    mutate(crossimm_seed = sprintf("Seed day %s\nCross-immunity %s%%", seed, crossimm*100)) %>% 
    mutate(crossimm_seed = fct_relevel(crossimm_seed, unique(str_sort(crossimm_seed, numeric = T)))) %>% 
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "More transmissable (r[init] matched to d[E] 0.5)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "More transmissable (r[init] matched to d[E] 1.5)",
                             preinf==0.5 ~ "Pre-infectious (d[E]) 0.5 days",
                             preinf==1.5 ~ "Pre-infectious (d[E]) 1.5 days",
                             preinf==2.5 & transmiss==1 ~ "Resident strain only")) %>% 
    mutate(label = fct_relevel(label,
                               c("Resident strain only",
                                 "Pre-infectious (d[E]) 0.5 days",
                                 "More transmissable (r[init] matched to d[E] 0.5)",
                                 "Pre-infectious (d[E]) 1.5 days",
                                 "More transmissable (r[init] matched to d[E] 1.5)"
                               ))) %>% 
    mutate(vacc_label = factor(sprintf("Vaccine coverage %s%%", vacc_start*100)),
           seed_label = factor(sprintf("Seed day\n%s", round(seed))),
           crossimm_perc = factor(sprintf("%s%%", crossimm*100)),
           crossimm_label = factor(sprintf("Cross-immunity\n%s%%", crossimm*100))) %>% 
    mutate(vacc_label = fct_relevel(vacc_label, unique(str_sort(vacc_label, numeric = T))),
           seed_label = fct_relevel(seed_label, unique(str_sort(seed_label, numeric = T))),
           crossimm_perc = fct_relevel(crossimm_perc, unique(str_sort(crossimm_perc, numeric = T))),
           crossimm_label = fct_relevel(crossimm_label, unique(str_sort(crossimm_label, numeric = T))))


## Seed day lines
# seed_line <- graph_table %>% 
#     group_by(crossimm_seed) %>% 
#     summarise(seed_line = mean(seed))

## Daily new cases, percentage population
# ggplot(graph_table, aes(x=day, y=proppop_new_I, colour=label)) +
#     geom_line() +
#     facet_grid(rows=vars(crossimm_seed),
#                cols=vars(vacc_label)) +
#     scale_colour_manual(values=legend_colour) +
#     scale_x_continuous(limits=c(0,graph_limit)) +
#     geom_vline(data=seed_line, 
#                aes(xintercept=seed_line), 
#                linetype="dashed", size=0.5) +
#     scale_y_continuous(labels = label_percent(accuracy=0.1)) +
#     labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
#          x = "Day",
#          y = "Overall daily new cases (% population)",
#          colour = "Variant characteristic") +
#     theme(plot.caption = element_text(hjust = 1),
#           strip.text.y.right = element_text(angle = 0))
#     
# graph_name <- file.path(".", output_folder, sprintf("newinf_r%s.png", graph_r0))
# 
# ggsave(graph_name, width = 28, height = 25, units = "cm", dpi = 300)

##
### Daily new cases, crossimm same graph, percentage population
seed_line <- graph_table %>% 
    group_by(seed_label) %>% 
    summarise(seed_line = mean(seed))

ggplot(graph_table, aes(x=day, y=proppop_new_I, colour=label, linetype=crossimm_perc)) +
    geom_line() +
    facet_grid(rows=vars(seed_label),
               cols=vars(vacc_label)) +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_continuous(labels = label_percent(accuracy=0.1)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day",
         y = "Daily new cases, as % population",
         colour = "Variant characteristic",
         linetype = "Cross-immunity") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0)) +
    scale_linetype_manual(values=legend_linetype)
    
graph_name <- file.path(".", output_folder, sprintf("newinfcrossimm_r%s.png", graph_r0))

ggsave(graph_name, width = 28, height = 25, units = "cm", dpi = 300)


## Daily new cases of resident and variant, percentage population
seed_line <- graph_table %>%
    group_by(crossimm_seed) %>%
    summarise(seed_line = mean(seed))

graph_table %>% 
    pivot_longer(c("proppop_new_I1", "proppop_new_I2"), names_to="strain", values_to="count") %>% 
    mutate(strain = fct_recode(strain, 
                               "Resident strain" = "proppop_new_I1",
                               "Novel variant" = "proppop_new_I2")) %>% 
    mutate(strain = fct_relevel(strain, "Novel variant")) %>% 
    ggplot(aes(x=day, y=count, colour=label, linetype=strain)) +
    geom_line() +
    facet_grid(rows=vars(crossimm_seed),
               cols=vars(vacc_label)) +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_continuous(labels = label_percent(accuracy=0.1)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day",
         y = "Daily new cases, as % population",
         colour = "Variant characteristic",
         linetype = "Strain") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0))

graph_name <- file.path(".", output_folder, sprintf("newinf1and2_r%s.png", graph_r0))

ggsave(graph_name, width = 28, height = 25, units = "cm", dpi = 300)

### Separated by vacc cov, new resident and new variant
graph_vacc_label <- levels(graph_table$vacc_label[graph_table$r0_strain1==graph_r0])

seed_line <- graph_table %>%
    group_by(seed_label) %>%
    summarise(seed_line = mean(seed))

i <- 3
graph_limit <- 350

# for (i in 1:length(graph_vacc_label)) {
    graph_table %>% 
        filter(vacc_label==graph_vacc_label[i]) %>% 
        pivot_longer(c("proppop_new_I1", "proppop_new_I2"), names_to="strain", values_to="count") %>% 
        mutate(strain = fct_recode(strain, 
                                   "Resident strain" = "proppop_new_I1",
                                   "Novel variant" = "proppop_new_I2")) %>% 
        mutate(strain = fct_relevel(strain, "Novel variant")) %>% 
        ggplot(aes(x=day, y=count, colour=label, linetype=strain)) +
        geom_line() +
        facet_grid(rows=vars(crossimm_label),
                   cols=vars(seed_label)) +
        scale_colour_manual(values=legend_colour,
                            labels=variant_label) +
        scale_x_continuous(limits=c(0,graph_limit)) +
        geom_vline(data=seed_line, 
                   aes(xintercept=seed_line), 
                   linetype="dashed", size=0.5) +
        scale_y_continuous(labels = label_percent(accuracy=0.1)) +
        labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
             x = "Day",
             y = "Daily new cases, as % population",
             colour = "Variant characteristic",
             title = graph_vacc_label[i],
             linetype = "Strain") +
        theme(plot.caption = element_text(hjust = 1),
              strip.text.y.right = element_text(angle = 0))
    
    graph_name <- file.path(".", output_folder, sprintf("newinf1and2_r%s_vacccov%s.png", graph_r0, str_extract(graph_vacc_label[i], "\\d+\\.?\\d*")))
    
    ggsave(graph_name, width = 40, height = 25, units = "cm", dpi = 300)
# }

### Separated by vacc cov, new infections overall
graph_vacc_label <- levels(graph_table$vacc_label[graph_table$r0_strain1==graph_r0])
    
seed_line <- graph_table %>%
        group_by(seed_label) %>%
        summarise(seed_line = mean(seed))
    
i <- 3
graph_limit <- 350
    
# for (i in 1:length(graph_vacc_label)) {
    graph_table %>% 
        filter(vacc_label==graph_vacc_label[i]) %>% 
        ggplot(aes(x=day, y=proppop_new_I, colour=label)) +
        geom_line() +
        facet_grid(rows=vars(crossimm_label),
                   cols=vars(seed_label)) +
        scale_colour_manual(values=legend_colour,
                            labels=variant_label) +
        scale_x_continuous(limits=c(0,graph_limit)) +
        geom_vline(data=seed_line, 
                   aes(xintercept=seed_line), 
                   linetype="dashed", size=0.5) +
        scale_y_continuous(labels = label_percent(accuracy=0.1)) +
        labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
             x = "Day",
             y = "Daily new cases, as % population",
             colour = "Variant characteristic",
             title = graph_vacc_label[i]) +
        theme(plot.caption = element_text(hjust = 1),
              strip.text.y.right = element_text(angle = 0))
    
    graph_name <- file.path(".", output_folder, sprintf("newinf_r%s_vacccov%s.png", graph_r0, str_extract(graph_vacc_label[i], "\\d+\\.?\\d*")))
    
    ggsave(graph_name, width = 40, height = 25, units = "cm", dpi = 300)
# }
    

## Daily new cases, log scale
seed_line <- graph_table %>% 
    group_by(seed_label) %>% 
    summarise(seed_line = mean(seed))

### Variant
variant_log <- graph_table %>% 
    filter(label!="Resident strain only") %>% 
    ggplot(aes(x=day, y=daily_new_I2, colour=label, linetype=crossimm_perc)) +
    geom_line() +
    facet_grid(rows=vars(seed_label),
               cols=vars(vacc_label)) +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_log10(
        limits = c(1e-14, 1e+05)
        ) +
    labs(
        # caption = bquote(atop("Resident strain" ~ R[0] == .(graph_r0), "Y-axis truncated to magnify initial dynamics")
        caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day",
         y = "Daily new variant cases",
         colour = "Variant characteristic",
         linetype = "Cross-immunity") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0)) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("newinf2_r%s_log.png", graph_r0))

ggsave(graph_name, width = 30, height = 20, units = "cm", dpi = 300)

## Resident
resident_log <- graph_table %>% 
    ggplot(aes(x=day, y=daily_new_I1, colour=label, linetype=crossimm_perc)) +
    geom_line() +
    facet_grid(rows=vars(seed_label),
               cols=vars(vacc_label)) +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_log10(
        limits = c(1e-14, 1e+05)
    ) +
    labs(
        # caption = bquote(atop("Resident strain" ~ R[0] == .(graph_r0), "Y-axis truncated to magnify initial dynamics")
        caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
        x = "Day",
        y = "Daily new resident cases",
        colour = "Variant characteristic",
        linetype = "Cross-immunity") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0)) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("newinf1_r%s_log.png", graph_r0))

ggsave(graph_name, width = 30, height = 20, units = "cm", dpi = 300)

variant_log / resident_log + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A')

graph_name <- file.path(".", output_folder, sprintf("newinf1_2_r%s_log_grid.png", graph_r0))

ggsave(graph_name, width = 30, height = 40, units = "cm", dpi = 300)


### Combined
graph_table %>% 
    ggplot(aes(x=day, y=daily_new_I, colour=label, linetype=crossimm_perc)) +
    geom_line() +
    facet_grid(rows=vars(seed_label),
               cols=vars(vacc_label)) +
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_log10(
        # limits = c(1e-8, 1e-01)
    ) +
    labs(
        # caption = bquote(atop("Resident strain" ~ R[0] == .(graph_r0), "Y-axis truncated to magnify initial dynamics")
        caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
        x = "Day",
        y = "Daily new cases",
        colour = "Variant characteristic",
        linetype = "Cross-immunity") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0)) +
    scale_linetype_manual(values=legend_linetype)

graph_name <- file.path(".", output_folder, sprintf("newinf_r%s_log.png", graph_r0))

ggsave(graph_name, width = 30, height = 20, units = "cm", dpi = 300)

#### Plotting susceptibles and R1, all seeds----
graph_r0 <- 4
graph_limit <- 300
    
plot_scenarios <- end %>% filter(r0_strain1==graph_r0)

# plot_scenarios <- rbind(plot_scenarios, end[end$seed==resident_seedproxy & end$r0_strain1==graph_r0,])

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

# S_R1 <- read.csv(file.path(".", output_folder, sprintf("S_R1_r%s.csv", graph_r0)))

graph_table <- S_R1 %>% 
    mutate(seed = if_else(seed==resident_seedproxy, 0, round(seed))) %>% 
    mutate(crossimm_seed = sprintf("Seed day %s\nCross-immunity %s%%", seed, crossimm*100)) %>% 
    mutate(crossimm_seed = fct_relevel(crossimm_seed, unique(str_sort(crossimm_seed, numeric = T)))) %>% 
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "More transmissable (r[init] matched to d[E] 0.5)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "More transmissable (r[init] matched to d[E] 1.5)",
                             preinf==0.5 ~ "Pre-infectious (d[E]) 0.5 days",
                             preinf==1.5 ~ "Pre-infectious (d[E]) 1.5 days",
                             preinf==2.5 & transmiss==1 ~ "Resident strain only")) %>% 
    mutate(label = fct_relevel(label,
                               c("Resident strain only",
                                 "Pre-infectious (d[E]) 0.5 days",
                                 "More transmissable (r[init] matched to d[E] 0.5)",
                                 "Pre-infectious (d[E]) 1.5 days",
                                 "More transmissable (r[init] matched to d[E] 1.5)"
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
    scale_colour_manual(values=legend_colour,
                        labels=variant_label) +
    scale_x_continuous(limits=c(0,graph_limit)) +
    geom_vline(data=seed_line, 
               aes(xintercept=seed_line), 
               linetype="dashed", size=0.5) +
    scale_y_continuous(labels = label_percent(accuracy=0.1)) +
    labs(caption = bquote("Resident strain" ~ R[0] == .(graph_r0)),
         x = "Day",
         y = "% population",
         colour = "Variant characteristic",
         linetype = "Compartment") +
    theme(plot.caption = element_text(hjust = 1),
          strip.text.y.right = element_text(angle = 0)) 


graph_name <- file.path(".", output_folder, sprintf("S_R1_r%s.png", graph_r0))

ggsave(graph_name, width = 28, height = 28, units = "cm", dpi = 300)
