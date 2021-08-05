################ Change vaccstart str extract

library(tidyverse)
library(ggrepel)
library(RColorBrewer)

theme_set(
    theme_gray(base_size = 12)
)

output_folder <- "seir_model_output"

end <- read.csv(file.path(".", output_folder, "end.csv"))

end %>% 
    count(match_set)

# Extract parameters
end <- end %>%
    mutate(scenario_pt = str_replace_all(scenario, "pt", "."),
           file = str_c(scenario, sys_time, sep="_")) %>% 
    mutate(
        r0_strain1 = str_extract(scenario_pt, "(?<=r)\\d*\\.?\\d*"),
        transmiss = str_extract(scenario_pt, "(?<=transmiss)\\d*\\.?\\d*"),
        crossimm = str_extract(scenario_pt, "(?<=crossimm)\\d*\\.?\\d*"),
        seed = str_extract(scenario_pt, "(?<=seed)\\d*"),
        vacc = str_extract(scenario_pt, "(?<=vacc)\\d*\\.?\\d*"),
        vacc_start = str_extract(scenario_pt, "(?<=start_)\\d*\\.?\\d*"),
        preinf = str_extract(scenario_pt, "(?<=preinf)\\d*\\.?\\d*")) %>% 
    mutate(across(c(r0_strain1, crossimm, seed, vacc, vacc_start, preinf, transmiss), as.numeric)) %>% 
    mutate(match_set_pt = str_replace_all(match_set, "pt", ".")) %>% 
    mutate(match_preinf = str_extract(match_set_pt, "(?<=preinf)\\d*\\.?\\d*")) %>% 
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "more transmiss (preinf 0.5)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "more transmiss (preinf 1.5)",
                             preinf==0.5 ~ "preinf 0.5",
                             preinf==1.5 ~ "preinf 1.5",
                             preinf==2.5 & transmiss==1 ~ "resident only")) 

# junk ------
# Combos
# combos <- end %>%
#     filter(transmiss!=1)
#     select(r0_strain1, transmiss=higher_r0_rel, crossimm, seed, preinf=higher_r0_preinf, match_set=scenario) %>%
#     mutate(vacc=1,
#            vacc_start30 = 0.3,
#            vacc_start60 = 0.6) %>%
#     pivot_longer(starts_with("vacc_start"), names_to = "name", values_to = "vacc_start") %>%
#     select(!name) %>%
#     rename_with(~ paste(.x, "range", sep="_"), !match_set) %>%
#     mutate(across(!match_set, as.numeric))
#     
# write.csv(combos, file.path(".", output_folder, "vacc_combos.csv"), row.names=F)

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

end %>%
    filter(r0_strain1 == 4, vacc_start==0, seed!=1200) %>% 
    ggplot(aes(x=scenario, y=allinf, fill=preinf)) +
        geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#### Total size ------
graph_r0 <- 4

allinf_strain1only_vacc0 <- end$allinf[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == 0]
allinf_strain1only_vacc0pt3 <- end$allinf[end$r0_strain1== graph_r0 & end$seed==1200 & end$vacc_start == 0.3]
allinf_strain1only_vacc0pt6 <- end$allinf[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == 0.6]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_allinf = case_when(vacc_start==0 ~ allinf/allinf_strain1only_vacc0,
                                      vacc_start==0.3 ~ allinf/allinf_strain1only_vacc0pt3,
                                      vacc_start==0.6 ~ allinf/allinf_strain1only_vacc0pt6))

# Relative, total cases
graph_table %>% 
    filter(seed!=1200) %>% 
    ggplot(aes(x=as.factor(seed), y=rel_allinf, colour=label, group=match_set, label=transmiss)) +
    geom_point(position=position_dodge(width=0.5)) +
    facet_grid(rows=vars(crossimm),
               cols=vars(vacc_start),
               labeller=label_both) +
    # geom_text_repel(data=subset(graph_table, transmiss > 1),
    #                 position=position_dodge(width=0.5),
    #                 size=3,
    #                 min.segment.length = 0) +
    labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
    scale_colour_manual(values=c(
        "#1F78B4",
        "#E31A1C",
        "#A6CEE3",
        "#FB9A99"))

# Total cases
graph_table %>%
    mutate(seed = factor(replace(seed, seed==1200, "NA"))) %>%
    ggplot(aes(x=seed, y=allinf, colour=label, group=match_set, label=transmiss)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start),
                   labeller=label_both) +
        # geom_text_repel(data=subset(graph_table, transmiss > 1),
        #                 position=position_dodge(width=0.5),
        #                 size=3,
        #                 min.segment.length = 0) +
        labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
        scale_colour_manual(values=c(
            "#1F78B4",
            "#E31A1C",
            "#A6CEE3",
            "#FB9A99",
            "#6A3D9A")) +
    scale_y_continuous(labels = scales::label_comma())

# Peak new cases --------
graph_r0 <- 4

peakinf_strain1only_vacc0 <- end$peak_new_I[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == 0]
peakinf_strain1only_vacc0pt3 <- end$peak_new_I[end$r0_strain1== graph_r0 & end$seed==1200 & end$vacc_start == 0.3]
peakinf_strain1only_vacc0pt6 <- end$peak_new_I[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == 0.6]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_peakinf = case_when(vacc_start==0 ~ peak_new_I/peakinf_strain1only_vacc0,
                                  vacc_start==0.3 ~ peak_new_I/peakinf_strain1only_vacc0pt3,
                                  vacc_start==0.6 ~ peak_new_I/peakinf_strain1only_vacc0pt6))

# Relative, peak cases
graph_table %>% 
    filter(seed!=1200) %>% 
    ggplot(aes(x=as.factor(seed), y=rel_peakinf, colour=label, group=match_set, label=transmiss)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start),
                   labeller=label_both) +
        # geom_text_repel(data=subset(graph_table, transmiss > 1),
        #                 position=position_dodge(width=0.5),
        #                 size=3,
        #                 min.segment.length = 0) +
        labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
        scale_colour_manual(values=c(
            "#1F78B4",
            "#E31A1C",
            "#A6CEE3",
            "#FB9A99"))

# Peak cases
graph_table %>% 
    mutate(seed = factor(replace(seed, seed==1200, "NA"))) %>%
    ggplot(aes(x=as.factor(seed), y=peak_new_I, colour=label, group=match_set, label=transmiss)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start),
                   labeller=label_both) +
        # geom_text_repel(data=subset(graph_table, transmiss > 1),
        #                 position=position_dodge(width=0.5),
        #                 size=3,
        #                 min.segment.length = 0) +
        labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
    scale_colour_manual(values=c(
        "#1F78B4",
        "#E31A1C",
        "#A6CEE3",
        "#FB9A99",
        "#6A3D9A")) +
    scale_y_continuous(labels = scales::label_comma())
    

# Peak time -----
graph_r0 <- 4

peaktime_strain1only_vacc0 <- end$peak_time[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == 0]
peaktime_strain1only_vacc0pt3 <- end$peak_time[end$r0_strain1== graph_r0 & end$seed==1200 & end$vacc_start == 0.3]
peaktime_strain1only_vacc0pt6 <- end$peak_time[end$r0_strain1==graph_r0 & end$seed==1200 & end$vacc_start == 0.6]

graph_table <- end %>%
    filter(r0_strain1 == graph_r0) %>%
    mutate(rel_peaktime = case_when(vacc_start==0 ~ peak_time-peaktime_strain1only_vacc0,
                                   vacc_start==0.3 ~ peak_time-peaktime_strain1only_vacc0pt3,
                                   vacc_start==0.6 ~ peak_time-peaktime_strain1only_vacc0pt6))

# Relative peak time
graph_table %>% 
    filter(seed!=1200) %>% 
    ggplot(aes(x=as.factor(seed), y=rel_peaktime, colour=label, group=match_set, label=transmiss)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start),
                   labeller=label_both) +
        # geom_text_repel(data=subset(graph_table, transmiss > 1),
        #                 position=position_dodge(width=0.5),
        #                 size=3,
        #                 min.segment.length = 0) +
        labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
        scale_colour_manual(values=c(
            "#1F78B4",
            "#E31A1C",
            "#A6CEE3",
            "#FB9A99"))

# Peak time
graph_table %>% 
    mutate(seed = factor(replace(seed, seed==1200, "NA"))) %>%
    ggplot(aes(x=as.factor(seed), y=peak_time, colour=label, group=match_set, label=transmiss)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_grid(rows=vars(crossimm),
                   cols=vars(vacc_start),
                   labeller=label_both) +
        # geom_text_repel(data=subset(graph_table, transmiss > 1),
        #                 position=position_dodge(width=0.5),
        #                 size=3,
        #                 min.segment.length = 0) +
        labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
        scale_colour_manual(values=c(
            "#1F78B4",
            "#E31A1C",
            "#A6CEE3",
            "#FB9A99",
            "#6A3D9A")) +
        scale_y_continuous(labels = scales::label_comma())


###


###### Cumulative --------
graph_r0 <- 4
graph_vacc_start <- 0

graph_table <- end %>% 
    filter(r0_strain1==graph_r0) %>% 
    select(r0_strain1, transmiss, crossimm, seed, vacc, vacc_start, preinf, match_preinf, label, strain1_only, bothstrains, strain2_only) %>%
    mutate(crossimm_vacc = sprintf("crossimm %s\nvacc %s", crossimm, vacc_start),
           seed_label = case_when(seed==1200 ~ "none", 
                               seed!=1200 ~ sprintf("seed %s", seed))) %>%
    mutate(seed_label = factor(seed_label)) %>% 
    pivot_longer(c(strain1_only, bothstrains, strain2_only), names_to = "strain", values_to = "count")

ggplot(graph_table, aes(x=label, y=count, fill=strain)) +
    geom_bar(position="stack", stat="identity") +
    facet_grid(cols = vars(seed_label),
               rows = vars(crossimm_vacc),
               scale = "free_x", 
               space = "free_x") +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(caption=sprintf("Strain 1 R0 = %s", graph_r0)) +
    theme(plot.caption = element_text(hjust = 1),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          strip.text.y.right = element_text(angle = 0))
    

#### Plotting new infections ----
graph_r0 <- 4
graph_seed <- 0

plot_scenarios <- end %>% filter(r0_strain1==graph_r0, seed==graph_seed)

plot_scenarios <- rbind(plot_scenarios, end[end$seed==1200 & end$r0_strain1==graph_r0,])

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

daily_inf_table %>% 
    mutate(label = case_when(preinf==2.5 & match_preinf==0.5 & transmiss > 1 ~ "more transmiss (preinf 0.5)",
                             preinf==2.5 & match_preinf==1.5 & transmiss > 1 ~ "more transmiss (preinf 1.5)",
                             preinf==0.5 ~ "preinf 0.5",
                             preinf==1.5 ~ "preinf 1.5",
                             preinf==2.5 & transmiss==1 ~ "resident only")) %>% 
    mutate(label=as.factor(label)) %>% 
    ggplot(aes(x=day, y=daily_new_I, colour=label)) +
    geom_line(size=0.8) +
    facet_grid(cols=vars(crossimm),
               rows=vars(vacc_start)) +
    scale_colour_manual(values=c(
        "#1F78B4",
        "#E31A1C",
        "#A6CEE3",
        "#FB9A99",
        "#6A3D9A"
    )) +
    scale_x_continuous(limits=c(0,300)) +
    geom_vline(xintercept = graph_seed, linetype="dashed", size=0.5) +
    scale_y_continuous(labels = scales::label_comma())
    
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

