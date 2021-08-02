library(tidyverse)
library(RColorBrewer)

theme_set(
    theme_gray(base_size = 12)
)

# end <- read.csv(file.path(".", output_folder, "end.csv"))
# files <- end$file
# calc_file <- files[2]
# scenario <- end$scenario[end$file==calc_file]
# sys_time <- end$sys_time[end$file==calc_file]
# output_folder <- "seir_model_output"

# scenario <- "transmiss1-25_seed60_novacc"
# seedtime2 <- as.numeric(str_extract(scenario, "(?<=seed)[:digit:]*"))
# sys_time <- "12-07-2021--14-56"
# calc_file <- sprintf("%s_%s_calc.csv", scenario, sys_time)

table <- read.csv(file.path(".", output_folder, calc_file)) 

table_collapsed <- table %>% filter(time %in% seq(0, max(time), 1))
    
# Plot growth rate
table_collapsed %>%  
    select(c(time, growth_rate1, growth_rate2)) %>% 
    pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
    ggplot(aes(x=time, y=count, colour=compartment)) +
    geom_line() +
    scale_y_continuous(labels = scales::label_comma()) +
    geom_vline(xintercept = seedtime2,
               linetype="dashed",
               size=0.5) +
    geom_hline(yintercept=0,
               size=0.5)

graph_name <- file.path(".", output_folder, sprintf("%s_%s_growthrate.png", scenario, sys_time))

ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)

# Plot cumulative
table_collapsed %>%  
    select(c(time, strain1_only, strain2_only, bothstrains)) %>% 
    pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
    ggplot(aes(x=time, y=count, colour=compartment)) +
    geom_line() +
    scale_y_continuous(labels = scales::label_comma()) +
    geom_vline(xintercept = seedtime2,
               linetype="dashed",
               size=0.5)

graph_name <- file.path(".", output_folder, sprintf("%s_%s_cumulative.png", scenario, sys_time))

ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)

## Stacked area
table_collapsed %>%  
    select(c(time, strain1_only, strain2_only, bothstrains)) %>% 
    pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
    ggplot(aes(x=time, y=count, fill=compartment)) +
    geom_area() +
    scale_y_continuous(labels = scales::label_comma()) +
    geom_vline(xintercept = seedtime2,
               linetype="dashed",
               size=0.5)


graph_name <- file.path(".", output_folder, sprintf("%s_%s_cumstacked.png", scenario, sys_time))

ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)

## Stacked proportion
table_collapsed %>%  
    select(c(time, prop_strain1_only, prop_strain2_only, prop_bothstrains)) %>% 
    pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
    ggplot(aes(x=time, y=count, fill=compartment)) +
    geom_area() +
    scale_y_continuous(labels = scales::label_comma()) +
    geom_vline(xintercept = seedtime2,
               linetype="dashed",
               size=0.5)

graph_name <- file.path(".", output_folder, sprintf("%s_%s_cumstackedprop.png", scenario, sys_time))

ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)

## Checking cumulative calculations
# table_collapsed %>%  
#     select(c(time, strain1_only, strain2_only, bothstrains, allinf)) %>%
#     filter(time==500)



# Plot new infectious
new_inf_table <- table %>% 
    mutate(day = as.integer(time)) %>% 
    group_by(day) %>% 
    summarise(daily_new_I1 = sum(new_I1),
              daily_new_I2 = sum(new_I2),
              daily_new_I = sum(new_I)) %>% 
    mutate(prop_new_I1 = daily_new_I1/daily_new_I,
           prop_new_I2 = daily_new_I2/daily_new_I)

dailyinf_file <- sprintf("%s_%s_dailyinf.csv", scenario, sys_time)

write.csv(new_inf_table, file.path(".", output_folder, dailyinf_file), row.names=F)

## Daily
new_inf_table %>%
    select(c(day, daily_new_I1, daily_new_I2)) %>% 
    pivot_longer(cols=!day, names_to = "compartment", values_to = "count") %>%
    ggplot(aes(x=day, y=count, colour=compartment)) +
    geom_line() +
    scale_y_continuous(labels = scales::label_comma()) +
    geom_vline(xintercept = seedtime2,
               linetype="dashed",
               size=0.5)

graph_name <- file.path(".", output_folder, sprintf("%s_%s_newinf.png", scenario, sys_time))

ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)

## Stacked area
new_inf_table %>%  
    select(c(day, daily_new_I1, daily_new_I2)) %>% 
    pivot_longer(cols=!day, names_to = "compartment", values_to = "count") %>% 
    ggplot(aes(x=day, y=count, fill=compartment)) +
    geom_area() +
    scale_y_continuous(labels = scales::label_comma()) +
    geom_vline(xintercept = seedtime2,
               linetype="dashed",
               size=0.5)

graph_name <- file.path(".", output_folder, sprintf("%s_%s_newinf_stacked.png", scenario, sys_time))

ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)

## Stacked proportion
new_inf_table %>%  
    select(c(day, prop_new_I1, prop_new_I2)) %>% 
    pivot_longer(cols=!day, names_to = "compartment", values_to = "count") %>% 
    ggplot(aes(x=day, y=count, fill=compartment)) +
    geom_area() +
    geom_vline(xintercept = seedtime2,
               linetype="dashed",
               size=0.5)

graph_name <- file.path(".", output_folder, sprintf("%s_%s_newinf_prop.png", scenario, sys_time))

ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
