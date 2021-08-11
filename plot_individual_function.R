library(tidyverse)
# library(RColorBrewer)

theme_set(
    theme_gray(base_size = 12)
)

plot_graphs <- function(output_folder, scenario, seedtime2, sys_time, vacc_start, end_day){

    calc_file <- sprintf("%s_%s_calc.csv", scenario, sys_time)
        
    table <- read.csv(file.path(".", output_folder, calc_file)) %>% 
        filter(time <= end_day)
    
    table_collapsed <- table %>% filter(time %in% seq(0, max(time), 1))
        
    # Plot growth rate ----
    table_collapsed %>%  
        select(c(time, growth_rate1, growth_rate2)) %>% 
        pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
        ggplot(aes(x=time, y=count, colour=compartment)) +
        geom_line() +
        scale_y_continuous(labels = scales::label_comma()) +
        geom_vline(xintercept = seedtime2,
                   linetype="dashed",
                   size=0.5) +
        geom_vline(xintercept = vacc_start,
                   linetype="dashed",
                   color="#6A3D9A",
                   size=0.5) +
        geom_hline(yintercept=0,
                   size=0.5)
    
    graph_name <- file.path(".", output_folder, sprintf("%s_%s_growthrate.png", scenario, sys_time))
    
    ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
    
    # Plot cumulative infected ----
    # ## Stacked area
    # table_collapsed %>%  
    #     select(c(time, strain1_only, strain2_only, bothstrains)) %>% 
    #     pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
    #     ggplot(aes(x=time, y=count, fill=compartment)) +
    #     geom_area() +
    #     scale_y_continuous(labels = scales::label_comma()) +
    #     geom_vline(xintercept = seedtime2,
    #                linetype="dashed",
    #                size=0.5) +
    #     geom_vline(xintercept = vacc_start,
    #                linetype="dashed",
    #                color="#6A3D9A",
    #                size=0.5)
    # 
    # 
    # graph_name <- file.path(".", output_folder, sprintf("%s_%s_cumul_inf.png", scenario, sys_time))
    # 
    # ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
    
    # ## Stacked proportion
    # table_collapsed %>%  
    #     select(c(time, prop_strain1_only, prop_strain2_only, prop_bothstrains)) %>% 
    #     pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
    #     ggplot(aes(x=time, y=count, fill=compartment)) +
    #     geom_area() +
    #     scale_y_continuous(labels = scales::label_comma()) +
    #     geom_vline(xintercept = seedtime2,
    #                linetype="dashed",
    #                size=0.5) +
    #     geom_vline(xintercept = vacc_start,
    #                linetype="dashed",
    #                color="#6A3D9A",
    #                size=0.5)
    # 
    # graph_name <- file.path(".", output_folder, sprintf("%s_%s_cumul_inf_prop.png", scenario, sys_time))
    # 
    # ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
    
    ## Percentage population
    table_collapsed %>%  
        select(c(time, proppop_strain1_only, proppop_strain2_only, proppop_bothstrains)) %>% 
        pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>% 
        ggplot(aes(x=time, y=count, fill=compartment)) +
        geom_area() +
        scale_y_continuous(labels = scales::percent) +
        geom_vline(xintercept = seedtime2,
                   linetype="dashed",
                   size=0.5) +
        geom_vline(xintercept = vacc_start,
                   linetype="dashed",
                   color="#6A3D9A",
                   size=0.5)
    
    graph_name <- file.path(".", output_folder, sprintf("%s_%s_cumulinf_proppop.png", scenario, sys_time))
    
    ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
    
    ## Checking cumulative calculations
    # table_collapsed %>%  
    #     select(c(time, strain1_only, strain2_only, bothstrains, allinf)) %>%
    #     filter(time==500)
    
    # Plot percentage S, V, R1, R2, sus1, sus2 ----
    table_collapsed %>%  
        select(c(time, proppop_S, proppop_V, sus1, sus2, proppop_R1, proppop_R2, prevalence_I1, prevalence_I2, prevalence_I)) %>% 
        pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>%
        mutate(grid = case_when(
            compartment %in% c("proppop_S", "proppop_V", "sus1", "sus2") ~ "A",
            compartment %in% c("proppop_R1", "proppop_R2") ~ "B",
            compartment %in% c("prevalence_I1", "prevalence_I2", "prevalence_I") ~ "C")) %>% 
        ggplot(aes(x=time, y=count, color=compartment)) +
        geom_line() +
        scale_y_continuous(labels = scales::percent) +
        geom_vline(xintercept = seedtime2,
                   linetype="dashed",
                   size=0.5) +
        geom_vline(xintercept = vacc_start,
                   linetype="dashed",
                   color="#6A3D9A",
                   size=0.5) +
        facet_grid(grid ~.,
                   scale="free") +
        theme(strip.background = element_blank(),
              strip.text.y = element_blank()
        ) +
        scale_colour_manual(values = c(
            "proppop_S" = "#B15928",
            "proppop_V" = "#6A3D9A",
            "sus1" = "#A6CEE3",
            "sus2" = "#1F78B4",
            "proppop_R1" = "#FDBF6F",
            "proppop_R2" = "#FF7F00",
            "prevalence_I1" = "#FB9A99",
            "prevalence_I2" = "#E31A1C",
            "prevalence_I"= "#67000D")
        )
    
    graph_name <- file.path(".", output_folder, sprintf("%s_%s_proppop.png", scenario, sys_time))
    
    ggsave(graph_name, width = 15, height = 10, units = "cm", dpi = 300)
    
    
    # Plot new infectious ----
    dailyinf_file <- sprintf("%s_%s_dailyinf.csv", scenario, sys_time)
    
    new_inf_table <- read.csv(file.path(".", output_folder, dailyinf_file)) %>% 
        filter(day <= end_day)
    
    ## Daily
    new_inf_table %>%
        select(c(day, daily_new_I1, daily_new_I2, daily_new_I)) %>% 
        pivot_longer(cols=!day, names_to = "compartment", values_to = "count") %>%
        ggplot(aes(x=day, y=count, colour=compartment)) +
        geom_line() +
        scale_y_continuous(labels = scales::label_comma()) +
        geom_vline(xintercept = seedtime2,
                   linetype="dashed",
                   size=0.5) +
        geom_vline(xintercept = vacc_start,
                   linetype="dashed",
                   color="#6A3D9A",
                   size=0.5)
    
    graph_name <- file.path(".", output_folder, sprintf("%s_%s_newinf.png", scenario, sys_time))
    
    ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
    
    # ## Stacked area
    # new_inf_table %>%  
    #     select(c(day, daily_new_I1, daily_new_I2)) %>% 
    #     pivot_longer(cols=!day, names_to = "compartment", values_to = "count") %>% 
    #     ggplot(aes(x=day, y=count, fill=compartment)) +
    #     geom_area() +
    #     scale_y_continuous(labels = scales::label_comma()) +
    #     geom_vline(xintercept = seedtime2,
    #                linetype="dashed",
    #                size=0.5)
    # 
    # graph_name <- file.path(".", output_folder, sprintf("%s_%s_newinf_stacked.png", scenario, sys_time))
    # 
    # ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
    
    
    ## Stacked proportion
    new_inf_table %>%  
        select(c(day, prop_new_I1, prop_new_I2)) %>% 
        pivot_longer(cols=!day, names_to = "compartment", values_to = "count") %>% 
        ggplot(aes(x=day, y=count, fill=compartment)) +
        geom_area() +
        geom_vline(xintercept = seedtime2,
                   linetype="dashed",
                   size=0.5) +
        geom_vline(xintercept = vacc_start,
                   linetype="dashed",
                   color="#6A3D9A",
                   size=0.5) +
        scale_y_continuous(labels = scales::percent)
        
    
    graph_name <- file.path(".", output_folder, sprintf("%s_%s_newinf_prop.png", scenario, sys_time))

    ggsave(graph_name, width = 15, height = 8, units = "cm", dpi = 300)
}