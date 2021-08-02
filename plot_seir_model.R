plot_seir_model <- function(output, plot_time, seedtime2, vacc_start, output_folder, scenario, sys_time, plot_id){
    # Plot compartments
        plot_states <- output %>% 
            mutate(I1 = I_s1 + I_pc1 + I_c1,
                   I2 = I_s2 + I_pc2 + I_c2) %>% 
            select(c(time, S, V, E1, E_v1, E_v2, E2, I1, I2, R1, R2)) %>%
            pivot_longer(cols=!time, names_to = "compartment", values_to = "count") %>%
            mutate(grid = case_when(
                compartment %in% c("S", "V") ~ "A",
                compartment %in% c("R1", "R2") ~ "B",
                compartment %in% c("E1", "E2", "I1", "I2") ~ "C",
                compartment %in% c("E_v1", "E_v2") ~ "D")) %>%
            filter(time %in% seq(0, plot_time, 1)) %>% 
            ggplot(aes(x=time, y=count, colour=compartment)) +
            geom_line() +
            facet_grid(grid ~., 
                       scale="free") +
            theme(strip.background = element_blank(),
                  strip.text.y = element_blank()
            ) +
            scale_colour_manual(values = c(
                "S" = "#B15928",
                "V" = "#6A3D9A",
                "E1" = "#A6CEE3",
                "E2" = "#1F78B4",
                "E_v1" = "#B2DF8A",
                "E_v2" = "#33A02C",
                "I1" = "#FB9A99",
                "I2" = "#E31A1C",
                "R1" = "#FDBF6F",
                "R2" = "#FF7F00")
            ) +
            scale_y_continuous(labels = scales::label_comma()) +
            geom_vline(xintercept = seedtime2,
                       linetype="dashed",
                       size=0.5) +
            geom_vline(xintercept = vacc_start,
                       linetype="dashed",
                       color="#6A3D9A",
                       size=0.5)
        
        plot_name <- file.path(".", output_folder, sprintf("%s_%s_plot_%s.png", scenario, sys_time, plot_id))
        
        ggsave(plot_name, plot_states, width = 15, height = 12, units = "cm", dpi = 300)
}