library(deSolve)
library(tidyverse)
library(RColorBrewer)

source("seir_model.R")
source("match_growthrate.R")
source("plot_individual_function.R")

output_folder <- "seir_model_output"
if (!dir.exists(output_folder)) {
    dir.create(output_folder)
}

theme_set(
    theme_gray(base_size = 12)
)

transmiss_range <- seq(1, 2, 0.25)
crossimm_range <- seq(0.25, 1, 0.25)
seed_range <- c(0, 14, 28, 70)
vacc_range <- c(0, 0.05, 0.1, seq(0.2, 0.8, 0.2))
combos <- crossing(transmiss_range, crossimm_range, seed_range, vacc_range)

for (i in 1:nrow(combos)){
    
    # Unique scenario name
    scenario <- sprintf("transmiss%s_crossimm%s_seed%s_vacc%s", combos$transmiss_range[i], combos$crossimm_range[i], combos$seed_range[i], combos$vacc_range[i]) %>% str_replace_all("\\.", "pt")
    
    # PARAMETERS
    population <- c(
        total_pop = 60000000)
    
    initial_inf <- 1000
    
    ## Vaccination
    vaccination <- c(
        yes = 0,
        init = population[["total_pop"]]*combos$vacc_range[i],
        weekly_rate = 2000000
    )
    vaccination[["daily_rate"]] <- vaccination[["yes"]]*(vaccination[["weekly_rate"]]/7)
    
    
    ## Strain 1
    strain1 <- c(
        preinf_period = 2.5,
        prop_clin = 0.45,
        preclin_period = 2.5,
        clin_period = 2.5,
        subclin_inf = 0.5,
        r0 = 3,
        crossimm = combos$crossimm_range[i],
        vacceff_i = 0.8*vaccination[["yes"]],
        vacceff_d = 0.85*vaccination[["yes"]]
    )
    
    strain1[["inf_rate"]] <- 1/strain1[["preinf_period"]]
    strain1[["clin_rate"]] <- 1/strain1[["preclin_period"]]
    strain1[["clin_rec_rate"]] <- 1/strain1[["clin_period"]]
    strain1[["inf_duration"]] <- strain1[["preclin_period"]] + strain1[["clin_period"]]
    strain1[["subclin_period"]] <- strain1[["inf_duration"]]
    strain1[["subclin_rec_rate"]] <- 1/strain1[["subclin_period"]]
    strain1[["gen_time"]] <- strain1[["preinf_period"]] + strain1[["inf_duration"]]
    strain1[["beta"]] <- strain1[["r0"]]/(population[["total_pop"]]*strain1[["inf_duration"]])
    
    ## Strain 2
    seedtime2 <- combos$seed_range[i]
    transmiss = combos$transmiss_range[i]
    vacceff_reduc_i = 1
    vacceff_reduc_d = 1
    
    strain2 <- c(
        preinf_period = strain1[["preinf_period"]],
        # preinf_period = 2.4789320,
        prop_clin = strain1[["prop_clin"]],
        preclin_period = strain1[["preclin_period"]],
        clin_period = strain1[["clin_period"]],
        subclin_inf = strain1[["subclin_inf"]],
        r0 = strain1[["r0"]]*transmiss,
        vacceff_i = strain1[["vacceff_i"]]*vacceff_reduc_i,
        vacceff_d = strain1[["vacceff_d"]]*vacceff_reduc_d
    )
    
    strain2[["inf_rate"]] <- 1/strain2[["preinf_period"]]
    strain2[["clin_rate"]] <- 1/strain2[["preclin_period"]]
    strain2[["clin_rec_rate"]] <- 1/strain2[["clin_period"]]
    strain2[["inf_duration"]] <- strain2[["preclin_period"]] + strain2[["clin_period"]]
    strain2[["subclin_period"]] <- strain2[["inf_duration"]]
    strain2[["subclin_rec_rate"]] <- 1/strain2[["subclin_period"]]
    strain2[["gen_time"]] <- strain2[["preinf_period"]] + strain2[["inf_duration"]]
    strain2[["beta"]] <- strain2[["r0"]]/(population[["total_pop"]]*strain2[["inf_duration"]])
    
    seed2 <- data.frame(var=c("I_s2", "I_pc2", "I_c2", "V", "S"), time = rep(seedtime2, 5), value = c(rep(initial_inf, 3), vaccination[["init"]], -vaccination[["init"]]), method = rep("add", 5))
    
    # Strain 1 only
    # seed2 <- data.frame(var=c("I_s2", "I_pc2", "I_c2"), time = rep(seedtime2, 3), value = rep(0, 3), method = rep("add", 3))
    
    # STATES
    state <- c(S = NA,
               E1 = 0,
               I_s1 = initial_inf,
               I_pc1 = initial_inf,
               I_c1 = initial_inf,
               R1 = 0,
               V = 0,
               E_v1 = 0,
               E_v2 = 0,
               E2 = 0,
               I_s2 = 0,
               I_pc2 = 0,
               I_c2 = 0,
               R2 = 0,
               Cum_Iu_s1 = 0,
               Cum_Iu_pc1 = 0,
               Cum_Iv_s1 = 0,
               Cum_Iv_pc1 = 0,
               Cum_Iu_s2 = 0,
               Cum_Iu_pc2 = 0,
               Cum_Iv_s2 = 0,
               Cum_Iv_pc2 = 0,
               Cum_reinf = 0)
    
    state[["S"]] <- population[["total_pop"]] - sum(state[names(state) != "S"])
    
    # Create parameter list
    parameters <- list(strain1 = strain1, strain2 = strain2, vaccination = vaccination, population = population)
    
    # TIME
    run_time <- 500
    dt <- 0.02
    
    sys_time <- format(Sys.time(), "%d-%m-%Y--%H-%M")
    
    # RUN MODEL
    table <- run_model(dt, run_time, state, parameters, seed2, scenario, sys_time, output_folder)
    
    # Calculate growth rates
    table <- table %>% 
        mutate(sus1 = (S + (1-strain1[["vacceff_i"]])*V)/population[["total_pop"]],
               sus2 = (S + (1-strain2[["vacceff_i"]])*V + (1-strain1[["crossimm"]])*R1)/population[["total_pop"]]) %>% 
        mutate(growth_rate1 = log(strain1[["r0"]]*sus1)/strain1[["gen_time"]],
               growth_rate2 = log(strain2[["r0"]]*sus2)/strain2[["gen_time"]])
    
    ## Extract situation at Strain 2 introduction
    intro_state <- table %>% 
        filter(time==seedtime2)
    
    ## Calculate matched parameters
    match_growth <- match_growthrate(parameters, intro_state, population[["total_pop"]], scenario, sys_time)
    
    # Calculations
    table <- table %>%
        mutate(
            I1 = I_s1 + I_pc1 + I_c1,
            I2 = I_s2 + I_pc2 + I_c2,
            # Cumulative new infectious
            Cum_Inf1 = Cum_Iu_s1 + Cum_Iu_pc1 + Cum_Iv_s1 + Cum_Iv_pc1,
            Cum_Inf2 = Cum_Iu_s2 + Cum_Iu_pc2 + Cum_Iv_s2 + Cum_Iv_pc2,
            Cum_Infv = Cum_Iv_s1 + Cum_Iv_pc1 + Cum_Iv_s2 + Cum_Iv_pc2) %>%
        mutate(
            # New infectious 
            new_I1 = Cum_Inf1-lag(Cum_Inf1, default=0),
            new_I2 = Cum_Inf2-lag(Cum_Inf2, default=0),
            new_Iv = Cum_Infv-lag(Cum_Infv, default=0)) %>% 
        mutate(
            new_I = new_I1 + new_I2,
            Cum_Inf = Cum_Inf1 + Cum_Inf2) %>%
        mutate(
            # Proportion of new infectious of each strain
            prop_new_I1 = new_I1/new_I,
            prop_new_I2 = new_I2/new_I,
            prop_new_Iv = new_Iv/new_I,
            prop_Cum_Infv = Cum_Infv/Cum_Inf) %>% 
        mutate(
            # Cumulative 
            strain1_only = R1,
            bothstrains = Cum_reinf,
            strain2_only = R2 - Cum_reinf,
            # allinf_total = population[["total_pop"]]-S-V, 
            allinf = R1+R2 # Slightly smaller than allinf_total, still some infection of R1 by strain 2 so still individuals in E2, I2
            ) %>% 
        mutate(
            # Proportion cumulative
            prop_strain1_only = strain1_only/allinf,
            prop_bothstrains = bothstrains/allinf,
            prop_strain2_only = strain2_only/allinf
        )
    
    calc_file <- sprintf("%s_%s_calc.csv", scenario, sys_time)
    
    write.csv(table, file.path(".", output_folder, calc_file), row.names=F)
    
    # Save end situation
    ## Create end file
    end_file <- file.path(".", output_folder, "end.csv")
    
    if (!file.exists(end_file)) {
        
        end <- table[1,] %>%
            mutate(
                scenario = NA,
                sys_time = NA,
                .before = time)
        
        # Empty table
        end <- end[0,]
        
        write.csv(end, end_file, row.names=F)
    }
    
    # Capture end state
    end_state <- table %>% 
        filter(time==run_time) %>% 
        mutate(
            scenario = scenario,
            sys_time = sys_time,
            .before = time)
    
    # Append to end file
    write.table(end_state, file=end_file, append=T, quote=F, sep=",", row.names=F, col.names=F)
    
    # Make graphs
    plot_graphs(output_folder, calc_file, scenario, seedtime2, sys_time)
}