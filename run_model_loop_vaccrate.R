library(deSolve)
library(tidyverse)
# library(RColorBrewer)

source("seir_model.R")
source("match_growthrate.R")
source("plot_seir_model.R")
source("plot_individual_function.R")

output_folder <- "seir_model_output"
if (!dir.exists(output_folder)) {
    dir.create(output_folder)
}

theme_set(
    theme_gray(base_size = 12)
)

r0_strain1_range <- c(1.5, 3)
transmiss_range <- seq(1, 2, 0.25)
crossimm_range <- c(0.5, 0.75, 0.95, 1)
seed_range <- c(0, 14, 35, 70)
vacc_range <- c(0, 1)
vacc_start_range <- c(0, 35)

combos <- crossing(r0_strain1_range, transmiss_range, crossimm_range, seed_range, vacc_range, vacc_start_range) %>% 
    filter(!(vacc_range==0 & vacc_start_range==35))

#### Running missing models
# combos <- as.data.frame(missing) %>% 
#     mutate(
#         r0_strain1_range = str_extract(missing, "(?<=r)\\d*\\.?\\d*"),
#         transmiss_range = str_extract(missing, "(?<=transmiss)\\d*\\.?\\d*"),
#         crossimm_range = str_extract(missing, "(?<=crossimm)\\d*\\.?\\d*"),
#         seed_range = str_extract(missing, "(?<=seed)\\d*"),
#         vacc_range = str_extract(missing, "(?<=vacc)\\d*\\.?\\d*"),
#         vacc_start_range = str_extract(missing, "(?<=start_)\\d*\\.?\\d*")
#     ) %>% 
#     mutate(across(!missing, as.numeric))

for (i in 1:nrow(combos)) {
    
    # Unique scenario name
    scenario <- sprintf("r%s_transmiss%s_crossimm%s_seed%s_vacc%s_start_%s", combos$r0_strain1_range[i], combos$transmiss_range[i], combos$crossimm_range[i], combos$seed_range[i], combos$vacc_range[i], combos$vacc_start_range[i]) %>% str_replace_all("\\.", "pt")
    
    # PARAMETERS
    population <- c(total_pop = 60000000)
    
    initial_inf <- 1000
    
    ## Vaccination
    vaccination <- c(
        yes = combos$vacc_range[i],
        init = 0,
        weekly_rate = 1500000,
        start = combos$vacc_start_range[i]
    )
    
    vaccination[["daily_rate"]] <- vaccination[["yes"]]*(vaccination[["weekly_rate"]]/7)
    
    
    ## Strain 1
    strain1 <- c(
        preinf_period = 2.5,
        prop_clin = 0.45,
        preclin_period = 2.5,
        clin_period = 2.5,
        subclin_inf = 0.5,
        r0 = combos$r0_strain1_range[i],
        crossimm = combos$crossimm_range[i],
        vacceff_i = 0.85*vaccination[["yes"]],
        vacceff_d = 0.9*vaccination[["yes"]]
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
    
    strain2 <- c(
        preinf_period = strain1[["preinf_period"]],
        # preinf_period = 2.4789320,
        prop_clin = strain1[["prop_clin"]],
        preclin_period = strain1[["preclin_period"]],
        clin_period = strain1[["clin_period"]],
        subclin_inf = strain1[["subclin_inf"]],
        r0 = strain1[["r0"]]*transmiss,
        vacceff_i = 0.75*vaccination[["yes"]],
        vacceff_d = 0.8*vaccination[["yes"]]
    )
    
    strain2[["inf_rate"]] <- 1/strain2[["preinf_period"]]
    strain2[["clin_rate"]] <- 1/strain2[["preclin_period"]]
    strain2[["clin_rec_rate"]] <- 1/strain2[["clin_period"]]
    strain2[["inf_duration"]] <- strain2[["preclin_period"]] + strain2[["clin_period"]]
    strain2[["subclin_period"]] <- strain2[["inf_duration"]]
    strain2[["subclin_rec_rate"]] <- 1/strain2[["subclin_period"]]
    strain2[["gen_time"]] <- strain2[["preinf_period"]] + strain2[["inf_duration"]]
    strain2[["beta"]] <- strain2[["r0"]]/(population[["total_pop"]]*strain2[["inf_duration"]])
    
    seed2 <- data.frame(var=c("I_s2", "I_pc2", "I_c2"), time = rep(seedtime2, 3), value = rep(initial_inf, 3), method = rep("add", 3))
    
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
    run_time <- 1200
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
    
    # Calculate daily new infections
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
    
    # Capture peak
    peak_new_I <- max(new_inf_table$daily_new_I)
    peak_time <- new_inf_table$day[new_inf_table$daily_new_I == peak_new_I]
    
    # Save end situation
    ## Create end file
    end_file <- file.path(".", output_folder, "end.csv")
    
    if (!file.exists(end_file)) {
        
        end <- table[1,] %>%
            mutate(
                scenario = NA,
                sys_time = NA,
                .before = time) %>% 
            mutate(
                peak_new_I = NA,
                peak_time = NA,
                growth_rate1_intro = NA,
                growth_rate2_intro = NA,
                higher_r0 = NA,
                shorter_preinf = NA,
                shorter_preinf_gentime = NA,
                reduc_crossimm = NA,
                reduc_vacceff = NA)
        
        # Empty table
        end <- end[0,]
        
        write.csv(end, end_file, row.names=F)
    }
    
    end_day <- new_inf_table %>% 
        filter(day > peak_time & daily_new_I < 1) %>%
        slice(1) %>% 
        select(day) %>% 
        as.numeric()
    
    if (is.na(end_day)) {
        end_day <- floor(max(table$time))
    }
    
    end_state <- table %>% 
        filter(time==end_day) %>% 
        mutate(
            scenario = scenario,
            sys_time = sys_time,
            .before = time) %>% 
        mutate(
            peak_new_I = peak_new_I,
            peak_time = peak_time,
            growth_rate1_intro = intro_state$growth_rate1,
            growth_rate2_intro = intro_state$growth_rate2,
            higher_r0 = match_growth$higher_r0[match_growth$name=="r0"],
            shorter_preinf = match_growth$shorter_preinf[match_growth$name=="preinf_period"],
            shorter_preinf_gentime = match_growth$shorter_preinf[match_growth$name=="gen_time"],
            reduc_crossimm = match_growth$reduc_crossimm[match_growth$name=="crossimm"],
            reduc_vacceff = match_growth$reduc_vacceff[match_growth$name=="vacceff_i"])
    
    # Append to end file
    write.table(end_state, file=end_file, append=T, quote=F, sep=",", row.names=F, col.names=F)
    
    # Plot seir model
    plot_seir_model(table, end_day, seedtime2, vaccination[["start"]], output_folder, scenario, sys_time, "trim")
    
    # Make graphs
    plot_graphs(output_folder, scenario, seedtime2, sys_time, vaccination[["start"]], end_day)
}