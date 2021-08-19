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
    theme_gray(base_size = 11)
)

##### !!!! Check match_set!!!!!
r0_strain1_range <- c(4)
transmiss_range <- c(1)
crossimm_range <- c(0.25, 0.75, 1)
# seed_range <- c(0, 65.82, 87.18, 102.24)
# seed_range <- c(0, 33.42, 41.72, 47.6)
seed_range <- c(0, 21.44, 26.32, 29.78)
preinf_range <- c(0.5, 1.5)
vacc_range <- c(0)
vacc_start_range <- c(0)
# vacc_start_range <- c(0.1667, 0.3334)
# vacc_start_range <- c(0.3, 0.6)
# vacc_start_range <- c(0.375, 0.75)


combos_specify <- crossing(r0_strain1_range, transmiss_range, crossimm_range, seed_range, preinf_range, vacc_range, vacc_start_range) 

## Save combos
# combos <- combos_specify[0,]
# write.csv(combos, file=file.path(".", output_folder, "combos.csv"), row.names=F)

combos <- rbind(combos, combos_specify)

# Append new combos
write.table(combos, file=file.path(".", output_folder, "combos.csv"), append=T, quote=F, sep=",", row.names=F, col.names=F)

# Read combos
#  <- read.csv(file.path(".", output_folder, "combos.csv"))

#### Running missing models ----
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

#### Run ----
for (i in 1:nrow(combos)) {
    
    # Unique scenario name
    scenario <- sprintf("r%s_transmiss%s_crossimm%s_seed%s_preinf%s_vacc%s_start_%s", combos$r0_strain1_range[i], round(combos$transmiss_range[i],2), combos$crossimm_range[i], combos$seed_range[i], combos$preinf_range[i], combos$vacc_range[i], combos$vacc_start_range[i]) %>% str_replace_all("\\.", "pt")
    
    # match_set <- sprintf("r%s_transmiss%s_crossimm%s_seed%s_preinf%s_vacc0_start_0", combos$r0_strain1_range[i], round(combos$transmiss_range[i],2), combos$crossimm_range[i], combos$seed_range[i], combos$preinf_range[i]) %>% str_replace_all("\\.", "pt")
    match_set <- combos$match_set[i]
    # match_set <- scenario
    
    # TIME
    if (combos$r0_strain1_range[i]==1.5 & combos$vacc_start_range[i]==0.3334) {
        run_time <- 1500
    } else {
        run_time <- 800
    }
    
    dt <- 0.02
    
    sys_time <- format(Sys.time(), "%d-%m-%Y--%H-%M")
    
    seedtime2 <- combos$seed_range[i]
    
    if (seedtime2 > run_time) {
        seedtime2 <- run_time
    }
    
    # PARAMETERS
    population <- c(total_pop = 1000000)
    
    initial_inf <- 10 # Initial seeding of 20 (0.002% population)
    
    ## Vaccination
    vaccination <- c(
        yes = combos$vacc_range[i],
        init = population[["total_pop"]]*combos$vacc_start_range[i],
        weekly_rate = 0,
        start = 0
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
        vacceff_i = 0.8*vaccination[["yes"]],
        vacceff_d = 0.9*vaccination[["yes"]]
    )
    
    strain1[["inf_rate"]] <- 1/strain1[["preinf_period"]]
    strain1[["clin_rate"]] <- 1/strain1[["preclin_period"]]
    strain1[["clin_rec_rate"]] <- 1/strain1[["clin_period"]]
    strain1[["inf_duration"]] <- strain1[["preclin_period"]] + strain1[["clin_period"]]
    strain1[["subclin_period"]] <- strain1[["inf_duration"]]
    strain1[["subclin_rec_rate"]] <- 1/strain1[["subclin_period"]]
    strain1[["gen_time"]] <- strain1[["preinf_period"]] + strain1[["inf_duration"]]
    strain1[["r0_scalingfactor"]] <- ((1-strain1[["prop_clin"]])*strain1[["subclin_inf"]]) + strain1[["prop_clin"]]
    strain1[["beta"]] <- (strain1[["r0"]]/strain1[["r0_scalingfactor"]])/(population[["total_pop"]]*strain1[["inf_duration"]])
    strain1[["vaccrisk_i"]] <- 1 - strain1[["vacceff_i"]]
    strain1[["vaccrisk_d_i"]] <- (1 - strain1[["vacceff_d"]])/strain1[["vaccrisk_i"]]
    
    ## Strain 2
    transmiss = combos$transmiss_range[i]
    
    strain2 <- c(
        preinf_period = combos$preinf_range[i],
        prop_clin = strain1[["prop_clin"]],
        preclin_period = strain1[["preclin_period"]],
        clin_period = strain1[["clin_period"]],
        subclin_inf = strain1[["subclin_inf"]],
        r0 = strain1[["r0"]]*transmiss,
        vacceff_i = 0.7*vaccination[["yes"]],
        vacceff_d = 0.8*vaccination[["yes"]]
    )
    
    strain2[["inf_rate"]] <- 1/strain2[["preinf_period"]]
    strain2[["clin_rate"]] <- 1/strain2[["preclin_period"]]
    strain2[["clin_rec_rate"]] <- 1/strain2[["clin_period"]]
    strain2[["inf_duration"]] <- strain2[["preclin_period"]] + strain2[["clin_period"]]
    strain2[["subclin_period"]] <- strain2[["inf_duration"]]
    strain2[["subclin_rec_rate"]] <- 1/strain2[["subclin_period"]]
    strain2[["gen_time"]] <- strain2[["preinf_period"]] + strain2[["inf_duration"]]
    strain2[["r0_scalingfactor"]] <- ((1-strain2[["prop_clin"]])*strain2[["subclin_inf"]]) + strain2[["prop_clin"]]
    strain2[["beta"]] <- (strain2[["r0"]]/strain2[["r0_scalingfactor"]])/(population[["total_pop"]]*strain2[["inf_duration"]])
    strain2[["vaccrisk_i"]] <- 1 - strain2[["vacceff_i"]]
    strain2[["vaccrisk_d_i"]] <- (1 - strain2[["vacceff_d"]])/strain2[["vaccrisk_i"]]
    
    seed2 <- data.frame(var=c("I_s2", "I_pc2"), time = rep(seedtime2, 2), value = rep(initial_inf, 2), method = rep("add", 2))
    
    # STATES
    state <- c(S = NA,
               E1 = 0,
               I_s1 = initial_inf,
               I_pc1 = initial_inf,
               I_c1 = 0,
               R1 = 0,
               V = vaccination[["init"]],
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
            I = I1 + I2,
            # Cumulative new infectious
            Cum_Inf1 = Cum_Iu_s1 + Cum_Iu_pc1 + Cum_Iv_s1 + Cum_Iv_pc1,
            Cum_Inf2 = Cum_Iu_s2 + Cum_Iu_pc2 + Cum_Iv_s2 + Cum_Iv_pc2,
            Cum_Infv = Cum_Iv_s1 + Cum_Iv_pc1 + Cum_Iv_s2 + Cum_Iv_pc2) %>%
        mutate(
            # New infectious per dt time step
            new_I1 = Cum_Inf1-lag(Cum_Inf1, default=0),
            new_I2 = Cum_Inf2-lag(Cum_Inf2, default=0),
            new_Iv = Cum_Infv-lag(Cum_Infv, default=0)) %>% 
        mutate(
            new_I = new_I1 + new_I2,
            Cum_Inf = Cum_Inf1 + Cum_Inf2) %>%
        mutate(
            # Cumulative percentage of cases
            Cumprop_Inf = Cum_Inf/sum(new_I)) %>% 
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
            prop_strain2_only = strain2_only/allinf,
            # Proportion of population
            proppop_strain1_only = strain1_only/population[["total_pop"]],
            proppop_bothstrains = bothstrains/population[["total_pop"]],
            proppop_strain2_only = strain2_only/population[["total_pop"]],
            proppop_allinf = allinf/population[["total_pop"]]
        ) %>% 
        mutate(
            cumprop_allinf = allinf/last(allinf)) %>% 
        mutate(
            # Prevalence of cases (infectious)
            prevalence_I1 = I1/population[["total_pop"]],
            prevalence_I2 = I2/population[["total_pop"]],
            prevalence_I = (I1+I2)/population[["total_pop"]],
            # Proportion of population
            proppop_S = S/population[["total_pop"]],
            proppop_V = V/population[["total_pop"]],
            proppop_R1 = R1/population[["total_pop"]],
            proppop_R2 = R2/population[["total_pop"]]
        )

    
    # Calculate daily new infections
    new_inf_table <- table %>% 
        mutate(day = as.integer(time)) %>% 
        group_by(day) %>% 
        summarise(daily_new_I1 = sum(new_I1),
                  daily_new_I2 = sum(new_I2),
                  daily_new_I = sum(new_I)) %>% 
        mutate(prop_new_I1 = daily_new_I1/daily_new_I,
               prop_new_I2 = daily_new_I2/daily_new_I,
               # Proportion population
               proppop_new_I1 = daily_new_I1/population[["total_pop"]],
               proppop_new_I2 = daily_new_I2/population[["total_pop"]],
               proppop_new_I = daily_new_I/population[["total_pop"]])
    

    # Capture peak
    peak_new_I <- max(new_inf_table$daily_new_I)
    peak_time <- new_inf_table$day[new_inf_table$daily_new_I == peak_new_I]
    
    peak_new_I1 <- max(new_inf_table$daily_new_I1)
    peak_time_new_I1 <- new_inf_table$day[new_inf_table$daily_new_I1 == peak_new_I1]
    
    peak_new_I2 <- max(new_inf_table$daily_new_I2)
    peak_time_new_I2 <- new_inf_table$day[new_inf_table$daily_new_I2 == peak_new_I2][1]
    
    # Proportion of peak
    new_inf_table <- new_inf_table %>% 
        mutate(
            prop_peak_new_I1 = daily_new_I1/peak_new_I1,
            prop_peak_new_I2 = daily_new_I2/peak_new_I2,
            prop_peak_new_I = daily_new_I/peak_new_I)
    
    # Save daily inf file
    dailyinf_file <- sprintf("%s_%s_dailyinf.csv", scenario, sys_time)
    
    write.csv(new_inf_table, file.path(".", output_folder, dailyinf_file), row.names=F)
    
    
    ## Peak prevalence of infectious cases
    peakprev_I <- max(table$I)
    proppop_peak_new_I <- peak_new_I/population[["total_pop"]]
    peakprev_I_time <- table$time[table$I == peakprev_I]
    
    peakprev_I1 <- max(table$I1)
    proppop_peak_new_I1 <- peak_new_I1/population[["total_pop"]]
    peakprev_I1_time <- table$time[table$I1 == peakprev_I1]
    
    peakprev_I2 <- max(table$I2)
    proppop_peak_new_I2 <- peak_new_I2/population[["total_pop"]]
    peakprev_I2_time <- table$time[table$I2 == peakprev_I2][1]
    
    ## Proportion of peak prev
    table <- table %>% 
        mutate(
            prop_peakprev_I = I/peakprev_I,
            prop_peakprev_I1 = I1/peakprev_I1,
            prop_peakprev_I2 = I2/peakprev_I2)
    
    # Save calculations
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
                .before = time) %>% 
            mutate(
                peak_new_I = NA,
                proppop_peak_new_I = NA,
                peak_time = NA,
                peak_new_I1 = NA,
                proppop_peak_new_I1 = NA,
                peak_time_new_I1 = NA,
                peak_new_I2 = NA,
                proppop_peak_new_I2 = NA,
                peak_time_new_I2 = NA,
                proppop_Cum_Inf = NA,
                peakprev_I = NA,
                peakprev_I_time = NA,
                peakprev_I1 = NA,
                peakprev_I1_time = NA,
                peakprev_I2 = NA,
                peakprev_I2_time = NA,
                growth_rate1_intro = NA,
                growth_rate2_intro = NA,
                sus2_intro = NA,
                higher_r0 = NA,
                higher_r0_rel = NA,
                higher_r0_preinf = NA,
                shorter_preinf = NA,
                shorter_preinf_gentime = NA,
                match_set = NA)
        
        # Empty table
        end <- end[0,]
        
        write.csv(end, end_file, row.names=F)
    }
    
    # end_day <- new_inf_table %>% 
    #     filter(day > peak_time & daily_new_I < 1) %>%
    #     slice(1) %>% 
    #     select(day) %>% 
    #     as.numeric()

    end_day <- table %>%
        filter(time > peak_time & I < 1) %>%
        slice(1) %>%
        select(time) %>%
        as.numeric() %>% 
        floor()
    
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
            proppop_peak_new_I = proppop_peak_new_I,
            peak_time = peak_time,
            peak_new_I1 = peak_new_I1,
            proppop_peak_new_I1 = proppop_peak_new_I1,
            peak_time_new_I1 = peak_time_new_I1,
            peak_new_I2 = peak_new_I2,
            proppop_peak_new_I2 = proppop_peak_new_I2,
            peak_time_new_I2 = peak_time_new_I2,
            proppop_Cum_Inf = Cum_Inf/population[["total_pop"]],
            peakprev_I = peakprev_I,
            peakprev_I_time = peakprev_I_time,
            peakprev_I1 = peakprev_I1,
            peakprev_I1_time = peakprev_I1_time,
            peakprev_I2 = peakprev_I2,
            peakprev_I2_time = peakprev_I2_time,
            growth_rate1_intro = intro_state$growth_rate1,
            growth_rate2_intro = intro_state$growth_rate2,
            sus2_intro = intro_state$sus2,
            higher_r0 = match_growth$higher_r0[match_growth$name=="r0"],
            higher_r0_rel = match_growth$higher_r0[match_growth$name=="r0"]/strain2[["r0"]],
            higher_r0_preinf = match_growth$higher_r0[match_growth$name=="preinf_period"],
            shorter_preinf = match_growth$shorter_preinf[match_growth$name=="preinf_period"],
            shorter_preinf_gentime = match_growth$shorter_preinf[match_growth$name=="gen_time"],
            match_set = match_set)
    
    # Append to end file
    write.table(end_state, file=end_file, append=T, quote=F, sep=",", row.names=F, col.names=F)
    
    # Plot seir model
    plot_seir_model(table, end_day, seedtime2, vaccination[["start"]], output_folder, scenario, sys_time, "trim")
    
    # Make graphs
    plot_graphs(output_folder, scenario, seedtime2, sys_time, vaccination[["start"]], end_day)
}