
run_model <- function(dt, run_time, state_init, parameters, seed2, scenario, sys_time, output_folder) {
    
    source("plot_seir_model.R")
    
    # Create output folder if doesn't exist
    if (!dir.exists(output_folder)) {
        dir.create(output_folder)
    }
    
    # Ensure compartments defined in correct order
    state <- c(
        S = state_init[["S"]],
        E1 = state_init[["E1"]],
        I_s1 = state_init[["I_s1"]],
        I_pc1 = state_init[["I_pc1"]],
        I_c1 = state_init[["I_c1"]],
        R1 = state_init[["R1"]],
        V = state_init[["V"]],
        E_v1 = state_init[["E_v1"]],
        E_v2 = state_init[["E_v2"]],
        E2 = state_init[["E2"]],
        I_s2 = state_init[["I_s2"]],
        I_pc2 = state_init[["I_pc2"]],
        I_c2 = state_init[["I_c2"]],
        R2 = state_init[["R2"]],
        Cum_Iu_s1 = state_init[["Cum_Iu_s1"]],
        Cum_Iu_pc1 = state_init[["Cum_Iu_pc1"]],
        Cum_Iv_s1 = state_init[["Cum_Iv_s1"]],
        Cum_Iv_pc1 = state_init[["Cum_Iv_pc1"]],
        Cum_Iu_s2 = state_init[["Cum_Iu_s2"]],
        Cum_Iu_pc2 = state_init[["Cum_Iu_pc2"]],
        Cum_Iv_s2 = state_init[["Cum_Iv_s2"]],
        Cum_Iv_pc2 = state_init[["Cum_Iv_pc2"]],
        Cum_reinf = state_init[["Cum_reinf"]]
    )
    
    # Define ODEs
    seir_ode <- function(times, state, parameters) {
        
        # Define compartments
        
        S <- state[["S"]]
        
        E1 <- state[["E1"]]
        I_s1 <- state[["I_s1"]]
        I_pc1 <- state[["I_pc1"]]
        I_c1 <- state[["I_c1"]]
        R1 <- state[["R1"]]
        
        V <- state[["V"]]
        E_v1 <- state[["E_v1"]]
        E_v2 <- state[["E_v2"]]
        
        E2 <- state[["E2"]]
        I_s2 <- state[["I_s2"]]
        I_pc2 <- state[["I_pc2"]]
        I_c2 <- state[["I_c2"]]
        R2 <- state[["R2"]]
        
        Cum_Iu_s1 <- state[["Cum_Iu_s1"]]
        Cum_Iu_pc1 <- state[["Cum_Iu_pc1"]]
        Cum_Iv_s1 <- state[["Cum_Iv_s1"]]
        Cum_Iv_pc1 <- state[["Cum_Iv_pc1"]]
        Cum_Iu_s2 <- state[["Cum_Iu_s2"]]
        Cum_Iu_pc2 <- state[["Cum_Iu_pc2"]]
        Cum_Iv_s2 <- state[["Cum_Iv_s2"]]
        Cum_Iv_pc2 <- state[["Cum_Iv_pc2"]]
        Cum_reinf <- state[["Cum_reinf"]]
        
        # Extract parameters
        strain1 <- parameters[["strain1"]]
        strain2 <- parameters[["strain2"]]
        vaccination <- parameters[["vaccination"]]
        
        if (vaccination[["yes"]] == 1 & times >= vaccination[["start"]] & S>0) {
            vacc_rate <- vaccination[["daily_rate"]]
        } else {
            vacc_rate <- 0
        }
        
        # if (S>0) {
        #     vacc_rate <- vaccination[["daily_rate"]]
        # } else {
        #     
        # }
        
        ## Strain 1
        beta1 <- strain1[["beta"]]
        
        inf_rate1 <- strain1[["inf_rate"]] 
        
        prop_clin1 <- strain1[["prop_clin"]] 
        clin_rate1 <- strain1[["clin_rate"]]
        clin_rec_rate1 <- strain1[["clin_rec_rate"]]
        
        subclin_inf1 <- strain1[["subclin_inf"]]
        subclin_rec_rate1 <- strain1[["subclin_rec_rate"]]
        
        crossimm <- strain1[["crossimm"]]
        
        vacceff_i1 <- strain1[["vacceff_i"]]
        vacceff_d1 <- strain1[["vacceff_d"]]
        
        foi1 <-  beta1*((subclin_inf1*I_s1) + I_pc1 + I_c1)
        
        ## Strain 2
        beta2 <- strain2[["beta"]]
        
        inf_rate2 <- strain2[["inf_rate"]] 
        
        prop_clin2 <- strain2[["prop_clin"]] 
        clin_rate2 <- strain2[["clin_rate"]]
        clin_rec_rate2 <- strain2[["clin_rec_rate"]]
        
        subclin_inf2 <- strain2[["subclin_inf"]]
        subclin_rec_rate2 <- strain2[["subclin_rec_rate"]]
        
        vacceff_i2 <- strain2[["vacceff_i"]]
        vacceff_d2 <- strain2[["vacceff_d"]]
        
        foi2 <-  beta2*((subclin_inf2*I_s2) + I_pc2 + I_c2)
        
        # Transition equations
        dS <- -S*(foi1 + foi2) - vacc_rate
        
        dE1 <- S*foi1 - inf_rate1*E1 
        dI_s1 <- inf_rate1*(1-prop_clin1)*E1 + inf_rate1*(1-prop_clin1)*vacceff_d1*E_v1 - subclin_rec_rate1*I_s1
        dI_pc1 <- inf_rate1*prop_clin1*E1 + inf_rate1*prop_clin1*(1-vacceff_d1)*E_v1 - clin_rate1*I_pc1
        dI_c1 <- clin_rate1*I_pc1 - clin_rec_rate1*I_c1
        dR1 <- subclin_rec_rate1*I_s1 + clin_rec_rate1*I_c1 - (1-crossimm)*foi2*R1
        
        dV <- vacc_rate - (1-vacceff_i1)*foi1*V - (1-vacceff_i2)*foi2*V
        dE_v1 <- (1-vacceff_i1)*foi1*V - inf_rate1*E_v1
        dE_v2 <- (1-vacceff_i2)*foi2*V - inf_rate2*E_v2
        
        dE2 <- S*foi2 + (1-crossimm)*foi2*R1 - inf_rate2*E2 
        dI_s2 <- inf_rate2*(1-prop_clin2)*E2 + inf_rate2*(1-prop_clin2)*vacceff_d2*E_v2 - subclin_rec_rate2*I_s2
        dI_pc2 <- inf_rate2*prop_clin2*E2 + inf_rate2*prop_clin2*(1-vacceff_d2)*E_v2 - clin_rate2*I_pc2
        dI_c2 <- clin_rate2*I_pc2 - clin_rec_rate2*I_c2
        dR2 <- subclin_rec_rate2*I_s2 + clin_rec_rate2*I_c2
        
        # Cumulative measures
        dCum_Iu_s1 <- inf_rate1*(1-prop_clin1)*E1
        dCum_Iu_pc1 <- inf_rate1*prop_clin1*E1 
        dCum_Iv_s1 <- inf_rate1*(1-prop_clin1)*vacceff_d1*E_v1
        dCum_Iv_pc1 <- inf_rate1*prop_clin1*(1-vacceff_d1)*E_v1
        dCum_Iu_s2 <- inf_rate2*(1-prop_clin2)*E2
        dCum_Iu_pc2 <- inf_rate2*prop_clin2*E2 
        dCum_Iv_s2 <- inf_rate2*(1-prop_clin2)*vacceff_d2*E_v2
        dCum_Iv_pc2 <- inf_rate2*prop_clin2*(1-vacceff_d2)*E_v2
        dCum_reinf <- (1-crossimm)*foi2*R1
        
        return(list(c(dS, dE1, dI_s1, dI_pc1, dI_c1, dR1, dV, dE_v1, dE_v2, dE2, dI_s2, dI_pc2, dI_c2, dR2, dCum_Iu_s1, dCum_Iu_pc1, dCum_Iv_s1, dCum_Iv_pc1, dCum_Iu_s2, dCum_Iu_pc2, dCum_Iv_s2, dCum_Iv_pc2, dCum_reinf)))
    }
    
    # Specify run time
    times <- seq(0, run_time, dt)
    
    # Run model
    output <- as.data.frame(ode(y = state, times = times, func = seir_ode, parms = parameters, events = list(data = seed2)))
    
    # Save output
    
    output_file <- sprintf("%s_%s_output.csv", scenario, sys_time)
    
    write.csv(output, file.path(".", output_folder, output_file), row.names=F)
    
    # Save model input
    strain1_parameters <- data.frame(parameter = names(parameters[["strain1"]]), value = parameters[["strain1"]], set = rep("strain1", length(parameters[["strain1"]])), row.names=NULL)
    
    strain2_parameters <- data.frame(parameter = names(parameters[["strain2"]]), value = parameters[["strain2"]], set = rep("strain2", length(parameters[["strain2"]])), row.names=NULL)
    
    vacc_parameters <- data.frame(parameter = names(parameters[["vaccination"]]), value = parameters[["vaccination"]], set = rep("vaccination", length(parameters[["vaccination"]])), row.names=NULL)
    
    state_parameters <- data.frame(parameter = names(state_init), value = state_init, set = rep("state", length(state_init)), row.names=NULL)
    
    population_parameters <- data.frame(parameter = names(parameters[["population"]]), value = parameters[["population"]], set = rep("population", length(parameters[["population"]])), row.names=NULL)
    
    seed_parameters <- seed2 %>% 
        select(var, value) %>% 
        rbind(c("seed_time", max(seed2$time))) %>%
        mutate(value=as.numeric(value)) %>% 
        mutate(set = "seed")
    
    names(seed_parameters) <- c("parameter", "value", "set")
    
    time_parameters <- data.frame(parameter = c("dt", "run_time"), value = c(dt, run_time), set = rep("time", 2))
    
    inputs <- bind_rows(strain1_parameters, strain2_parameters, vacc_parameters, population_parameters, seed_parameters, time_parameters) 
    
    input_file <- sprintf("%s_%s_input.csv", scenario, sys_time)
    
    write.csv(inputs, file.path(".", output_folder, input_file), row.names=F)
    
    # Plot model
    plot_seir_model(output, run_time, seed2$time[1], vaccination[["start"]], output_folder, scenario, sys_time, "full")
    
    # Return table 
    return(output)
}