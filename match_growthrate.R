match_growthrate <- function(parameters, intro_state, total_pop, scenario, sys_time) {
    # More transmissible
    match_r0 <- function(parameters, intro_state, total_pop) {
        S <- intro_state$S
        V <- intro_state$V
        R1 <- intro_state$R1
        growth_rate <- intro_state$growth_rate2
        
        strain1 <- parameters[["strain1"]]
        strain2 <- parameters[["strain2"]]
        
        vacceff_i <- strain2[["vacceff_i"]]
        crossimm <- strain1[["crossimm"]]
        inf_duration <- strain1[["inf_duration"]]
        preinf_period <- strain1[["preinf_period"]]
        gen_time <- preinf_period + inf_duration
        
        sus2 <- (S + (1-vacceff_i)*V + (1-crossimm)*R1)/total_pop
        
        r0 <- exp(growth_rate*gen_time)/sus2
        
        return(enframe(c("preinf_period" = preinf_period, 
                         "growth_rate" = growth_rate, 
                         "r0" = r0, 
                         "sus2" = sus2, 
                         "inf_duration" = inf_duration,
                         "gen_time" = gen_time,
                         "vacceff_i" = vacceff_i, 
                         "crossimm" = crossimm), value="higher_r0"))
    }
    
    # Shorter pre-infectious period, shorter generation time
    match_preinf <- function(parameters, intro_state, total_pop) {
        S <- intro_state$S
        V <- intro_state$V
        R1 <- intro_state$R1
        growth_rate <- intro_state$growth_rate2
        
        strain1 <- parameters[["strain1"]]
        strain2 <- parameters[["strain2"]]
        
        vacceff_i <- strain2[["vacceff_i"]]
        r0 <- strain1[["r0"]]
        crossimm <- strain1[["crossimm"]]
        inf_duration <- strain1[["inf_duration"]]
        
        sus2 <- (S + (1-vacceff_i)*V + (1-crossimm)*R1)/total_pop
        
        preinf_period <- log(r0*sus2)/growth_rate - inf_duration
        
        gen_time <- preinf_period + inf_duration
        
        return(enframe(c("preinf_period" = preinf_period, 
                         "growth_rate" = growth_rate, 
                         "r0" = r0, 
                         "sus2" = sus2, 
                         "inf_duration" = inf_duration,
                         "gen_time" = gen_time,
                         "vacceff_i" = vacceff_i, 
                         "crossimm" = crossimm), value="shorter_preinf"))
    }
    
    # Reduce cross-immunity
    match_crossimm <- function(parameters, intro_state, total_pop) {
        S <- intro_state$S
        V <- intro_state$V
        R1 <- intro_state$R1
        growth_rate <- intro_state$growth_rate2
        
        strain1 <- parameters[["strain1"]]
        strain2 <- parameters[["strain2"]]
        
        vacceff_i <- strain2[["vacceff_i"]]
        r0 <- strain1[["r0"]]
        preinf_period <- strain1[["preinf_period"]]
        inf_duration <- strain1[["inf_duration"]]
        gen_time <- preinf_period + inf_duration
        
        sus2 <- exp(growth_rate*gen_time)/r0
        
        # crossimm <- 1 - ((((total_pop*exp(growth_rate*gen_time))/r0) - S - (1-vacceff_i)*V)/R1)
        crossimm <- 1 - ((total_pop*sus2) - S - (1-vacceff_i)*V)/R1
        
        return(enframe(c("preinf_period" = preinf_period, 
                         "growth_rate" = growth_rate, 
                         "r0" = r0, 
                         "sus2" = sus2, 
                         "inf_duration" = inf_duration,
                         "gen_time" = gen_time,
                         "vacceff_i" = vacceff_i, 
                         "crossimm" = crossimm), value="reduc_crossimm"))
    }
    
    # Reduce vaccine efficacy for infection
    match_vacceff_i <- function(parameters, intro_state, total_pop) {
        S <- intro_state$S
        V <- intro_state$V
        R1 <- intro_state$R1
        growth_rate <- intro_state$growth_rate2
        
        strain1 <- parameters[["strain1"]]
        
        r0 <- strain1[["r0"]]
        crossimm <- strain1[["crossimm"]]
        
        preinf_period <- strain1[["preinf_period"]]
        inf_duration <- strain1[["inf_duration"]]
        gen_time <- preinf_period + inf_duration
        
        sus2 <- exp(growth_rate*gen_time)/r0
        
        vacceff_i <- 1 - ((total_pop*sus2) - S - (1-crossimm)*R1)/V
        
        return(enframe(c("preinf_period" = preinf_period, 
                         "growth_rate" = growth_rate, 
                         "r0" = r0, 
                         "sus2" = sus2, 
                         "inf_duration" = inf_duration,
                         "gen_time" = gen_time,
                         "vacceff_i" = vacceff_i, 
                         "crossimm" = crossimm), value="reduc_vacceff"))
    }
    
    matched_r0 <- match_r0(parameters, intro_state, total_pop)
    matched_preinf <- match_preinf(parameters, intro_state, total_pop)
    matched_crossimm <- match_crossimm(parameters, intro_state, total_pop)
    matched_vacceff_i <- match_vacceff_i(parameters, intro_state, total_pop) # check
    
    matched <- list(matched_r0, matched_preinf, matched_crossimm, matched_vacceff_i) %>% reduce(left_join, by = "name")
    
    matched_file <- sprintf("%s_%s_matched.csv", scenario, sys_time)
    
    write.csv(matched, file.path(".", "seir_model_output", matched_file), row.names=F)
    
    return(matched)
}
    
    # # Longer infectious period
    # match_inf <- function(parameters, intro_state, total_pop, matched_preinf) {
    #     S <- intro_state$S
    #     V <- intro_state$V
    #     R1 <- intro_state$R1
    #     growth_rate <- intro_state$growth_rate2
    #     
    #     strain1 <- parameters[["strain1"]]
    #     
    #     r0 <- strain1[["r0"]]
    #     vacceff_i <- strain1[["vacceff_i"]]
    #     crossimm <- strain1[["crossimm"]]
    # 
    #     preinf_period <- matched_preinf$short_preinf[matched_preinf$name=="preinf_period"]
    #     
    #     sus2 <- (S + (1-vacceff_i)*V + (1-crossimm)*R1)/total_pop
    #     
    #     inf_duration <- log(r0*sus2)/growth_rate - preinf_period
    #     
    #     gen_time <- preinf_period + inf_duration
    #     
    #     return(enframe(c("preinf_period" = preinf_period, 
    #                      "growth_rate" = growth_rate, 
    #                      "r0" = r0, 
    #                      "sus2" = sus2, 
    #                      "inf_duration" = inf_duration,
    #                      "gen_time" = gen_time,
    #                      "vacceff_i" = vacceff_i, 
    #                      "crossimm" = crossimm), value="long_inf"))
    # }