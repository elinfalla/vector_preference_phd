############################################################################
########### FUNCTIONS USED TO COMPARE DONNELLY AND MADDEN MODELS ###########
############################################################################

vary_param <- function(init_states, times, parms, varied_parm_name, varied_parm_vals, func, model) {
  
  # function to run a model ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of final incidence, I, of all runs, along with values of varied parameter
  
  # initialise output data frame
  output <- data.frame(parm_name = rep(varied_parm_name, length(varied_parm_vals)),
                       parm_val = varied_parm_vals,
                       final_I = rep(NA, length(varied_parm_vals)),
                       model = rep(model, length(varied_parm_vals))
  )
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states, 
                                 times = times, 
                                 parms = parms, 
                                 func = func))
    
    if ("I" %in% names(trajectory)) {
      names(trajectory)[which(names(trajectory) == "I")] <- "i"
    }
    
    output[run, "final_I"] <- round(trajectory[nrow(trajectory), "i"], 3)
  }
  
  return(output)
}

sensitivity_analysis <- function(parms_mad, parms_don, init_states_don, init_states_mad, madden_func, 
                                 madden_cunniffe_func, donnelly_func, times, default_parms) {
  
  sen_analysis_df <- data.frame(parm_name = character(),
                                parm_val = numeric(),
                                final_I = numeric(),
                                model = character())
  
  print("Starting Madden sensitivity analysis")
  
  for (i in 1:length(parms_mad)) {
    
    output_mad <- vary_param(init_states_mad, times, default_parms,
                             varied_parm_name = names(parms_mad[i]),
                             varied_parm_vals = parms_mad[[i]],
                             func = madden_func,
                             model = "Madden") 
    
    output_mad_cun <- vary_param(init_states_mad, times, default_parms,
                                 varied_parm_name = names(parms_mad[i]),
                                 varied_parm_vals = parms_mad[[i]],
                                 func = madden_cunniffe_func,
                                 model = "Madden_Cunniffe") 
    
    sen_analysis_df <- rbind(sen_analysis_df, output_mad, output_mad_cun)
    
  }
  
  print("Starting Donnelly sensitivity analysis")
  
  for (i in 1:length(parms_don)) {
    
    output_don <- vary_param(init_states_don, times, default_parms,
                             varied_parm_name = names(parms_don[i]),
                             varied_parm_vals = parms_don[[i]],
                             func = donnelly_func,
                             model = "Donnelly") 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output_don)
    
  }
  
  return(sen_analysis_df) 
}

find_theta_phi_vals <- function(theta_phi_data, init_states_don, init_states_mad, madden_func, madden_cunniffe_func,
                                donnelly_func, times, parms) {
  
  # function to find point at which phi and theta value the final proportion of infected plants is 0.5. PHI AND THETA
  # DATA MUST BE IN FORMAT OF PROPOTION OF INFECTED PLANTS (rather than number)
  
  phi_data <- theta_phi_data %>% filter(parm_name == "phi")
  theta_data <- theta_phi_data %>% filter(parm_name == "theta")
  
  phi_index <- match(TRUE, phi_data$final_I > 0.5) # finds first instance of final_I > 0.5
  theta_index <- match(TRUE, theta_data$final_I > 0.5)
  
  if (is.na(phi_index)) {
    if (is.na(theta_index)) {
      stop("PHI AND THETA: No threshold found - final incidence always 0")
    } else {
      stop("PHI: No threshold found - final incidence always 0")
    }
  } 
  else if (is.na(theta_index)) {
    stop("THETA: No threshold found - final incidence always 0")
  } 
  else if (phi_index == 1 | theta_index == 1) {
    stop("PHI/THETA: No threshold found - final incidence never 0")
  }
  # extract rows of data around the threshold
  phi_subset <- phi_data[(phi_index-1):phi_index,]
  theta_subset <- theta_data[(theta_index-1):theta_index,]
  
  # if the first instance where final disease incidence > 0 is smaller than 0.01, consider
  # the estimate for the threshold accurate enough and return it
  if (phi_subset[2, "final_I"] < 0.502 | theta_subset[2, "final_I"] < 0.502) {
    
    phi_thresh <- mean(c(phi_subset[1, "parm_val"], phi_subset[2, "parm_val"]))
    theta_thresh <- mean(c(theta_subset[1, "parm_val"], theta_subset[2, "parm_val"]))
    
    return(c(phi_threshold = phi_thresh,
             theta_threshold = theta_thresh))
  }
  
  # else run find_epidemic_threshold() with new smaller parameter limits recursively until 
  # above condition is reached
  
  phi_lower_lim <- phi_subset[1, "parm_val"]
  phi_upper_lim <- phi_subset[2, "parm_val"]
  
  theta_lower_lim <- theta_subset[1, "parm_val"]
  theta_upper_lim <- theta_subset[2, "parm_val"]
  
  num_parm_runs <- 100
  theta_vals_condensed <- list(theta = seq(theta_lower_lim, theta_upper_lim, length.out = num_parm_runs))
  
  phi_vals_condensed <- list(phi = seq(phi_lower_lim, phi_upper_lim, length.out = num_parm_runs))
  
  # re-run sensitivity analysis
  new_theta_phi_vals <- sensitivity_analysis(parms_mad = phi_vals_condensed,
                                             parms_don = theta_vals_condensed,
                                             init_states_don, 
                                             init_states_mad, 
                                             madden_func,
                                             madden_cunniffe_func,
                                             donnelly_func,
                                             times, 
                                             parms)
  
  new_theta_phi_vals$final_I <- new_theta_phi_vals$final_I / parms[["H"]]
  
  # feed back into function recursively
  find_theta_phi_vals(new_theta_phi_vals, 
                      init_states_don, 
                      init_states_mad, 
                      times, 
                      parms)
}

