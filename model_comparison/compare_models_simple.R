#############################################################################
##### COMPARISON OF SIMPLIFIED MADDEN (2000) AND DONNELLY (2019) MODELS #####
#############################################################################

### MADDEN MODEL:
### Simplified version of model by Madden et al. 2000 to make it more comparable to model by Donnelly et al. 2019.
### No vector dynamics (immigration, emigration birth, death), plants are represented by an SI rather than an SEIR model, vectors 
### represented by an SI rather than SEI model. Vectors cannot be born infective (as this model represents 
### NPT viruses).

### DONNELLY MODEL:
### Simplification of deterministic model by Donnelly et al. 2019 to make it more comparable to Madden et al. 2000.
### No vector dynamics - constant vector population size, no vector preference (v and e always set to 1), no vector
### emigration (p always set to 0). Also model records number of infected plants, not proportion.

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gridtext)

#### DEFINE MADDEN ODE
madden_simple_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  become_infected <- par[["phi"]] * par[["b"]] * y[["Z"]] * (par[["H"]] - y[["I"]]) / par[["H"]]
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death

  # VECTOR EQUATIONS
  acquisition <- par[["phi"]] * par[["a"]] * y[["I"]] * y[["X"]] / par[["H"]]
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dX <- - acquisition + stop_being_infective
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dX, dZ)))
}

#### DEFINE DONNELLY ODE
donnelly_simple_ode <- function(times, states, parms) {
  
  # STATES
  I <- states[["I"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  theta <- parms[["theta"]]
  H <- parms[["H"]]
  p <- parms[["p"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w"]]
  Pacq <- parms[["Pacq"]]
  Pinoc <- parms[["Pinoc"]]
  A <- parms[["A"]]
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  # q = prob of surviving flight between plants (i.e. doesn't emigrate/die)
  #q <- 1 - p  
  
  # i_hat = weighted frequency of infected plants, accounting for attraction towards infected plants
  # i.e. adapted version of i / (1 - i + i) = i / 1, to include v
  #i_hat <- v*i / ((1-i) + v*i) 
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  #xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  xI <- Pacq*Pinoc*(1 - w)*I*(H - I)/(w*H)
  
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xI - gamma*I
  
  return(list(di))
}

# define timeframe
times <- seq(0, 40, by = 0.4)

# FIXED parameters (from Donnelly model)
p <- 0 # aphid emigration/death rate per journey between plants
v <- 1 # infected plant attractiveness
e <- 1 # infected plant acceptability 

# parameters
parms <- c(
  
  ## DONNELLY PARMS
  gamma = 20/3, # rate of recovery/death of I plants
  theta = 1, # aphid dispersal rate per day
  p = p, # aphid emigration/death rate per journey between plants
  v = v, # infected plant attractiveness
  e = e, # infected plant acceptability 
  w = 0.2, # feeding rate on healthy plant
  Pacq = 1, # chance of virus acquisition by vector from infected plant
  Pinoc = 1, # chance of inoculating healthy plant
  
  ## MADDEN PARMS
  k1 = 1/0.021, # inoculation rate by insect (set to NPT virus, equiv to 0.5hr)
  lamda = 1/0.021, # rate of acquisition of virus by vector (set to NPT virus, equiv to 0.5hr)
  T = 0.5, # # time feeding/probing per plant visit (set to 0.5/phi in Madden model)
  tau = 1/0.25, # rate of moving through infectious state in vector (set to NPT virus, equiv to 6hr)
  
  ## PARMS FOR BOTH
  A = 1200, # total (constant) number of vectors
  H = 400 # number host plants # vector population size
)


### SET MADDEN PARMS THAT ARE DETERMINED BY DONNELLY PARMS
# probability of vector acquisition of virus from an infected plant per visit
a <- parms[["Pacq"]] #1 - exp(-parms[["lamda"]] * parms[["T"]])

# probability of plant inoculation per infective insect visit 
b <- parms[["Pinoc"]] #1 - exp(-parms[["k1"]] * parms[["T"]])

# plants visited per day by an insect
phi <- sqrt(parms[["theta"]] * parms[["tau"]] * (1 - parms[["w"]]) / 
              (parms[["w"]] #* parms[["a"]] * parms[["b"]]
               )) 

# natural plant death rate (equivalent to beta in original Madden model)
c <- parms[["gamma"]] / 2

# plant death due to infection
d <- parms[["gamma"]] / 2

# add to parms
parms <- c(parms, c(a = a, 
                    b = b, 
                    phi = phi,
                    c = c, 
                    d = d))

# STATES
I <- 1 # initial number of infected plants

init_states_don <- c(
  I = I # number of infected plants
)

init_states_mad <- c(
  I = I, # number of infected plants
  X = parms[["A"]], # number healthy insects
  Z = 0 # number infective insects
)

# ### RUN DONNELLY EPIDEMIC
# trajectory_don <- data.frame(ode(y = init_states_don, 
#                              times = times, 
#                              parms = parms, 
#                              func = donnelly_simple_ode))
# 
# # plot trajectory of infected plants and number of aphids
# plant_trajec <- ggplot(data = trajectory_don, aes(x = time, y = I)) +
#   geom_line() +
#   labs(x = "Time (days)", 
#        y = "Number of infected plants, I")
# 
# ### RUN MADDEN EPIDEMIC
# trajectory_mad <- data.frame(ode(y = init_states_mad,
#                       times = times,
#                       func = madden_simple_ode,
#                       parms = parms))
# 
# trajectory_mad_long <- reshape2::melt(trajectory_mad, id.vars = "time")
# names(trajectory_mad_long) <- c("time", "compartment", "number")
# 
# # plot trajectory of infected plants and vector compartments
# plant_trajec <- plant_trajec +
#   geom_line(data = trajectory_mad_long %>% filter(compartment == "I"), 
#             aes(x = time, y = number), color = "red") 
# plant_trajec


# param_vs_final_I <- function(init_states_don, init_states_mad, times, parms, varied_parms_vals) {
#   
#   ### describe function here
#   
#   # initialise output data frame
#   params <- varied_parms_vals %>% slice(rep(1:n(), times = 2))
#   
#   output <- data.frame(params,
#                        final_I = rep(NA, length.out = nrow(params)),
#                        model = rep(c("Donnelly", "Madden"), each = nrow(params)/2))
#   
#   final_I_don <- c()
#   final_I_mad <- c()
#   
#   parm_names_grep <- paste("^", names(varied_parms_vals), sep = "", collapse = "|")
#   
#   for (run in 1:length(varied_parms_vals)) {
#     
#     parms[grepl(parm_names_grep, names(parms))] <- varied_parms_vals[run,]
#     
#     trajectory_don <- data.frame(ode(y = init_states_don, 
#                                  times = times, 
#                                  parms = parms, 
#                                  func = donnelly_simple_ode))
#     final_I_don <- c(final_I_don, round(trajectory_don[nrow(trajectory_don), "I"], 3))
#     
#     trajectory_mad <- data.frame(ode(y = init_states_mad, 
#                                  times = times, 
#                                  parms = parms, 
#                                  func = madden_simple_ode))
#     final_I_mad <- c(final_I_mad, round(trajectory_mad[nrow(trajectory_mad), "I"], 3))
#     
#   }
#   output[,"final_I"] <- c(final_I_don, final_I_mad)
#   
#   return(output)
# }



vary_param_don <- function(init_states_don, times, parms, varied_parm_name, varied_parm_vals) {
  
  # function to run the Donnelly ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of final incidence, I, of all runs, along with values of varied parameter
  
  # initialise output data frame
  output <- data.frame(parm_name = rep(varied_parm_name, length(varied_parm_vals)),
                             parm_val = varied_parm_vals,
                             final_I = rep(NA, length(varied_parm_vals)),
                             model = rep("Donnelly", length(varied_parm_vals))
                             )

  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states_don, 
                                 times = times, 
                                 parms = parms, 
                                 func = donnelly_simple_ode))
    output[run, "final_I"] <- round(trajectory[nrow(trajectory), "I"], 3)
    
  }
  
  return(output)
}

vary_param_mad <- function(init_states_mad, times, parms, varied_parm_name, varied_parm_vals) {
  
  # function to run the Madden ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of final incidence, I, of all runs, along with values of varied parameter
  
  # initialise output data frame
  output <- data.frame(parm_name = rep(varied_parm_name, length(varied_parm_vals)),
                       parm_val = varied_parm_vals,
                       final_I = rep(NA, length(varied_parm_vals)),
                       model = rep("Madden", length(varied_parm_vals))
  )
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states_mad, 
                                 times = times, 
                                 parms = parms, 
                                 func = madden_simple_ode))
    output[run, "final_I"] <- round(trajectory[nrow(trajectory), "I"], 3)
    
  }
  
  return(output)
}

sensitivity_analysis <- function(parms_mad, parms_don, init_states_don, init_states_mad, times, default_parms) {
  
  sen_analysis_df <- data.frame(parm_name = character(),
                                 parm_val = numeric(),
                                 final_I = numeric(),
                                 model = character())
  
  print("Starting Madden sensitivity analysis")
  
  for (i in 1:length(parms_mad)) {
    
    output <- vary_param_mad(init_states_mad, times, parms,
                             varied_parm_name = names(parms_mad[i]),
                             varied_parm_vals = parms_mad[[i]]) 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output)

  }
  
  print("Starting Donnelly sensitivity analysis")
  for (i in 1:length(parms_don)) {
    
    output <- vary_param_don(init_states_don, times, parms,
                             varied_parm_name = names(parms_don[i]),
                             varied_parm_vals = parms_don[[i]]) 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output)
    
  }
  
  return(sen_analysis_df) 
}
  

analysis_parms_don <- list(theta = seq(0, 20, length.out = 30),
                           w = seq(0.01, 1, length.out = 30),
                           Pacq = seq(0, 1, length.out = 30),
                           Pinoc = seq(0, 1, length.out = 30),
                           gamma = seq(0, 20, length.out = 30))

analysis_parms_mad <- list(phi = seq(0, 20, length.out = 30),
                           tau = seq(0, 20, length.out = 30),
                           a = seq(0, 1, length.out = 30),
                           b = seq(0, 1, length.out = 30),
                           c = seq(0, 20, length.out = 30),
                           d = seq(0, 20, length.out = 30))


sens_analysis_res <- sensitivity_analysis(parms_mad = analysis_parms_mad,
                                          parms_don = analysis_parms_don,
                                          init_states_don, 
                                          init_states_mad, 
                                          times, 
                                          parms)

###  CREATE DONNELLY AND MADDEN MODEL SENSITIVITY ANALYSIS PLOTS
don_res_only <- sens_analysis_res %>%
  filter(model == "Donnelly")
mad_res_only <- sens_analysis_res %>%
  filter(model == "Madden")

don_model_plots <- lapply(unique(don_res_only$parm_name), function(p) ggplot(data = don_res_only[don_res_only$parm_name == p,],
                                                                                 aes(x = parm_val, 
                                                                                     y = final_I)) +
                            geom_line() +
                            labs(x = p))

mad_model_plots <- lapply(unique(mad_res_only$parm_name), function(p) ggplot(data = mad_res_only[mad_res_only$parm_name == p,],
                                                                             aes(x = parm_val, 
                                                                                 y = final_I)) +
                            geom_line() +
                            labs(x = p))
  
# create table to be given alongside graphs giving default parameter values
parms_table <- round(data.frame(parms[names(parms) != "k1" &
                                     names(parms) != "lamda" &
                                     names(parms) != "T"]), 2)
parms_grob <- tableGrob(parms_table, cols = c("Value")) # turn into grob (gtable) for plotting

# ARRANGE PLOTS (AND PARAMETER TABLE) AND SAVE TO PDF
don_plots <- gridExtra::arrangeGrob(don_model_plots[[1]], 
                        don_model_plots[[2]], 
                        don_model_plots[[3]],
                        don_model_plots[[4]],
                        don_model_plots[[5]], 
                        grid.rect(gp=gpar(col="white")), # empty space
                        nrow = 3, ncol = 2)

mad_plots <- gridExtra::arrangeGrob(mad_model_plots[[1]], 
                        mad_model_plots[[2]], 
                        mad_model_plots[[3]],
                        mad_model_plots[[4]],
                        mad_model_plots[[5]], 
                        mad_model_plots[[6]],
                        nrow = 3, ncol = 2)

layout <- rbind(c(1,1,3),
                c(1,1,4),
                c(2,2,4),
                c(2,2,4))
title <- textbox_grob("Simple models comparison - no vector dynamics or preference", 
                      gp = gpar(fontface = "bold",
                                fontsize = 13),
                      padding = unit(c(0, 1, 0, 1), "cm"))

pdf(file = "sens_analysis_simple_models.pdf")
all_plots <- gridExtra::grid.arrange(don_plots, 
                                     mad_plots,
                                     title,
                                     parms_grob,
                                     layout_matrix = layout)
all_plots
dev.off()
