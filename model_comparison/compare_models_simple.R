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
    
    output <- vary_param_mad(init_states_mad, times, default_parms,
                             varied_parm_name = names(parms_mad[i]),
                             varied_parm_vals = parms_mad[[i]]) 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output)

  }
  
  print("Starting Donnelly sensitivity analysis")
  for (i in 1:length(parms_don)) {
    
    output <- vary_param_don(init_states_don, times, default_parms,
                             varied_parm_name = names(parms_don[i]),
                             varied_parm_vals = parms_don[[i]]) 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output)
    
  }
  
  return(sen_analysis_df) 
}
  
num_parm_runs <- 100

analysis_parms_don <- list(theta = seq(0, 20, length.out = num_parm_runs),
                           w = seq(0.01, 1, length.out = num_parm_runs),
                           Pacq = seq(0, 1, length.out = num_parm_runs),
                           Pinoc = seq(0, 1, length.out = num_parm_runs),
                           gamma = seq(0, 20, length.out = num_parm_runs))

analysis_parms_mad <- list(phi = seq(0, 20, length.out = num_parm_runs),
                           tau = seq(0, 20, length.out = num_parm_runs),
                           a = seq(0, 1, length.out = num_parm_runs),
                           b = seq(0, 1, length.out = num_parm_runs),
                           c = seq(0, 20, length.out = num_parm_runs),
                           d = seq(0, 20, length.out = num_parm_runs))


# sens_analysis_res <- sensitivity_analysis(parms_mad = analysis_parms_mad,
#                                           parms_don = analysis_parms_don,
#                                           init_states_don, 
#                                           init_states_mad, 
#                                           times, 
#                                           parms)
# 
# ###  CREATE DONNELLY AND MADDEN MODEL SENSITIVITY ANALYSIS PLOTS
# don_res_only <- sens_analysis_res %>%
#   filter(model == "Donnelly")
# mad_res_only <- sens_analysis_res %>%
#   filter(model == "Madden")
# 
# don_model_plots <- lapply(unique(don_res_only$parm_name), function(p) ggplot(data = don_res_only[don_res_only$parm_name == p,],
#                                                                                  aes(x = parm_val, 
#                                                                                      y = final_I)) +
#                             geom_line() +
#                             labs(x = p))
# 
# mad_model_plots <- lapply(unique(mad_res_only$parm_name), function(p) ggplot(data = mad_res_only[mad_res_only$parm_name == p,],
#                                                                              aes(x = parm_val, 
#                                                                                  y = final_I)) +
#                             geom_line() +
#                             labs(x = p))
#   
# # create table to be given alongside graphs giving default parameter values
# parms_table <- round(data.frame(parms[names(parms) != "k1" &
#                                      names(parms) != "lamda" &
#                                      names(parms) != "T"]), 2)
# parms_grob <- tableGrob(parms_table, cols = c("Value")) # turn into grob (gtable) for plotting
# 
# # ARRANGE PLOTS (AND PARAMETER TABLE) AND SAVE TO PDF
# don_plots <- gridExtra::arrangeGrob(don_model_plots[[1]], 
#                         don_model_plots[[2]], 
#                         don_model_plots[[3]],
#                         don_model_plots[[4]],
#                         don_model_plots[[5]], 
#                         grid.rect(gp=gpar(col="white")), # empty space
#                         nrow = 3, ncol = 2)
# 
# mad_plots <- gridExtra::arrangeGrob(mad_model_plots[[1]], 
#                         mad_model_plots[[2]], 
#                         mad_model_plots[[3]],
#                         mad_model_plots[[4]],
#                         mad_model_plots[[5]], 
#                         mad_model_plots[[6]],
#                         nrow = 3, ncol = 2)
# 
# layout <- rbind(c(1,1,3),
#                 c(1,1,4),
#                 c(2,2,4),
#                 c(2,2,4))
# title <- textbox_grob("Simple models comparison - no vector dynamics or preference", 
#                       gp = gpar(fontface = "bold",
#                                 fontsize = 13),
#                       padding = unit(c(0, 1, 0, 1), "cm"))
# 
# # create pdf file to print plot to
# pdf(file = "sens_analysis_simple_models.pdf")
# all_plots <- gridExtra::grid.arrange(don_plots, 
#                                      mad_plots,
#                                      title,
#                                      parms_grob,
#                                      layout_matrix = layout)
# all_plots
# dev.off()


###################################################
### INVESTIGATE RELATIONSHIP BETWEEN THETA AND PHI
###################################################
## theta = Donnelly model, aphid dispersal rate per day
## phi = Madden model, number plants visited per day by vector

## rerun sensitivity analysis for just phi and theta
upper_lim <- 50
theta_vals <- list(theta = seq(0, upper_lim, length.out = num_parm_runs))

phi_vals <- list(phi = seq(0, upper_lim, length.out = num_parm_runs))

theta_phi_data <- sensitivity_analysis(parms_mad = phi_vals,
                                          parms_don = theta_vals,
                                          init_states_don, 
                                          init_states_mad, 
                                          times, 
                                          parms)

plot_reverse_final_I <- ggplot(data = theta_phi_data, aes(x = final_I, y = parm_val)) +
  geom_line() +
  facet_wrap(~parm_name)
plot_reverse_final_I

plot_final_I <- ggplot(data = theta_phi_data, aes(x = parm_val, y = final_I)) +
  geom_line() +
  facet_wrap(~parm_name) +
  annotate("text", label = paste("w =", parms[["w"]]), x = 35, y = 1)

plot_final_I

find_epidemic_threshold <- function(theta_phi_data, init_states_don, init_states_mad, times, parms) {
  
  # function to find point at which epidemic takes off (final_I > 0) for phi and theta
  
  phi_data <- theta_phi_data %>% filter(parm_name == "phi")
  theta_data <- theta_phi_data %>% filter(parm_name == "theta")
  
  phi_index <- match(TRUE, phi_data$final_I > 0) # finds first instance of final_I > 0
  theta_index <- match(TRUE, theta_data$final_I > 0)

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
  if (phi_subset[2, "final_I"] < 0.002 | theta_subset[2, "final_I"] < 0.002) {
    
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
  epidemic_threshold_res <- sensitivity_analysis(parms_mad = phi_vals_condensed,
                                                 parms_don = theta_vals_condensed,
                                                 init_states_don, 
                                                 init_states_mad, 
                                                 times, 
                                                 parms)
  # feed back into function recursively
  find_epidemic_threshold(epidemic_threshold_res, 
                          init_states_don, 
                          init_states_mad, 
                          times, 
                          parms)
}

out <- find_epidemic_threshold(theta_phi_data,
                               init_states_don, 
                               init_states_mad, 
                               times, 
                               parms)
out["phi_threshold"] # 2.88
out["theta_threshold"] # 0.54

############
### relationship between phi and theta depends on w (feeding rate)
# w = 0.5 gives expected value of 1 probe per dispersal, making phi and theta equivalent
# see lab book for explanation

# set w to 1
parms_new <- parms
parms_new["w"] <- 0.5

# re-run sensitivity analysis for phi and theta when w = 1
new_w_res <- sensitivity_analysis(parms_mad = phi_vals,
                                      parms_don = theta_vals,
                                      init_states_don, 
                                      init_states_mad, 
                                      times, 
                                      parms_new)

w_plot_final_I <- ggplot(data = new_w_res, aes(x = parm_val, y = final_I)) +
  geom_line() +
  facet_wrap(~parm_name) +
  annotate("text", label = paste("w =", parms_new[["w"]]), x = 35, y = 1)
w_plot_final_I

# same plot zoomed in on epidemic threshold
zoom_phi_vals <- list(phi = seq(2.5, 3, length.out = num_parm_runs))
zoom_theta_vals <- list(theta = seq(2, 2.5, length.out = num_parm_runs))

w_res_zoom <- sensitivity_analysis(parms_mad = zoom_phi_vals,
                                   parms_don = zoom_theta_vals,
                                   init_states_don, 
                                   init_states_mad, 
                                   times, 
                                   parms_new)
zoom_plot_final_I <- ggplot(data = w_res_zoom, aes(x = parm_val, y = final_I)) +
  geom_line() +
  facet_wrap(~parm_name, scales = "free")
zoom_plot_final_I

# compare graphs when w = 0.2 and 0.5 - theta looks much more like phi when w=0.5
grid.arrange(plot_final_I, w_plot_final_I)

# compare parameter values at which epidemic takes off for phi and theta when w=0.5
threshold <- find_epidemic_threshold(new_w_res,
                                     init_states_don,
                                     init_states_mad,
                                     times,
                                     parms_new)
threshold["phi_threshold"] # 2.882359 
threshold["theta_threshold"] # 2.157943

## plot w (feeding rate) against 1-w/w (expected number of probes per dispersal)
# see lab book for derivation
w <- seq(0,1,by=0.02)
plot(w, (1-w)/w, 
     type = "l",
     xlab = "w, feeding rate",
     ylab = "(1-w)/w, expected number of plant probes per dispersal")


#####
### try varying phi and tau simultaneously
# making phi = tau may help equate the models (see lab book)

vary_2_param_mad <- function(init_states_mad, times, parms, varied_parm_name1, varied_parm_name2, varied_parm_vals) {
  
  # function to run the Madden ODE multiple times with *two* parameters varying, as specified by user.
  # both varying parameters are the same at any given point.
  # returns data frame of final incidence, I, of all runs, along with values of varied parameter
  
  # initialise output data frame
  parm_names <- paste(varied_parm_name1, varied_parm_name2, sep = ",")
  
  output <- data.frame(parm_names = rep(parm_names, length(varied_parm_vals)),
                       parm_val = varied_parm_vals,
                       final_I = rep(NA, length(varied_parm_vals)),
                       model = rep("Madden", length(varied_parm_vals))
  )
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name1]] <- varied_parm_vals[run]
    parms[[varied_parm_name2]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states_mad, 
                                 times = times, 
                                 parms = parms, 
                                 func = madden_simple_ode))
    output[run, "final_I"] <- round(trajectory[nrow(trajectory), "I"], 3)
    
  }
  
  return(output)
}

phi_tau_df <- vary_2_param_mad(init_states_mad, 
                              times, 
                              parms_new, 
                              "phi", 
                              "tau", 
                              unlist(unname(phi_vals)))

tau_phi_plot_final_I <- ggplot(data = phi_tau_df, aes(x = parm_val, y = final_I)) +
  geom_line() +
  annotate("text", label = unique(phi_tau_df$parm_names), x = 0, y = 350) +
  ylim(0, 400)
tau_phi_plot_final_I

w_plot_final_I_theta_only <- ggplot(data = new_w_res %>% filter(parm_name == "theta"), 
                                    aes(x = parm_val, y = final_I)) +
  geom_line() +
  annotate("text", label = paste("w =", parms_new[["w"]]), x = 35, y = 1) +
  annotate("text", label = "theta", x = 0, y = 350) +
  ylim(0, 400)

grid.arrange(tau_phi_plot_final_I, w_plot_final_I_theta_only)

### zoom in to see epidemic threshold
phi_tau_zoom_df <- vary_2_param_mad(init_states_mad, 
                               times, 
                               parms_new, 
                               "phi", 
                               "tau", 
                               unlist(unname(zoom_phi_vals)))

zoom_tau_phi_plot_final_I <- ggplot(data = phi_tau_zoom_df, aes(x = parm_val, y = final_I)) +
  geom_line() +
  annotate("text", label = unique(phi_tau_df$parm_names), x = 0, y = 350)

zoom_plot_final_I_theta_only <- ggplot(data = w_res_zoom %>% filter(parm_name == "theta"),
                                       aes(x = parm_val, y = final_I)) +
  geom_line() +
  annotate("text", label = paste("w =", parms_new[["w"]]), x = 4.5, y = 1) +
  annotate("text", label = "theta", x = 0, y = 200)

grid.arrange(zoom_tau_phi_plot_final_I, zoom_plot_final_I_theta_only)

## look at similarity between phi and theta vs final_I when tau = 2, not 4
parms_new["tau"] <- 2

# re-run sensitivity analysis for phi and theta when tau = 2
new_tau_res <- sensitivity_analysis(parms_mad = phi_vals,
                                    parms_don = theta_vals,
                                    init_states_don, 
                                    init_states_mad, 
                                    times, 
                                    parms_new)
tau_plot_final_I <- ggplot(data = new_tau_res, aes(x = parm_val, y = final_I)) +
  geom_line() +
  facet_wrap(~parm_name) +
  annotate("text", label = paste("tau =", parms_new[["tau"]]), x = 35, y = 1)
tau_plot_final_I

grid.arrange(tau_plot_final_I, w_plot_final_I)

# compare disease thresholds
threshold <- find_epidemic_threshold(new_tau_res,
                                     init_states_don,
                                     init_states_mad,
                                     times,
                                     parms_new)
threshold["phi_threshold"] 
threshold["theta_threshold"] 

## test idea that epidemic trajectories are similar when phi = theta, w=0.5, tau=2
vary_trajec <- function(varied_parm_vals, H_vals, 
                        init_states_don, init_states_mad, times, parms) {
  
  ### function that runs the Madden and Donnelly ODEs for combinations of parameter values specified by
  ### user for different values of H and plots I over time. Returns list of all plots.
  

  # create df of all combinations of H with the other parameters
  all_combos <- expand.grid(parm_1 = varied_parm_vals[,1], H = H_vals)
  all_combos <- cbind(all_combos, varied_parm_vals[,-1])
  names(all_combos)[1] <- names(varied_parm_vals)[1]
  
  # initialise plot list
  plots <- list()
  
  for (i in 1:nrow(all_combos)) {
    
    for (col in 1:ncol(all_combos)) {
      parm_name <- names(all_combos)[col]
      
      parms[[parm_name]] <- all_combos[i, parm_name]
    }
    
    # run donnelly epidemic
    trajectory_don <- data.frame(ode(y = init_states_don,
                                     times = times,
                                     parms = parms,
                                     func = donnelly_simple_ode))
    
    ### run madden epidemic
    trajectory_mad <- data.frame(ode(y = init_states_mad,
                                     times = times,
                                     func = madden_simple_ode,
                                     parms = parms))
    
    # turn madden output into long format
    trajectory_mad_long <- reshape2::melt(trajectory_mad, id.vars = "time")
    names(trajectory_mad_long) <- c("time", "compartment", "number")
    
    # plot donnelly trajectory 
    plant_trajec <- ggplot(data = trajectory_don, aes(x = time, y = I)) +
      geom_line()  +
      ylim(0, parms[["H"]])
    
    # plot madden trajectory on same plot
    plant_trajec <- plant_trajec +
      geom_line(data = trajectory_mad_long %>% filter(compartment == "I"),
                aes(x = time, y = number), color = "red") +
      annotate("text", label = paste("red = Madden\n", all_combos[i,1], all_combos[i,2]), 
               x = times[length(times)/3*2], y = parms[["H"]]*2/3, size = 3)
    
    # store plot in list
    plots[[i]] <- plant_trajec
    
  }
  
  return(plots)
}

### TEST 1: tau = phi = theta, w = 0.5
parms_test1 <- parms
parms_test1[["w"]] <- 0.5

varied_parm_vals <- data.frame(theta = 1:5, 
                         phi = 1:5, 
                         tau = 1:5)
H_vals <- c(200, 400, 600, 800)

tau_equal_phi_plots <- vary_trajec(varied_parm_vals = varied_parm_vals, 
                     H_vals = H_vals, 
                     parms = parms_test1,
                     init_states_don = init_states_don, 
                     init_states_mad = init_states_mad, 
                     times = times)

# create pdf of results
pdf("tau_equal_phi_test.pdf", height = 12, width = 12#, paper = "a4"
)
do.call("grid.arrange", c(tau_equal_phi_plots, ncol = nrow(varied_parm_vals)))

dev.off()

### TEST 2: tau = 0, phi=theta, w = 0.5
parms_test2 <- parms
parms_test2[["w"]] <- 0.5
parms_test2[["tau"]] <- 0.01

varied_parm_vals <- data.frame(theta = 1:5, 
                               phi = 1:5)

tau_equal_zero_plots <- vary_trajec(varied_parm_vals = varied_parm_vals, 
                                   H_vals = H_vals, 
                                   parms = parms_test2,
                                   init_states_don = init_states_don, 
                                   init_states_mad = init_states_mad, 
                                   times = times)

# create pdf of results
pdf("tau_equal_zero_test.pdf", height = 12, width = 12)
do.call("grid.arrange", c(tau_equal_zero_plots, ncol = nrow(varied_parm_vals)))
dev.off()

### TEST 3: tau = 2, phi=theta, w = 0.5
parms_test3 <- parms
parms_test3[["w"]] <- 0.5
parms_test3[["tau"]] <- 2

varied_parm_vals <- data.frame(theta = 1:5, 
                               phi = 1:5)

tau_equal_one_plots <- vary_trajec(varied_parm_vals = varied_parm_vals, 
                                    H_vals = H_vals, 
                                    parms = parms_test3,
                                    init_states_don = init_states_don, 
                                    init_states_mad = init_states_mad, 
                                    times = times)

# create pdf of results
pdf("tau_equal_one_test.pdf", height = 12, width = 12)
do.call("grid.arrange", c(tau_equal_one_plots, ncol = nrow(varied_parm_vals)))
dev.off()

### TEST 4: tau = H/20, phi=theta, w = 0.5
parms_test4 <- parms
parms_test4[["w"]] <- 0.5

varied_parm_vals <- data.frame(theta = 1:5, 
                               phi = 1:5,
                               tau = H_vals/20)

tau_prop_to_H_plots <- vary_trajec(varied_parm_vals = varied_parm_vals, 
                                   H_vals = H_vals, 
                                   parms = parms_test4,
                                   init_states_don = init_states_don, 
                                   init_states_mad = init_states_mad, 
                                   times = times)

# create pdf of results
pdf("tau_prop_to_H_test.pdf", height = 12, width = 12)
do.call("grid.arrange", c(tau_prop_to_H_plots, ncol = nrow(varied_parm_vals)))
dev.off()

# TEST 5: same as test 1 (tau=phi=theta) but A is bigger
parms_test5 <- parms
parms_test5[["w"]] <- 0.5
parms_test5[["A"]] <- 1600

varied_parm_vals <- data.frame(theta = 1:5, 
                               phi = 1:5, 
                               tau = 1:5)
H_vals <- c(200, 400, 600, 800)

tau_prop_to_phi_plots <- vary_trajec(varied_parm_vals = varied_parm_vals, 
                                   H_vals = H_vals, 
                                   parms = parms_test5,
                                   init_states_don = init_states_don, 
                                   init_states_mad = init_states_mad, 
                                   times = times)

# create pdf of results
pdf("tau_prop_to_phi_test.pdf", height = 12, width = 12)
do.call("grid.arrange", c(tau_prop_to_phi_plots, ncol = nrow(varied_parm_vals)))
dev.off()


#########
### RUN SENSITIVITY ANALYSIS OF PHI AND THETA FOR VARYING TAU AND W VALUES

multiple_sensitivity_analysis <- function(parms_mad, parms_don, 
                                          varied_parm_name_mad, varied_parm_val_mad, 
                                          varied_parm_name_don, varied_parm_val_don, 
                                          init_states_don, init_states_mad, times, parms) {
  
  ## function that runs a sensitivity analysis multiple times for varying values of one madden and one donnelly
  ## parameter. returns dataframe of output and final graph displaying result in the form of a grob. use grid.arrange()
  ## to display the plot.
  
  if (length(varied_parm_val_don) != length(varied_parm_val_mad)) {
    stop("varied_parm_val_mad and varied_parm_val_don must be the same length")
  }
  
  # initialise output dataframe
  all_sens_analysis_df <- data.frame(
                          analysis_parm_name = character(),
                          analysis_parm_val = numeric(),
                          final_I = numeric(),
                          model = character(),
                          varied_parm_name = character(),
                          varied_parm_val = numeric()
                          )
  
  # calculate how many rows of the output df are for each model
  num_mad_rows <- sum(lengths(parms_mad))
  num_don_rows <- sum(lengths(parms_don))
  
  # for loop to run sensitivity analysis for all values of varied_parms
  for (run in 1:length(varied_parm_val_don)) {
  
    # specify varying parameter values for this run
    parms[varied_parm_name_mad] <- varied_parm_val_mad[run]
    parms[varied_parm_name_don] <- varied_parm_val_don[run]
    
    # run sensitivity analysis with these parameter values
    sens_analysis <- sensitivity_analysis(parms_mad,
                           parms_don,
                           init_states_don, 
                           init_states_mad, 
                           times, 
                           parms)
    
    # create df for information about the varying parameters in this run
    varying_parm_data <- data.frame(c(rep(varied_parm_name_mad, num_mad_rows), 
                                      rep(varied_parm_name_don, num_don_rows)), # varied_parm_name column
                                    c(rep(varied_parm_val_mad[run], num_mad_rows),
                                      rep(varied_parm_val_don[run], num_don_rows)) # varied_parm_val column
                                    )
    
    # add varying parameters data to output of sensitivity analysis
    full_sens_analysis_df <- cbind(sens_analysis, varying_parm_data)
    names(full_sens_analysis_df) <- c("analysis_parm_name", "analysis_parm_val", "final_I", "model", "varied_parm_name", "varied_parm_val")
    
    # add rows of this sensitivity analysis to overall dataframe
    all_sens_analysis_df <- rbind(all_sens_analysis_df, full_sens_analysis_df)
  }
  
  # split data by model
  split_data <- split(all_sens_analysis_df, f = all_sens_analysis_df$model)
  
  # plot madden model data
  madden_plot <- ggplot(data = split_data$Madden,
                         aes(x = analysis_parm_val, y = final_I, col = as.factor(round(varied_parm_val,2)))) +
    facet_wrap(~model) +
    geom_line() +
    labs(color = unique(split_data$Madden$varied_parm_name),
         x = unique(split_data$Madden$analysis_parm_name), 
         y = "Final disease incidence")
  
  # plot donnelly
  donnelly_plot <- madden_plot %+% split_data$Donnelly + # %+% means same plot but using Donnelly data
    labs(color = unique(split_data$Donnelly$varied_parm_name),
         x = unique(split_data$Donnelly$analysis_parm_name))
  
  final_plot <- arrangeGrob(madden_plot, donnelly_plot)
  
  # return final dataframe and grob of final plot (use grid.arrange to plot it)
  return(list(all_sens_analysis_df, final_plot))
}

w_vals <- seq(0.1, 0.99999, length.out = 10)
tau_vals_1 <- 1/w_vals

mult_sens_analysis_phi_theta <- multiple_sensitivity_analysis(parms_mad = phi_vals, 
                              parms_don = theta_vals, 
                              varied_parm_name_mad = "tau", 
                              varied_parm_val_mad = tau_vals_1, 
                              varied_parm_name_don = "w",
                              varied_parm_val_don = w_vals, 
                              init_states_don, init_states_mad, times, parms)

# plot result
grid.arrange(mult_sens_analysis_phi_theta[[2]])

tau_vals_2 <- seq(0.1, 400, by = 40)

mult_sens_analysis_phi_theta_2 <- multiple_sensitivity_analysis(parms_mad = phi_vals, 
                                                        parms_don = theta_vals, 
                                                        varied_parm_name_mad = "tau", 
                                                        varied_parm_val_mad = tau_vals_2, 
                                                        varied_parm_name_don = "w",
                                                        varied_parm_val_don = w_vals, 
                                                        init_states_don, init_states_mad, times, parms)
grid.arrange(mult_sens_analysis_phi_theta_2[[2]])

#########
### RUN SENSITIVITY ANALYSIS OF W AND TAU FOR VARYING PHI AND THETA VALUES

w_vals_list <- list(w = seq(0.001, 0.999, length.out = 100))
tau_vals_list <- list(tau = seq(0, 50, length.out = 100))

phi_vals_vec <- seq(0, 40, length.out = 15)
theta_vals_vec <- seq(0, 40, length.out = 15)

mult_sens_analysis_w_tau <- multiple_sensitivity_analysis(parms_mad = tau_vals_list,
                                                              parms_don = w_vals_list,
                                                              varied_parm_name_mad = "phi",
                                                              varied_parm_val_mad = phi_vals_vec,
                                                              varied_parm_name_don = "theta",
                                                              varied_parm_val_don = theta_vals_vec,
                                                              init_states_don, init_states_mad, times, parms)
grid.arrange(mult_sens_analysis_w_tau[[2]])


#################
## INVESTIGATING RELATIONSHIP BETWEEN TAU AND W
#################

## rerun sensitivity analysis for just phi and theta
w_vals <- list(w = seq(0.0001, 1, length.out = num_parm_runs))

tau_vals <- list(tau = seq(0, 20, length.out = num_parm_runs))

w_tau_data <- sensitivity_analysis(parms_mad = tau_vals,
                                       parms_don = w_vals,
                                       init_states_don, 
                                       init_states_mad, 
                                       times, 
                                       parms)

plot_final_I <- ggplot(data = w_tau_data, aes(x = parm_val, y = final_I)) +
  geom_line() +
  facet_wrap(~parm_name, scales = "free")

plot_final_I
