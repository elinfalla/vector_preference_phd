############################################################################################
##### COMPARISON OF VECTOR PREFERENCE IN MADDEN (2000) AND DONNELLY (2019) TYPE MODELS #####
############################################################################################

### MADDEN MODEL:
### Modified version of model by Madden et al. 2000 to include vector preference - both of plant attractiveness and
### plant acceptability. No vector population dynamics (immigration, emigration birth, death), plants are represented 
### by an SI rather than an SEIR model, vectors represented by an SI rather than SEI model. Vectors cannot be born 
### infective (as this model represents NPT viruses).

### DONNELLY MODEL:
### Modification of deterministic model by Donnelly et al. 2019 to make it more comparable to this version of the 
### Madden et al. 2000 model. Vector preference included but no vector population dynamics - constant vector population 
### size, no vector emigration (p always set to 0). Also model records number of infected plants, not proportion.

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
madden_vpref_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  become_infected <- par[["phi"]] * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  
  # VECTOR EQUATIONS
  X <- par[["A"]] - y[["Z"]]
  
  acquisition <- par[["phi"]] * par[["a"]] * (1 - par[["e"]]*par[["w"]]) * X * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
}

#### DEFINE DONNELLY ODE
donnelly_vpref_ode <- function(times, states, parms) {
  
  # STATES
  I <- states[["I"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  theta <- parms[["theta"]]
  H <- parms[["H"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w"]]
  Pacq <- parms[["Pacq"]]
  Pinoc <- parms[["Pinoc"]]
  A <- parms[["A"]]
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  
  # I_hat = weighted number of infected plants, accounting for attraction towards infected plants
  i_hat <- v*I / (H - I + v*I)
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  #xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  xI <- (Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(w*(1 - i_hat*(1 - e)))
  
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A)*xI - gamma*I
  
  return(list(di))
}

# define timeframe
times <- seq(0, 40, by = 0.4)

# parameters
parms <- c(
  
  ## DONNELLY PARMS
  gamma = 20/3, # rate of recovery/death of I plants
  theta = 1, # aphid dispersal rate per day
  Pacq = 1, # chance of virus acquisition by vector from infected plant
  Pinoc = 1, # chance of inoculating healthy plant
  
  ## MADDEN PARMS
  k1 = 1/0.021, # inoculation rate by insect (set to NPT virus, equiv to 0.5hr)
  lamda = 1/0.021, # rate of acquisition of virus by vector (set to NPT virus, equiv to 0.5hr)
  T = 0.2, # # time feeding/probing per plant visit (set to 0.5/phi in Madden model)
  tau = 1/0.25, # rate of moving through infectious state in vector
  
  ## PARMS FOR BOTH
  A = 1200, # total (constant) number of vectors
  H = 400, # number host plants 
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  w = 0.2 # feeding rate on healthy plant
)


### SET MADDEN PARMS THAT ARE DETERMINED BY DONNELLY PARMS
# probability of vector acquisition of virus from an infected plant per visit
a <- parms[["Pacq"]] #1 - exp(-parms[["lamda"]] * parms[["T"]])

# probability of plant inoculation per infective insect visit 
b <- parms[["Pinoc"]] #1 - exp(-parms[["k1"]] * parms[["T"]])

# plants visited per day by an insect
phi <- parms[["theta"]] #sqrt(parms[["theta"]] * parms[["tau"]] * (1 - parms[["w"]]) / 
#               (parms[["w"]] #* parms[["a"]] * parms[["b"]]
#               )) 

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
  Z = 0 # number infective insects
)


### RUN DONNELLY EPIDEMIC
trajectory_don <- data.frame(ode(y = init_states_don,
                             times = times,
                             parms = parms,
                             func = donnelly_vpref_ode))

# plot trajectory of infected plants and number of aphids
plant_trajec <- ggplot(data = trajectory_don, aes(x = time, y = I)) +
  geom_line() +
  labs(x = "Time (days)",
       y = "Number of infected plants, I")

### RUN MADDEN EPIDEMIC
trajectory_mad <- data.frame(ode(y = init_states_mad,
                      times = times,
                      func = madden_vpref_ode,
                      parms = parms))

trajectory_mad_long <- reshape2::melt(trajectory_mad, id.vars = "time")
names(trajectory_mad_long) <- c("time", "compartment", "number")

# plot trajectory of infected plants and vector compartments
plant_trajec <- plant_trajec +
  geom_line(data = trajectory_mad_long %>% filter(compartment == "I"),
            aes(x = time, y = number), color = "red") +
  annotate("text", label = "red = Madden\nblack = Donnelly", x = times[length(times)/3*2], y = 0)
plant_trajec


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
                                 func = donnelly_vpref_ode))
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
                                 func = madden_vpref_ode))
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

num_parm_runs <- 50
analysis_parms_don <- list(theta = seq(0, 20, length.out = num_parm_runs),
                           w = seq(0.01, 1, length.out = num_parm_runs),
                           Pacq = seq(0, 1, length.out = num_parm_runs),
                           Pinoc = seq(0, 1, length.out = num_parm_runs),
                           gamma = seq(0, 20, length.out = num_parm_runs),
                           v = seq(0, 20, length.out = num_parm_runs),
                           e = seq(0, 20, length.out = num_parm_runs))

analysis_parms_mad <- list(phi = seq(0, 20, length.out = num_parm_runs),
                           tau = seq(0, 20, length.out = num_parm_runs),
                           a = seq(0, 1, length.out = num_parm_runs),
                           b = seq(0, 1, length.out = num_parm_runs),
                           c = seq(0, 20, length.out = num_parm_runs),
                           d = seq(0, 20, length.out = num_parm_runs),
                           v = seq(0, 20, length.out = num_parm_runs),
                           e = seq(0, 20, length.out = num_parm_runs))


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
                                    don_model_plots[[6]],
                                    don_model_plots[[7]],
                                    ncol = 1)

mad_plots <- gridExtra::arrangeGrob(mad_model_plots[[1]], 
                                    mad_model_plots[[2]], 
                                    mad_model_plots[[3]],
                                    mad_model_plots[[4]],
                                    mad_model_plots[[5]], 
                                    mad_model_plots[[6]],
                                    mad_model_plots[[7]],
                                    mad_model_plots[[8]],
                                    ncol = 1)

layout <- rbind(c(1,2,3),
                c(1,2,4),
                c(1,2,4),
                c(1,2,4))
title <- textbox_grob("Model comparison with vector preference - constant vector population", 
                      gp = gpar(fontface = "bold",
                                fontsize = 13),
                      padding = unit(c(0, 1, 0, 1), "cm"))

# create pdf file to print plot to
pdf(file = "sens_analysis_vec_pref_models.pdf")
all_plots <- gridExtra::grid.arrange(don_plots,
                                     mad_plots,
                                     title,
                                     parms_grob,
                                     layout_matrix = layout,
                                     widths = c(2,2,1))
all_plots
dev.off()

##### CREATE HEATMAPS OF V VS E for final incidence of disease

find_I_equilibrium <- function(parms, model) {
  
  ## function that for given parameter values, finds the non-zero equilibrium value of I,
  ## number of infected plants, for either the Donnelly model or the Madden model
  
  if (model == "Donnelly") {
    gamma <- parms[["gamma"]]
    theta <- parms[["theta"]]
    H <- parms[["H"]]
    v <- parms[["v"]]
    e <- parms[["e"]]
    w <- parms[["w"]]
    Pacq <- parms[["Pacq"]]
    Pinoc <- parms[["Pinoc"]]
    A <- parms[["A"]]
    
    # coefficients worked out by hand (see lab book)
    coeff_0 <- theta*A*v*(1 - e*w)*H*Pacq*Pinoc - gamma*w*H^2
    coeff_1 <- gamma*w*v*(1 - e)*H + 2*gamma*w*H - 2*gamma*w*v*H - theta*A*v*(1 - e*w)*Pacq*Pinoc
    coeff_2 <- gamma*w*(v^2*(1 - e) - v*(1 - e) - v^2 + 2*v - 1)
  }
  else if (model == "Madden") {
    c <- parms[["c"]]
    d <- parms[["d"]]
    phi <- parms[["phi"]]
    H <- parms[["H"]]
    v <- parms[["v"]]
    e <- parms[["e"]]
    w <- parms[["w"]]
    a <- parms[["a"]]
    b <- parms[["b"]]
    A <- parms[["A"]]
    tau <- parms[["tau"]]
    
    coeff_0 <- (c + d)*tau*H^2 - phi^2*b*H*a*(1 - e*w)*A*v
    coeff_1 <- (c + d)*(2*tau*H*v - 2*tau*H + phi*a*(1 - e*w)*v*H) + phi^2*b*a*(1 - e*w)*A*v
    coeff_2 <- (c + d)*(tau + tau*v^2 - 2*tau*v - phi*a*(1 - e*w)*v + phi*a*(1 - e*w)*v^2)
  }
  else {
    stop("model must be either 'Donnelly' or 'Madden'")
  }
  
  # solve for I using coefficients to get all equilibrimus
  equilibriums <- polyroot(c(coeff_0, coeff_1, coeff_2))
  
  # subset to those that have no imaginary part, then remove imaginary part (as it is 0)
  equilibriums <- equilibriums[Im(zapsmall(equilibriums)) == 0]
  equilibriums <- Re(equilibriums)
  
  # subset to biologically possible equilibrium i.e. I > 0, I <= H (number host plants)
  equilibriums <- equilibriums[equilibriums > 0 & equilibriums <= parms[["H"]]]
  
  if (length(equilibriums) == 0) { # there is no stable positive equilibrium for I
    return(0)
  } 
  else if (length(equilibriums) > 1) {
      #stop("Multiple positive equilibriums for I")
    return(NA)
  }
  
  # turn equilibrium into proportion of infected plants, i
  equilibriums <- equilibriums / parms[["H"]]
  
  return(equilibriums)

}

equilibrium_apply_func <- function(e_v_vals, parms, col_names, model) {
  
  ## function to be passed to an apply function - for given e and v vals, it changes
  ## e and v to those values and calculates the equilibrium value of I, number of 
  ## plants infected, for the Donnelly model
  
  parms[[ col_names[1] ]] <- e_v_vals[1]
  parms[[ col_names[2] ]] <- e_v_vals[2]
  
  equilibrium <- find_I_equilibrium(parms, model)
  
  return(equilibrium)
}

# define e and v vals across which the equilibrim of I will be found
v_vals <- seq(0.001, 5, length.out = 100)
e_vals <- seq(0.001, 1/parms[["w"]], length.out = 100) # in madden e*w must be < 1


# create df with all combinations of v and e
v_e_df <- expand.grid(v = v_vals, e = e_vals) 

# create new I_eq column for each model by calculating the I equilibrium for each row of v_e_df (i.e. each combination
# of v and e values) for first the Donnelly model, then the Madden model
I_eq_values <- lapply(c("Donnelly", "Madden"), 
                      function(model) apply(v_e_df, 1, equilibrium_apply_func, 
                                            parms = parms, 
                                            col_names = names(v_e_df), 
                                            model = model))
v_e_df$I_eq_don <- I_eq_values[[1]]
v_e_df$I_eq_mad <- I_eq_values[[2]]


# turn any bistable equilibriums into a high number so they're a different colour on the heatmap
# v_e_df$I_eq <- apply(v_e_df, 1, function(row) if(length(row$I_eq) > 1) {row$I_eq = 1.2}
#                                               else {row$I_eq = row$I_eq})

heatmap_don <- ggplot(data = v_e_df, aes(x = v, y = e, fill = I_eq_don)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "yellow", name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Donnelly")


heatmap_mad <- ggplot(data = v_e_df, aes(x = v, y = e, fill = I_eq_mad)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "yellow", name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Madden")

grid.arrange(heatmap_don, heatmap_mad)

c <- parms[["c"]]
d <- parms[["d"]]
phi <- parms[["phi"]]
H <- parms[["H"]]
v <- parms[["v"]]
e <- parms[["e"]]
w <- parms[["w"]]
a <- parms[["a"]]
b <- parms[["b"]]
A <- parms[["A"]]
tau <- parms[["tau"]]

coeff_0 <- (c + d)*tau*H^2 - phi^2*b*H*a*(1 - e*w)*A*v
coeff_1 <- (c + d)*(2*tau*H*v - 2*tau*H + phi*a*(1 - e*w)*v*H) + phi^2*b*a*(1 - e*w)*A*v
coeff_2 <- (c + d)*(tau + tau*v^2 - 2*tau*v - phi*a*(1 - e*w)*v + phi*a*(1 - e*w)*v^2)

polyroot(c(coeff_0, coeff_1, coeff_2))

