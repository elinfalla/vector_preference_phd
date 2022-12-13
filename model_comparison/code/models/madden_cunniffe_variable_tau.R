#############################################################################
##### CREATION OF MADDEN-CUNNIFFE MODEL WHERE TAU VARIES ALONG WITH PHI #####
#############################################################################

rm(list = ls())
library(deSolve)
library(dplyr)

# source model comparison functions - vary_param(), sensitivity_analysis(), find_theta_phi_vals()
source("./code/comparison/model_comparison_FUNCTIONS.R")

# source other versions of madden, donnelly and madden-cunniffe models
source("./code/comparison/model_ode_versions_FUNCTIONS.R")

### PARMS
# parameters
parms <- c(
  
  ## DONNELLY PARMS
  gamma = 4, # rate of recovery/death of I plants
  theta = 1.05, # aphid dispersal rate per day
  Pacq = 0.8, # chance of virus acquisition by vector from infected plant
  Pinoc = 0.8, # chance of inoculating healthy plant
  w_don = 0.2, # feeding probability on healthy plant

  ## MADDEN PARMS

  tau = 4, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 3.5, # plants visited per day by an insect
  
  ## PARMS FOR BOTH
  A = 1200, # total (constant) number of vectors
  H = 400, # number host plants 
  v = 1, # infected plant attractiveness
  e = 1 # infected plant acceptability 
  
)


### SET MADDEN PARMS THAT ARE DETERMINED BY DONNELLY PARMS
# probability of vector acquisition of virus from an infected plant per visit
parms[["a"]] <- parms[["Pacq"]] #1 - exp(-parms[["lamda"]] * parms[["T"]])

# probability of plant inoculation per infective insect visit 
parms[["b"]] <- parms[["Pinoc"]] #1 - exp(-parms[["k1"]] * parms[["T"]])


# natural plant death rate (equivalent to beta in original Madden model)
parms[["c"]] <- parms[["gamma"]] / 2

# plant death due to infection
parms[["d"]] <- parms[["gamma"]] / 2

parms[["eta"]] <- 1/(parms[["phi"]]*parms[["w_mad"]])

# STATES
I <- 1 # initial number of infected plants

init_states_don <- c(
  I = I # number of infected plants
)

init_states_mad <- c(
  I = I, # number of infected plants
  Z = 0 # number infective insects
)

times <- seq(0, 40, by = 0.4)


madden_cunniffe_vpref_varytau_ode <- function(times, y, par) {
  
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  # define variable phi
  phi <- (S + y[["I"]]*par[["v"]]) / (par[["w_mad"]]*par[["eta"]]*(S + y[["I"]]*par[["v"]]*par[["e"]]))
  
  # make tau equivalent to phi
  tau <- phi
  # PLANT EQUATIONS
  
  become_infected <- phi * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  
  # VECTOR EQUATIONS
  X <- par[["A"]] - y[["Z"]]
  
  acquisition <- phi * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * X * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- tau * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
}


#### FIND NEW PHI VAL FOR 0.5 EQUILIBRIUM INCIDENCE ####
num_parm_runs <- 100
parms_mad <- list(phi = seq(0, 20, length.out = num_parm_runs))
parms_don <- list(theta = seq(0, 20, length.out = num_parm_runs))

phi_theta_sense_anal <- sensitivity_analysis(parms_mad, 
                                             parms_don, 
                                             init_states_don, 
                                             init_states_mad, 
                                             madden_func = madden_vpref_ode, 
                                             madden_cunniffe_func = madden_cunniffe_vpref_varytau_ode,
                                             donnelly_func = donnelly_vpref_ode,
                                             times = times, 
                                             default_parms = parms)
phi_theta_sense_anal <- phi_theta_sense_anal %>% # make final incidence into proportion of plants rather than number
  mutate(final_I = final_I/parms[["H"]])

phi_theta_thresholds <- find_theta_phi_vals(theta_phi_data = phi_theta_sense_anal, 
                                            init_states_don, 
                                            init_states_mad, 
                                            madden_func = madden_vpref_ode, 
                                            madden_cunniffe_func = madden_cunniffe_vpref_varytau_ode,
                                            donnelly_func = donnelly_vpref_ode,
                                            times = times, 
                                            parms = parms)

parms[["eta"]] <- 1/(phi_theta_thresholds[["phi_threshold"]]*parms[["w_mad"]]) # set eta using phi threshold found
parms[["phi"]] <- phi_theta_thresholds[["phi_threshold"]]
parms[["theta"]] <- phi_theta_thresholds[["theta_threshold"]]

trajec_mad_cun <- data.frame(ode(y = init_states_don,
                                 times = times,
                                 parms = parms,
                                 func = donnelly_vpref_ode))
trajec_mad_cun$I <- trajec_mad_cun$I/(parms[["H"]]) # make into a proportion

trajec_mad_cun_long <- reshape2::melt(trajec_mad_cun, id.vars = "time")
names(trajec_mad_cun_long) <- c("time", "compartment", "number")

mad_cun_plant_trajec <- ggplot(data = trajec_mad_cun_long %>% filter(compartment == "I"),
                           aes(x = time, y = number)) +
  geom_line() +
  ggtitle("Madden model (fixed phi)") +
  labs(y = "Proportion of infected plants", x = "Time (days)")


mad_cun_vec_trajec <- ggplot(data = trajec_mad_cun_long %>% filter(!compartment == "I"),
                         aes(x = time, y = number, col = compartment)) +
  geom_line() +
  labs(col = "Vector state", y = "Number of vectors", x = "Time (days)")

grid.arrange(mad_cun_plant_trajec, mad_cun_vec_trajec)

