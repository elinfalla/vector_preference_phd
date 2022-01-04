
#######################################################################
###### RECREATION OF DETERMINISTIC MODEL BY DONNELLY ET AL. 2019 ######
#######################################################################

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

# parameters
parms <- c(
  gamma = 20/3, # rate of recovery/death of I plants
  b = 0.1, # aphid mortality rate per day
  theta = 1, # aphid dispersal rate per day
  a = 2, # reproduction rate per aphid per day
  K = 10, # aphid reproduction limit - maximum aphids per plant (diff to donnelly_stoch_sim.R)
  H = 400, # number host plants
  p = 0.2, # aphid emigration/death rate per journey between plants
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  w = 0.2, # feeding rate on healthy plant
  Pacq = 1, # chance of virus acqusition from infected plant
  Pinoc = 1 # chance of inoculating healthy plant
)

#states
init_states <- c(
i = 1/parms[["H"]], # frequency of infected plants
As = 1200 / parms[["H"]] / 2, # per plant aphid density on average healthy plant
Ai = 1200 / parms[["H"]] / 2 # per plant aphid density on average infected plant
)

# timeframe
times <- seq(0, 10, by = 0.2)

# define ode function
donnelly_ode_func <- function(times, states, parms) {
  
  # STATES
  i <- states[["i"]]
  As <- states[["As"]]
  Ai <- states[["Ai"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  b <- parms[["b"]]
  theta <- parms[["theta"]]
  a <- parms[["a"]]
  K <- parms[["K"]]
  H <- parms[["H"]]
  p <- parms[["p"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w"]]
  Pacq <- parms[["Pacq"]]
  Pinoc <- parms[["Pinoc"]]
  
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  
  # q = prob of surviving flight between plants (i.e. doesn't emigrate/die)
  q <- 1 - p  
  
  # i_hat = weighted frequency of infected plants, accounting for attraction towards infected plants
  # i.e. adapted version of i / (1 - i + i) = i / 1, to include v
  i_hat <- v*i / ((1-i) + v*i) 
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  # A = total number of aphids i.e. aphids per plant type * number of that plant type
  # i.e. S*As + I*Ai where S is number of susceptible plants and I is number of infected plants
  A <- As*(H - i*H) + Ai*(i*H)  
  
  
  # STATE EQUATIONS
  # aphids - change in per plant aphid density on susceptible (As)  and infected (Ai) plants
  # = aphid births - aphid deaths - aphid settling on other plant type + aphids moving from other plant type
  dAs <- a*As*(1 - As/K) - b*As - theta*As*(1 - Fs) + theta*Ai*Fs*i/(1 - i)
  dAi <- a*Ai*(1 - Ai/K) - b*Ai - theta*Ai*(1 - Fi) + theta*As*Fi*(1 - i)/i

  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(c(di, dAs, dAi)))
}


vary_param_trajec <- function(init_states, times, parms, varied_parm_name, varied_parm_vals) {
  
  # function to run the ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of all runs
  
  # initialise output data frame
  all_trajects <- data.frame(time = numeric(),
                             compartment = character(),
                             value = numeric(),
                             param_val = numeric())
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states, 
                                 times = times, 
                                 parms = parms, 
                                 func = donnelly_ode_func))
    
    trajectory_long <- reshape2::melt(trajectory, id.vars="time")
    
    # add the value of the varying parameter as an extra column
    trajectory_long$param_val <- rep(varied_parm_vals[run], times=nrow(trajectory_long))
    
    # add to dataframe of all trajectories
    all_trajects <- rbind(all_trajects, trajectory_long)
  }
  
  return(all_trajects)
}

########
# run epidemic with default parameters (as donnelly et al. 2019)
########
trajectory <- data.frame(ode(y = init_states, 
                             times = times, 
                             parms = parms, 
                             func = donnelly_ode_func))

trajectory_long <- reshape2::melt(trajectory, id.vars="time")
names(trajectory_long) <- c("time", "compartment", "number")

# plot trajectory of infected plants and number of aphids
plant_trajec <- ggplot(data=trajectory_long %>% filter(compartment == "i"), 
                       aes(x = time, y = number)) +
  geom_line() +
  labs(x = "Time (days)", y = "Frequency of infected plants, i")
plant_trajec

aphid_trajec <- ggplot(data=trajectory_long %>% filter(!(compartment == "i")),
                       aes(x = time, y = number, col = compartment)) +
  geom_line()
aphid_trajec

##########
# run epidemic with varying values of e and v
##########

#VARY E (plant acceptability)
e_vals <- c(0.25, 0.5, 1, 1.5, 2)

e_trajecs <- vary_param_trajec(init_states, times, parms, "e", e_vals)

# plot output
e_vals_plot <- ggplot(data=e_trajecs %>% filter(variable == "i"), aes(x = time, y = value, col = as.factor(param_val))) +
  geom_line() +
  labs(col = "e value")
e_vals_plot

# VARY V (plant attractiveness)
v_vals <- c(1, 3, 0.6, 0.55)

long_times <- seq(0, 20, by=0.5)
v_trajecs <- vary_param_trajec(init_states, long_times, parms, "v", v_vals)

# plot output
v_vals_plot <- ggplot(data=v_trajecs %>% filter(variable == "i"), aes(x = time, y = value, col = as.factor(param_val))) +
  geom_line() +
  labs(col = "v value")
v_vals_plot
