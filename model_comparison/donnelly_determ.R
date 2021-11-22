
#######################################################################
##### RECREATION OF STOCHASTIC SIMULATION BY DONNELLY ET AL. 2019 #####
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
i = 1/parms[["H"]],# frequency of infected plants
As = 1200 / parms[["H"]] / 2, # per plant aphid density on average healthy plant
Ai = 1200 / parms[["H"]] / 2 # per plant aphid density on average infected plant
)

# timeframe
times <- seq(0, 10, by = 0.2)

# define ode function
donnelly_ode_func <- function(times, states, parms) {
  
  #states
  i <- states[["i"]]
  As <- states[["As"]]
  Ai <- states[["Ai"]]
  
  # parameters
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
  
  # set up equations for state equations
  q <- 1 - p  # q = prob of surviving
  
  i_hat <- v*i / ((1-i) + v*i)
  
  xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e))
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  A <- As*(H - i*H) + Ai*(i*H)  # A = total number of aphids i.e. aphids per plant type * number of that plant type
  
  # STATE EQUATIONS
  # aphids
  dAs <- a*As*(1 - As/K) - b*As - theta*As*(1 - Fs) + theta*Ai*Fs*i/(1 - i)
  dAi <- a*Ai*(1 - Ai/K) - b*Ai - theta*Ai*(1 - Fi) + theta*As*Fi*(1 - i)/i

  # infected plants
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(c(di, dAs, dAi)))
}

vary_param_trajec <- function(init_states, times, parms, varied_parm_name, varied_parm_vals) {
  
  list_of_trajects <- vector(mode = "list", length = length(varied_parm_vals))
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states, 
                                 times = times, 
                                 parms = parms, 
                                 func = donnelly_ode_func))
    
    trajectory_long <- reshape2::melt(trajectory, id.vars="time")
    names(trajectory_long) <- c("time", "compartment", "number")
    
    list_of_trajects[[run]] <- trajectory_long
  }
  
  return(list_of_trajects)
}

# run epidemic with default parameters (as donnelly et al. 2019)
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

# run epidemic with varying values of e and w
e_vals <- c(0.25, 0.5, 1, 1.5, 2)

e_trajecs <- vary_param_trajec(init_states, times, parms, "e", e_vals)
e_trajecs_i_only <- list(e_trajecs[[1]] %>% filter(compartment == "i"),
                         e_trajecs[[2]] %>% filter(compartment == "i"),
                         e_trajecs[[3]] %>% filter(compartment == "i"),
                         e_trajecs[[4]] %>% filter(compartment == "i"),
                         e_trajecs[[5]] %>% filter(compartment == "i"))

vary_e_plot <- ggplot(data=e_trajecs_i_only[[1]], aes(x = time, y = number)) +
  geom_line(col="red") +
  geom_line(data=e_trajecs_i_only[[2]], col="orange") +
  geom_line(data=e_trajecs_i_only[[3]], col="yellow") +
  geom_line(data=e_trajecs_i_only[[4]], col="green") +
  geom_line(data=e_trajecs_i_only[[5]], col="blue")
vary_e_plot
