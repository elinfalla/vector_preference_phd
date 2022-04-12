
########################################################
#### CHECK I AND i VERSIONS OF DONNELLY MODEL ALIGN ####
########################################################

## This is a version of the Donnelly model with vector preference but no aphid 
## popultion dynamics - constant aphid population size

rm(list = ls())

# packages
library(deSolve)
library(ggplot2)
library(gridExtra)


donnelly_vpref_ode_I <- function(times, states, parms) {
  
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
  dI <- (theta*A)*xI - gamma*I
  
  return(list(dI))
}

donnelly_vpref_ode_i <- function(times, states, parms) {
  
  # STATES
  i <- states[["i"]]
  
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
  i_hat <- v*i / ((1-i) + v*i) 
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  #xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  xi <- (Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(w*(1 - i_hat*(1 - e)))
  
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(di))
}

parms <- c(
  gamma = 3/20, # rate of recovery/death of I plants
  b = 0.1, # aphid mortality rate per day
  theta = 1, # aphid dispersal rate per day
  a = 2, # reproduction rate per aphid per day
  K = 10, # aphid reproduction limit - maximum aphids per plant (diff to donnelly_stoch_sim.R)
  H = 100, # number host plants
  p = 0.2, # aphid emigration/death rate per journey between plants
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  w = 0.2, # feeding rate on healthy plant
  Pacq = 0.5, # chance of virus acquisition from infected plant
  Pinoc = 0.5, # chance of inoculating healthy plant
  A = 1200 # number aphids
)

#states
num_I <- 1

init_states_i <- c(
  i = num_I/parms[["H"]] # frequency of infected plants
)

init_states_I <- c(
  I = num_I # frequency of infected plants
)

# timeframe
times <- seq(0, 8, by = 0.2)

########
# run epidemic
########

run_epidemic <- function(init_states_i, init_states_I, times, parms) {
  
  ### runs ode for both versions of the Donnelly model and plots graph of 
  ### number of infected over time for each, on the same plot
  
  trajectory_I <- data.frame(ode(y = init_states_I, 
                                 times = times, 
                                 parms = parms, 
                                 func = donnelly_vpref_ode_I))
  
  trajectory_long_I <- reshape2::melt(trajectory_I, id.vars="time")
  names(trajectory_long_I) <- c("time", "compartment", "number")
  
  trajectory_i <- data.frame(ode(y = init_states_i, 
                                 times = times, 
                                 parms = parms, 
                                 func = donnelly_vpref_ode_i))
  
  trajectory_long_i <- reshape2::melt(trajectory_i, id.vars="time")
  names(trajectory_long_i) <- c("time", "compartment", "number")
  trajectory_long_i$number <- trajectory_long_i$number * parms[["H"]]
  
  unname(trajectory_long_i$number) == unname(trajectory_long_I$number)
  
  plot <- ggplot(data = trajectory_long_i, aes(x = time, y = number)) +
    geom_line(color = "red") +
    geom_line(data = trajectory_long_I, linetype = "dashed") + 
    annotate("text", label = paste("dashed = I\nH =", parms[["H"]]), x = 6, y = 50)
  
  return(plot)
}

## run epidemic for different values of H (number of plants) to verify the models give the same
# result

H_vals <- c(200, 400, 600, 800, 1000, 1200)
out <- list()

for (i in 1:length(H_vals)) {
  parms[["H"]] <- H_vals[i]
  init_states_i <- c(i = num_I/parms[["H"]])
  out[[i]] <- run_epidemic(init_states_i, init_states_I, times, parms)
}

# plot result
grid.arrange(out[[1]], out[[2]], out[[3]], out[[4]], out[[5]], out[[6]])

