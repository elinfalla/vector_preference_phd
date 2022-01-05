
###########################################################################
###### SIMPLIFICATION OF DETERMINISTIC MODEL BY DONNELLY ET AL. 2019 ######
###########################################################################

### Simplification of deterministic model by Donnelly et al. 2019 to make it more comparable to Madden et al. 2000.
### No vector dynamics - constant vector population size, no vector preference (v and e always set to 1), no vector
### emigration (p always set to 0).

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

# FIXED parameters
p = 0 # aphid emigration/death rate per journey between plants
v = 1 # infected plant attractiveness
e = 1 # infected plant acceptability 

# parameters
parms <- c(
  gamma = 20/3, # rate of recovery/death of I plants
  theta = 1, # aphid dispersal rate per day
  H = 400, # number host plants
  p = p, # aphid emigration/death rate per journey between plants
  v = v, # infected plant attractiveness
  e = e, # infected plant acceptability 
  w = 0.2, # feeding rate on healthy plant
  Pacq = 1, # chance of virus acquisition from infected plant
  Pinoc = 1, # chance of inoculating healthy plant
  A = 1200 # vector population size
)

#states
init_states <- c(
  i = 1/parms[["H"]] # frequency of infected plants
)

# timeframe
times <- seq(0, 20, by = 0.2)

# define ode function
donnelly_simple_ode <- function(times, states, parms) {
  
  # STATES
  i <- states[["i"]]
  
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
  q <- 1 - p  
  
  # i_hat = weighted frequency of infected plants, accounting for attraction towards infected plants
  # i.e. adapted version of i / (1 - i + i) = i / 1, to include v
  i_hat <- v*i / ((1-i) + v*i) 
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))

  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(di))
}




########
# run epidemic with default parameters (as donnelly et al. 2019)
########
trajectory <- data.frame(ode(y = init_states, 
                             times = times, 
                             parms = parms, 
                             func = donnelly_simple_ode))

# plot trajectory of infected plants and number of aphids
plant_trajec <- ggplot(data=trajectory, aes(x = time, y = i)) +
  geom_line() +
  labs(x = "Time (days)", y = "Frequency of infected plants, i")
plant_trajec

