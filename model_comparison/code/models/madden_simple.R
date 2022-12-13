
#########################################################
##### SIMPLIFICATION OF MODEL BY MADDEN ET AL. 2000 #####
#########################################################

### Simplified version of model by Madden et al. 2000 to make it more comparable to model by Donnelly et al. 2019.
### Constant vector population size: no vector immigration, emigration, birth or death. Plants are represented by 
### an SI rather than an SEIR model, vectors represented by an SI rather than SEI model. Vectors cannot be born 
### infective (as this model represents NPT viruses).

rm(list=ls())

#packages
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)

# TIME FRAME
times <- seq(0, 20, by = 0.2)

# INITIAL STATES
P <- 1200 # total number of vectors
N <- 400 # number of plants
I <- 1 # number of infected plants

init_states <- c(
  I = I, # number of infected plants
  #S = N - I, # number of susceptible plants
  X = P, # number healthy insects
  Z = 0 # number infective insects
)

## DONNELLY (2019) MODEL PARAMETERS - for determination of phi
theta <- 1
w <- 0.2


# PARAMETERS
alpha <- 0.2 # mortality rate
k1 <- 1/0.021
T <- 0.5 # usually 0.5/phi
lamda = 1/0.021
tau <- 1/0.25
b <- 1 #- exp(-k1 * T)
a <- 1 #- exp(-lamda * T)

phi <- sqrt(theta*tau*(1-w)/(w#*a*b
                             )) # plants visited per day by an insect

parms <- c(
  N = N, # number of host plants (equivalent to K in original Madden model)
  P = P, # total number of vectors
  phi = 4, # plants visited per day by insect
  k1 = k1, # inoculation rate by insect (set to NPT virus, equiv to 0.5hr)
  lamda = lamda, # rate of acquisition of virus by vector (set to NPT virus, equiv to 0.5hr)
  alpha = alpha, # vector population mortality rate
  T = T, # # time feeding/probing per plant visit
  c = 0.5,#20/3/2, # natural plant death rate (equivalent to beta in original Madden model)
  d = 0.5,#20/3/2, # plant death due to infection
  v_t = alpha * (init_states[["X"]] + init_states[["Z"]]), # birth rate of vectors per day
  tau = 4, # rate of moving through infectious state in vector (set to NPT virus, equiv to 6hr)
  b = b, # prob of plant inoculation per infective insect visit 
  a = a # probability of vector acquisition of virus from an infected plant per visit
  
)


madden_simple_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  
  become_infected <- par[["phi"]] * par[["b"]] * y[["Z"]] * (par[["N"]] - y[["I"]]) / par[["N"]]
  
  natural_death_I <- par[["c"]] * y[["I"]]
  #natural_death_S <- par[["c"]] * y[["S"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  #birth <- par[["c"]] * par[["N"]] + par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  #dS <- birth - become_infected - natural_death_S
  
  # VECTOR EQUATIONS
  
  acquisition <- par[["phi"]] * par[["a"]] * y[["I"]] * y[["X"]] / par[["N"]]
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dX <- - acquisition + stop_being_infective
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dX, dZ)))
}

run <- data.frame(ode(y = init_states,
                         times = times,
                         func = madden_simple_ode,
                         parms = parms))

run_long <- reshape2::melt(run, id.vars = "time")
names(run_long) <- c("time", "compartment", "number")

# plot trajectory of infected plants and vector compartments
plant_trajec <- ggplot(data = run_long %>% filter(compartment == "I"), 
                       aes(x = time, y = number)) +
  geom_line() +
  labs(x = "Time (days)", y = "Number of infected plants, I")
plant_trajec

vector_trajec <- ggplot(data = run_long %>% filter(!(compartment == "I")),
                       aes(x = time, y = number, col = compartment)) +
  geom_line()
vector_trajec
