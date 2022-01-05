
#########################################################
##### SIMPLIFICATION OF MODEL BY MADDEN ET AL. 2000 #####
#########################################################

### Simplified version of model by Madden et al. 2000 to make it more comparable to model by Donnelly et al. 2019.
### No vector immigration, plants are represented by an SI rather than an SEIR model.

rm(list=ls())

#packages
library(deSolve)
library(reshape2)

# TIME FRAME
times <- seq(0, 20, by = 0.2)

# INITIAL STATES
P <- 1200 # total number of vectors
N <- 400 # number of plants
I <- 100

init_states <- c(
  I = I, # number of infected plants
  #S = N - I, # number of susceptible plants
  X = P, # number healthy insects
  Y = 0, # number latent insects
  Z = 0 # number infective insects
)

# PARAMETERS
phi <- 1 # plants visited per day by an insect
alpha <- 0.2 # mortality rate

parms <- c(
  N = N, # number of host plants (equivalent to K in original Madden model)
  P = P, # total number of vectors
  phi = phi, # plants visited per day by insect
  k1 = 1/0.021, # inoculation rate by insect (set to NPT virus, equiv to 0.5hr)
  T = 0.5/phi, # # time feeding/probing per plant visit
  c = 20/3/2, # natural plant death rate (equivalent to beta in original Madden model)
  d = 20/3/2, # plant death due to infection
  lamda = 1/0.021, # rate of acquisition of virus by vector (set to NPT virus, equiv to 0.5hr)
  v_t = alpha * (init_states[["X"]] + init_states[["Y"]] + init_states[["Z"]]), # birth rate of vectors per day
  q = 0, # probability of vector offspring being viliferous (set to NPT virus)
  alpha = alpha, # vector population mortality rate
  Ex = 0.2/3, # emigration rate for healthy insects
  Ey = 0.2/3, # emigration rate for latent state insects
  Ez = 0.2/3, # emigration rate for infective insects
  tau = 1/0.25, # rate of moving through infectious state in vector (set to NPT virus, equiv to 6hr)
  eta = 99999999999 # rate of moving through latent state in vector (set to NPT virus, no latent stage)
)


madden_simple <- function(times, y, par) {
  
  # PLANT EQUATIONS
  b <- 1 - exp(-par[["k1"]] * par[["T"]]) # prob of plant inoculation per infective insect visit 
  
  become_infected <- par[["phi"]] * b * y[["Z"]] * (par[["N"]] - y[["I"]]) / par[["N"]]
  
  natural_death_I <- par[["c"]] * y[["I"]]
  #natural_death_S <- par[["c"]] * y[["S"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  #birth <- par[["c"]] * par[["N"]] + par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  #dS <- birth - become_infected - natural_death_S
  
  # VECTOR EQUATIONS
  a <- 1 - exp(-par[["lamda"]] * par[["T"]]) # probability of vector acquisition of virus from an infected plant per visit
  
  # birth and death
  birth <- par[["v_t"]]
  born_infective <- par[["v_t"]] * par[["q"]] * (y[["Z"]] / (y[["X"]] + y[["Y"]] + y[["Z"]]))
  death_X <- par[["alpha"]] * y[["X"]]
  death_Y <- par[["alpha"]] * y[["Y"]]
  death_Z <- par[["alpha"]] * y[["Z"]]
  
  # emigration
  emigration_X <- par[["Ex"]] * y[["X"]]
  emigration_Y <- par[["Ey"]] * y[["Y"]]
  emigration_Z <- par[["Ez"]] * y[["Z"]]
  
  # virus
  acquisition <- par[["phi"]] * a * y[["I"]] * y[["X"]] / par[["N"]]
  become_infective <- par[["eta"]] * y[["Y"]]
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dX <- birth - born_infective - death_X - emigration_X - acquisition + stop_being_infective
  dY <- acquisition - become_infective - emigration_Y - death_Y
  dZ <- become_infective - stop_being_infective - emigration_Z + born_infective - death_Z
  
  return(list(c(dI, #dS, 
                dX, dY, dZ)))
}

run <- data.frame(ode(y = init_states,
                         times = times,
                         func = madden_simple,
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
