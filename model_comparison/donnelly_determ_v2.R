
#######################################################################
###### RECREATION OF DETERMINISTIC MODEL BY DONNELLY ET AL. 2019 ######
#######################################################################

### This version (v2) works out the equilibrium values of As and Ai to feed into
### the equation for di/dt, removing the need for As/dt and Ai/dt equations.

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

### FUNCTIONS ###

As_Ai_to_optimise <- function(aphid_pop_density, parms) {
  
  ### Function that contains the equations for As and Ai to be optimised. Will be passed to calc_A_equilibrium() 
  ### to find the equilibrium As and Ai values
  
  As <- aphid_pop_density[1]
  Ai <- aphid_pop_density[2]
  
  # define necessary parameters
  i <- parms[["i"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  a <- parms[["a"]]
  b <- parms[["b"]]
  K <- parms[["K"]]
  theta <- parms[["theta"]]
  
  # equations
  i_hat <- v*i / ((1-i) + v*i) 
  
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  # equations to be minimised
  dAs <- a*As*(1 - As/K) - b*As - theta*As*(1 - Fs) + theta*Ai*Fs*i/(1 - i)
  dAi <- a*Ai*(1 - Ai/K) - b*Ai - theta*Ai*(1 - Fi) + theta*As*Fi*(1 - i)/i
  
  return(dAs^2 + dAi^2) # return sum of squares. this way the minimum value is 0 which indicates a root (equilibrium)
}

calc_A_equilibrium_optim <- function(parms, i, n_tests) {
  
  ### Function to calculate As and Ai equilibrium by optimising As_Ai_to_optimise(). Returns vector of As
  ### then Ai value
  
  # add current value of i as a parameter
  parms["i"] <- i
  
  # run optim() many times from different starting values to find all possible equilibriums
  optim_equilibrium <- lapply(1:n_tests, function(x) round(optim(par = runif(2, min=0, max=100), 
                                                                 fn = As_Ai_to_optimise, 
                                                                 parms = parms,
                                                                 method = "L-BFGS-B")$par, 3))
  
  # filter down to unique equilibriums
  unique_equilibrium <- optim_equilibrium[!duplicated(optim_equilibrium)]
  unique_equilibrium # only equilibrium where both numbers (As and Ai) are positive = 9.5, 9.5
  
  # filter down to the positive (i.e. biologically possible) solution(s)
  find_positive_solution <- function(equil) {
    if (equil[1] > 0 && equil[2] > 0) {
      return (T)
    }
    else {
      return(F)
    }
  }
  
  positive_equilibrium <- unlist(unique_equilibrium[sapply(unique_equilibrium, find_positive_solution)==T])
  
  # check if multiple positive solutions
  if (length(positive_equilibrium) > 2) {
    # if the sum of the difference between equilibriums for As and Ai are more than 0.1 (i.e. not just a rounding)
    # issue, stop with an error
    if(sum(abs(diff(positive_equilibrium, lag = 2))) > 0.1) { 
      stop("More than 1 positive solution for As and Ai values")
    }
  }
  
  return(positive_equilibrium) # return equilibrium As and Ai values
}

calc_A_equilibrium <- function(coeff_1, coeff_2, coeff_3, coeff_4, parms, states, i_hat) {
  
  # define states
  i <- states[["i"]]
  
  # define parameters
  e <- parms[["e"]]
  a <- parms[["a"]]
  b <- parms[["b"]]
  K <- parms[["K"]]
  theta <- parms[["theta"]]
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  ## define coefficients 
  # (worked out by hand from dAs and dAi equations- see lab book). will give roots of Ai
  coeff_1 <- i*(a - b - theta*(1 - Fs))*(b + theta*(1 - Fi) - a)/(theta*Fi*(1 - i)) + theta*Fs*i/(1 - i)
  
  coeff_2 <- (a/K)*i*(a - b - theta*(1 - Fs))/(theta*Fi*(1 - i)) - 
    a*i^2*(b^2 + 2*b*theta*(1 - Fi) - 2*b*a + theta^2*(1 - Fi)^2 - 2*a*theta*(1 - Fi) + a^2)/(K*theta^2*(1 - i)^2*Fi^2)
  
  coeff_3 <- -(a*i^2/K)*(2*b*a + 2*a*theta*(1 - Fi) - 2*a^2)/(K*theta^2*(1 - i)^2*Fi^2)
  
  coeff_4 <- -a^3*i^2/(K^3*theta^2*(1 - i)^2*Fi^2)
  
  
  # get Ai equilibrium values
  Ai_equilibriums <- polyroot(c(0, coeff_1, coeff_2, coeff_3, coeff_4))
  
  # subset to those that have no imaginary part, then remove imaginary part (as it is 0)
  Ai_equilibriums <- Ai_equilibriums[Im(zapsmall(Ai_equilibriums)) == 0]
  Ai_equilibriums <- Re(Ai_equilibriums)
  
  # for each Ai equilibrium, calculate value of As (equation calculated manually - see lab book)
  As_equilibriums <- i*(b*Ai_equilibriums + theta*Ai_equilibriums*(1 - Fi) - a*Ai_equilibriums*(1-Ai_equilibriums/K)) / 
    (theta*Fi*(1-i))
  
  # find equilibrium where both As and Ai are positive i.e. biologically possible equilibrium
  Ai <- Ai_equilibriums[which(Ai_equilibriums > 0 & As_equilibriums > 0)]
  As <- As_equilibriums[which(Ai_equilibriums > 0 & As_equilibriums > 0)]
  
  # if multiple positive equilibriums, stop with error
  if (length(Ai) > 1 | length(As) > 1) {
    stop("multiple positive As, Ai equilibriums")
  }
  return(c(As, Ai))
}

donnelly_ode_func <- function(times, states, parms) {
  
  ### ODE function for deterministic Donnelly model, with As and Ai (aphids per healthy and infected plants 
  ### respectively) assumed to be always at equilibrium. Therefore calculates equilibrium As and Ai values 
  ### for each value of i (proportion infected plants).
  
  # STATES
  i <- states[["i"]]
  
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
  

  
  # calculate equilibrium level of As and Ai for current value of i
  # = per plant aphid density on susceptible (As)  and infected (Ai) plants
  

  # 2. plug into calc_A_equilibrium function
  As_Ai <- calc_A_equilibrium(coeff_1, coeff_2, coeff_3, coeff_4, parms, states, i_hat)
  As <- As_Ai[1]
  Ai <- As_Ai[2]
  
  ### OLD WAY USING OPTIM
  # As_Ai <- calc_A_equilibrium_optim(parms, i, n_tests = 150)
  # As <- As_Ai[1] 
  # Ai <- As_Ai[2]
  
  # A = total number of aphids i.e. aphids per plant type * number of that plant type
  # i.e. S*As + I*Ai where S is number of susceptible plants and I is number of infected plants
  A <- As*(H - i*H) + Ai*(i*H)  
  
  # STATE EQUATION
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(c(di)))
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

#### DEFINE ALL PARAMETERS AND INITIAL STATES

## donnelly param vals:
# gamma = 20/3, # rate of recovery/death of I plants
# b = , # aphid mortality rate per day
# theta = 1, # aphid dispersal rate per day
# a = , # reproduction rate per aphid per day
# K = , # aphid reproduction limit - maximum aphids per plant (diff to donnelly_stoch_sim.R)
# H = 400, # number host plants
# p = , # aphid emigration/death rate per journey between plants
# v = 1, # infected plant attractiveness
# e = 1, # infected plant acceptability 
# w = 0.2, # feeding rate on healthy plant
# Pacq = 0.5, # chance of virus acquisition from infected plant
# Pinoc = 0.5 # chance of inoculating healthy plant
# A = 1200 ????

# parameters
parms <- c(
  gamma = 3/20, # rate of recovery/death of I plants
  b = 0.1, # aphid mortality rate per day
  theta = 1, # aphid dispersal rate per day
  a = 2, # reproduction rate per aphid per day
  K = 10, # aphid reproduction limit - maximum aphids per plant (diff to donnelly_stoch_sim.R)
  H = 600, # number host plants
  p = 0.2, # aphid emigration/death rate per journey between plants
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  w = 0.2, # feeding rate on healthy plant
  Pacq = 0.5, # chance of virus acquisition from infected plant
  Pinoc = 0.5 # chance of inoculating healthy plant
)

#states
init_states <- c(
  i = 1/parms[["H"]] # frequency of infected plants
)

# timeframe
times <- seq(0, 8, by = 0.2)

########
# run epidemic with default parameters (as donnelly et al. 2019)
########
trajectory <- data.frame(ode(y = init_states, 
                             times = times, 
                             parms = parms, 
                             func = donnelly_ode_func))

trajectory_long <- reshape2::melt(trajectory, id.vars="time")
names(trajectory_long) <- c("time", "compartment", "number")

# plot trajectory of infected plants
plant_trajec <- ggplot(data=trajectory_long %>% filter(compartment == "i"), 
                       aes(x = time, y = number)) +
  geom_line() +
  labs(x = "Time (days)", y = "Frequency of infected plants, i")
plant_trajec

##########
# run epidemic with varying values of e and v
##########

#VARY E (plant acceptability)
e_vals <- c(0.25, 1, 2)

e_trajecs <- vary_param_trajec(init_states, times, parms, "e", e_vals)

# plot output
e_vals_plot <- ggplot(data=e_trajecs %>% filter(variable == "i"), aes(x = time, y = value, col = as.factor(param_val))) +
  geom_line() +
  labs(col = "e value")
e_vals_plot

# VARY V (plant attractiveness)
v_vals <- c(1, 3, 0.5)

v_trajecs <- vary_param_trajec(init_states, times, parms, "v", v_vals)

# plot output
v_vals_plot <- ggplot(data=v_trajecs %>% filter(variable == "i"), aes(x = time, y = value, col = as.factor(param_val))) +
  geom_line() +
  labs(col = "v value", y = "frequency of infected plants, i")
v_vals_plot
