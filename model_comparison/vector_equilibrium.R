
###################################################################################
###### DETERMINING APHID POPULATION EQUILIBRIUM -  FROM DONNELLY ET AL. 2019 ######
###################################################################################

### Script to test methods of finding an equilibrium value for dAs/dt and dAi/dt
### (rate of change of per plant aphid density on susceptible and infected plants
### respectively). Methods: ODE solve, optim(), polyroot(), fsolve().

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)


##############
### Aphid population ODE solve
##############

parms <- c( # IF THESE ARE CHANGED THE EQUILIBRIUM VALUES WILL CHANGE
  i = 0.1, # proportion of host plants that are infected
  v = 5, # infected plant attractiveness
  e = 3, # infected plant acceptability 
  a = 2, # reproduction rate per aphid per day
  b = 0.1, # aphid mortality rate per day
  K = 10, # aphid reproduction limit - maximum aphids per plant (diff to donnelly_stoch_sim.R)
  theta = 1 # aphid dispersal rate per day
)

init_states <- c(
  As = 12, # per plant aphid density on average healthy plant
  Ai = 12 # per plant aphid density on average infected plant
)

times <- seq(0, 10, by = 0.2)

aphid_pop_ode <- function(times, y, parms) {
  
  # define states
  As <- y[["As"]]
  Ai <- y[["Ai"]]
  
  # define parameters
  i <- parms[["i"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  a <- parms[["a"]]
  b <- parms[["b"]]
  K <- parms[["K"]]
  theta <- parms[["theta"]]
  
  # i_hat = weighted frequency of infected plants, accounting for attraction towards infected plants
  # i.e. adapted version of i / (1 - i + i) = i / 1, to include v
  i_hat <- v*i / ((1-i) + v*i) 
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  # STATE EQUATIONS
  # aphids - change in per plant aphid density on susceptible (As)  and infected (Ai) plants
  # = aphid births - aphid deaths - aphid settling on other plant type + aphids moving from other plant type
  dAs <- a*As*(1 - As/K) - b*As - theta*As*(1 - Fs) + theta*Ai*Fs*i/(1 - i)
  dAi <- a*Ai*(1 - Ai/K) - b*Ai - theta*Ai*(1 - Fi) + theta*As*Fi*(1 - i)/i
  
  return(list(c(dAs, dAi)))
}

vary_init_states <- function(init_states_list, times, parms) {
  
  # function to run the ODE multiple times with the initial states varying, as specified by user.
  # returns data frame of all runs. init_states_list is a list of vectors with the initial conditions
  # for each run
  
  # initialise output data frame
  all_trajects <- data.frame(time = numeric(),
                             compartment = character(),
                             value = numeric(),
                             run_num = numeric())
  
  for (run in 1:length(init_states_list)) {
    
    trajectory <- data.frame(ode(y = init_states_list[[run]], 
                                 times = times, 
                                 parms = parms, 
                                 func = aphid_pop_ode))
    
    trajectory_long <- reshape2::melt(trajectory, id.vars="time")

    trajectory_long$run_num <- rep(run, nrow(trajectory_long))
      
    # add to dataframe of all trajectories
    all_trajects <- rbind(all_trajects, trajectory_long)
  }
  
  return(all_trajects)
}

init_states_list <- list(c(As = 1, Ai = 0),
                        c(As = 3, Ai = 2),
                        c(As = 5, Ai = 4),
                        c(As = 7, Ai = 6),
                        c(As = 9, Ai = 8),
                        c(As = 11, Ai = 10),
                        c(As = 13, Ai = 12),
                        c(As = 15, Ai = 14),
                        c(As = 17, Ai = 16))


aphid_trajecs <- vary_init_states(init_states_list, times, parms)



aphid_plot <- ggplot(data = aphid_trajecs %>% filter(run_num == 1), aes(x = time, y = value, col = variable)) +
  geom_line() +
  labs(title = "Aphid population ODE solve")

num_runs <- length(init_states_list)

for (run in 1:(num_runs - 1)) {
  aphid_plot <- aphid_plot +
    geom_line(data = aphid_trajecs %>% filter(run_num == run))
    
}

aphid_plot

## Check that the equilibrium doesn't change for varying i
i_vals <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.98)

vary_i_vals <- function(parms, init_states, times, i_vals) {
  
  ### Function that runs ODE with varying values of i (proportion infected plants) and checks the
  ### equilibrium values of As and Ai don't change
  
  end_vals <- vector(mode = "list", length = length(i_vals))
  
  for (run in 1:length(i_vals)) {
    
    parms[["i"]] <- i_vals[run]

    trajectory <- data.frame(ode(y = init_states, 
                                 times = times, 
                                 parms = parms, 
                                 func = aphid_pop_ode))
    end_As <- round(trajectory[nrow(trajectory), "As"], 3)
    end_Ai <- round(trajectory[nrow(trajectory), "Ai"], 3)
    
    end_vals[[run]] <- c(end_As, end_Ai)
  
    trajectory_long <- reshape2::melt(trajectory, id.vars = "time")
    names(trajectory_long) <- c("time", "compartment", "value")
    
    # create iterative plot that plots trajectories for each i value
    trajec_colours <- rainbow(length(i_vals))
    names(trajec_colours) <- as.character(i_vals)

    if (run == 1) {
      trajecs_plot <- ggplot(data = trajectory_long, aes(x = time, y = value, group = compartment)) +
        geom_line(color = trajec_colours[run]#as.character(i_vals[run])
                  )
    } else {
      trajecs_plot <- trajecs_plot +
        geom_line(data = trajectory_long, aes(x = time, y = value, group = compartment), color = trajec_colours[run]
                  #as.character(i_vals[run])
                  )
    }
   }
  # # give plot a legend
  # trajecs_plot <- trajecs_plot +
  #   scale_color_manual(name = 'i value',
  #                      breaks = as.character(i_vals),
  #                      values = trajec_colours)
  
  # narrow down end vals (equilibriums) to only unique values
  unique_equilibriums <- end_vals[!duplicated(end_vals)]
  
  # print message about number of equilibriums found
  if (length(unique_equilibriums) == 1) {
    print("Equilibrium doesn't change")
  } else{
    print(paste0("Equilibrium changes: number i vals tested = ", length(i_vals), 
                 ", number equilibriums = ", length(unique_equilibriums)))
  }
  
  return(list(unique_equilibriums, trajecs_plot))
}

vary_i_equilibriums <- vary_i_vals(parms, init_states, times, i_vals)
vary_i_equilibriums[[1]] # equilibrium values
vary_i_equilibriums[[2]] # plot of all trajectories

##############
### Find aphid population equilibrium using optim()
##############

RHS_to_optimise <- function(aphid_pop_density, parms) {
  
  ### Function that contains the equations for As and Ai to be optimised. Will be passed to optim()
  
  As <- aphid_pop_density[1]
  Ai <- aphid_pop_density[2]
  
  # define parameters
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

# run optim() many times from different starting values to find all possible equilibriums
n_tests <- 100
optim_equilibrium <- lapply(1:n_tests, function(x) round(optim(par = runif(2, min=0, max=100), 
                                                fn = RHS_to_optimise, 
                                                parms = parms,
                                                method = "L-BFGS-B")$par, 4))

# filter down to unique equilibriums
unique_equilibrium <- optim_equilibrium[!duplicated(optim_equilibrium)]
unique_equilibrium # only equilibrium where both numbers (As and Ai) are positive = 9.5, 9.5

# filter down to the positive (i.e. biologically possible) solution(s)
find_positive_solution <- function(eq) {
  if (eq[1] > 0 && eq[2] > 0) {
    return (T)
  }
  else {
    return(F)
  }
}

positive_equilibrium <- unlist(unique_equilibrium[sapply(unique_equilibrium, find_positive_solution)==T])
positive_equilibrium


##############
### Find aphid population equilibrium using polyroot()
##############

# define parameters
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

## coefficients (worked out by hand from dAs and dAi equations- see lab book). will give roots of Ai
coeff_1 <- i*(a - b - theta*(1 - Fs))*(b + theta*(1 - Fi) - a)/(theta*Fi*(1 - i)) + theta*Fs*i/(1 - i)

coeff_2 <- (a/K)*i*(a - b - theta*(1 - Fs))/(theta*Fi*(1 - i)) - 
  a*i^2*(b^2 + 2*b*theta*(1 - Fi) - 2*b*a + theta^2*(1 - Fi)^2 - 2*a*theta*(1 - Fi) + a^2)/(K*theta^2*(1 - i)^2*Fi^2)

coeff_3 <- -(a*i^2/K)*(2*b*a + 2*a*theta*(1 - Fi) - 2*a^2)/(K*theta^2*(1 - i)^2*Fi^2)

coeff_4 <- -a^3*i^2/(K^3*theta^2*(1 - i)^2*Fi^2)

proot_equilibrium <- polyroot(c(0, coeff_1, coeff_2, coeff_3, coeff_4))
proot_equilibrium # same as equilibrium Ai values for optim() method

##############
### Find aphid population equilibrium using fsolve()
##############

library(pracma) # contains fsolve() function

fsolve(RHS_to_optimise, c(2,2), parms = parms)

# run fsolve() many times from different starting values to find all possible equilibriums
n_tests <- 100
fsolve_equilibrium <- lapply(1:n_tests, function(x) round(fsolve(RHS_to_optimise, 
                                                                 runif(2, min=0, max=20), 
                                                                 parms = parms,
                                                                 maxiter = 300)$x, 4))

# filter down to unique equilibriums
fsolve_unique_equilibrium <- fsolve_equilibrium[!duplicated(fsolve_equilibrium)]
fsolve_unique_equilibrium # only equilibrium where both numbers (As and Ai) are positive = 9.5, 9.5

