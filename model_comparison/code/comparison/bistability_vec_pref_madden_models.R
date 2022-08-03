###################################################################
#### DEMONSTRATION OF BISTABLE BEHAVIOUR IN MADDEN-TYPE MODELS ####
###################################################################

# Simple code demonstrating the bistable behaviour seen in the Madden and 
# Madden-Cunniffe models. V and e are set so the model is in the region of
# bistability. These models have vector preference dynamics but no vector 
# population dynamics.

rm(list=ls())

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gridtext)
library(scales)

#### DEFINE MODEL ODES ####
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
  
  acquisition <- par[["phi"]] * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * X * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
}

madden_cunniffe_vpref_ode <- function(times, y, par) {
  
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  # define variable phi
  phi <- (S + y[["I"]]*par[["v"]]) / (par[["w_mad"]]*par[["eta"]]*(S + y[["I"]]*par[["v"]]*par[["e"]]))
  
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
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
}

#### DEFINE TIMESCALES, PARAMETERS, INITAL CONDITIONS ####

# define timeframes
times_short <- seq(0, 1, length.out = 200)
times_med <- seq(0, 15, length.out = 200)
times_long <- seq(0, 1000, length.out = 2000)

# parameters
parms <- c(
  tau = 4, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 5.4730983, # plants visited per day by an insect
  A = 1200, # total (constant) number of vectors
  H = 400, # number host plants 
  v = 0.3039697, # infected plant attractiveness - value at which bistability is seen
  e = 0.4554545, # infected plant acceptability  - value at which bistability is seen
  a = 0.8, # virus acquisition rate
  b = 0.8, # virus inoculation rate
  c = 2, # natural plant death rate
  d = 2 # virus-induced plant death rate
)
parms[["eta"]] <- 1/(parms[["phi"]]*parms[["w_mad"]]) # average length of visit to a plant


# initial I states
all_init_states_mad <- seq(1, 0.8*parms[["H"]], length.out = 30)


#### DEFINE FUNCTION TO PLOT MANY TRAJECTORIES WITH DIFFERENT INITIAL CONDITIONS #### 
plot_many_trajec_madden <- function(init_I_states, times, parms, func) {
  
  ## function that for given times, parms, and madden-type ode function, runs
  ## the ode for many initial starting I states (initial Z = 0) and plots 
  ## proportion I vs time, proportion Z vs time, and proportion I vs proportion Z
  
  # use different colours for I vs Z plot to more easily differentiate between lines
  colours <- c("red", "blue", "darkgreen", "black") 
  
  # use for loop to loop over initial I values and run ode to plot
  for (state in init_I_states) {
    
    # run ode
    trajec <- data.frame(ode(y = c(I = state, Z = 0),
                             times = times,
                             parms = parms,
                             func = func))
    trajec$I <- trajec$I/parms[["H"]] # make proportion
    trajec$Z <- trajec$Z/parms[["A"]]
    
    
    if (state == init_I_states[1]) {
      
      # plot trajectory of infected plants vs time
      I_trajec_plot <- ggplot(data = trajec, aes(x = time, y = I)) +
        geom_line() +
        labs(x = "Time (arbitrary units)",
             y = "Proportion of infected plants")
      
      # plot trajectory of infective aphids vs time
      Z_trajec_plot <- ggplot(data = trajec, aes(x = time, y = Z)) +
        geom_line() +
        labs(x = "Time (arbitrary units)",
             y = "Proportion of infective aphids")
      
      # plot trajectory of infected plants vs infective aphids
      I_vs_Z_plot <- ggplot(data = trajec, aes(x = I, y = Z)) +
        geom_line(colour = sample(colours, 1)) +
        labs(x = "Proportion of infected plants",
             y = "Proportion of infective aphids")
    } 
    else { # once plot created, add subsequent trajectories (different initial I) to the plots
      I_trajec_plot <- I_trajec_plot +
        geom_line(data = trajec, aes(x = time, y = I))
      
      Z_trajec_plot <- Z_trajec_plot +
        geom_line(data = trajec, aes(x = time, y = Z))
      
      I_vs_Z_plot <- I_vs_Z_plot +
        geom_line(data = trajec, aes(x = I, y = Z), colour = sample(colours, 1))
    }
    
  }
  return(list(I_trajec_plot, Z_trajec_plot, I_vs_Z_plot)) # return plots
}

#### DEMONSTRATE BISTABLE BEHAVIOUR - tmax = 1 ####

# madden plots
trajec_plot_mad_zoom <- plot_many_trajec_madden(all_init_states_mad,
                                                times = times_short, 
                                                parms = parms,
                                                func = madden_vpref_ode)
trajec_plot_mad_zoom[[1]] <- trajec_plot_mad_zoom[[1]] +
  labs(title = "Madden - tmax = 1") # add title

# madden-cunniffe plots
trajec_plot_mad_cun_zoom <- plot_many_trajec_madden(all_init_states_mad,
                                                    times = times_short,
                                                    parms = parms,
                                                    func = madden_cunniffe_vpref_ode)
trajec_plot_mad_cun_zoom[[1]] <- trajec_plot_mad_cun_zoom[[1]] +
  labs(title = "Madden-Cunniffe - zoomed (tmax=1)") # add title

# display all together
lay <- matrix(c(1,2,3,4,5,6), nrow = 3, byrow = F)
grid.arrange(trajec_plot_mad_cun_zoom[[1]],
             trajec_plot_mad_cun_zoom[[2]],
             trajec_plot_mad_cun_zoom[[3]],
             trajec_plot_mad_zoom[[1]],
             trajec_plot_mad_zoom[[2]],
             trajec_plot_mad_zoom[[3]],
             ncol = 2, layout_matrix = lay)

# pdf("results/bistable_trajecs_madden-type_zoom.pdf")
# grid.arrange(trajec_plot_mad_cun_zoom[[1]],
#              trajec_plot_mad_cun_zoom[[2]],
#              trajec_plot_mad_cun_zoom[[3]],
#              trajec_plot_mad_zoom[[1]],
#              trajec_plot_mad_zoom[[2]],
#              trajec_plot_mad_zoom[[3]],
#              ncol = 2, layout_matrix = lay)
# dev.off()

#### DEMONSTRATE BISTABLE BEHAVIOUR - tmax = 1000 ####

# madden plots
trajec_plot_mad_long <- plot_many_trajec_madden(all_init_states_mad,
                                                times = times_long,
                                                parms = parms,
                                                func = madden_vpref_ode)
trajec_plot_mad_long[[1]] <- trajec_plot_mad_long[[1]] + labs(title = "Madden - tmax = 1000")

# madden-cunniffe plots
trajec_plot_mad_cun_long <- plot_many_trajec_madden(all_init_states_mad,
                                                    times = times_long,
                                                    parms = parms,
                                                    func = madden_cunniffe_vpref_ode)
trajec_plot_mad_cun_long[[1]] <- trajec_plot_mad_cun_long[[1]] + labs(title = "Madden-Cunniffe - tmax = 1000")

# display together
grid.arrange(trajec_plot_mad_cun_long[[1]],
             trajec_plot_mad_cun_long[[2]],
             trajec_plot_mad_cun_long[[3]],
             trajec_plot_mad_long[[1]],
             trajec_plot_mad_long[[2]],
             trajec_plot_mad_long[[3]],
             ncol = 2, layout_matrix = lay)

# pdf("results/bistable_trajecs_madden-type_long.pdf")
# grid.arrange(trajec_plot_mad_cun_long[[1]],
#              trajec_plot_mad_cun_long[[2]],
#              trajec_plot_mad_cun_long[[3]],
#              trajec_plot_mad_long[[1]],
#              trajec_plot_mad_long[[2]],
#              trajec_plot_mad_long[[3]],
#              ncol = 2, layout_matrix = lay)
# dev.off()

#### DEMONSTRATE BISTABLE BEHAVIOUR - tmax = 15 ####

# madden plots
trajec_plot_mad <- plot_many_trajec_madden(all_init_states_mad,
                                           times = times_med,
                                           parms = parms,
                                           func = madden_vpref_ode)
trajec_plot_mad[[1]] <- trajec_plot_mad[[1]] + labs(title = "Madden - tmax = 15")

# madden-cunniffe plots
trajec_plot_mad_cun <- plot_many_trajec_madden(all_init_states_mad,
                                               times = times_med,
                                               parms = parms,
                                               func = madden_cunniffe_vpref_ode)
trajec_plot_mad_cun[[1]] <- trajec_plot_mad_cun[[1]] + labs(title = "Madden-Cunniffe - tmax = 15")

# display together
grid.arrange(trajec_plot_mad_cun[[1]],
             trajec_plot_mad_cun[[2]],
             trajec_plot_mad_cun[[3]],
             trajec_plot_mad[[1]],
             trajec_plot_mad[[2]],
             trajec_plot_mad[[3]],
             ncol = 2, layout_matrix = lay)

# pdf("results/bistable_trajecs_madden-type_short.pdf")
# grid.arrange(trajec_plot_mad_cun[[1]],
#              trajec_plot_mad_cun[[2]],
#              trajec_plot_mad_cun[[3]],
#              trajec_plot_mad[[1]],
#              trajec_plot_mad[[2]],
#              trajec_plot_mad[[3]],
#              ncol = 2, layout_matrix = lay)
# dev.off()


# compare shorter and longer timeframes of Z vs time for madden
#grid.arrange(trajec_plot_mad_long[[2]], trajec_plot_mad[[2]])

# compare mad and mad-cun plots of Z vs time
#grid.arrange(trajec_plot_mad_long[[2]], trajec_plot_mad[[2]])
