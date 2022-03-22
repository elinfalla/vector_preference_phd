
rm(list = ls())

library(deSolve)
library(ggplot2)
library(gridExtra)

madden_rplant_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  Ss <- par[["Hs"]] - y[["Is"]] # number of susceptible plants
  Sr <- par[["Hr"]] - y[["Ir"]]
  
  become_infected_s <- par[["q"]] * par[["phi"]] * par[["b"]] * y[["Z"]] * Ss/
    (Ss + par[["v"]]*y[["Is"]] + Sr + par[["v"]]*y[["Ir"]])
  become_infected_r <- par[["r"]] * par[["q"]] * par[["phi"]] * par[["b"]] * y[["Z"]] * Sr/
    (Ss + par[["v"]]*y[["Is"]] + Sr + par[["v"]]*y[["Ir"]])
  
  natural_death_s <- par[["c"]] * y[["Is"]]
  natural_death_r <- par[["c"]] * y[["Ir"]]
  
  virus_induced_death_s <- par[["d"]] * y[["Is"]]
  virus_induced_death_r <- par[["d"]] * y[["Ir"]]
  
  # state equations
  dIs <- become_infected_s - natural_death_s - virus_induced_death_s
  dIr <- become_infected_r - natural_death_r - virus_induced_death_r
  
  # VECTOR EQUATIONS
  A <- y[["X"]] + y[["Z"]]
  H <- par[["Hs"]] + par[["Hr"]]
  
  acquisition <- par[["q"]] * par[["phi"]] * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * y[["X"]] * par[["v"]]*
    (y[["Is"]] + y[["Ir"]] * par[["r"]] - y[["Ir"]]*(1 - par[["r"]])) /
    (Ss + par[["v"]]*y[["Is"]] + Sr + par[["v"]]*y[["Ir"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  birth <- par[["lamda"]] * A * (1 - ((A/H)/par[["k"]]))
  death_X <- par[["alpha"]] * y[["X"]]
  death_Z <- par[["alpha"]] * y[["Z"]]
  
  
  # state equations
  dX <- birth - acquisition + stop_being_infective - death_X
  dZ <- acquisition - stop_being_infective - death_Z
  
  return(list(c(dIs, dIr, dX, dZ)))
}

times <- seq(0, 40, by = 0.4)

# parameters
parms <- c(

  tau = 4, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 5.5, # plants visited per day by an insect
  Hs = 400, # number susceptible host plants 
  Hr = 200, # number resistant host plants
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  q = 0.8, # probability of not emigrating per flight
  r = 0.2, # resistance parameter, large = less resistant
  lamda = 2, # birth rate of vectors
  k = 10, # carrying capacity of vectors per plant
  a = 0.8, # virus acquisition probability by vectors
  b = 0.8, # plant innoculation probability by vectors
  c = 2, # natural plant death rate
  d = 2, # plant death rate due to infection
  alpha = 0.2 # vector death rate
  
)

init_states <- c(
  Is = 1, # number of infected plants
  Ir = 0,
  X = 1200,
  Z = 0 # number infective insects
)

plot_ode <- function(init_states, times, parms) {
  
  trajectory_mad <- data.frame(ode(y = init_states,
                                   times = times,
                                   func = madden_rplant_ode,
                                   parms = parms))
  trajectory_mad$Ir <- trajectory_mad$Ir/(parms[["Hr"]]) # make proportion
  trajectory_mad$Is <- trajectory_mad$Is/(parms[["Hs"]])

  trajectory_mad_long <- reshape2::melt(trajectory_mad, id.vars = "time")
  names(trajectory_mad_long) <- c("time", "compartment", "number")
  
  # plot trajectory of infected plants and vector compartments
  plant_trajec <- ggplot(data = trajectory_mad_long %>% filter(compartment == "Ir")) +
    geom_line(aes(x = time, y = number, color = "Ir")) +
    geom_line(data = trajectory_mad_long %>% filter(compartment == "Is"),
              aes(x = time, y = number, color = "Is")) +
    labs(x = "Time (days)",
         y = "Proportion of infected plants per class") +
    scale_color_manual(breaks = c("Ir", "Is"),
                       values = c(Ir = "blue", Is = "red")) +
    theme(legend.title = element_blank()) +
    annotate("text", label = paste("r =", parms[["r"]]), x = times[length(times)], y = 0)
  
  vector_trajec <- ggplot(data = trajectory_mad_long %>% filter(compartment == "X")) +
    geom_line(aes(x = time, y = number, color = "X_healthy")) +
    geom_line(data = trajectory_mad_long %>% filter(compartment == "Z"), 
              aes(x = time, y = number, color = "Z_infective")) +
    labs(x = "Time (days)",
         y = "Number of vectors") +
    scale_color_manual(breaks = c("X_healthy", "Z_infective"),
                       values = c(X_healthy = "black", Z_infective = "green")) +
    theme(legend.title = element_blank()) +
    annotate("text", label = paste("r =", parms[["r"]]), x = times[length(times)], y = 0)
  
  
  return(list(plant_trajec, vector_trajec))
}

run <- plot_ode(init_states, times, parms) 
parms_no_resistance <- parms
parms_no_resistance[["r"]] <- 1

parms_full_resistance <- parms
parms_full_resistance[["r"]] <- 0

run_no_resistance <- plot_ode(init_states, times, parms_no_resistance)
run_full_resistance <- plot_ode(init_states, times, parms_full_resistance)

grid.arrange(run[[1]], run_no_resistance[[1]], run_full_resistance[[1]])
grid.arrange(run[[2]], run_no_resistance[[2]], run_full_resistance[[2]])

# plot trajectory of infective and healthy vectors


grid.arrange(plant_trajec, vector_trajec)

