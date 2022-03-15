
rm(list = ls())

library(deSolve)
library(ggplot2)

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
    (y[["Is"]] + y[["Ir"]] * par[["r"]]) /
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

  tau = 1.25, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 3.5, # plants visited per day by an insect
  Hs = 400, # number susceptible host plants 
  Hr = 200, # number resistant host plants
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  q = 0.8, # probability of not emigrating per flight
  r = 0.2, # resistance paramter, large = less resistant
  lamda = 0, # birth rate of vectors
  k = 10, # carrying capacity of vectors per plant
  a = 0.8,
  b = 0.8,
  c = 2,
  d = 2,
  alpha = 0.2 # death rate
  
)

init_states <- c(
  Is = 1, # number of infected plants
  Ir = 0,
  X = 1200,
  Z = 0 # number infective insects
)

trajectory_mad <- data.frame(ode(y = init_states,
                                 times = times,
                                 func = madden_rplant_ode,
                                 parms = parms))
#trajectory_mad$I <- trajectory_mad$I/parms[["H"]] # make proportion


trajectory_mad_long <- reshape2::melt(trajectory_mad, id.vars = "time")
names(trajectory_mad_long) <- c("time", "compartment", "number")

# plot trajectory of infected plants and vector compartments
plant_trajec <- ggplot(data = trajectory_mad_long %>% filter(compartment == "Ir"),
                       aes(x = time, y = number)) +
  geom_line() +
  geom_line(data = trajectory_mad_long %>% filter(compartment == "Is"), color = "red") +
  labs(x = "Time (days)",
       y = "Number of infected plants, I")

plant_trajec


