
#################################################################################
##### EXTENSION OF MODEL BY MADDEN ET AL. 2000  TO INCLUDE RESISTANT PLANTS #####
#################################################################################

### Based on simplified version of model by Madden et al. 2000. Vector preference added in: virus-
### induced plant attractiveness, v, and acceptability, e. No vector immigration or emigration. 
### Plants are represented by an SI rather than an SEIR model, vectors represented by an SI rather 
### than SEI model. Vectors cannot be born infective (as this model represents NPT viruses).

### Extension to contain resistant plants: 2 versions. 1 = resistant plants have full resistance.
### 2 = resistant plants have partial resistance, level specified by parameter r.

rm(list = ls())

library(deSolve)
library(ggplot2)
library(gridExtra)
library(dplyr)

######## FULL RESISTANCE ###########

madden_fullr_rplant_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  S <- par[["H"]] - y[["I"]] - par[["R"]]
  
  become_infected <- par[["phi"]] * par[["b"]] * 
    y[["Z"]] * S / (S + par[["v_i"]]*y[["I"]] + par[["v_r"]]*par[["R"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  
  # VECTOR EQUATIONS
  X <- par[["A"]] - y[["Z"]]
  
  acquisition <- par[["phi"]] * par[["a"]] * 
    (1 - par[["e"]]*par[["w_mad"]]) * X * par[["v_i"]]*y[["I"]] / 
    (S + par[["v_i"]]*y[["I"]] + par[["v_r"]] * par[["R"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
}

madden_cunniffe_fullr_rplant_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  I <-  y[["I"]]
  R <- par[["R"]]
  S <- par[["H"]] - I - R
  
  # define variable phi
  phi <- (S + par[["v_i"]]*I + par[["v_r"]]*R) / 
    (par[["w_mad"]]*par[["eta"]] * (S + par[["v_i"]]*par[["e"]]*I + par[["v_r"]]*R))
  
  become_infected <- phi * par[["b"]] * 
  y[["Z"]] * S / (S + par[["v_i"]]*I + par[["v_r"]]*R)
  
  natural_death_I <- par[["c"]] * I
  virus_induced_death <- par[["d"]] * I
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  
  # VECTOR EQUATIONS
  X <- par[["A"]] - y[["Z"]]
  
  acquisition <- phi * par[["a"]] * 
    (1 - par[["e"]]*par[["w_mad"]]) * X * par[["v_i"]]*I / 
    (S + par[["v_i"]]*I + par[["v_r"]] * R)
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
}

donnelly_fullr_rplant_ode <- function(times, y, parms) {
  
  # STATES
  I <- y[["I"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  theta <- parms[["theta"]]
  H <- parms[["H"]]
  v_i <- parms[["v_i"]]
  v_r <- parms[["v_r"]]
  e <- parms[["e"]]
  w <- parms[["w_don"]]
  Pacq <- parms[["Pacq"]]
  Pinoc <- parms[["Pinoc"]]
  A <- parms[["A"]]
  R <- parms[["R"]]
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  
  # I_hat = weighted number of infected plants, accounting for attraction towards infected plants
  S <- H - I - R
  i_hat <- v_i*I / (S + v_r*R + v_i*I)
  r_hat <- v_r*R / (S + v_r*R + v_i*I)
  s_hat <- 1 - r_hat - i_hat
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  #xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  xI <- Pacq*Pinoc*
           i_hat*(1-e*w)*s_hat / ((s_hat+r_hat)*w + i_hat*e*w)

  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  dI <- (theta*A)*xI - gamma*I
  
  return(list(dI))
}


vary_param_AUDPC <- function(init_states, times, parms, varied_parm_name, varied_parm_vals, func, model) {
  
  # function to run a model ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of Area Under Disease Progress Curve (AUDPC) of all runs (proxy for epidemic 
  # severity), along with values of varied parameter

  # initialise output data frame
  output <- data.frame(parm_name = rep(varied_parm_name, length(varied_parm_vals)),
                       parm_val = varied_parm_vals,
                       AUDPC = rep(NA, length(varied_parm_vals)),
                       model = rep(model, length(varied_parm_vals))
  )
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states, 
                                 times = times, 
                                 parms = parms, 
                                 func = func))
    
    output[run, "AUDPC"] <- flux::auc(trajectory[["time"]], trajectory[["I"]])
  }
  
  return(output)
}

vary_param_finalI <- function(init_states, times, parms, varied_parm_name, varied_parm_vals, func, model) {
  
  # function to run a model ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of final disease incidence (I) of all runs, along with values of varied parameter
  
  # initialise output data frame
  output <- data.frame(parm_name = rep(varied_parm_name, length(varied_parm_vals)),
                       parm_val = varied_parm_vals,
                       final_I = rep(NA, length(varied_parm_vals)),
                       model = rep(model, length(varied_parm_vals))
  )
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states, 
                                 times = times, 
                                 parms = parms, 
                                 func = func))
    
    output[run, "final_I"] <- trajectory[nrow(trajectory),"I"]
  }
  
  return(output)
}

## parameters
parms_fullr <- c(
  
  tau = 4, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 5.5, # plants visited per day by an insect
  H = 400, # total number host plants 
  R = 100, # number resistant host plants
  v_i = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  v_r = 1.2, # resistant plant attractiveness
  #q = 0.8, # probability of not emigrating per flight
  #lamda = 2, # birth rate of vectors
  #k = 10, # carrying capacity of vectors per plant
  a = 1, # virus acquisition probability by vectors
  b = 1, # plant innoculation probability by vectors
  c = 2, # natural plant death rate
  d = 2, # plant death rate due to infection
  #alpha = 0.2 # vector death rate
  A = 1200, # total number of aphids
  
  # DONNELLY PARMS
  gamma = 4, # rate of recovery/death of I plants
  theta = 1.05, # aphid dispersal rate per day
  Pacq = 1, # chance of virus acquisition by vector from infected plant
  Pinoc = 1, # chance of inoculating healthy plant
  w_don = 0.2 # feeding probability on healthy plant
)
parms_fullr[["eta"]] <- 1/(parms_fullr[["phi"]]*parms_fullr[["w_mad"]]) # length of average feed

init_states_fullr_mad <- c(
  I = 1, # number of infected plants
  Z = 0 # number infective insects
)

init_states_fullr_don <- c(
  I = 1 # number of infected plants
)

times <- seq(0, 20, by = 0.4)

create_rplant_fig_data <- function(init_states_don, init_states_mad, times, parms,
                               varied_parm_name, varied_parm_vals, v_e_vals) {
  
  # function to create plots of a parameter (e.g. R, v_r) vs final I for pairs of v_i and e vals, for
  # all 3 models (don, mad, mad_cun). v_e_vals must be dataframe with columns named v_i and e, and rows
  # corresponding to pairs of values to be plotted
  
  all_plot_data <- data.frame(matrix(ncol = 6, nrow = 0))
  
  for (row in 1:nrow(v_e_vals)) {
    
    parms[["v_i"]] <- v_e_vals[row, "v_i"]
    parms[["e"]] <- v_e_vals[row, "e"]
    
    out_don <- vary_param_finalI(init_states_don, times, parms, 
                                varied_parm_name = varied_parm_name, 
                                varied_parm_vals = varied_parm_vals,
                                func = donnelly_fullr_rplant_ode,
                                model = "Donnelly")
    
    out_mad <- vary_param_finalI(init_states_mad, times, parms, 
                                varied_parm_name = varied_parm_name, 
                                varied_parm_vals = varied_parm_vals,
                                func = madden_fullr_rplant_ode,
                                model = "Madden")
    
    out_mad_cun <- vary_param_finalI(init_states_mad, times, parms, 
                                    varied_parm_name = varied_parm_name, 
                                    varied_parm_vals = varied_parm_vals,
                                    func = madden_cunniffe_fullr_rplant_ode,
                                    model = "Madden-Cunniffe")
    
    out_all_models <- rbind(out_don, out_mad, out_mad_cun)
    out_all_models <- cbind(out_all_models,
                            data.frame(v_i = rep(parms[["v_i"]], nrow(out_all_models)),
                                       e = rep(parms[["e"]], nrow(out_all_models))
                                       ))
    all_plot_data <- rbind(all_plot_data, out_all_models)
    
  }
  all_plot_data <- all_plot_data %>%
    mutate(v_i_e = paste0("v[i]*' = ", .$v_i, ", '*epsilon*' = ", .$e, "'"))
  return(all_plot_data)
}

R_vals <- seq(0, parms_fullr[["H"]], by = parms_fullr[["H"]]*0.05)
v_e_vals <- data.frame(v_i = c(2, 1, 3, 1.2),
                       e = c(0.5, 1, 3, 0.8))

R_plot_data <- create_rplant_fig_data(init_states_don = init_states_fullr_don,
                                  init_states_mad = init_states_fullr_mad,
                                  times = c(0, 40),
                                  parms = parms_fullr,
                                  varied_parm_name = "R",
                                  varied_parm_vals = R_vals,
                                  v_e_vals = v_e_vals)

R_plots <- ggplot(data = R_plot_data, aes(x = parm_val, y = final_I, color = model)) +
  facet_wrap(~v_i_e, label = "label_parsed") +
  geom_line() +
  scale_colour_manual(values = c("blue", "red", "green3")) + #don, mad, mad-cun
  theme_bw() +
  theme(text = element_text(size = 15),
        strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        legend.position = "none") +
  labs(x = "Number of resistant plants (R)",
       y = "Final disease incidence")

R_plots

v_r_vals <- seq(0, 5, by = 0.2)

v_r_plot_data <- create_rplant_fig_data(init_states_don = init_states_fullr_don,
                                        init_states_mad = init_states_fullr_mad,
                                        times = c(0,40),
                                        parms = parms_fullr,
                                        varied_parm_name = "v_r",
                                        varied_parm_vals = v_r_vals,
                                        v_e_vals = v_e_vals)

v_r_plots_w_legend <- ggplot(data = v_r_plot_data, aes(x = parm_val, y = final_I, color = model)) +
  facet_wrap(~v_i_e, label = "label_parsed") +
  geom_line() +
  scale_colour_manual(values = c("blue", "red", "green3")) + #don, mad, mad-cun
  theme_bw() +
  theme(text = element_text(size = 15),
        strip.background = element_blank(), 
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = expression("Attractiveness of resistant plants ("*v[r]*")"),
       y = "Final disease incidence")
#v_r_plots_w_legend

get_legend <- function(a.gplot){ # function to get legend from plot
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}

r_plots_legend <- get_legend(v_r_plots_w_legend)

v_r_plots <- v_r_plots_w_legend +
  theme(legend.position = "none")

layout_matrix <- matrix(c(1,2,3),
                        ncol = 1)

pdf("results/resistant_plants_plots.pdf", height = 10, width = 8)
grid.arrange(R_plots, v_r_plots, r_plots_legend,
             layout_matrix = layout_matrix,
             heights = c(1,1,0.2))
dev.off()

parms_new <- parms_fullr
parms_new[["v_r"]] <- 3
parms_new[["v_i"]] <- 1
parms_new[["e"]] <- 1

t <- data.frame(deSolve::ode(y = init_states_fullr_mad,
                             times = times,
                             func = madden_cunniffe_fullr_rplant_ode,
                             parms = parms_new))
t_d <- data.frame(deSolve::ode(y = init_states_fullr_don,
                               times = times,
                               func = donnelly_fullr_rplant_ode,
                               parms = parms_new))
par(mfrow=c(1,2))
plot(t$time, t$I, type = "l", ylim = c(0,350))
plot(t_d$time, t_d$I, type = "l", ylim = c(0,350))
par(mfrow=c(1,1))

######## PARTIAL RESISTANCE ###########

madden_partial_rplant_ode <- function(times, y, par) {
  
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

# parameters
parms_partialr <- c(

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

init_states_partialr <- c(
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

# run <- plot_ode(init_states, times, parms) 
# parms_no_resistance <- parms
# parms_no_resistance[["r"]] <- 1
# 
# parms_full_resistance <- parms
# parms_full_resistance[["r"]] <- 0
# 
# run_no_resistance <- plot_ode(init_states, times, parms_no_resistance)
# run_full_resistance <- plot_ode(init_states, times, parms_full_resistance)
# 
# grid.arrange(run[[1]], run_no_resistance[[1]], run_full_resistance[[1]])
# grid.arrange(run[[2]], run_no_resistance[[2]], run_full_resistance[[2]])
# 
# # plot trajectory of infective and healthy vectors
# 
# 
# grid.arrange(plant_trajec, vector_trajec)
# 
