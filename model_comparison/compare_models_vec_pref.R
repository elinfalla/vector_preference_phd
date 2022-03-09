############################################################################################
##### COMPARISON OF VECTOR PREFERENCE IN MADDEN (2000) AND DONNELLY (2019) TYPE MODELS #####
############################################################################################

### MADDEN MODEL:
### Modified version of model by Madden et al. 2000 to include vector preference - both of plant attractiveness and
### plant acceptability. No vector population dynamics (immigration, emigration birth, death), plants are represented 
### by an SI rather than an SEIR model, vectors represented by an SI rather than SEI model. Vectors cannot be born 
### infective (as this model represents NPT viruses).

### DONNELLY MODEL:
### Modification of deterministic model by Donnelly et al. 2019 to make it more comparable to this version of the 
### Madden et al. 2000 model. Vector preference included but no vector population dynamics - constant vector population 
### size, no vector emigration (p always set to 0). Also model records number of infected plants, not proportion.

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gridtext)

#### DEFINE MADDEN ODE
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

#### DEFINE DONNELLY ODE
donnelly_vpref_ode <- function(times, states, parms) {
  
  # STATES
  I <- states[["I"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  theta <- parms[["theta"]]
  H <- parms[["H"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w_don"]]
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
  di <- (theta*A)*xI - gamma*I
  
  return(list(di))
}

# define timeframe
times <- seq(0, 40, by = 0.4)

# parameters
parms <- c(
  
  ## DONNELLY PARMS
  gamma = 4, # rate of recovery/death of I plants
  theta = 1.05, # aphid dispersal rate per day
  Pacq = 0.8, # chance of virus acquisition by vector from infected plant
  Pinoc = 0.8, # chance of inoculating healthy plant
  w_don = 0.2, # feeding probability on healthy plant
  
  ## MADDEN PARMS
  # k1 = 1/0.021, # inoculation rate by insect (set to NPT virus, equiv to 0.5hr)
  # lamda = 1/0.021, # rate of acquisition of virus by vector (set to NPT virus, equiv to 0.5hr)
  # T = 0.2, # time feeding/probing per plant visit (set to 0.5/phi in Madden model)
  tau = 1.25, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 3.5, # plants visited per day by an insect
  
  ## PARMS FOR BOTH
  A = 1200, # total (constant) number of vectors
  H = 400, # number host plants 
  v = 1, # infected plant attractiveness
  e = 1 # infected plant acceptability 

)


### SET MADDEN PARMS THAT ARE DETERMINED BY DONNELLY PARMS
# probability of vector acquisition of virus from an infected plant per visit
parms[["a"]] <- parms[["Pacq"]] #1 - exp(-parms[["lamda"]] * parms[["T"]])

# probability of plant inoculation per infective insect visit 
parms[["b"]] <- parms[["Pinoc"]] #1 - exp(-parms[["k1"]] * parms[["T"]])

# plants visited per day by an insect
# parms[["phi"]] <- #parms[["theta"]] 
#   sqrt(parms[["theta"]] * parms[["tau"]] * (1 - parms[["w_don"]]) /
#               (parms[["w_don"]] #* parms[["a"]] * parms[["b"]]
#               ))

# rate of moving through infectious state in vector
# parms[["tau"]] <- parms[["theta"]]

# natural plant death rate (equivalent to beta in original Madden model)
parms[["c"]] <- parms[["gamma"]] / 2

# plant death due to infection
parms[["d"]] <- parms[["gamma"]] / 2

# feeding probability on healthy plant
# parms[["w_mad"]] <- (parms[["w_don"]]*(1 + parms[["e"]]) - 1) / (parms[["e"]]*parms[["w_don"]]) 
# 
# parms[["w_mad"]] <- ((1 - parms[["e"]]*parms[["w_don"]])*(parms[["tau"]] + parms[["phi"]]) - 
#                     parms[["phi"]]*(parms[["w_don"]] - parms[["w_don"]]*(1 - parms[["e"]]))) /
#   ((1 - parms[["e"]]*parms[["w_don"]])*parms[["phi"]]*parms[["e"]] - 
#   parms[["e"]]*(parms[["w_don"]] - parms[["w_don"]]*(1 - parms[["e"]])))


# STATES
I <- 1 # initial number of infected plants

init_states_don <- c(
  I = I # number of infected plants
)

init_states_mad <- c(
  I = I, # number of infected plants
  Z = 0 # number infective insects
)


### RUN DONNELLY EPIDEMIC
trajectory_don <- data.frame(ode(y = init_states_don,
                             times = times,
                             parms = parms,
                             func = donnelly_vpref_ode))
trajectory_don$I <- trajectory_don$I/parms[["H"]] # make proportion

# plot trajectory of infected plants and number of aphids
plant_trajec <- ggplot(data = trajectory_don, aes(x = time, y = I)) +
  geom_line() +
  labs(x = "Time (days)",
       y = "Number of infected plants, I")

### RUN MADDEN EPIDEMIC
trajectory_mad <- data.frame(ode(y = init_states_mad,
                      times = times,
                      func = madden_vpref_ode,
                      parms = parms))
trajectory_mad$I <- trajectory_mad$I/parms[["H"]] # make proportion


trajectory_mad_long <- reshape2::melt(trajectory_mad, id.vars = "time")
names(trajectory_mad_long) <- c("time", "compartment", "number")

# plot trajectory of infected plants and vector compartments
plant_trajec <- plant_trajec +
  geom_line(data = trajectory_mad_long %>% filter(compartment == "I"),
            aes(x = time, y = number), color = "red") +
  annotate("text", label = "red = Madden\nblack = Donnelly", x = times[length(times)/3*2], y = 0)
plant_trajec


vary_param_don <- function(init_states_don, times, parms, varied_parm_name, varied_parm_vals) {
  
  # function to run the Donnelly ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of final incidence, I, of all runs, along with values of varied parameter
  
  # initialise output data frame
  output <- data.frame(parm_name = rep(varied_parm_name, length(varied_parm_vals)),
                       parm_val = varied_parm_vals,
                       final_I = rep(NA, length(varied_parm_vals)),
                       model = rep("Donnelly", length(varied_parm_vals))
  )
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states_don, 
                                 times = times, 
                                 parms = parms, 
                                 func = donnelly_vpref_ode))
    output[run, "final_I"] <- round(trajectory[nrow(trajectory), "I"], 3)
    
  }
  
  return(output)
}

vary_param_mad <- function(init_states_mad, times, parms, varied_parm_name, varied_parm_vals) {
  
  # function to run the Madden ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of final incidence, I, of all runs, along with values of varied parameter
  
  # initialise output data frame
  output <- data.frame(parm_name = rep(varied_parm_name, length(varied_parm_vals)),
                       parm_val = varied_parm_vals,
                       final_I = rep(NA, length(varied_parm_vals)),
                       model = rep("Madden", length(varied_parm_vals))
  )
  
  for (run in 1:length(varied_parm_vals)) {
    
    parms[[varied_parm_name]] <- varied_parm_vals[run]
    
    trajectory <- data.frame(ode(y = init_states_mad, 
                                 times = times, 
                                 parms = parms, 
                                 func = madden_vpref_ode))
    output[run, "final_I"] <- round(trajectory[nrow(trajectory), "I"], 3)
    
  }
  
  return(output)
}

sensitivity_analysis <- function(parms_mad, parms_don, init_states_don, init_states_mad, times, default_parms) {
  
  sen_analysis_df <- data.frame(parm_name = character(),
                                parm_val = numeric(),
                                final_I = numeric(),
                                model = character())
  
  print("Starting Madden sensitivity analysis")
  
  for (i in 1:length(parms_mad)) {
    
    output <- vary_param_mad(init_states_mad, times, default_parms,
                             varied_parm_name = names(parms_mad[i]),
                             varied_parm_vals = parms_mad[[i]]) 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output)
    
  }
  
  print("Starting Donnelly sensitivity analysis")
  for (i in 1:length(parms_don)) {
    
    output <- vary_param_don(init_states_don, times, default_parms,
                             varied_parm_name = names(parms_don[i]),
                             varied_parm_vals = parms_don[[i]]) 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output)
    
  }
  
  return(sen_analysis_df) 
}

num_parm_runs <- 50
analysis_parms_don <- list(theta = seq(0, 20, length.out = num_parm_runs),
                           w_don = seq(0.01, 1, length.out = num_parm_runs),
                           Pacq = seq(0, 1, length.out = num_parm_runs),
                           Pinoc = seq(0, 1, length.out = num_parm_runs),
                           gamma = seq(0, 20, length.out = num_parm_runs),
                           v = seq(0, 20, length.out = num_parm_runs),
                           e = seq(0, 20, length.out = num_parm_runs))

analysis_parms_mad <- list(phi = seq(0, 20, length.out = num_parm_runs),
                           tau = seq(0, 20, length.out = num_parm_runs),
                           w_mad = seq(0.01, 1, length.out = num_parm_runs),
                           a = seq(0, 1, length.out = num_parm_runs),
                           b = seq(0, 1, length.out = num_parm_runs),
                           c = seq(0, 20, length.out = num_parm_runs),
                           d = seq(0, 20, length.out = num_parm_runs),
                           v = seq(0, 20, length.out = num_parm_runs),
                           e = seq(0, 20, length.out = num_parm_runs))


sens_analysis_res <- sensitivity_analysis(parms_mad = analysis_parms_mad,
                                          parms_don = analysis_parms_don,
                                          init_states_don, 
                                          init_states_mad, 
                                          times, 
                                          parms)

###  CREATE DONNELLY AND MADDEN MODEL SENSITIVITY ANALYSIS PLOTS
don_res_only <- sens_analysis_res %>%
  filter(model == "Donnelly")
mad_res_only <- sens_analysis_res %>%
  filter(model == "Madden")

don_model_plots <- lapply(unique(don_res_only$parm_name), function(p) ggplot(data = don_res_only[don_res_only$parm_name == p,],
                                                                             aes(x = parm_val, 
                                                                                 y = final_I)) +
                            geom_line() +
                            labs(x = p))

mad_model_plots <- lapply(unique(mad_res_only$parm_name), function(p) ggplot(data = mad_res_only[mad_res_only$parm_name == p,],
                                                                             aes(x = parm_val, 
                                                                                 y = final_I)) +
                            geom_line() +
                            labs(x = p))

# create table to be given alongside graphs giving default parameter values
parms_table <- round(data.frame(parms[names(parms) != "k1" &
                                        names(parms) != "lamda" &
                                        names(parms) != "T"]), 2)
parms_grob <- tableGrob(parms_table, cols = c("Value")) # turn into grob (gtable) for plotting

# ARRANGE PLOTS (AND PARAMETER TABLE) AND SAVE TO PDF
don_plots <- gridExtra::arrangeGrob(don_model_plots[[1]], 
                                    don_model_plots[[2]], 
                                    grid.rect(gp=gpar(col="white")), # empty space
                                    don_model_plots[[3]],
                                    don_model_plots[[4]],
                                    don_model_plots[[5]],
                                    grid.rect(gp=gpar(col="white")), # empty space
                                    don_model_plots[[6]],
                                    don_model_plots[[7]],
                                    ncol = 1)

mad_plots <- gridExtra::arrangeGrob(mad_model_plots[[1]], 
                                    mad_model_plots[[2]], 
                                    mad_model_plots[[3]],
                                    mad_model_plots[[4]],
                                    mad_model_plots[[5]], 
                                    mad_model_plots[[6]],
                                    mad_model_plots[[7]],
                                    mad_model_plots[[8]],
                                    mad_model_plots[[9]],
                                    ncol = 1)

layout <- rbind(c(1,2,3),
                c(1,2,4),
                c(1,2,4),
                c(1,2,4))
title <- textbox_grob("Model comparison with vector preference - constant vector population", 
                      gp = gpar(fontface = "bold",
                                fontsize = 13),
                      padding = unit(c(0, 1, 0, 1), "cm"))

# create pdf file to print plot to
pdf(file = "sens_analysis_vec_pref_models.pdf", width = 9, height = 12)
all_plots <- gridExtra::grid.arrange(don_plots,
                                     mad_plots,
                                     title,
                                     parms_grob,
                                     layout_matrix = layout,
                                     widths = c(2,2,1))
all_plots
dev.off()

##### CREATE HEATMAPS OF V VS E for final incidence of disease

find_I_equilibrium <- function(parms, model, mult_eq_return_NA=T) {
  
  ## function that for given parameter values, finds the non-zero equilibrium value of I,
  ## number of infected plants, for either the Donnelly model or the Madden model
  
  if (model == "Donnelly") {
    gamma <- parms[["gamma"]]
    theta <- parms[["theta"]]
    H <- parms[["H"]]
    v <- parms[["v"]]
    e <- parms[["e"]]
    w <- parms[["w_don"]]
    Pacq <- parms[["Pacq"]]
    Pinoc <- parms[["Pinoc"]]
    A <- parms[["A"]]
    
    # coefficients worked out by hand (see lab book)
    coeff_0 <- theta*A*v*(1 - e*w)*H*Pacq*Pinoc - gamma*w*H^2
    coeff_1 <- gamma*w*v*(1 - e)*H + 2*gamma*w*H - 2*gamma*w*v*H - theta*A*v*(1 - e*w)*Pacq*Pinoc
    coeff_2 <- gamma*w*(v^2*(1 - e) - v*(1 - e) - v^2 + 2*v - 1)
  }
  else if (model == "Madden") {
    c <- parms[["c"]]
    d <- parms[["d"]]
    phi <- parms[["phi"]]
    H <- parms[["H"]]
    v <- parms[["v"]]
    e <- parms[["e"]]
    w <- parms[["w_mad"]]
    a <- parms[["a"]]
    b <- parms[["b"]]
    A <- parms[["A"]]
    tau <- parms[["tau"]]
    
    coeff_0 <- (c + d)*tau*H^2 - phi^2*b*H*a*(1 - e*w)*A*v
    coeff_1 <- (c + d)*(2*tau*H*v - 2*tau*H + phi*a*(1 - e*w)*v*H) + phi^2*b*a*(1 - e*w)*A*v
    coeff_2 <- (c + d)*(tau + tau*v^2 - 2*tau*v - phi*a*(1 - e*w)*v + phi*a*(1 - e*w)*v^2)
  }
  else if (model == "Madden_Cunniffe") {
    c <- parms[["c"]]
    d <- parms[["d"]]
    eta <- parms[["eta"]]
    H <- parms[["H"]]
    v <- parms[["v"]]
    e <- parms[["e"]]
    w <- parms[["w_mad"]]
    a <- parms[["a"]]
    b <- parms[["b"]]
    A <- parms[["A"]]
    tau <- parms[["tau"]]
    
    coeff_0 <- (c + d)*tau*w^2*eta^2*H^2 - a*(1 - e*w)*A*v*b*H
    coeff_1 <- (c + d)*(tau*w^2*eta^2*(2*H*v*e - 2*H) + a*(1 - e*w)*v*w*eta*H) + a*(1 - e*w)*A*v*b
    coeff_2 <- (c + d)*(tau*w^2*eta^2*(1 + v^2*e^2 - 2*v*e) + a*(1 - e*w)*v*w*eta*(v*e - 1))
  }
  else {
    stop("model must be either 'Donnelly', 'Madden' or 'Madden_Cunniffe")
  }
  #if (model == "Madden" & v > 4 & e < 1) {browser()}
  # solve for I using coefficients to get all equilibrimus
  equilibriums <- polyroot(c(coeff_0, coeff_1, coeff_2))
  
  # subset to those that have no imaginary part, then remove imaginary part (as it is 0)
  equilibriums <- equilibriums[Im(zapsmall(equilibriums)) == 0]
  equilibriums <- Re(equilibriums)
  
  # subset to biologically possible equilibrium i.e. I > 0, I <= H (number host plants)
  equilibriums <- equilibriums[equilibriums > 0 & equilibriums <= parms[["H"]]]
  
  if (length(equilibriums) == 0) { # there is no stable positive equilibrium for I
    return(0)
  } 
  else if (length(equilibriums) > 1) {
      if (mult_eq_return_NA == T) { # if input has mult_eq_return_NA = T, return NA
        return(NA)
      }
  }
  
  # turn equilibrium into proportion of infected plants, i
  equilibriums <- equilibriums / parms[["H"]]
  
  return(equilibriums)

}

equilibrium_apply_func <- function(e_v_vals, parms, col_names, model) {
  
  ## function to be passed to an apply function - for given e and v vals, it changes
  ## e and v to those values and calculates the equilibrium value of I, number of 
  ## plants infected, for the Donnelly model
  
  parms[[ col_names[1] ]] <- e_v_vals[1]
  parms[[ col_names[2] ]] <- e_v_vals[2]
  
  equilibrium <- find_I_equilibrium(parms, model)
  
  return(equilibrium)
}

create_heatmap_data <- function(e_vals, v_vals, models, parms) {
  
  ### function that returns dataframe containing equilibrium values of I for a range of e 
  ### and v values, specified by user, for models specified by user. models parameter is a
  ### vector containing any of "Donnelly", "Madden" and "Madden_Cunniffe". output is ideal 
  ### for creating heatmap plot.
  
  # create df with all combinations of v and e
  v_e_df <- expand.grid(v = v_vals, e = e_vals) 
  
  # create new I_eq column for each model by calculating the I equilibrium for each row of v_e_df (i.e. each combination
  # of v and e values) for first the Donnelly model, then the Madden model
  I_eq_values <- lapply(models, 
                        function(model) apply(v_e_df, 1, equilibrium_apply_func, 
                                              parms = parms, 
                                              col_names = names(v_e_df), 
                                              model = model))
  col_names <- paste0(models, "_I_eq")
  
  for (i in 1:length(col_names)) {
    v_e_df[, col_names[i]] <- I_eq_values[[i]]
  }
  
  return(v_e_df)
}

# turn any bistable equilibriums into a high number so they're a different colour on the heatmap
# v_e_df$I_eq <- apply(v_e_df, 1, function(row) if(length(row$I_eq) > 1) {row$I_eq = 1.2}
#                                               else {row$I_eq = row$I_eq})

# define e and v vals across which the equilibrium of I will be found
v_vals <- seq(0.001, 5, length.out = 100)
e_vals <- seq(0.001, 1/parms[["w_mad"]], length.out = 100) # in madden e*w must be < 1

# create heatmap data
v_e_df <- create_heatmap_data(e_vals, v_vals, models = c("Donnelly", "Madden"), parms)

heatmap_don <- ggplot(data = v_e_df, aes(x = v, y = e, fill = Donnelly_I_eq)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                      name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Donnelly")


heatmap_mad <- ggplot(data = v_e_df, aes(x = v, y = e, fill = Madden_I_eq)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                       name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Madden")

pdf("heatmaps_v_e_final_I.pdf")
grid.arrange(heatmap_don, heatmap_mad)
dev.off()

##### INVESTIGATE BISTABILITY
parms_new <- parms
parms_new[["v"]] <- 0.40495960
parms_new[["e"]] <- 0.001

init_states <- c(I = 0.15*parms[["H"]])

trajectory_don <- data.frame(ode(y = init_states,
                                 times = times,
                                 parms = parms_new,
                                 func = donnelly_vpref_ode))

# plot trajectory of infected plants and number of aphids
plant_trajec <- ggplot(data = trajectory_don, aes(x = time, y = I)) +
  geom_line() +
  labs(x = "Time (days)",
       y = "Number of infected plants, I")
plant_trajec

## plot polynomial for I against I for both models
gamma <- parms[["gamma"]]
theta <- parms[["theta"]]
Pacq <- parms[["Pacq"]]
Pinoc <- parms[["Pinoc"]]
c <- parms[["c"]]
d <- parms[["d"]]
phi <- parms[["phi"]]
H <- parms[["H"]]
v <- parms[["v"]]
e <- parms[["e"]]
w_don <- parms[["w_don"]]
w_mad <- parms[["w_mad"]]
a <- parms[["a"]]
b <- parms[["b"]]
A <- parms[["A"]]
tau <- parms[["tau"]]


coeff_0_don <- theta*A*v*(1 - e*w_don)*H*Pacq*Pinoc - gamma*w_don*H^2
coeff_1_don <- gamma*w_don*v*(1 - e)*H + 2*gamma*w_don*H - 2*gamma*w_don*v*H - theta*A*v*(1 - e*w_don)*Pacq*Pinoc
coeff_2_don <- gamma*w_don*(v^2*(1 - e) - v*(1 - e) - v^2 + 2*v - 1)

coeff_0_mad <- (c + d)*tau*H^2 - phi^2*b*H*a*(1 - e*w_mad)*A*v
coeff_1_mad <- (c + d)*(2*tau*H*v - 2*tau*H + phi*a*(1 - e*w_mad)*v*H) + phi^2*b*a*(1 - e*w_mad)*A*v
coeff_2_mad <- (c + d)*(tau + tau*v^2 - 2*tau*v - phi*a*(1 - e*w_mad)*v + phi*a*(1 - e*w_mad)*v^2)

I <- seq(150, 600, length.out = 200)
polynom_don <- coeff_0_don + coeff_1_don*I + coeff_2_don*I^2
polynom_mad <- coeff_0_mad + coeff_1_mad*I + coeff_2_mad*I^2
plot(x = I, y = polynom_don, type = "l")
abline(h = 0)
plot(x = I, y = polynom_mad, type = "l")
abline(h = 0)


## find stability of fixed points
v <- 0.3
e <- 0.2

coeff_0_don <- theta*A*v*(1 - e*w_don)*H*Pacq*Pinoc - gamma*w_don*H^2
coeff_1_don <- gamma*w_don*v*(1 - e)*H + 2*gamma*w_don*H - 2*gamma*w_don*v*H - theta*A*v*(1 - e*w_don)*Pacq*Pinoc
coeff_2_don <- gamma*w_don*(v^2*(1 - e) - v*(1 - e) - v^2 + 2*v - 1)

I <- Re(polyroot(c(coeff_0_don, coeff_1_don, coeff_2_don)))

alpha <- Pacq*Pinoc*theta*A*v*(1 - e*w_don)
dI2 <- (alpha*(H - 2*I)*(w_don*(H + (v-1)*I)^2 - w_don*v*(1 - e)*(H*I + (v - 1)*I^2)) -
          alpha*(H*I - I^2)*(2*w_don*(H + (v - 1)*I)*(v - 1) - w_don*v*(1 - e)*(I + 2*(v - 1)*I))) /
  ((w_don*(H + (v - 1)*I)^2 - w_don*v*(1 - e)*(H*I + (v - 1)*I^2))^2) - gamma
dI2 # one stable, one unstable equilibrium


#### INVESTIGATE RELATIONSHIP BETWEEN W_DON AND W_MAD
# feeding probability for the 2 respective models

multiple_sensitivity_analysis <- function(parms_mad, parms_don, 
                                          varied_parm_name_mad, varied_parm_val_mad, 
                                          varied_parm_name_don, varied_parm_val_don, 
                                          init_states_don, init_states_mad, times, parms) {
  
  ## function that runs a sensitivity analysis multiple times for varying values of one madden and one donnelly
  ## parameter. returns dataframe of output and final graph displaying result in the form of a grob. use grid.arrange()
  ## to display the plot.
  
  if (length(varied_parm_val_don) != length(varied_parm_val_mad)) {
    stop("varied_parm_val_mad and varied_parm_val_don must be the same length")
  }
  
  # initialise output dataframe
  all_sens_analysis_df <- data.frame(
    analysis_parm_name = character(),
    analysis_parm_val = numeric(),
    final_I = numeric(),
    model = character(),
    varied_parm_name = character(),
    varied_parm_val = numeric()
  )
  
  # calculate how many rows of the output df are for each model
  num_mad_rows <- sum(lengths(parms_mad))
  num_don_rows <- sum(lengths(parms_don))
  
  # for loop to run sensitivity analysis for all values of varied_parms
  for (run in 1:length(varied_parm_val_don)) {
    
    # specify varying parameter values for this run
    parms[varied_parm_name_mad] <- varied_parm_val_mad[run]
    parms[varied_parm_name_don] <- varied_parm_val_don[run]
    
    # run sensitivity analysis with these parameter values
    sens_analysis <- sensitivity_analysis(parms_mad,
                                          parms_don,
                                          init_states_don, 
                                          init_states_mad, 
                                          times, 
                                          parms)
    
    # create df for information about the varying parameters in this run
    varying_parm_data <- data.frame(c(rep(varied_parm_name_mad, num_mad_rows), 
                                      rep(varied_parm_name_don, num_don_rows)), # varied_parm_name column
                                    c(rep(varied_parm_val_mad[run], num_mad_rows),
                                      rep(varied_parm_val_don[run], num_don_rows)) # varied_parm_val column
    )
    
    # add varying parameters data to output of sensitivity analysis
    full_sens_analysis_df <- cbind(sens_analysis, varying_parm_data)
    names(full_sens_analysis_df) <- c("analysis_parm_name", "analysis_parm_val", "final_I", "model", "varied_parm_name", "varied_parm_val")
    
    # add rows of this sensitivity analysis to overall dataframe
    all_sens_analysis_df <- rbind(all_sens_analysis_df, full_sens_analysis_df)
  }
  
  # split data by model
  split_data <- split(all_sens_analysis_df, f = all_sens_analysis_df$model)
  
  # plot madden model data
  madden_plot <- ggplot(data = split_data$Madden,
                        aes(x = analysis_parm_val, y = final_I, col = as.factor(round(varied_parm_val,2)))) +
    facet_wrap(~model) +
    geom_line() +
    labs(color = unique(split_data$Madden$varied_parm_name),
         x = unique(split_data$Madden$analysis_parm_name), 
         y = "Final disease incidence")
  
  # plot donnelly
  donnelly_plot <- madden_plot %+% split_data$Donnelly + # %+% means same plot but using Donnelly data
    labs(color = unique(split_data$Donnelly$varied_parm_name),
         x = unique(split_data$Donnelly$analysis_parm_name))
  
  final_plot <- arrangeGrob(madden_plot, donnelly_plot)
  
  # return final dataframe and grob of final plot (use grid.arrange to plot it)
  return(list(all_sens_analysis_df, final_plot))
}

e_vals <- seq(0.001, 5, length.out = 11)

w_don_vals <- list(w_don = seq(0.0001, 1, length.out = num_parm_runs))

w_mad_vals <- list(w_mad = seq(0, 1, length.out = num_parm_runs))


mult_sens_analysis_w_epsilon <- multiple_sensitivity_analysis(parms_mad = w_mad_vals, 
                                                              parms_don = w_don_vals, 
                                                              varied_parm_name_mad = "e", 
                                                              varied_parm_val_mad = e_vals, 
                                                              varied_parm_name_don = "e",
                                                              varied_parm_val_don = e_vals, 
                                                              init_states_don, init_states_mad, times, parms)
# plot result
grid.arrange(mult_sens_analysis_w_epsilon[[2]])


########################################
## ADDING VARIABLE PHI INTO MADDEN MODEL
########################################

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

## add eta (average duration of a plant feed)
parms[["eta"]] <- 2

### RUN DONNELLY EPIDEMIC
trajectory_don <- data.frame(ode(y = init_states_don,
                                 times = times,
                                 parms = parms,
                                 func = donnelly_vpref_ode))
trajectory_don$I <- trajectory_don$I/parms[["H"]] # make proportion

# plot trajectory of infected plants and number of aphids
new_plant_trajec <- ggplot(data = trajectory_don, aes(x = time, y = I)) +
  geom_line() +
  labs(x = "Time (days)",
       y = "Number of infected plants, I")

### RUN MADDEN EPIDEMIC
trajectory_mad <- data.frame(ode(y = init_states_mad,
                                 times = times,
                                 func = madden_cunniffe_vpref_ode,
                                 parms = parms))
trajectory_mad$I <- trajectory_mad$I/parms[["H"]] # make proportion


trajectory_mad_long <- reshape2::melt(trajectory_mad, id.vars = "time")
names(trajectory_mad_long) <- c("time", "compartment", "number")

# plot trajectory of infected plants and vector compartments
new_plant_trajec <- new_plant_trajec +
  geom_line(data = trajectory_mad_long %>% filter(compartment == "I"),
            aes(x = time, y = number), color = "red") +
  annotate("text", label = "red = Madden\nblack = Donnelly", x = times[length(times)/3*2], y = 0)

grid.arrange(plant_trajec, new_plant_trajec)

#### RE-RUN HEATMAPS FOR NEW MADDEN MODEL

v_vals <- seq(0.001, 5, length.out = 100)
e_vals <- seq(0.001, 1/parms[["w_mad"]], length.out = 100) # in madden e*w must be < 1

# create heatmap data
v_e_df2 <- create_heatmap_data(e_vals, v_vals, models = c("Donnelly", "Madden", "Madden_Cunniffe"), parms)

heatmap_don_new <- ggplot(data = v_e_df2, aes(x = v, y = e, fill = Donnelly_I_eq)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                       name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Donnelly")


heatmap_mad_new <- ggplot(data = v_e_df2, aes(x = v, y = e, fill = Madden_I_eq)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                       name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Madden")

heatmap_mad_cun <- ggplot(data = v_e_df2, aes(x = v, y = e, fill = Madden_Cunniffe_I_eq)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                       name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Madden (variable phi)")

pdf("heatmaps_3_models.pdf", width = 7, height = 10)
grid.arrange(heatmap_don_new, heatmap_mad_new, heatmap_mad_cun)
dev.off()
