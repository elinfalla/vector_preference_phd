################################################################################################################
##### COMPARISON OF VECTOR PREFERENCE AND VECTOR DYNAMICS IN MADDEN (2000) AND DONNELLY (2019) TYPE MODELS #####
################################################################################################################

### MADDEN MODEL:
### Modified version of model by Madden et al. 2000 to include vector preference - both of plant attractiveness and
### plant acceptability. Vector population dynamics (emigration birth, death) are included to emulate the Donnelly model,
### plants are represented by an SI rather than an SEIR model, vectors represented by an SI rather than SEI model. 
### Vectors cannot be born infective (as this model represents NPT viruses).

### One version of the model has a fixed value of phi (number of plants visited pre day), the other has a variable value
### of phi to account for the fact that feeding takes longer than probing, determined as defined by Cunniffe et al. 2021.

### DONNELLY MODEL:
### Modification of deterministic model by Donnelly et al. 2019 to so it tracks I, number of infected plants, rather than 
### i, proportion of infected plants. Otherwise this model is equivalent to the one in Donnelly et al. 2019.

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gridtext)
library(flux)

#### DEFINE MADDEN (fixed phi) ODE
madden_vdynamic_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  become_infected <- par[["phi"]] * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  #di <- dI/parms[["H"]]
  
  # VECTOR EQUATIONS
  A <- y[["X"]] + y[["Z"]]
  
  acquisition <- par[["phi"]] * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * y[["X"]] * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  birth <- par[["lamda"]] * A * (1 - ((A/par[["H"]])/par[["K"]]))
  death_X <- par[["alpha"]] * y[["X"]]
  death_Z <- par[["alpha"]] * y[["Z"]]
  
  emigration_X <- par[["phi"]] * par[["p"]] * y[["X"]]
  emigration_Z <- par[["phi"]] * par[["p"]] * y[["Z"]]
  
  # state equations
  dX <- stop_being_infective - acquisition + birth - death_X #- emigration_X
  dZ <- acquisition - stop_being_infective - death_Z #- emigration_Z
  

  return(list(c(dI, dX, dZ)))
}

#### DEFINE MADDEN (variable phi) ODE
madden_cunniffe_vdynamic_ode <- function(times, y, par) {

  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  # define variable phi
  phi <- (S + y[["I"]]*par[["v"]]) / (par[["w_mad"]]*par[["eta"]]*(S + y[["I"]]*par[["v"]]*par[["e"]]))
  
  # PLANT EQUATIONS
  become_infected <- phi * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  #di <- dI/parms[["H"]]
  
  # VECTOR EQUATIONS
  A <- y[["X"]] + y[["Z"]]
  
  acquisition <- phi * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * y[["X"]] * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  birth <- par[["lamda"]] * A * (1 - ((A/par[["H"]])/par[["K"]]))
  death_X <- par[["alpha"]] * y[["X"]]
  death_Z <- par[["alpha"]] * y[["Z"]]
  
  emigration_X <- phi * par[["p"]] * y[["X"]]
  emigration_Z <- phi * par[["p"]] * y[["Z"]]
  
  # state equations
  dX <- stop_being_infective - acquisition + birth - death_X - emigration_X
  dZ <- acquisition - stop_being_infective - death_Z - emigration_Z
  
  
  return(list(c(dI, dX, dZ)))
}

### DEFINE DONNELLY ODE (full model - 3 state equations)
donnelly_full_vdynamic_ode <- function(times, states, parms) {
  
  ### ODE function for deterministic Donnelly model
  
  # STATES
  i <- states[["i"]]
  As <- states[["As"]]
  Ai <- states[["Ai"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  alpha <- parms[["alpha"]]
  theta <- parms[["theta"]]
  lamda <- parms[["lamda"]]
  K <- parms[["K"]]
  H <- parms[["H"]]
  p <- parms[["p"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w_don"]]
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
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  # A = total number of aphids i.e. aphids per plant type * number of that plant type
  # i.e. S*As + I*Ai where S is number of susceptible plants and I is number of infected plants
  A <- As*(H - i*H) + Ai*(i*H)  
  
  
  # STATE EQUATIONS
  # aphids - change in per plant aphid density on susceptible (As)  and infected (Ai) plants
  # = aphid births - aphid deaths - aphid settling on other plant type + aphids moving from other plant type
  dAs <- lamda*As*(1 - As/K) - alpha*As - theta*As*(1 - Fs) + theta*Ai*Fs*i/(1 - i)
  dAi <- lamda*Ai*(1 - Ai/K) - alpha*Ai - theta*Ai*(1 - Fi) + theta*As*Fi*(1 - i)/i
  
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(c(di, dAs, dAi)))
}

### DEFINE DONNELLY ODE (simplified - 1 state equation) 
# define function to calculate As and Ai equilibriums for a given value of i.
calc_A_equilibrium <- function(parms, states, i_hat, return_all=F) {
  
  # define states
  i <- states[["i"]]
  
  # define parameters
  e <- parms[["e"]]
  lamda <- parms[["lamda"]]
  alpha <- parms[["alpha"]]
  K <- parms[["K"]]
  theta <- parms[["theta"]]
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  ## define coefficients 
  # (worked out by hand from dAs and dAi equations- see lab book). will give roots of Ai
  coeff_1 <- i*(lamda - alpha - theta*(1 - Fs))*(alpha + theta*(1 - Fi) - lamda)/(theta*Fi*(1 - i)) + theta*Fs*i/(1 - i)
  
  coeff_2 <- (lamda/K)*i*(lamda - alpha - theta*(1 - Fs))/(theta*Fi*(1 - i)) - 
    lamda*i^2*(alpha^2 + 2*alpha*theta*(1 - Fi) - 2*alpha*lamda + theta^2*(1 - Fi)^2 - 2*lamda*theta*(1 - Fi) + lamda^2)/(K*theta^2*(1 - i)^2*Fi^2)
  
  coeff_3 <- -(lamda*i^2/K)*(2*alpha*lamda + 2*lamda*theta*(1 - Fi) - 2*lamda^2)/(K*theta^2*(1 - i)^2*Fi^2)
  
  coeff_4 <- -lamda^3*i^2/(K^3*theta^2*(1 - i)^2*Fi^2)
  
  # get Ai equilibrium values
  Ai_equilibriums <- polyroot(c(0, coeff_1, coeff_2, coeff_3, coeff_4))
  
  # subset to those that have no imaginary part, then remove imaginary part (as it is 0)
  Ai_equilibriums <- Ai_equilibriums[Im(zapsmall(Ai_equilibriums)) == 0]
  Ai_equilibriums <- Re(Ai_equilibriums)
  
  # for each Ai equilibrium, calculate value of As (equation calculated manually - see lab book)
  As_equilibriums <- i*(alpha*Ai_equilibriums + theta*Ai_equilibriums*(1 - Fi) - lamda*Ai_equilibriums*(1-Ai_equilibriums/K)) / 
    (theta*Fi*(1-i))
  
  # find equilibrium where both As and Ai are positive i.e. biologically possible equilibrium
  Ai <- Ai_equilibriums[which(Ai_equilibriums > 0 & As_equilibriums > 0)]
  As <- As_equilibriums[which(Ai_equilibriums > 0 & As_equilibriums > 0)]
  
  # if multiple positive equilibriums and return_all = F, arbitrarily select the first one. Give warning
  # if return_all = T, then return all positive equilibriums
  if (length(Ai) > 1 | length(As) > 1) {
    if (return_all == F) {
      Ai <- Ai[1]
      As <- As[1]
    }
    warning("Multiple positive As/Ai equilibriums")
    
  }
  return(c(As, Ai))
}

donnelly_simp_vdynamic_ode <- function(times, states, parms) {
  
  ### ODE function for deterministic Donnelly model, with As and Ai (aphids per healthy and infected plants 
  ### respectively) assumed to be always at equilibrium. Therefore calculates equilibrium As and Ai values 
  ### for each value of i (proportion infected plants).
  
  # STATES
  i <- states[["i"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  alpha <- parms[["alpha"]]
  theta <- parms[["theta"]]
  lamda <- parms[["lamda"]]
  K <- parms[["K"]]
  H <- parms[["H"]]
  p <- parms[["p"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w_don"]]
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
  As_Ai <- calc_A_equilibrium(parms, states, i_hat)
  As <- As_Ai[1]
  Ai <- As_Ai[2]
  
  # A = total number of aphids i.e. aphids per plant type * number of that plant type
  # i.e. S*As + I*Ai where S is number of susceptible plants and I is number of infected plants
  A <- As*(H - i*H) + Ai*(i*H)  
  
  # STATE EQUATION
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(c(di)))
}


# parameters
parms <- c(
  
  ## MADDEN PARMS
  tau = 4, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 5.473098, # plants visited per day by an insect
  a = 0.8, # virus acquisition probability by vectors
  b = 0.8, # plant innoculation probability by vectors
  c = 2, # natural plant death rate
  d = 2, # plant death rate due to infection

  # DONNELLY PARMS
  gamma = 4, # rate of recovery/death of I plants
  theta = 1.041022, # aphid dispersal rate per day
  Pacq = 0.8, # chance of virus acquisition by vector from infected plant
  Pinoc = 0.8, # chance of inoculating healthy plant
  w_don = 0.2, # feeding probability on healthy plant
  
  # PARMS FOR BOTH
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  #q = 0.8, # probability of not emigrating per flight
  lamda = 2.2, # birth rate of vectors
  K = 10, # carrying capacity of vectors per plant
  alpha = 0.2, # vector death rate
  H = 400, # number of host plants
  p = 0.2 # chance of vector emigration per flight
  
)

parms[["eta"]] <- 1/(parms[["phi"]]*parms[["w_mad"]]) # average length of vector feed
parms[["a"]] <- parms[["Pacq"]]  # virus acquisition probability by vectors
parms[["b"]] <- parms[["Pinoc"]] # plant innoculation probability by vectors
parms[["c"]] <- parms[["gamma"]] / 2 # natural plant death rate
parms[["d"]] <- parms[["gamma"]] / 2 # plant death rate due to infection

num_vectors <- 1200

init_states_mad <- c(
  I = 1,
  X = num_vectors,
  Z = 0
)

init_states_don_simp <- c(
  i = 1/parms[["H"]]
)

init_states_don_full <- c(
  i = 1/parms[["H"]],
  As = num_vectors / parms[["H"]] / 2,
  Ai = num_vectors / parms[["H"]] / 2
)

times <- seq(0, 30, by = 0.2)


###### run epidemic - madden fixed phi ######
trajec_mad <- data.frame(ode(y = init_states_mad,
                             times = times,
                             parms = parms,
                             func = madden_vdynamic_ode))
trajec_mad$I <- trajec_mad$I/(parms[["H"]]) # make into a proportion

trajec_mad_long <- reshape2::melt(trajec_mad, id.vars = "time")
names(trajec_mad_long) <- c("time", "compartment", "number")

mad_plant_trajec <- ggplot(data = trajec_mad_long %>% filter(compartment == "I"),
                           aes(x = time, y = number)) +
  geom_line() +
  ggtitle("Madden model (fixed phi)") +
  labs(y = "Proportion of infected plants", x= "Time (days)")
  

mad_vec_trajec <- ggplot(data = trajec_mad_long %>% filter(!compartment == "I"),
                         aes(x = time, y = number, col = compartment)) +
  geom_line() +
  labs(col = "Vector state", y = "Number of vectors", x = "Time (days)")

grid.arrange(mad_plant_trajec, mad_vec_trajec)


###### run epidemic - madden variable phi ######
trajec_mad_cun <- data.frame(ode(y = init_states_mad,
                             times = times,
                             parms = parms,
                             func = madden_cunniffe_vdynamic_ode))
trajec_mad_cun$I <- trajec_mad_cun$I/(parms[["H"]]) # make into a proportion

trajec_mad_cun_long <- reshape2::melt(trajec_mad_cun, id.vars = "time")
names(trajec_mad_cun_long) <- c("time", "compartment", "number")

mad_cun_plant_trajec <- ggplot(data = trajec_mad_cun_long %>% filter(compartment == "I"),
                           aes(x = time, y = number)) +
  geom_line() +
  ggtitle("Madden model (variable phi)") +
  labs(y = "Proportion of infected plants", x = "Time (days)")


mad_cun_vec_trajec <- ggplot(data = trajec_mad_cun_long %>% filter(!compartment == "I"),
                         aes(x = time, y = number, col = compartment)) +
  geom_line() +
  labs(col = "Vector state", y = "Number of vectors", x = "Time (days)")

grid.arrange(mad_cun_plant_trajec, mad_cun_vec_trajec)

###### run epidemic - donnelly simple ######
trajec_don_simp <- data.frame(ode(y = init_states_don_simp,
                                 times = times,
                                 parms = parms,
                                 func = donnelly_simp_vdynamic_ode))

don_simp_plant_trajec <- ggplot(data = trajec_don_simp,
                               aes(x = time, y = i)) +
  geom_line() +
  ggtitle("Simplified Donnelly model") +
  labs(y = "Proportion of infected plants", x = "Time (days)")
don_simp_plant_trajec


###### run epidemic - donnelly full ######
trajec_don_full <- data.frame(ode(y = init_states_don_full,
                                  times = times,
                                  parms = parms,
                                  func = donnelly_full_vdynamic_ode))

don_full_plant_trajec <- ggplot(data = trajec_don_full,
                                aes(x = time, y = i)) +
  geom_line() +
  ggtitle("Full Donnelly model") +
  labs(y = "Proportion of infected plants", x = "Time (days)")
don_full_plant_trajec



###### ANALYSE WHETHER SIMPLIFIED DONNELLY MODEL IS A GOOD APPROX FOR FULL MODEL ######

calc_AUDPC_diff <- function(parm_vals, x_parm_names, y_parm_names, 
                            func1, func2, parms, init_states_func1, init_states_func2, times) {
  
  ### apply function for get_AUDPC_diff_data(). For 2 given parameter values, it calculates the
  ### AUDPC with those values for 2 models, and returns the absolute difference in AUDPC between 
  ### the models.
  
  for (parm in x_parm_names) {
    parms[[parm]] <- parm_vals[1]
  }
  for (parm in y_parm_names) {
    parms[[parm]] <- parm_vals[2]
  }
  
  func1_trajec <- data.frame(ode(y = init_states_func1,
                                 times = times,
                                 func = func1,
                                 parms = parms))
  func2_trajec <- data.frame(ode(y = init_states_func2,
                                 times = times,
                                 func = func2,
                                 parms = parms))
  
  # return absolute difference between the areas under the curves of the 2 models 
  return(abs(flux::auc(func1_trajec[,1], func1_trajec[,2]) - flux::auc(func2_trajec[,1], func2_trajec[,2])))
}
  

get_AUDPC_diff_data <- function(x_parm_names, y_parm_names, x_parm_vals, y_parm_vals, 
                            func1, func2, parms, init_states_func1, init_states_func2, times) {
  
  ### function that for given parameters (and values) and 2 given ODE functions, calculates
  ### the absolute difference in the Area Under Disease Progress Curve (AUDPC) between runs of the 2 ODEs
  ### for each combination of parameter values. Returns a dataframe of the parameter values and 
  ### their respective AUDPC difference.
  
  #expand.grid for parm1 and parm2 vals
  all_parm_combos <- expand.grid(x_parm_vals, y_parm_vals)
  names(all_parm_combos) <- c(paste(x_parm_names, collapse = "_"), paste(y_parm_names, collapse = "_"))
  
  # for each row of resulting df calculate difference in AUDPC
  AUDPC_diffs <- apply(all_parm_combos, 1, calc_AUDPC_diff, x_parm_names, y_parm_names, 
                      func1, func2, parms, init_states_func1, init_states_func2, times)
  
  all_parm_combos$AUDPC_diff <- AUDPC_diffs
  
  return(all_parm_combos)
}

AUDPC_Pacq_Pinoc_theta <- get_AUDPC_diff_data(x_parm_name = "theta",
                                              y_parm_name = c("Pacq", "Pinoc"),
                                              x_parm_vals = seq(0.00001, 5, length.out = 30),
                                              y_parm_vals = seq(0, 1, length.out = 30),
                                              func1 = donnelly_full_vdynamic_ode,
                                              func2 = donnelly_simp_vdynamic_ode,
                                              parms = parms,
                                              init_states_func1 = init_states_don_full,
                                              init_states_func2 = init_states_don_simp,
                                              times = times)

Pacq_Pinoc_theta_don_comparison <- ggplot(data = AUDPC_Pacq_Pinoc_theta,
                                          aes(x = theta, y = Pacq_Pinoc, fill = AUDPC_diff)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "blue",
                       name = "Absolute\ndifference\nin AUDPC", na.value = "white") +
  labs(y = "Pacq/Pinoc")
Pacq_Pinoc_theta_don_comparison

AUDPC_v_e <- get_AUDPC_diff_data(x_parm_name = "v",
                                 y_parm_name = "e",
                                 x_parm_vals = seq(0.001, 1/parms[["w_don"]], length.out = 100),
                                 y_parm_vals = seq(0.001, 1/parms[["w_mad"]], length.out = 100),
                                 func1 = donnelly_full_vdynamic_ode,
                                 func2 = donnelly_simp_vdynamic_ode,
                                 parms = parms,
                                 init_states_func1 = init_states_don_full,
                                 init_states_func2 = init_states_don_simp,
                                 times = times)

v_e_don_comparison <- ggplot(data = AUDPC_v_e,
                             aes(x = v, y = e, fill = AUDPC_diff)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("darkred", "red", "blue", "navyblue"),
                       values = scales::rescale(c(0, min(AUDPC_v_e$AUDPC_diff), 1, max(AUDPC_v_e$AUDPC_diff))),
                       limits = c(0, max(AUDPC_v_e$AUDPC_diff)),
                       name = "Absolute difference\nin AUDPC",
                       na.value = "white",
                       guide = guide_colorbar(barheight = unit(8, "cm"), barwidth = unit(1.5, "cm")))
v_e_don_comparison

#############

### DO SENSITIVITY ANALYSIS 

vary_param <- function(init_states, times, parms, varied_parm_name, varied_parm_vals, func, model) {
  
  # function to run the Madden ODE multiple times with one parameter varying, as specified by user.
  # returns data frame of final incidence, I, of all runs, along with values of varied parameter
  
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
    
    if ("I" %in% names(trajectory)) {
      names(trajectory)[which(names(trajectory) == "I")] <- "i"
    }
    
    output[run, "final_I"] <- round(trajectory[nrow(trajectory), "i"], 3)
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
    
    output_mad <- vary_param(init_states_mad, times, default_parms,
                             varied_parm_name = names(parms_mad[i]),
                             varied_parm_vals = parms_mad[[i]],
                             func = madden_vdynamic_ode,
                             model = "Madden") 
    
    output_mad_cun <- vary_param(init_states_mad, times, default_parms,
                             varied_parm_name = names(parms_mad[i]),
                             varied_parm_vals = parms_mad[[i]],
                             func = madden_cunniffe_vdynamic_ode,
                             model = "Madden_Cunniffe") 
    
    sen_analysis_df <- rbind(sen_analysis_df, output_mad, output_mad_cun)
    
  }
  
  print("Starting Donnelly sensitivity analysis")
  
  for (i in 1:length(parms_don)) {
    
    output_don <- vary_param(init_states_don, times, default_parms,
                             varied_parm_name = names(parms_don[i]),
                             varied_parm_vals = parms_don[[i]],
                             func = donnelly_full_vdynamic_ode,
                             model = "Donnelly") 
    
    
    sen_analysis_df <- rbind(sen_analysis_df, output_don)
    
  }
  
  return(sen_analysis_df) 
}

num_parm_runs <- 100
analysis_parms_don <- list(theta = seq(0.0001, 20, length.out = num_parm_runs),
                           w_don = seq(0.01, 1, length.out = num_parm_runs),
                           Pacq = seq(0, 1, length.out = num_parm_runs),
                           Pinoc = seq(0, 1, length.out = num_parm_runs),
                           gamma = seq(0, 20, length.out = num_parm_runs),
                           v = seq(0, 20, length.out = num_parm_runs),
                           e = seq(0, 20, length.out = num_parm_runs),
                           lamda = seq(0.0001, 20, length.out = num_parm_runs),
                           K = seq(0.0001, 20, length.out = num_parm_runs),
                           alpha = seq(0.0001, 20, length.out = num_parm_runs),
                           p = seq(0.0001, 1, length.out = num_parm_runs))

analysis_parms_mad <- list(phi = seq(0, 20, length.out = num_parm_runs),
                           eta = seq(0.0001, 20, length.out = num_parm_runs),
                           tau = seq(0, 20, length.out = num_parm_runs),
                           w_mad = seq(0.01, 1, length.out = num_parm_runs),
                           a = seq(0, 1, length.out = num_parm_runs),
                           b = seq(0, 1, length.out = num_parm_runs),
                           c = seq(0, 20, length.out = num_parm_runs),
                           d = seq(0, 20, length.out = num_parm_runs),
                           v = seq(0, 20, length.out = num_parm_runs),
                           e = seq(0, 20, length.out = num_parm_runs),
                           lamda = seq(0.0001, 20, length.out = num_parm_runs),
                           K = seq(0.0001, 20, length.out = num_parm_runs),
                           alpha = seq(0.0001, 20, length.out = num_parm_runs),
                           p = seq(0.0001, 1, length.out = num_parm_runs))

# sensitivity analysis for phi and theta
sens_analysis_res <- sensitivity_analysis(parms_mad = analysis_parms_mad,
                                          parms_don = analysis_parms_don,
                                          init_states_don_full, 
                                          init_states_mad, 
                                          times, 
                                          parms)

# plot (to pdf)
sens_anal_plots <- lapply(unique(sens_analysis_res$model), function(m){
  apply(unique(sens_analysis_res %>% filter(model == m) %>% select(parm_name)), 1, function(p) {
    ggplot(data = sens_analysis_res %>% filter(model == m & parm_name == p),
           aes(x = parm_val,
               y = final_I)) +
    geom_line() +
    labs(x = p) }
    ) }
)


pdf("sens_analysis_vec_dynamics_models.pdf")
do.call("grid.arrange", sens_anal_plots[[1]])
do.call("grid.arrange", sens_anal_plots[[2]])
do.call("grid.arrange", sens_anal_plots[[3]])
dev.off()

### FIND OPTIMUM PHI AND THETA VALS
find_theta_phi_vals <- function(theta_phi_data, init_states_don, init_states_mad, times, parms) {
  
  # function to find point at which epidemic takes off (final_I > 0) for phi and theta
  
  phi_data <- theta_phi_data %>% filter(parm_name == "phi")
  theta_data <- theta_phi_data %>% filter(parm_name == "theta")
  
  phi_index <- match(TRUE, phi_data$final_I > 0.5) # finds first instance of final_I > 0.5
  theta_index <- match(TRUE, theta_data$final_I > 0.5)
  
  if (is.na(phi_index)) {
    if (is.na(theta_index)) {
      stop("PHI AND THETA: No threshold found - final incidence always >0.5")
    } else {
      stop("PHI: No threshold found - final incidence always >0.5")
    }
  } 
  else if (is.na(theta_index)) {
    stop("THETA: No threshold found - final incidence always >0.5")
  } 
  else if (phi_index == 1 | theta_index == 1) {
    stop("PHI/THETA: No threshold found - final incidence never >0.5")
  }
  # extract rows of data around the threshold
  phi_subset <- phi_data[(phi_index-1):phi_index,]
  theta_subset <- theta_data[(theta_index-1):theta_index,]
  
  # if the first instance where final disease incidence > 0 is smaller than 0.01, consider
  # the estimate for the threshold accurate enough and return it
  if (phi_subset[2, "final_I"] < 0.502 | theta_subset[2, "final_I"] < 0.502) {
    
    phi_thresh <- mean(c(phi_subset[1, "parm_val"], phi_subset[2, "parm_val"]))
    theta_thresh <- mean(c(theta_subset[1, "parm_val"], theta_subset[2, "parm_val"]))
    
    return(c(phi_threshold = phi_thresh,
             theta_threshold = theta_thresh))
  }
  
  # else run find_epidemic_threshold() with new smaller parameter limits recursively until 
  # above condition is reached
  
  phi_lower_lim <- phi_subset[1, "parm_val"]
  phi_upper_lim <- phi_subset[2, "parm_val"]
  
  theta_lower_lim <- theta_subset[1, "parm_val"]
  theta_upper_lim <- theta_subset[2, "parm_val"]
  
  num_parm_runs <- 100
  theta_vals_condensed <- list(theta = seq(theta_lower_lim, theta_upper_lim, length.out = num_parm_runs))
  
  phi_vals_condensed <- list(phi = seq(phi_lower_lim, phi_upper_lim, length.out = num_parm_runs))
  
  # re-run sensitivity analysis
  new_theta_phi_vals <- sensitivity_analysis(parms_mad = phi_vals_condensed,
                                             parms_don = theta_vals_condensed,
                                             init_states_don, 
                                             init_states_mad, 
                                             times, 
                                             parms)
  new_theta_phi_vals <- new_theta_phi_vals %>% 
    mutate(final_I = case_when(model == "Madden" ~ final_I/parms[["H"]],
                               model == "Donnelly" ~ final_I))
  
  # feed back into function recursively
  find_theta_phi_vals(new_theta_phi_vals, 
                      init_states_don, 
                      init_states_mad, 
                      times, 
                      parms)
}

theta_phi_data_3model <- sens_analysis_res %>%
  filter(parm_name == "theta" | parm_name == "phi") %>%
  mutate(final_I = case_when(model == "Madden" | model == "Madden_Cunniffe" ~ final_I/parms[["H"]],
                             model == "Donnelly" ~ final_I))

theta_phi_data_2model <- theta_phi_data_3model %>%
  filter(!model == "Madden_Cunniffe")

out <- find_theta_phi_vals(theta_phi_data_2model,
                           init_states_don_full, 
                           init_states_mad, 
                           times, 
                           parms)

out[["phi_threshold"]] # 4.5515764 
out[["theta_threshold"]] # 0.9683637

parms[["phi"]] <- out[["phi_threshold"]]
parms[["theta"]] <- out[["theta_threshold"]]


# CALCULATE EQUILIBRIUM DISEASE INCIDENCE FOR DONNELLY

find_i_equilibrium <- function(v_e_vals, init_states, parms, func) {
  
  parms[["v"]] <- v_e_vals[1]
  parms[["e"]] <- v_e_vals[2]
  
  if ("i" %in% names(init_states)) {
    low_i <- 0.0001
    high_i <- 0.9999
  } else {
    low_i <- 0.0001*parms[["H"]]
    high_i <- 0.9999*parms[["H"]]
  }
  
  init_states[1] <- low_i
  times <- seq(0, 10000, length.out = 2)
  equilibrium_val_low_start <- round(data.frame(ode(y = init_states,
                                     times = times,
                                     func = func,
                                     parms = parms))[length(times), 2], 5)
 
  init_states[1] <- high_i
  equilibrium_val_high_start <-  round(data.frame(ode(y = init_states,
                                                times = times,
                                                func = func,
                                                parms = parms))[length(times), 2], 5)
  
  if ("I" %in% names(init_states)) {
    equilibrium_val_low_start <- equilibrium_val_low_start/parms[["H"]]
    equilibrium_val_high_start <- equilibrium_val_high_start/parms[["H"]]
  }
  
  if (equilibrium_val_low_start != equilibrium_val_high_start) {
    return(NA)
  }
  return(equilibrium_val_low_start)

}
v_vals <- seq(0.001, 5, length.out = 100)
e_vals <- seq(0.001, 1/parms[["w_mad"]], length.out = 100) # in madden e*w must be < 1
v_e_vals <- expand.grid(v = v_vals, e = e_vals)

v_e_vals$I_eq_don <- apply(X = v_e_vals, MARGIN = 1, FUN = find_i_equilibrium, init_states_don_full, parms, donnelly_full_vdynamic_ode)
v_e_vals$I_eq_mad <- apply(X = v_e_vals, MARGIN = 1, FUN = find_i_equilibrium, init_states_mad, parms, madden_vdynamic_ode)
v_e_vals$I_eq_mad_cun <- apply(X = v_e_vals, MARGIN = 1, FUN = find_i_equilibrium, init_states_mad, parms, madden_cunniffe_vdynamic_ode)

# v_e_vals <- v_e_vals %>%
#   mutate(I_eq_mad = I_eq_mad / parms[["H"]],
#          I_eq_mad_cun = I_eq_mad_cun / parms[["H"]])

heatmap_don <- ggplot(data = v_e_vals, aes(x = v, y = e, fill = I_eq_don)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                       name = "Equilibrium\ni value", 
                       na.value = "white",) +
  labs(title = "Donnelly")

heatmap_mad <- ggplot(data = v_e_vals, aes(x = v, y = e, fill = I_eq_mad)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                       name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Madden")

heatmap_mad_cun <- ggplot(data = v_e_vals, aes(x = v, y = e, fill = I_eq_mad_cun)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("darkred", "red", "blue"), values = c(0, 0.00001, 1),
                       name = "Equilibrium\ni value", na.value = "white") +
  labs(title = "Madden_Cunniffe")

grid.arrange(heatmap_mad_cun, heatmap_mad, heatmap_don)

