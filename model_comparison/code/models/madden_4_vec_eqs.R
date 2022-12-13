##########################################################################################
##### CREATION OF MADDEN-CUNNIFFE MODEL THAT TRACKS THE TYPE OF PLANT VECTORS ARE ON #####
##########################################################################################

### Version of Madden-Cunniffe model with 4 vector state equations, allowing tracking of
### vectors on infected versus healthy plants. In addition, infective vectors can only
### infect one plant before losing infectivity, represented by making tau (infectivity
### loss rate) = phi (varible flight rate). Flight rate is variable according to Cunniffe
### 2021.

rm(list = ls())

library(deSolve)

# import default parms, times, initial states
source("./code/comparison/DEFAULT_PARMS_TIMES_INIT_STATES.R")

# change initial states to incorporate new vector states
init_states_mad <- init_states_mad[-which(names(init_states_mad) == "Z")] # remove Z state

parms[["H"]] <- 9
parms[["A"]] <- 27
aphid_per_plant <- parms[["A"]]/parms[["H"]]
S <- parms[["H"]] - I # I imported from default parms doc

init_states_mad <- append(init_states_mad, c(Xs = aphid_per_plant*S,
                                             Xi = aphid_per_plant*I,
                                             Zs = 0,
                                             Zi = 0))

madden_cunniffe_4_vec_ode <- function(times, y, par) {
  
  # define states
  I <- y[["I"]]
  S <- par[["H"]] - I
  
  Xs <- y[["Xs"]]
  Xi <- y[["Xi"]]
  Zs <- y[["Zs"]]
  Zi <- y[["Zi"]]
  
  # define parameters
  w <- par[["w_mad"]]
  A <- par[["A"]]
  H <- par[["H"]]
  v <- par[["v"]]
  e <- par[["e"]]
  eta <- par[["eta"]]
  c <- par[["c"]]
  d <- par[["d"]]
  a <- par[["a"]]
  b <- par[["b"]]
  
  # define phi - vector flight rate, and tau, vector infectivity loss rate
  phi <- (S + v*I) / (w*eta*(S + v*e*I))
  tau <- phi
  
  ## PLANT STATE EQUATION
  dI <- phi*b*(Zi)*S / (S + v*I) - (c + d)*I
  
  ## VECTOR STATE EQUATIONS
  S_infectivity_loss <- tau*Zs
  I_infectivity_loss <- tau*Zi
  acquisition <- phi*a*(1 - e*w)*v*I / (S + v*I) # needs to be multiplied by relevant vector compartment
  flight_Xs_to_Xi <- phi*(1 - a + a*e*w)*Xs*v*I / (S + v*I)
  flight_Xi_to_Xs <- phi*Xi*S / (S + v*I)
  flight_Zi_to_Zs <- phi*(1 - b)*(1 - w)*Zi*S / (S + v*I)
  
  dXs <- S_infectivity_loss - acquisition*Xs - flight_Xs_to_Xi + flight_Xi_to_Xs
  dXi <- I_infectivity_loss - acquisition*Xi + flight_Xs_to_Xi - flight_Xi_to_Xs
  dZs <- flight_Zi_to_Zs - S_infectivity_loss - acquisition*Zs
  dZi <- acquisition*(Zs + Xi + Xs) - I_infectivity_loss - flight_Zi_to_Zs
  
  return(list(c(dI, dXs, dXi, dZs, dZi)))
}

parms["eta"] <- 1
parms["b"] <- 1
parms["a"] <- 1
parms["c"] <- 2
parms["d"] <- 2

run_ode <- data.frame(deSolve::ode(y = init_states_mad,#[-4], # no Zs
                                   times = times[1:50],
                                   parms = parms,
                                   func = madden_cunniffe_4_vec_ode))
#plot(run_ode$time, run_ode$I)
plot(run_ode$time, run_ode$Xs, type = "l", xlab = "Time", ylab = "Num vectors",
     ylim = c(0, parms[["A"]]), col = "green")
lines(run_ode$time, run_ode$Xi, col = "orange")
lines(run_ode$time, run_ode$Zs, col = "purple")
lines(run_ode$time, run_ode$Zi, col = "red")
lines(run_ode$time, run_ode$I, col = "blue")

legend("topright", legend = c("I", "Xi", "Xs", "Zi", "Zs"),
       col = c("blue", "orange", "green", "red", "purple"),
       lty = 1)

parms
init_states_mad

