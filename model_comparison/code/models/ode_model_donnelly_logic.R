##############################################################
######### ODE MODEL WITH DONNELLY TRANSMISSION LOGIC #########
##############################################################

## ODE model where you lose infectivity after probing/feeding on
## one plant, rather than having a loss of infectivity rate 
## parameter as in the Madden model.

rm(list = ls())

library(deSolve)

model_ode <- function(times, y, parms) {
  
  omega <- parms[["omega"]]
  eta <- parms[["eta"]]
  H <- parms[["H"]]
  A <- parms[["A"]]
  b <- parms[["b"]]
  a <- parms[["a"]]
  c <- parms[["c"]]
  d <- parms[["d"]]
  v <- parms[["v"]]
  epsilon <- parms[["epsilon"]]
  
  I <- y[["I"]]
  Z <- y[["Z"]]
  
  S <- H - I
  X <- A - Z
  #X <- y[["X"]]
  
  phi <- (S + v*I) / (omega*eta*(S + v*epsilon*I))

  dI <- phi*Z*b*S/(S + v*I) - (c + d)*I
  dZ <- phi*X*a*(1 - epsilon*omega)*v*I/(S + v*I) - phi*Z*(1 - a*(1 - epsilon*omega))*v*I/(S + v*I) - phi*Z*S/(S + v*I)
  #dX <- phi*Z*(1 - a*(1 - epsilon*omega))*v*I/(S + v*I) + phi*Z*S/(S + v*I) - phi*(X + Z)*a*(1 - epsilon*omega)*v*I/(S + v*I)
  
  return(list(c(dI, dZ#, dX
                )))
}

parms <- c(
  omega = 0.2,
  eta = 1,
  H = 25,
  A = 25*3,
  b = 1,
  a = 1,
  c = 1.5,
  d = 1.5,
  v = 1,
  epsilon = 1
)

init_I <- 5
init_Z <- parms[["A"]]*parms[["a"]]*(1-parms[["epsilon"]]*parms[["omega"]])*parms[["v"]]*init_I/
  (parms[["v"]]*init_I + parms[["H"]] - init_I)

init_states <- c(I = init_I, Z = init_Z#, X = 0
)


times <- seq(0, 10, by = 0.05)


run <- data.frame(ode(y = init_states,
                      times = times,
                      parms = parms,
                      func = model_ode))

plot(run$time, run$I, ylim=c(0,max(run)), xlab = "time", ylab = "number", type = "l")
lines(run$time, run$Z, col="blue")
#lines(run$time, run$X, col="red")


