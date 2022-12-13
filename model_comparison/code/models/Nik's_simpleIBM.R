#
# Nik Cunniffe
# 9th Dec 2022
#

#
# Clear R environment
#
rm(list = ls())

# DEFINE ODE MODEL
compartmental_don_logic_ode <- function(times, y, parms) {
  
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

#
# Set up parameters
#
#nFact <- 1
N_base <- 70 
A_base <- 140
c_plus_d <- 3  # plant net death rate (=c+d in Elin model)
feedLength <- 1.3
v <- 1
epsilon <- 1
omega <- 0.2

parms <- c(
  H = N_base,
  A = A_base,
  c = c_plus_d/2,
  d = c_plus_d/2,
  eta = feedLength,
  v = v,
  epsilon = epsilon,
  omega = omega,
  a = 1,
  b = 1
)

updateIHat <- T

i0 <- N_base
tMax <- 15
dumpStep <- 0.1

nRuns <- 500

scale <- c(0.25, 0.5, 1, 2, 4, 8)
#num_plants_vec <- c(20, 40, 60, 80, 100, 120, 150)
nik_I_eq_means <- c()
nik_I_eq_medians <- c()
ode_I_eqs <- c()

for (sim in 1:length(scale)) {
  print(paste0("param set ", sim, "/", length(scale)))
  
  # set up parameters (num plants and aphids)
  N <- N_base * scale[sim]
  A <- A_base * scale[sim]
  parms["H"] <- N
  parms["A"] <- A
  
  i0 <- N
  init_I <- i0
  init_Z <- parms[["A"]]*parms[["a"]]*(1-parms[["epsilon"]]*parms[["omega"]])*parms[["v"]]*init_I/
    (parms[["v"]]*init_I + parms[["H"]] - init_I)
  
  init_states <- c(I = init_I, Z = init_Z)
  
  ## RUN ODE and save final I
  trajectory <- data.frame(ode(y = init_states,
                               times = c(0, tMax),
                               parms = parms,
                               func = compartmental_don_logic_ode))
  ode_I_eqs[sim] <- trajectory[nrow(trajectory), "I"]
  
  ### RUN NIK SIMULATION
  allFinalI <- NULL
  for(run in 1:nRuns)
  {
    
    #
    # Set up list
    #
    plantInf <- rep(FALSE,N)
    toInf <- sample(1:N,i0)
    plantInf[toInf] <- TRUE
    I <- i0
    S <- N - i0
    
    #
    # Loop around doing simulation
    #
    t <- 0
    tNextDump <- 0
    while(t < tMax)
    {
      if(t >= tNextDump)
      {
        while(tNextDump < t)
        {
          #s <- sprintf("%.3f %d %d", tNextDump, S, I)
         # print(s)
          tNextDump <- tNextDump + dumpStep 
        }
      }
      rPlant <- N * c_plus_d
      rVect <- A * (1/ feedLength)
      rTot <- rPlant + rVect
      # update time to that of next event
      tDelay <- rexp(1,rTot)
      t <- t + tDelay
      # find which type of event
      if(runif(1) <= rPlant/(rPlant + rVect))
      {
        # next event death of a plant
        # note way this is coded, healthy plants "recover" too, but not important
        plantID <- sample(1:N,1)
        if(plantInf[plantID]) # only update state for those that are infected
        {
          I <- I - 1
          S <- S + 1
        }
        plantInf[plantID] <- FALSE
      }else{
        # next event flight of an aphid
        bLanded <- FALSE
        isInf <- FALSE
  
        iHat <- v * I / ( S + v * I)
        while(!bLanded)
        {
          if(runif(1) <= iHat)
          {
            # probe I
            if(runif(1) <= epsilon * omega)
            {
              # landed on I
              bLanded <- TRUE
            }else{
              isInf <- TRUE
            }
          }else{
            # probe S
            if(runif(1) <= omega)
            {
              bLanded <- TRUE
            } #else{
              if(isInf)
              {
                toInf <- sample(which(plantInf == FALSE),1)
                plantInf[toInf] <- TRUE
                I <- I + 1
                S <- S - 1
                if(updateIHat)
                {
                  iHat <- v * I / ( S + v * I)
                }
              }
              isInf <- FALSE
            #}
          }
        }
      }
    }
    #s <- sprintf("    %d %.3f %d %d", run, t, S, I)
    allFinalI <- c(allFinalI,I)
    #print(s)
  }
  
  #print(mean(allFinalI))
  nik_I_eq_means[sim] <- mean(allFinalI)
  nik_I_eq_medians[sim] <- median(allFinalI)


}

plot(scale, ode_I_eqs, 
     xlab = "population scale (base: H=70, A=140)",
     ylab = paste("final I (after", tMax, "timesteps)"),
     ylim = c(0, max(c(nik_I_eq_means, nik_I_eq_medians, ode_I_eqs)))
)
lines(scale, nik_I_eq_means, col = "red")
lines(scale, nik_I_eq_medians, col = "blue")
legend("topleft", 
       legend = c("Elin's ode", "Nik sim (means)", "Nik sim (medians)"),
       col = c("black", "red", "blue"),
       lwd = 1)

#expected <- c(1.721429e+01,  2.885714e+01,  3.492857e+01,  3.542857e+01,  3.035714e+01,  1.971429e+01, -9.502546e-10)


