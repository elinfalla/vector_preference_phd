
#########################################################
##### SIMPLIFICATION OF MODEL BY MADDEN ET AL. 2000 #####
#########################################################

### Simplified version of model by Madden et al. 2000 to make it more comparable to model by Donnelly et al. 2019.
### No vector immigration, plants are represented by an SI rather than an SEIR model.

library(deSolve)

# PARAMETERS
parms <- c(
  N = .., # number of host plants
  P = .., # total number of vectors
  phi = .., # plants visited per day by insect
  k1 = .., # inoculation rate by insect (set to NPT virus, equiv to 0.5hr)
  T = .., # # time feeding/probing per plant visit
  K = .., # plant pop size
  c = .., # natural plant death rate
  d = .., # plant death due to infection
  Ex = .., # healthy vector emigration rate
  Ey = .., # latent vector emigration rate
  Ez = .., # infective vector emigration rate
)

# INITIAL STATES
init_states <- c(
  I = .., # number of infected plants
  X = .., # number healthy insects
  Y = .., # number latent insects
  Z = .., # number infective insects
)


madden_simp <- function(times, y, parms) {
  
  # plant equations
  b <- 1 - exp(-par[["k1"]] * par[["T"]]) # prob of plant inoculation per infective insect visit 
  
  become_infected <- par[["phi"]] * b * y[["Z"]] * y[["H"]] / par[["K"]]
  
  dI <- become_infected - par[["c"]]*I - par[["d"]]*I
  
  # vector equations
  
}