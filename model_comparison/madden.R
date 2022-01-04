
#####################################################
##### RECREATION OF MODEL BY MADDEN ET AL. 2000 #####
#####################################################

rm(list=ls())

#packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

# define K
K = 1000 # plant population size

# INITIAL STATES
#PLANTS
L <- 0.002 * K # number of latently infected plants
S <- 0 # number of infectious plants
R <- 0 # number of removed plants
H <- K - L - S - R # number of healthy plants

#VECTORS
P <- K * 4 # total number of vectors
Y <- 0 # number of insects in the latent stage
Z <- 0.005 * P # number of infective insects
X <- P - Y - Z # number of virus-free insects

init_states <- c(H = H,
                 L = L,
                 S = S,
                 R = R,
                 X = X,
                 Y = Y,
                 Z = Z)


# PARAMETERS
phi <- 2 # plants visited per day by an insect
alpha <- 0.2

init_parms <- c(
K = K, # plant population size
phi = phi, # plants visited per day by an insect
T = 0.5/phi, # time feeding/probing per plant visit
k1 = 1/0.021, # inoculation rate by vector (set to NPT virus, equiv to 0.5hr)
k2 = 1/5, # rate of moving through latent period of virus in plant (equiv of 5 day latent period)
k3 = 1/12.5, # rate of moving through infectious period (equiv of 12.5 day infectious period)
beta = 0.01, # plant mortality + replanting rate
lamda = 1/0.021, # rate of acquisition of virus by vector (set to NPT virus, equiv to 0.5hr)
eta = 99999999999, # rate of moving through latent state in vector (set to NPT virus, no latent stage)
tau = 1/0.25, # rate of moving through infectious state in vector (set to NPT virus, equiv to 6hr)
v_t = alpha * (X+Y+Z), # birth rate of vectors per day within field
Ex = 0, # emigration rate for healthy insects
Ey = 0, # emigration rate for latent state insects
Ez = 0, # emigration rate for infective insects
Ix = 0, # immigration rate for healthy insects
Iy = 0, # immigration rate for latent state insects
Iz = 0, # immigration rate for infective insects
q = 0, # probability of vector offspring being viliferous (set to NPT virus)
alpha = alpha # vector population mortality rate
)


#time frame
times <- 1:250

madden_ode_func <- function(times, y, par) {
  
  a <- 1 - exp(-par[["lamda"]] * par[["T"]]) # probability of vector acquisition of virus from an infected plant per visit
  b <- 1 - exp(-par[["k1"]] * par[["T"]]) # probability of inoculation of a plant by an infective vector per visit
  
  # plant equations
  replanting <- par[["beta"]] * (par[["K"]] - y[["H"]])
  L_death <- par[["beta"]] * y[["L"]]
  S_death <- par[["beta"]] * y[["S"]]
  R_death <- par[["beta"]] * y[["R"]]
  
  become_infected <- par[["phi"]] * b * y[["Z"]] * y[["H"]] / par[["K"]]
  become_infectious <- par[["k2"]] * y[["L"]]
  recover <- par[["k3"]] * y[["S"]]
  
  dH <- replanting - become_infected
  dL <- become_infected - L_death - become_infectious
  dS <- become_infectious - S_death - recover
  dR <- recover - R_death

  # dH <- par[["beta"]] * (par[["K"]] - y[["H"]]) - par[["phi"]] * b * y[["Z"]] * y[["H"]] / par[["K"]]
  # dL <- par[["phi"]] * b * y[["Z"]] * y[["H"]] / par[["K"]] - (par[["k2"]] + par[["beta"]]) * y[["L"]]
  # dS <- par[["k2"]] * y[["L"]] - (par[["k3"]] + par[["beta"]]) * y[["S"]]
  # dR <- par[["k3"]] * y[["S"]] - par[["beta"]] * y[["R"]]

  # vector equations
  # birth and death
  birth <- par[["v_t"]]
  born_infective <- par[["v_t"]] * par[["q"]] * (y[["Z"]] / (y[["X"]] + y[["Y"]] + y[["Z"]]))
  death_X <- par[["alpha"]] * y[["X"]]
  death_Y <- par[["alpha"]] * y[["Y"]]
  death_Z <- par[["alpha"]] * y[["Z"]]
  
  # migration
  immigration_X <- par[["Ix"]]
  immigration_Y <- par[["Iy"]]
  immigration_Z <- par[["Iz"]]
  emigration_X <- par[["Ex"]] * y[["X"]]
  emigration_Y <- par[["Ey"]] * y[["Y"]]
  emigration_Z <- par[["Ez"]] * y[["Z"]]
  
  # virus
  acquisition <- par[["phi"]] * a * y[["S"]] * y[["X"]] / par[["K"]]
  become_infective <- par[["eta"]] * y[["Y"]]
  stop_being_infective <- par[["tau"]] * y[["Z"]]

  dX <- birth - born_infective - death_X + immigration_X - emigration_X - acquisition + stop_being_infective
  dY <- acquisition - become_infective + immigration_Y - emigration_Y - death_Y
  dZ <- become_infective - stop_being_infective + immigration_Z - emigration_Z + born_infective - death_Z

  # dX <- par[["v_t"]] - par[["v_t"]] * par[["q"]] * (y[["Z"]] / (y[["X"]] + y[["Y"]] + y[["Z"]])) + par[["tau"]] * y[["Z"]] +
  #   par[["Ix"]] - par[["Ex"]] * y[["X"]] - par[["alpha"]] * y[["X"]] - par[["phi"]] * a * y[["S"]] * y[["X"]] / par[["K"]]
  # dY <- par[["phi"]] * a * y[["S"]] * y[["X"]] / par[["K"]] + par[["Iy"]] - par[["Ey"]] * y[["Y"]] - par[["alpha"]] * y[["Y"]] -
  #   par[["eta"]] * y[["Y"]]
  # dZ <- par[["eta"]] * y[["Y"]] - par[["tau"]] * y[["Z"]] + par[["Iz"]] - par[["Ez"]] * y[["Z"]] +
  #   par[["v_t"]] * par[["q"]] * (y[["Z"]] / (y[["X"]] + y[["Y"]] + y[["Z"]])) - par[["alpha"]] * y[["Z"]]

  
  return(list(c(dH, dL, dS, dR, dX, dY, dZ)))
} 

trajectory <- data.frame(ode(y = init_states, 
                             times = times, 
                             parms = init_parms, 
                             func = madden_ode_func))

trajectory_long <- reshape2::melt(trajectory, id.vars="time")
names(trajectory_long) <- c("time", "compartment", "number")
  
plant_trajec <- ggplot(data=trajectory_long %>% filter(compartment == "H" |
                                              compartment == "L" |
                                              compartment == "S" |
                                              compartment == "R"), 
                       aes(x = time, y = number, col = compartment)) +
  geom_line() +
  labs(y = "number of plants")
plant_trajec

vector_trajectory <- ggplot(data=trajectory_long %>% filter(compartment == "X" |
                                                             compartment == "Y" |
                                                             compartment == "Z"), 
                            aes(x = time, y = number, col = compartment)) +
  geom_line() +
  labs(y = "Nnumber of vectors")
vector_trajectory

