
###################################################################################
###### DETERMINING APHID POPULATION EQUILIBRIUM -  FROM DONNELLY ET AL. 2019 ######
###################################################################################

### Aphid population ODE solve

parms <- c(
  i = 0.5, # proportion of host plants that are infected
  v = 1, # infected plant attractiveness
  e = 1, # infected plant acceptability 
  a = 2, # reproduction rate per aphid per day
  b = 0.1, # aphid mortality rate per day
  K = 10, # aphid reproduction limit - maximum aphids per plant (diff to donnelly_stoch_sim.R)
  theta = 1, # aphid dispersal rate per day
  
)



aphid_pop_ODE <- function(times, y, parms) {
  
  
  
  # i_hat = weighted frequency of infected plants, accounting for attraction towards infected plants
  # i.e. adapted version of i / (1 - i + i) = i / 1, to include v
  i_hat <- v*i / ((1-i) + v*i) 
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  # STATE EQUATIONS
  # aphids - change in per plant aphid density on susceptible (As)  and infected (Ai) plants
  # = aphid births - aphid deaths - aphid settling on other plant type + aphids moving from other plant type
  dAs <- a*As*(1 - As/K) - b*As - theta*As*(1 - Fs) + theta*Ai*Fs*i/(1 - i)
  dAi <- a*Ai*(1 - Ai/K) - b*Ai - theta*Ai*(1 - Fi) + theta*As*Fi*(1 - i)/i
  
}

