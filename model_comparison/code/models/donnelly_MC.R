
##############################################################
##### RECREATION OF MARCOV CHAIN BY DONNELLY ET AL. 2019 #####
##############################################################

rm(list = ls())
library(ggplot2)

## function to simulate feeding dispersal of one aphid, starting with flight and ending with feeding
feeding_dispersal <- function(I, v, w, e) {
  
  S <- 1 - I
  i_hat <- (v*I)/(S + v*I)
  
  states <- c("probeS", "probeI", "feedS", "feedI", "flight")
  
  transition_matrix <- matrix(c(((1-w)*(1-i_hat)), ((1-w)*i_hat), w, 0,   0,
                                (1-e*w)*(1-i_hat), (1-e*w)*i_hat, 0, e*w, 0,
                                0,                 0,             1, 0,   0,
                                0,                 0,             0, 1,   0,
                                1-i_hat,           i_hat,         0, 0,   0),
                              nrow = 5, byrow=T,
                              dimnames = list(states))
  
  state <- "flight"
  state_changes <- "flight"
  transmissions <- 0
  feeding <- F
  
  ## MARKOV CHAIN ##
  while (feeding == F) {
    
    
    if (state == "flight") {
      
      # use transition matrix to work out new state of markov chain
      state <- sample(states, 1, replace=T, prob = transition_matrix["flight",])
    }
    
    else if (state == "probeI") {
      
      state <- sample(states, 1, replace=T, prob = transition_matrix["probeI",])
      
      if (state == "probeS") {
        transmissions <- transmissions + 1 # probeI -> probeS means transmission has taken place
      }
      else if (state == "feedI") {
        feeding <- T
      }
    }
    
    else if (state == "probeS") {
      
      state <- sample(states, 1, replace=T, prob = transition_matrix["probeS",])
      
      if (state == "feedS") {
        feeding <- T
      }
    }
    
    # add state to vector of all states
    state_changes <- c(state_changes, state)
  }
  #print(state_changes)
  #length(state_changes)
  
  return(transmissions)
}



## runs num_runs feeding dispersals and returns number of transmissions per dispersal
run_simulations <- function(num_runs, I, v, w ,e) {
  
  transmission_results <- data.frame(num_transmissions = vector(length=num_runs))
  
  for (run in 1:num_runs) {
    transmission_results[run,] <- feeding_dispersal(I, v, w, e)
  }
  
  return(transmission_results)
}



## runs num_runs feeding dispersals where the num transmissions feeds back into number infected plants
## and returns number of transmissions per dispersal + num infected at each step
run_simulations_feedback <- function(num_runs, init_I, num_plants, v, w, e) {
  
  # initialise results dataframe
  transmission_results <- data.frame(time = 1:num_runs,
                                     num_transmissions = vector(length = num_runs),
                                     freq_I = vector(length = num_runs))
  I <- init_I
  
  for (run in 1:num_runs) {
    
    num_transmissions <- feeding_dispersal(I, v, w, e)
    transmission_results[run, 2:3] <- c(num_transmissions, I)
    I <- (I*num_plants + num_transmissions)/num_plants # update I: convert to number from freq, add transmissions from dispersal, convert back to freq
    if (I > 1) {
      transmission_results <- transmission_results[-(run+1):num_runs,] # remove empty rows from end of dataframe
      break
    }
  }
  
  #hist(transmission_results, xlab = "Number of transmissions per feeding dispersal")
  return(transmission_results)
}


v <- 1.2
I <- 0.5
w <- 0.2
e <- 0.8
i_hat <- (v*I)/(1-I + v*I)

result <- run_simulations(10000, I=I, v=v, w=w, e=e)

# plot expected geometric distribution (model is zero-deflated), p = probability of success = probability of feeding = 1-prob of probing = 1-P_Sk
s_hat <- 1-i_hat
P_Sk <- ((1-w)*(1-e*w)*i_hat*s_hat)/((1-(1-e*w)*i_hat)*(1-(1-w)*s_hat)) # from Donnelly 2019, Appendix 1, Eq S8

geom_dist_points <- data.frame(rgeom(n = 10000, p = 1-P_Sk))
names(geom_dist_points) <- "num_transmissions"

ggplot(data.frame(result), aes(x = num_transmissions)) + 
  geom_histogram(binwidth = 1, fill="blue", color="blue") +
  geom_histogram(data = geom_dist_points, aes(x = num_transmissions), binwidth = 1, fill=NA, color="red") +
  labs(x="n, number of transmissions per dispersal")
  
### now run simulations with feedback on I
num_plants <- 3000
I <- 2/num_plants

feedback_result <- run_simulations_feedback(10000, init_I=I, num_plants=num_plants, v=v, w=w, e=e)

ggplot(data=feedback_result, aes(x = num_transmissions)) + 
  geom_histogram(binwidth = 1, fill="green", color="green") +
  geom_histogram(data = geom_dist_points, aes(x = num_transmissions), binwidth = 1, fill=NA, color="red") +
  labs(x="n, number of transmissions per dispersal")
ggplot(feedback_result, aes(x = time, y = freq_I)) +
  geom_point()
