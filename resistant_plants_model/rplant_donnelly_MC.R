
rm(list = ls())

library(dplyr)

## function to simulate feeding dispersal of one aphid, starting with flight and ending with feeding
feeding_dispersal <- function(I, R, v_i, v_r, w, e, H) {
  
  S <- H - I - R
  r_hat <- v_r*R / (S + v_i*I + v_r*R)
  i_hat <- v_i*I / (S + v_i*I + v_r*R)
  s_hat <- 1 - i_hat - r_hat
  
  states <- c("probeS", "probeI", "probeR", "feedS", "feedI", "feedR", "flight")
  
  transition_matrix <- matrix(c((1-w)*s_hat,   (1-w)*i_hat,   (1-w)*r_hat,   w, 0,   0, 0,
                                (1-e*w)*s_hat, (1-e*w)*i_hat, (1-e*w)*r_hat, 0, e*w, 0, 0,
                                (1-w)*s_hat,   (1-w)*i_hat,   (1-w)*r_hat,   0, 0,   w, 0,
                                0,             0,             0,             1, 0,   0, 0,
                                0,             0,             0,             0, 1,   0, 0,
                                0,             0,             0,             0, 0,   1, 0,
                                s_hat,         i_hat,         r_hat,         0, 0,   0, 0),
                              nrow = 7, byrow=T,
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
    
    else if (state == "probeR") {
      
      state <- sample(states, 1, replace=T, prob = transition_matrix["probeR",])
      
      if (state == "feedR") {
        feeding <- T
      }
    }
    
    # add state to vector of all states
    state_changes <- c(state_changes, state)
  }
 
  return(transmissions)
}



## runs num_runs feeding dispersals and returns number of transmissions per dispersal
run_simulations <- function(num_runs, I, R, v_i, v_r, w, e, H) {
  
  transmission_results <- data.frame(num_transmissions = vector(length=num_runs))
  
  for (run in 1:num_runs) {
    transmission_results[run,] <- feeding_dispersal(I, R, v_i, v_r, w, e, H)
  }
  
  return(transmission_results)
}

I <- 30
R <- 60
v_i <- 1.4
v_r <- 1 
w <- 0.2
e <- 1 
H <- 100
num_runs <- 100000
S <- H - I - R

output <- run_simulations(num_runs, I, R, v_i, v_r, w, e, H)

i_hat <- v_i*I / (S + v_i*I + v_r*R)
r_hat <- v_r*R / (S + v_i*I + v_r*R)
s_hat <- 1 - i_hat - r_hat

x_i <- i_hat*(1-e*w)*s_hat / ((s_hat+r_hat)*w + i_hat*e*w)

exp_transmissions <- output %>% 
  count(num_transmissions) %>%
  mutate(sum = num_transmissions*n)


sum(exp_transmissions$sum)/num_runs
x_i_old
x_i_new

