
#######################################################################
##### RECREATION OF STOCHASTIC SIMULATION BY DONNELLY ET AL. 2019 #####
#######################################################################
# winged aphids only

rm(list=ls())
library(matrixStats)
  
run_simulation <- function(t, timeframe, I, A, parms) {

  # initialise vectors of A, I and time vals
  I_vals <- I
  A_vals <- A
  time <- t
  
  # define vector preference parameters
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w"]]
  
  while (t < timeframe) {
    
    i_plant_recovery <- parms[["gamma"]] * I
    aphid_death <- A * parms[["b"]]
    aphid_dispersal <- A * parms[["theta"]]
    aphid_reproduction <- parms[["a"]] * A * (1 - A/parms[["K"]])

    total_rate <- sum(i_plant_recovery, aphid_death, aphid_dispersal, aphid_reproduction)

    # generate time step according to gillespie (exponential distribution)
    dt <- -log(runif(1, min=0, max=1)) / total_rate

    t <- t + dt
    
    events <- c("i_plant_recovery", "aphid_death", "aphid_dispersal", "aphid_reproduction")
    event_rates <- c(i_plant_recovery, aphid_death, aphid_dispersal, aphid_reproduction)
    
    current_event <- sample(events, 1, replace=T, prob=event_rates)
    
    if (current_event == "i_plant_recovery") {
      I <- I - 1
      #I_vals <- c(I_vals, I)
      
    }
    else if (current_event == "aphid_death") {
      A <- A - 1
      #A_vals <- c(A_vals, A)
      
    }
    else if (current_event == "aphid_reproduction") {
      A <- A + 1
      #A_vals <- c(A_vals, A)
      
    }
    else { # current_event == "aphid_dispersal"
      S <- parms["H"] - I
      i_hat <- (v*I)/(S + v*I)
      
      states <- c("probeS", "probeI", "feedS", "feedI", "flight")
      
      transition_matrix <- matrix(c(((1-w)*(1-i_hat)), ((1-w)*i_hat), w, 0,   0,
                                    (1-e*w)*(1-i_hat), (1-e*w)*i_hat, 0, e*w, 0,
                                    0,                 0,             1, 0,   0,
                                    0,                 0,             0, 1,   0,
                                    1-i_hat,           i_hat,         0, 0,   0),
                                  nrow = 5, byrow=T,
                                  dimnames = list(states))
      
      state <- "flight" # dispersal event always starts with flight from plant
      #state_changes <- "flight"
      #transmissions <- 0
      feeding <- F
      
      ## MARKOV CHAIN ##
      while (feeding == F) {
        
        # if uniform rv is less than p, aphid dies/emigrates. break out of feeding dispersal
        if (runif(1, min=0, max=1) < parms["p"]) {
          A <- A - 1
         # A_vals <- c(A_vals, A)
          break
        }
        
        if (state == "flight") {
          
          # use transition matrix to work out new state of markov chain
          state <- sample(states, 1, replace=T, prob = transition_matrix["flight",])
        }
        
        else if (state == "probeI") { # aphid has now acquired the virus
          
          state <- sample(states, 1, replace=T, prob = transition_matrix["probeI",])
          
          if (state == "probeS") {
            I <- I + 1 # inoculation event
            #I_vals <- c(I_vals, I)
            #transmissions <- transmissions + 1 # probeI -> probeS means transmission has taken place
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
        #state_changes <- c(state_changes, state)
      }
      
      
    }
    

    #}
    
    #S_vals <- c(S_vals, S)
    I_vals <- c(I_vals, I)
    A_vals <- c(A_vals, A)
    time <- c(time, t)
    
  }
  
  run_df <- data.frame(time, I_vals, A_vals)

  #find S, I and R levels at set times and put in data frame to return from function
  run.I.times <- approx(x = run_df$time,
                        y = run_df$I_vals,
                        xout = time_points,
                        method = "constant")
  run.I.times.df <- data.frame(run.I.times[1], run.I.times[2])
  colnames(run.I.times.df) <- c("time", "num_infected")

  run.A.times <- approx(x = run_df$time,
                        y = run_df$A_vals,
                        xout = time_points,
                        method = "constant")
  run.A.times.df <- data.frame(run.A.times[1], run.A.times[2])
  colnames(run.A.times.df) <- c("time", "num_aphids")



  return(list(run.I.times.df, run.A.times.df))
  #return(list(I_vals, A_vals))
}

# define time frame
t <- 0
timeframe <- 8
time_points <- seq(t, timeframe, by = 0.5)

# define parameters
parms <- c(
  gamma = 20/3, # rate of recovery/death of I plants
  b = 0.1, # aphid mortality rate per day
  theta = 1, # aphid dispersal rate per day
  a = 2, # reproduction rate per aphid per day
  K = 4000, # aphid reproduction limit - maximum aphids PER FIELD (diff to donnelly)
  H = 400, # number host plants
  p = 0.2, # aphid emigration/death rate per journey between plants
  v = 1.2, # infected plant attractiveness
  e = 0.8, # infected plant acceptability 
  w = 0.2 # feeding rate on healthy plant
)

# define states
A <- 1200 # num aphids
I <- 40 # num infected plants

# plot a single run
# out <- run_simulation(t, timeframe, I, A, parms)
#   
# plot(out[[1]]$time, out[[1]]$num_infected)
# plot(out[[2]]$time, out[[2]]$num_aphids)

# run 50 times in order to find median run
for (run in 1:50) {
  print(run)
  if (run == 1) {
    sims <- run_simulation(t, timeframe, I, A, parms)[[1]]
  } else {
    sims[,run+1] <- run_simulation(t, timeframe, I, A, parms)[[1]]$num_infected
  }
}

# plot number of infected plants over time
sims$median <- rowMedians(as.matrix(sims[,-1], ncol = ncol(sims)-1))
plot(sims$time, sims$median)

# plot frequency of infected plants (I/H) over time
sims$median_freq <- sims$median / parms[["H"]]
plot(sims$time, sims$median_freq)


