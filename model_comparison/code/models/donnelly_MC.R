
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
  probes <- 0
  feeding <- F
  
  ## MARKOV CHAIN ##
  while (feeding == F) {
    
    
    if (state == "flight") {
      
      # use transition matrix to work out new state of markov chain
      state <- sample(states, 1, replace=T, prob = transition_matrix["flight",])
    }
    
    else if (state == "probeI") {
      
      state <- sample(states, 1, replace=T, prob = transition_matrix["probeI",])
      probes <- probes + 1
      
      if (state == "probeS") {
        transmissions <- transmissions + 1 # probeI -> probeS means transmission has taken place
        probes <- probes + 1
      }
      else if (state == "feedI") {
        feeding <- T
      }
    }
    
    else if (state == "probeS") {
      
      state <- sample(states, 1, replace=T, prob = transition_matrix["probeS",])
      probes <- probes + 1
      
      if (state == "feedS") {
        feeding <- T
      }
    }
    
    # add state to vector of all states
    state_changes <- c(state_changes, state)
  }
  #print(state_changes)
  #length(state_changesx)
  
  return(c(transmissions, probes))
}



## runs num_runs feeding dispersals and returns number of transmissions per dispersal
run_simulations <- function(num_runs, I, v, w ,e) {
  
  transmission_results <- data.frame(num_transmissions = vector(length=num_runs),
                                     num_probes = vector(length=num_runs))
  
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
    
    num_transmissions <- feeding_dispersal(I, v, w, e)[1]
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


v <- 1
I <- 0.5
w <- 0.2
e <- 1
i_hat <- (v*I)/(1-I + v*I)
set.seed(43)
result <- run_simulations(10000, I=I, v=v, w=w, e=e)

transmissions_count <- as.data.frame(table(result$num_transmissions))
names(transmissions_count)[1] <- "num_transmissions"
transmissions_count$num_transmissions <- as.numeric(transmissions_count$num_transmissions) - 1

transmissions_barplot <- ggplot(data.frame(result), aes(x = num_transmissions)) +
  geom_bar(fill="black") +
  #geom_line(data = result_count, aes(x = result, y = Freq), colour = "red", size = 0.8) +
  labs(x = "Number of transmissions per aphid feeding dispersal",
       y = "Number of runs",
       title = "c)") +
  theme_bw() +
  theme(text = element_text(size = 15),
    strip.background = element_blank(), 
    )

probes_count <- as.data.frame(table(result$num_probes))
names(probes_count)[1] <- "num_probes"
probes_count$num_probes <- as.numeric(probes_count$num_probes) - 1

probes_barplot <- ggplot(data.frame(result), aes(x = num_probes)) +
  geom_bar(fill="black", width=0.4) +
  labs(x = "Number of plant probes per aphid feeding dispersal",
       y = "Number of runs",
       title = "b)") +
  theme_bw() +
  theme(text = element_text(size = 15),
        strip.background = element_blank(), 
  )

pdf("results/MC_simulation_barplot.pdf")
grid.arrange(probes_barplot, transmissions_barplot)
dev.off()

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
