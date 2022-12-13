rm(list=ls())

library(deSolve)
library(ggplot2)
library(dplyr)

# source ode functions - I am using donnelly_vpref_ode() and compartmental_don_logic_ode()
source("./code/comparison/model_ode_versions_FUNCTIONS.R")

times <- seq(0, 10, length.out = 200)

num_plants <- 70
num_aphids <- num_plants * 2
plant_recovery_rate <- 3
acquisition_rate <- 1
inoculation_rate <- 1
flight_rate <- 1/1.3
plant_attractivenes <- 1
plant_acceptability <- 1
feeding_probability <- 0.2

parms <- c(
  v = plant_attractivenes,
  epsilon = plant_acceptability,
  e = plant_acceptability, # same as epsilon
  H = num_plants,
  A = num_aphids,
  w = feeding_probability,
  w_don = feeding_probability,
  omega = feeding_probability, # same as w
  
  # donnelly parms
  gamma = plant_recovery_rate,
  theta = flight_rate,
  Pacq = acquisition_rate,
  Pinoc = inoculation_rate,
  
  # new model parms
  eta = 1/flight_rate,
  b = inoculation_rate,
  a = acquisition_rate,
  c = plant_recovery_rate/2,
  d = plant_recovery_rate/2
)

init_I <- num_plants
init_Z <- num_aphids*acquisition_rate*(1 - parms[["epsilon"]]*parms[["omega"]])*
  parms[["v"]]*init_I / (parms[["v"]]*init_I + num_plants - init_I)

don_init_states <- c(I = init_I)
new_init_states <- c(I = init_I, Z = init_Z)


don_trajec <- data.frame(ode(y = don_init_states,
                             times = times,
                             parms = parms,
                             func = donnelly_vpref_ode))
new_mod_trajec <- data.frame(ode(y = new_init_states,
                                 times = times,
                                 parms = parms,
                                 func = compartmental_don_logic_ode))

combined_I <- data.frame(time = c(don_trajec$time, new_mod_trajec$time),
                       I = c(don_trajec$I, new_mod_trajec$I),
                       model = rep(c("Donnelly", "New"), each = nrow(don_trajec)))

trajec_plot <- ggplot(data = combined_I, aes(x = time, y = I, col = model)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.2),
        legend.title = element_blank()) +
  labs(x = "Time (arbitrary units)",
       y = "Number infected plants")
trajec_plot

## PLOT NUMBER INFECTIVE APHIDS OVER TIME
don_trajec <- don_trajec %>%
  mutate(Z = (I/parms[["H"]]*parms[["a"]]*(1-parms[["omega"]])*parms[["A"]]))

combined_Z <- data.frame(time = c(don_trajec$time, new_mod_trajec$time),
                         Z = c(don_trajec$Z, new_mod_trajec$Z),
                         model = rep(c("Donnelly", "New"), each = nrow(don_trajec)))

aphid_plot <- ggplot(data = combined_Z, aes(x = time, y = Z, col = model)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.2),
        legend.title = element_blank()) +
  labs(x = "Time (arbitrary units)",
       y = "Number infective aphids")
aphid_plot

###### GET I EQUILIBRIUM FOR VARYING NUM PLANTS
num_plants_vec <- c(20, 40, 60, 80, 100, 120, 150)
eq_I <- rep(NA, length(num_plants_vec))

for (run in 1:length(num_plants_vec)) {
  if (run == 7) {browser()}
  don_init_states <- c(I = num_plants_vec[run])
  timeframe <- c(0, 400)
  parms["H"] <- num_plants_vec[run]
  
  trajec <- data.frame(ode(y = don_init_states,
                               times = timeframe,
                               parms = parms,
                               func = donnelly_vpref_ode))
  
  eq_I[run] <- trajec[nrow(trajec), "I"]

}
print(eq_I)


# init_Z_vals <- c(12, 25, 35, 45, 55)
# all_trajecs <- data.frame(time = don_trajec$time,
#                           I = don_trajec$I,
#                           model = "Donnelly model") # add don model to all_trajecs
# 
# 
# for (Z in init_Z_vals) {
#   init_states <- c(I = init_I, Z = Z)
#   new_mod_trajec <- data.frame(ode(y = init_states,
#                                    times = times,
#                                    parms = parms,
#                                    func = compartmental_don_logic_ode))
#   
#   all_trajecs <- rbind(all_trajecs, data.frame(time = new_mod_trajec$time, 
#                                                I = new_mod_trajec$I, 
#                                                model = rep(paste0("New model, Z=", Z), nrow(new_mod_trajec))))
#   
# }
# 
# plot <- ggplot(data = all_trajecs, aes(x = time, y = I, col = model)) +
#   geom_line(aes(alpha = model)) + 
#   theme_bw() +
#   theme(legend.position = c(0.75, 0.2),
#         legend.title = element_blank()) +
#   labs(x = "Time (arbitrary units)",
#        y = "Number infected plants") +
#   scale_color_manual(values = c("black", rainbow(length(init_Z_vals)))) +
#   scale_alpha_manual(values = c(1, rep(0.7, length(init_Z_vals))))
# plot

# pdf("./results/don_vs_don_logic_compartment_ode.pdf")
# plot
# dev.off()
