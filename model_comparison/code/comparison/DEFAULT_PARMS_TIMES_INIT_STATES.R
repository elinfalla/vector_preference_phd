# TIMEFRAME
times <- seq(0, 40, by = 0.4)

# INITIAL STATES
I <- 1 # initial number of infected plants

init_states_don <- c(
  I = I # number of infected plants
)

init_states_mad <- c(
  I = I, # number of infected plants
  Z = 0 # number infective insects
)

# PARAMETERS
parms <- c(
  
  ## DONNELLY PARMS
  gamma = 4, # rate of recovery/death of I plants
  theta = 1.041, # aphid dispersal rate per day
  Pacq = 0.8, # chance of virus acquisition by vector from infected plant
  Pinoc = 0.8, # chance of inoculating healthy plant
  w_don = 0.2, # feeding probability on healthy plant
  
  ## MADDEN PARMS
  # k1 = 1/0.021, # inoculation rate by insect (set to NPT virus, equiv to 0.5hr)
  # lambda = 1/0.021, # rate of acquisition of virus by vector (set to NPT virus, equiv to 0.5hr)
  # T = 0.2, # time feeding/probing per plant visit (set to 0.5/phi in Madden model)
  tau = 4, # rate of moving through infectious state in vector
  w_mad = 0.2, # feeding probability on healthy plant
  phi = 5.473, # plants visited per day by an insect
  
  ## PARMS FOR BOTH
  A = 1200, # total (constant) number of vectors
  H = 400, # number host plants 
  v = 1, # infected plant attractiveness
  e = 1 # infected plant acceptability 
  
)

### SET ETA SO IT'S THE SAME AS MADDEN UNDER NO VEC PREF
parms[["eta"]] <- 1/(parms[["phi"]]*parms[["w_mad"]])

### SET MADDEN PARMS THAT ARE DETERMINED BY DONNELLY PARMS
# probability of vector acquisition of virus from an infected plant per visit
parms[["a"]] <- parms[["Pacq"]] #1 - exp(-parms[["lamda"]] * parms[["T"]])

# probability of plant inoculation per infective insect visit 
parms[["b"]] <- parms[["Pinoc"]] #1 - exp(-parms[["k1"]] * parms[["T"]])

# natural plant death rate (equivalent to beta in original Madden model)
parms[["c"]] <- parms[["gamma"]] / 2

# plant death due to infection
parms[["d"]] <- parms[["gamma"]] / 2


