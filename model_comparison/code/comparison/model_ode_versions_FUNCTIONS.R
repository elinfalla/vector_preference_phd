madden_vpref_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  become_infected <- par[["phi"]] * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  
  # VECTOR EQUATIONS
  X <- par[["A"]] - y[["Z"]]
  
  acquisition <- par[["phi"]] * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * X * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
  
}

donnelly_vpref_ode <- function(times, states, parms) {
  
  # STATES
  I <- states[["I"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  theta <- parms[["theta"]]
  H <- parms[["H"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w_don"]]
  Pacq <- parms[["Pacq"]]
  Pinoc <- parms[["Pinoc"]]
  A <- parms[["A"]]
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  
  # I_hat = weighted number of infected plants, accounting for attraction towards infected plants
  i_hat <- v*I / (H - I + v*I)
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  #xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  xI <- (Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(w*(1 - i_hat*(1 - e)))
  
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A)*xI - gamma*I
  
  return(list(di))
}

madden_cunniffe_vpref_ode <- function(times, y, par) {
  
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  # define variable phi
  phi <- (S + y[["I"]]*par[["v"]]) / (par[["w_mad"]]*par[["eta"]]*(S + y[["I"]]*par[["v"]]*par[["e"]]))
  
  # PLANT EQUATIONS
  
  become_infected <- phi * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  
  # VECTOR EQUATIONS
  X <- par[["A"]] - y[["Z"]]
  
  acquisition <- phi * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * X * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  # state equations
  dZ <- acquisition - stop_being_infective
  
  return(list(c(dI, dZ)))
}


madden_vdynamic_ode <- function(times, y, par) {
  
  # PLANT EQUATIONS
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  become_infected <- par[["phi"]] * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  
  # VECTOR EQUATIONS
  A <- y[["X"]] + y[["Z"]]
  
  acquisition <- par[["phi"]] * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * y[["X"]] * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  birth <- par[["lamda"]] * A * (1 - ((A/par[["H"]])/par[["K"]]))
  death_X <- par[["alpha"]] * y[["X"]]
  death_Z <- par[["alpha"]] * y[["Z"]]
  
  emigration_X <- par[["phi"]] * par[["p"]] * y[["X"]]
  emigration_Z <- par[["phi"]] * par[["p"]] * y[["Z"]]
  
  # state equations
  dX <- stop_being_infective - acquisition + birth - death_X - emigration_X
  dZ <- acquisition - stop_being_infective - death_Z - emigration_Z
  
  
  return(list(c(dI, dX, dZ)))
}


madden_cunniffe_vdynamic_ode <- function(times, y, par) {
  
  S <- par[["H"]] - y[["I"]] # number of susceptible plants
  
  # define variable phi
  phi <- (S + y[["I"]]*par[["v"]]) / (par[["w_mad"]]*par[["eta"]]*(S + y[["I"]]*par[["v"]]*par[["e"]]))
  
  # PLANT EQUATIONS
  become_infected <- phi * par[["b"]] * y[["Z"]] * S/(S + par[["v"]]*y[["I"]])
  
  natural_death_I <- par[["c"]] * y[["I"]]
  virus_induced_death <- par[["d"]] * y[["I"]]
  
  # state equation
  dI <- become_infected - natural_death_I - virus_induced_death
  #di <- dI/parms[["H"]]
  
  # VECTOR EQUATIONS
  A <- y[["X"]] + y[["Z"]]
  
  acquisition <- phi * par[["a"]] * (1 - par[["e"]]*par[["w_mad"]]) * y[["X"]] * par[["v"]]*y[["I"]] /
    (S + par[["v"]]*y[["I"]])
  stop_being_infective <- par[["tau"]] * y[["Z"]]
  
  birth <- par[["lamda"]] * A * (1 - ((A/par[["H"]])/par[["K"]]))
  death_X <- par[["alpha"]] * y[["X"]]
  death_Z <- par[["alpha"]] * y[["Z"]]
  
  emigration_X <- phi * par[["p"]] * y[["X"]]
  emigration_Z <- phi * par[["p"]] * y[["Z"]]
  
  # state equations
  dX <- stop_being_infective - acquisition + birth - death_X - emigration_X
  dZ <- acquisition - stop_being_infective - death_Z - emigration_Z
  
  
  return(list(c(dI, dX, dZ)))
}


donnelly_full_vdynamic_ode <- function(times, states, parms) {
  
  ### ODE function for deterministic Donnelly model
  
  # STATES
  i <- states[["i"]]
  As <- states[["As"]]
  Ai <- states[["Ai"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  alpha <- parms[["alpha"]]
  theta <- parms[["theta"]]
  lamda <- parms[["lamda"]]
  K <- parms[["K"]]
  H <- parms[["H"]]
  p <- parms[["p"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w_don"]]
  Pacq <- parms[["Pacq"]]
  Pinoc <- parms[["Pinoc"]]
  
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  
  # q = prob of surviving flight between plants (i.e. doesn't emigrate/die)
  q <- 1 - p  
  
  # i_hat = weighted frequency of infected plants, accounting for attraction towards infected plants
  # i.e. adapted version of i / (1 - i + i) = i / 1, to include v
  i_hat <- v*i / ((1-i) + v*i) 
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  # A = total number of aphids i.e. aphids per plant type * number of that plant type
  # i.e. S*As + I*Ai where S is number of susceptible plants and I is number of infected plants
  A <- As*(H - i*H) + Ai*(i*H)  
  
  
  # STATE EQUATIONS
  # aphids - change in per plant aphid density on susceptible (As)  and infected (Ai) plants
  # = aphid births - aphid deaths - aphid settling on other plant type + aphids moving from other plant type
  dAs <- lamda*As*(1 - As/K) - alpha*As - theta*As*(1 - Fs) + theta*Ai*Fs*i/(1 - i)
  dAi <- lamda*Ai*(1 - Ai/K) - alpha*Ai - theta*Ai*(1 - Fi) + theta*As*Fi*(1 - i)/i
  
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(c(di, dAs, dAi)))
}

calc_A_equilibrium <- function(parms, states, i_hat, return_all=F) {
  
  # define states
  i <- states[["i"]]
  
  # define parameters
  e <- parms[["e"]]
  lamda <- parms[["lamda"]]
  alpha <- parms[["alpha"]]
  K <- parms[["K"]]
  theta <- parms[["theta"]]
  
  # probability of settling on susceptible (Fs) and infected (Fi) plants,
  # derived from analysis of Markov chain, see Donnelly et al. 2019, Appendix S4
  Fs <- (1 - i_hat) / (1 - i_hat*(1 - e)) 
  Fi <- e*i_hat / (1 - i_hat*(1 - e))
  
  ## define coefficients 
  # (worked out by hand from dAs and dAi equations- see lab book). will give roots of Ai
  coeff_1 <- i*(lamda - alpha - theta*(1 - Fs))*(alpha + theta*(1 - Fi) - lamda)/(theta*Fi*(1 - i)) + theta*Fs*i/(1 - i)
  
  coeff_2 <- (lamda/K)*i*(lamda - alpha - theta*(1 - Fs))/(theta*Fi*(1 - i)) - 
    lamda*i^2*(alpha^2 + 2*alpha*theta*(1 - Fi) - 2*alpha*lamda + theta^2*(1 - Fi)^2 - 2*lamda*theta*(1 - Fi) + lamda^2)/(K*theta^2*(1 - i)^2*Fi^2)
  
  coeff_3 <- -(lamda*i^2/K)*(2*alpha*lamda + 2*lamda*theta*(1 - Fi) - 2*lamda^2)/(K*theta^2*(1 - i)^2*Fi^2)
  
  coeff_4 <- -lamda^3*i^2/(K^3*theta^2*(1 - i)^2*Fi^2)
  
  # get Ai equilibrium values
  Ai_equilibriums <- polyroot(c(0, coeff_1, coeff_2, coeff_3, coeff_4))
  
  # subset to those that have no imaginary part, then remove imaginary part (as it is 0)
  Ai_equilibriums <- Ai_equilibriums[Im(zapsmall(Ai_equilibriums)) == 0]
  Ai_equilibriums <- Re(Ai_equilibriums)
  
  # for each Ai equilibrium, calculate value of As (equation calculated manually - see lab book)
  As_equilibriums <- i*(alpha*Ai_equilibriums + theta*Ai_equilibriums*(1 - Fi) - lamda*Ai_equilibriums*(1-Ai_equilibriums/K)) / 
    (theta*Fi*(1-i))
  
  # find equilibrium where both As and Ai are positive i.e. biologically possible equilibrium
  Ai <- Ai_equilibriums[which(Ai_equilibriums > 0 & As_equilibriums > 0)]
  As <- As_equilibriums[which(Ai_equilibriums > 0 & As_equilibriums > 0)]
  
  # if multiple positive equilibriums and return_all = F, arbitrarily select the first one. Give warning
  # if return_all = T, then return all positive equilibriums
  if (length(Ai) > 1 | length(As) > 1) {
    if (return_all == F) {
      Ai <- Ai[1]
      As <- As[1]
    }
    warning("Multiple positive As/Ai equilibriums")
    
  }
  return(c(As, Ai))
}

donnelly_simp_vdynamic_ode <- function(times, states, parms) {
  
  ### ODE function for deterministic Donnelly model, with As and Ai (aphids per healthy and infected plants 
  ### respectively) assumed to be always at equilibrium. Therefore calculates equilibrium As and Ai values 
  ### for each value of i (proportion infected plants).
  
  # STATES
  i <- states[["i"]]
  
  # PARAMETERS
  gamma <- parms[["gamma"]]
  alpha <- parms[["alpha"]]
  theta <- parms[["theta"]]
  lamda <- parms[["lamda"]]
  K <- parms[["K"]]
  H <- parms[["H"]]
  p <- parms[["p"]]
  v <- parms[["v"]]
  e <- parms[["e"]]
  w <- parms[["w_don"]]
  Pacq <- parms[["Pacq"]]
  Pinoc <- parms[["Pinoc"]]
  
  # DEFINE PARAMETERS NEEDED FOR STATE EQUATIONS
  
  # q = prob of surviving flight between plants (i.e. doesn't emigrate/die)
  q <- 1 - p  
  
  # i_hat = weighted frequency of infected plants, accounting for attraction towards infected plants
  # i.e. adapted version of i / (1 - i + i) = i / 1, to include v
  i_hat <- v*i / ((1-i) + v*i) 
  
  # expected number of transmissions per dispersal, 
  # derived from analysis of markov chain (see Donnelly 2019, Appendix S1)
  xi <- (q^2*Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(p + q*w*(1 - i_hat*(1 - e)))
  
  # calculate equilibrium level of As and Ai for current value of i
  # = per plant aphid density on susceptible (As)  and infected (Ai) plants
  As_Ai <- calc_A_equilibrium(parms, states, i_hat)
  As <- As_Ai[1]
  Ai <- As_Ai[2]
  
  # A = total number of aphids i.e. aphids per plant type * number of that plant type
  # i.e. S*As + I*Ai where S is number of susceptible plants and I is number of infected plants
  A <- As*(H - i*H) + Ai*(i*H)  
  
  # STATE EQUATION
  # infected plants - change in incidence of infected plants
  # = rate of aphid dispersal * mean number of transmissions - plant death
  di <- (theta*A/H)*xi - gamma*i
  
  return(list(c(di)))
}


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