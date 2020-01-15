#'
#' Author of the code: Julien Chevallier (julien.chevallier1 at univ-grenoble-alpes.fr)
#'


library(expm)

simul_Hawkes_Erlang<-function(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec = rep(1,nb_pop), eta_vec, X_init, nb_points = 10000, bound_method = "exact", stop_time)
{
  # The function simulates the multiclass Hawkes process model of Ditlevsen and Locherbach (2016)
  # "nb_pop" is the number of sub-populations
  # "nb_neuron" is a vector of length nb_pop. Coordinate k is the number of neurons in population k.
  # "intensity_function" is a list of length nb_pop. Each coordinate is the function that maps the potential to the intensity. We assume that the intensity functions are non decreasing
  # c, nu and eta refer to the paper Ditlevsen and Locherbach (2016). CAREFUL : eta corresponds to eta+1 in the paper.
  # "c_vec" is a vector of length nb_pop. Coordinate k is the multiplicative constant for the interaction function from population k+1 to population k. Can be positive or negative
  # "nu_vec" is a vector of length nb_pop. Coordinate k is the time constant in the exponential for the interaction function from population k+1 to population k.
  # "eta_vec" is a vector of length nb_pop. Coordinate k is the memory order for the interaction function from population k+1 to population k.
  # "X_init" is a vector of length eta_1 x eta_2 x ... x eta_n. It gives the initial condition for the potential and memory variables.
  # The simulation stops when the dominating processes have produced "nb_points" spike.
  # Output : a list of length 3
  # Let N denote the total number of spikes simulated
  # $spike_train is a vector of length N with all the spiking times
  # $type is a vector of length N with the indices of the spiking neurons
  # $intensity is a matrix "nb_pop" x N with the intensity values for each population at each spiking time
  t = 0
  X = X_init
  A = linear_ODE_matrix(nu_vec, eta_vec)
  lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, A, bound_method)
  points = rep(0, nb_points)
  type = rep(0, nb_points)
  intensity = matrix(0, nrow = nb_pop, ncol = nb_points)
  X_in_time = matrix(0, nrow = length(X_init), ncol = nb_points)
  for (i in 1:nb_points)
  {
    dominant_ISI = rexp(nb_pop, nb_neuron*lambda_bound) # For each population we have the next possible spike
    spiking_pop = which.min(dominant_ISI)
    possible_ISI = dominant_ISI[spiking_pop]
    if ( possible_ISI==0 ){stop("blow-up")} # To get out of a possible blow-up
    t = t + possible_ISI
    new_X = expAtv(A, X, possible_ISI)$eAtv # Update of X via the flow of the ODE
    # new_X = X + possible_ISI*(A%*%X)        Take only the linear term in the exponential (this is not so bad because our time steps between two spikes are very short)
    new_lambda = intensity(new_X, nb_pop, eta_vec, intensity_function)
    if ( new_lambda[spiking_pop]==Inf ){stop("blow-up")} # To get out of a possible blow-up
    if ( any(new_lambda > lambda_bound) ){print("The upperbound is not good")}
    test1 = rbinom(1,1,prob = new_lambda[spiking_pop]/lambda_bound[spiking_pop]) # If test = 1 then one neuron of the population emits a spike
    if (test1)
    {
      neuron_number = sample(1:(nb_neuron[spiking_pop]), 1) # the spiking neuron is uniform over the population
      if (spiking_pop == 1) {influenced_pop = nb_pop} # Get the index of the population influenced by the spike
      else {influenced_pop = spiking_pop - 1}
      X = new_X
      X[ sum(eta_vec[0:influenced_pop]) ] = X[ sum(eta_vec[0:influenced_pop]) ] + c_vec[influenced_pop]/nb_neuron[spiking_pop]
      lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, A, bound_method)
      # Store the values
      points[i] = t
      type[i] = sum(nb_neuron[0:(spiking_pop-1)]) + neuron_number # the index of the spiking neuron
      intensity[,i] = new_lambda
      X_in_time[,i] = X
    }
    else
    {
      X = new_X
      lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, A, bound_method)
    }
    if (t>stop_time) {break}
  }
  good_index=which(points!=0)
  percent_reject = 100*(nb_points - length(good_index))/nb_points
  print(paste("The algorithm rejects",percent_reject,"% of every simulated points."))
  return(list(spike_train = points[good_index], type = type[good_index], intensity = intensity[,good_index], X = X_in_time[,good_index]))
}

simul_Hawkes_Erlang_Anna <-function(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec = rep(1,nb_pop), eta_vec, X_init, stop_time, nb_points = 10000, bound_method = "exact")
{
  # The function simulates the multiclass Hawkes process model of Ditlevsen and Locherbach (2016)
  # "nb_pop" is the number of sub-populations
  # "nb_neuron" is a vector of length nb_pop. Coordinate k is the number of neurons in population k.
  # "intensity_function" is a list of length nb_pop. Each coordinate is the function that maps the potential to the intensity. We assume that the intensity functions are non decreasing
  # c, nu and eta refer to the paper Ditlevsen and Locherbach (2016). CAREFUL : eta corresponds to eta+1 in the paper.
  # "c_vec" is a vector of length nb_pop. Coordinate k is the multiplicative constant for the interaction function from population k+1 to population k. Can be positive or negative
  # "nu_vec" is a vector of length nb_pop. Coordinate k is the time constant in the exponential for the interaction function from population k+1 to population k.
  # "eta_vec" is a vector of length nb_pop. Coordinate k is the memory order for the interaction function from population k+1 to population k.
  # "X_init" is a vector of length eta_1 x eta_2 x ... x eta_n. It gives the initial condition for the potential and memory variables.
  # The simulation stops when the dominating processes have produced "nb_points" spike.
  # Output : a list of length 3
  # Let N denote the total number of spikes simulated
  # $spike_train is a vector of length N with all the spiking times
  # $type is a vector of length N with the indices of the spiking neurons
  # $intensity is a matrix "nb_pop" x N with the intensity values for each population at each spiking time
  t = 0
  i = 0
  j = 0  
  type = 0
  points = 0
  X = X_init
  A = linear_ODE_matrix(nu_vec, eta_vec)
  lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, A, bound_method)
  #points = rep(0, nb_points)
  #type = rep(0, nb_points)
  intensity = rep(0, nb_pop) # matrix(0, nrow = nb_pop, ncol = nb_points)
  X_in_time = rep(0, length(X_init)) # matrix(0, nrow = length(X_init), ncol = nb_points
  while (t<stop_time)
  {
    i = i+1
    dominant_ISI = rexp(nb_pop, nb_neuron*lambda_bound) # For each population we have the next possible spike
    spiking_pop = which.min(dominant_ISI)
    possible_ISI = dominant_ISI[spiking_pop]
    if ( possible_ISI==0 ){stop("blow-up")} # To get out of a possible blow-up
    t = t + possible_ISI
    new_X = expAtv(A, X, possible_ISI)$eAtv # Update of X via the flow of the ODE
    # new_X = X + possible_ISI*(A%*%X)        Take only the linear term in the exponential (this is not so bad because our time steps between two spikes are very short)
    new_lambda = intensity(new_X, nb_pop, eta_vec, intensity_function)
    if ( new_lambda[spiking_pop]==Inf ){stop("blow-up")} # To get out of a possible blow-up
    if ( any(new_lambda > lambda_bound) ){print("The upperbound is not good")}
    test1 = rbinom(1,1,prob = new_lambda[spiking_pop]/lambda_bound[spiking_pop]) # If test = 1 then one neuron of the population emits a spike
    if (test1)
    {
      j = j+1
      neuron_number = sample(1:(nb_neuron[spiking_pop]), 1) # the spiking neuron is uniform over the population
      if (spiking_pop == 1) {influenced_pop = nb_pop} # Get the index of the population influenced by the spike
      else {influenced_pop = spiking_pop - 1}
      X = new_X
      X[ sum(eta_vec[0:influenced_pop]) ] = X[ sum(eta_vec[0:influenced_pop]) ] + c_vec[influenced_pop]/nb_neuron[spiking_pop]
      lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, A, bound_method)
      # Store the values
      points = c(points, t)
      type = c(type, sum(nb_neuron[0:(spiking_pop-1)]) + neuron_number) # the index of the spiking neuron
      # intensity[,i] = new_lambda
      # X_in_time[,i] = X
      intensity = rbind(intensity, new_lambda)
      X_in_time = rbind(X_in_time,X)
    }
    else
    {
      X = new_X
      lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, A, bound_method)
    }
    if (i%%1e4==0) {
    percent = (1-j/i)*100
    print(paste("We have passed i=",i,"! Current time is t=",t))
    print(paste("We had ",percent ,"% percent of rejections on the total interval!"))
    }
  }
  # good_index=which(points!=0)
  # percent_reject = 100*(nb_points - length(good_index))/nb_points
  # print(paste("The algorithm rejects",percent_reject,"% of every simulated points."))
  # return(list(spike_train = points[good_index], type = type[good_index], intensity = intensity[,good_index], X = X_in_time[,good_index]))

  return(list(spike_train = points, type = type, intensity = intensity, X = X_in_time))
}

intensity <- function(X, nb_pop, eta_vec, intensity_function)
  # Output : a vector of length "nb_pop". Coordinate k is the upper-bound of the intensity in population k in case of no spike until +Inf.
  # If "potential" is positive then the bound is f("potential")
  # If "potential" is negative then the bound is f(0)
{
  value = rep(0, nb_pop)
  for (k in 1:nb_pop)
  {
    potential = X[sum(eta_vec[0:(k-1)]) + 1]
    value[k] = intensity_function[[k]](potential)
  }
  return(value)
}


intensity_bound <- function(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, A, bound_method)
  # Output : a vector of length "nb_pop". Coordinate k is the upper-bound of the intensity in population k in case of no spike.
  # There are three methods :
  # "old" is the first one which is FALSE
  # "exact" gives a good upperbound. No chance to get an error in the simulation algorithm.
  # "approx" gives an approximative upperbound for a short time interval. Small chance to get an error in the simulation algorithm.
  # CAREFUL : An error occurs if the upperbound is wrong over the set of spiking times.
  # Even if we get no error in the simulation algorithm, there is small chance that the upperboung is wrong over the whole real line (but this we cannot check)
{
  if (bound_method == "old")
  {
    bound = rep(0, nb_pop)
    for (k in 1:nb_pop)
    {
      potential = X[sum(eta_vec[0:(k-1)]) + 1]
      bound[k] = intensity_function[[k]]( max(potential, 0) )
    }
    return(bound) 
  }
  if (bound_method == "exact")
  {
    bound = rep(0, nb_pop)
    for (k in 1:nb_pop)
    {
      i_min = sum(eta_vec[1:(k-1)]) + 1
      i_max = sum(eta_vec[1:k])
      bound[k] = intensity_function[[k]]( X_bound(X[i_min:i_max], nu_vec[k]) )
    }
    return(bound)
  }
  if (bound_method == "approx")
  {
    bound = rep(0, nb_pop)
    apriori = intensity(X, nb_pop, eta_vec, intensity_function)
    total_firing_rate = sum(nb_neuron*apriori)
    expected_time = 1/total_firing_rate
    new_X = expAtv(A, X, expected_time)$eAtv
    for (k in 1:nb_pop)
    {
      new_potential = new_X[sum(eta_vec[0:(k-1)]) + 1]
      potential_bound = (1 + 5*sign(new_potential)/sqrt(sum(nb_neuron)))*new_potential 
      # the coeficient 5 is chosen to prevent the fluctuations. A way to scale with the number of points could be to replace 5 by a normal distribution quantile.
      bound[k] = intensity_function[[k]]( max(X[sum(eta_vec[0:(k-1)]) + 1], potential_bound) )
    }
    return(bound)
  }
}


X_bound <- function(X, nu)
{
  # Output : a real number. It is a bound for max_{t\in R} e^{-tA}X where A is the Jordan matrix with -nu on the diagonal.
  eta = length(X)
  power = (1:eta) - 1
  Y = X/(nu^power)
  return( max(Y, 0) )
}



linear_ODE_matrix <- function(nu_vec, eta_vec)
  # Constructs the matrix of the linear ODE from nu_vec and eta_vec
{
  dimension = sum(eta_vec)
  A = matrix(0, nrow = dimension, ncol = dimension)
  nb_pop = length(nu_vec)
  D = c()
  for (k in 1:nb_pop)
  {
    D = c(D, rep(-nu_vec[k], eta_vec[k]))
  }
  diag(A) <- D
  upper = c()
  for (k in 1:(nb_pop-1))
  {
    upper = c(upper, rep(1, eta_vec[k]-1), 0)
  }
  upper = c(upper, rep(1, eta_vec[nb_pop]-1))
  indx <- 1:(dimension-1)
  A[cbind(indx,indx+1)] <- upper
  return(A)
}


### USAGE ###
## To reproduce Figure 1 in Eva and Susane's paper
### USAGE ###
## To reproduce Figure 1 in Eva and Susane's paper
nb_pop = 2
nb_neuron = rep(100, 2)
f1 <- function(x)
{
  if (x<log(20)) {return(10*exp(x))}
  else {return( 400/(1+400*exp(-2*x)) )}
}
f2 <- function(x)
{
  if (x<log(20)) {return(exp(x))}
  else {return( 40/(1+400*exp(-2*x)) )}
}
intensity_function = list(f1,f2)
c_vec = c(-1, 1)
nu_vec = rep(1, 2)
eta_vec = c(3, 2) + 1   # In comparison with the paper of Eva and Susanne, we take eta+1 instead of eta
X_init = rep(0, sum(eta_vec))

# test = simul_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec, X_init, nb_points = 100000, bound_method = "exact")
# 
# plot(test$spike_train, test$intensity[1,], type='l', ylim = c(0,40), col='blue')
# lines(test$spike_train, test$intensity[2,],col='black')
# 
# plot(test$spike_train, test$X[1,], type = 'n', ylim = c(-40,10))
# for (i in 2:4)
# {
#   lines(test$spike_train, test$X[i,], col='grey')
# }
# lines(test$spike_train, test$X[1,])
# for (i in 6:7)
# {
#   lines(test$spike_train, test$X[i,], col='grey')
# }
# lines(test$spike_train, test$X[5,])
