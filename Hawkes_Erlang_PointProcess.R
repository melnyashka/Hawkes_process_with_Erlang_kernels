#'
#' Author of the code: Julien Chevallier (julien.chevallier1 at univ-grenoble-alpes.fr)
#'


library(expm)

simul_Hawkes_Erlang<-function(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec = rep(1,nb_pop), eta_vec, X_init, percent_small_ISI = 0.1, nb_points = 10000, bound_method = "exact")
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
  # "time_grid_scale" is a positive real number. The a priori time grid used for the simulation is multiplied by this coeficient.
  # The simulation stops when the dominating processes have produced "nb_points" spike.
  # Output : a list of length 3
  # Let N denote the total number of spikes simulated
  # $spike_train is a vector of length N with all the spiking times
  # $type is a vector of length N with the indices of the spiking neurons
  # $intensity is a matrix "nb_pop" x N with the intensity values for each population at each spiking time
  
  t = 0
  X = X_init
  A = linear_ODE_matrix(nu_vec, eta_vec)
  # Initialize the intensity bound
  lambda_bound = intensity(X, nb_pop, eta_vec, intensity_function)
  apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
  lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
  # Initialize the storage
  points = rep(0, nb_points)
  type = rep(0, nb_points)
  intensity = matrix(0, nrow = nb_pop, ncol = nb_points)
  X_in_time = matrix(0, nrow = length(X_init), ncol = nb_points)
  for (i in 1:nb_points)
  {
    dominant_ISI = rexp(nb_pop, nb_neuron*lambda_bound) # For each population we have the next possible spike
    if (min(dominant_ISI) > apriori_ISI)
    {
      t = t + apriori_ISI
      new_X = expAtv(A, X, apriori_ISI)$eAtv # Update of X via the flow of the ODE
      X = new_X
      # Store the values
      type[i] = Inf # It means that there is no point in that time window because it was too short
      # Update lambda_bound
      apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
      lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
    }
    else
    {
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
        # Store the values
        points[i] = t
        type[i] = sum(nb_neuron[0:(spiking_pop-1)]) + neuron_number # the index of the spiking neuron
        intensity[,i] = new_lambda
        X_in_time[,i] = X
        # Update lambda_bound
        apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
        lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
      }
      else
      {
        X = new_X
        # No value to store
        # Update lambda_bound
        apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
        lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
      }
    }
  }
  # good_index = which(type!=0)
  reject_index = which(type==0)
  smalltime_index = which(type==Inf)
  good_index = which( (type > 0) & (type < (sum(nb_neuron)+1)) )
  # percent_reject = 100*(nb_points - length(good_index))/nb_points
  percent_reject = 100*(length(reject_index))/nb_points
  percent_smalltime = 100*length(smalltime_index)/nb_points
  print(paste("The algorithm rejects",percent_reject,"% of every simulated points."))
  print(paste("The algorithm takes ISI too small",percent_smalltime,"% of the time."))
  return(list(spike_train = points[good_index], type = type[good_index], intensity = intensity[,good_index], X = X_in_time[,good_index]))
}

simul_Hawkes_Erlang_Anna <-function(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec = rep(1,nb_pop), eta_vec, X_init, percent_small_ISI = 0.1, stop_time, bound_method = "polyroot")
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
  
  # Initial bound
  lambda_bound = intensity(X, nb_pop, eta_vec, intensity_function)
  apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
  lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
  
  # Initialize the vectors
  intensity = rep(0, nb_pop) # matrix(0, nrow = nb_pop, ncol = nb_points)
  X_in_time = rep(0, length(X_init)) # matrix(0, nrow = length(X_init), ncol = nb_points
  while (t<stop_time)
  {
    i = i+1
    dominant_ISI = rexp(nb_pop, nb_neuron*lambda_bound) # For each population we have the next possible spike
    if (min(dominant_ISI)>apriori_ISI){
      t = t+apriori_ISI
      new_X = expAtv(A, X, apriori_ISI)$eAtv
      X = new_X
      type = c(type, Inf)
      apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
      lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
    } else {
      spiking_pop = which.min(dominant_ISI)
      possible_ISI = dominant_ISI[spiking_pop]
      if ( possible_ISI==0 ){stop("blow-up")} # To get out of a possible blow-up
      t = t + possible_ISI
      new_X = expAtv(A, X, possible_ISI)$eAtv # Update of X via the flow of the ODE
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
        # lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
        # Store the values
        points = c(points, t)
        type = c(type, sum(nb_neuron[0:(spiking_pop-1)]) + neuron_number) # the index of the spiking neuron
        # intensity[,i] = new_lambda
        # X_in_time[,i] = X
        intensity = rbind(intensity, new_lambda)
        X_in_time = rbind(X_in_time,X)
        apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
        lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
      }
      else
      {
        X = new_X
        apriori_ISI = -log(percent_small_ISI)/sum(nb_neuron*lambda_bound)
        lambda_bound = intensity_bound(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
      }
      
      }
    
    if (i%%1e4==0) {
    percent = (1-j/i)*100
    print(paste("We have passed i=",i,"! Current time is t=",t))
    print(paste("We had ",percent ,"% percent of rejections on the total interval!"))
    }
  }
  
  print(paste("We had ",round((sum(type == Inf)/i)*100, digits = 3),"% percents of too small time steps!"))
  print(paste("We had ",round((1-j/i-sum(type==Inf)/i)*100, digits = 3),"% percent of rejections (not going in condition (2))"))
  # write.table(cbind(X,points), file = name, sep = ",", append = FALSE, col.names = FALSE, row.names = FALSE)
  return(list(spike_train = points, type = type[type!=Inf], intensity = intensity, X = X_in_time))
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


intensity_bound <- function(X, nb_pop, eta_vec, nu_vec, intensity_function, nb_neuron, apriori_ISI, bound_method)
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
      bound[k] = intensity_function[[k]]( X_bound(X[i_min:i_max], nu_vec[k], apriori_ISI) )
    }
    return(bound)
  }
  if (bound_method == "polyroot")
  {
    bound = rep(0, nb_pop)
    apriori = intensity(X, nb_pop, eta_vec, intensity_function)
    total_firing_rate = sum(nb_neuron*apriori)
    expected_time = 1/total_firing_rate
    for (k in 1:nb_pop)
    {
      i_min = sum(eta_vec[0:(k-1)]) + 1
      i_max = sum(eta_vec[1:k])
      # bound[k] = intensity_function[[k]]( X_bound_root(X[i_min:i_max], nu_vec[k], Tmax = 10*expected_time) ) + 10^(-3)
      # CHANGE BELOW
      bound[k] = (1.0001)*intensity_function[[k]]( X_bound_root(X[i_min:i_max], nu_vec[k], Tmax = apriori_ISI) )
    }
    return(bound)
  }

}


X_bound <- function(X, nu, Tmax)
{
  # Output : a real number. It is a bound for max_{t \leq Tmax} (e^{-tA}X)_1 where A is the Jordan matrix with -nu on the diagonal.
  eta = length(X)
  power = (1:eta) - 1
  Y = X/(nu^power)
  alltime_bound = max(Y, 0)
  smalltime_bound = max(X)*exp(-(nu-1)*Tmax)
  bound_type = which.min(c(alltime_bound, smalltime_bound))
  # print(paste("type ",bound_type))
  return( min(alltime_bound, smalltime_bound) )
}

X_bound_root <- function(X, nu, Tmax)
  # Output : a real number. It is a bound for max_{t \leq Tmax} (e^{-tA}X)_1 where A is the Jordan matrix with -nu on the diagonal.
  # Takes the maximal value over all the critical points computed as roots of the derivative
{
  n = length(X)-1
  fact = factorial(0:n)
  roots = Re(polyroot((-nu*X+c(X[-1],0))/fact)) # Polynomial coeficients of the time derivative
  index = (roots<Tmax)&(roots>0)
  # if (any(index)){print(paste("Obtained time: ", unique(roots[index]), "Tmax = ", Tmax))}
  interest_times = c(unique(roots[index]), Tmax)
  k = length(interest_times)
  values = rep(0,k)
  for (i in 1:k)
  {
    values[i] = exp(-nu*interest_times[i])*sum(X*interest_times[i]^(0:n)/fact)
  }
  # if (any(index)){print(paste("First value:", X[1],"Computed value:", values))}
  return( max(X[1], values) )
}

time_evaluation <- function(time, n)
{
  sum(time^(0:n))
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
