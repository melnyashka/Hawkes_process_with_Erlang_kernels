#'
#' Author of the code: Anna Melnykova (anna.melnykova at univ-grenoble-alpes.fr)
#'
library(expm)

rate_function <- function(x, const=1){
  # here for simplicity we use the same family of functions
  if (x<log(20)){
    value <- const*exp(x)
  } else {
    value <- 40*const/(1+400*exp(-2*x))
  }
  return(value)
}

linear_ODE_matrix <- function(nu_vec, eta_vec)
  # Constructs the matrix of the linear ODE from nu_vec and eta_vec
  # we use this guy also for splitting
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

hawkes_splitting <- function(N, delta, nb_pop, nb_neur, eta_vec, nu_vec, c_vec, K){
  nb_total = nb_pop + sum(eta_vec)
  Z = matrix(nrow = N, ncol = nb_total)
  ind_rough = numeric(nb_pop) 
  A = linear_ODE_matrix(nu_vec, eta_vec+1)
  for (i in 1:nb_pop){ind_rough[i] = sum(eta_vec[1:i])+i}
  Z[1,] = rep(0,nb_total)
  p_vec = sqrt(nb_neur*sum(nb_neur)) # proportions of the neurons in each population
  for (i in 1:(N-1)){
    z = expAtv(A, Z[i,], delta/2)$eAtv # we are doing the first step of the approximation scheme
    for (ind in ind_rough){ # doing a step with the constant solution
      l = min(which(ind_rough >= ind)) # give a population
      if (ind+1 <= nb_total) {j = ind + 1} else {j = 1} # keep track on the next variable
      if (l+1 <= nb_pop) {l_ = l + 1} else {l_ = 1} # keep track on the index of population
      z[ind] = z[ind] + delta*c_vec[l]*rate_function(x = z[j], const = K[l_]) + c_vec[l]*sqrt(delta)*rnorm(1)*sqrt(rate_function(x = z[j], const = K[l_])/p_vec[l])
    } 
    Z[i+1,] = expAtv(A, z, delta/2)$eAtv #doing the last step
  }
  return(Z)
}

hawkes_approximation <- function(N, delta, n_pop, n_neur, eta, nu, c_rate, K){
  n_total = n_pop+sum(eta)
  Z = matrix(nrow = n_total, ncol = N) # slot for the process
  ind_rough = numeric(n_pop)
  # now we are checking the indexes of the rough variables
  for (i in 1:n_pop){ind_rough[i] = sum(eta[1:i])+i}
  Z[,1] = rep(0,n_total)
  for (n in 1:(N-1)){
    for (k in 1:n_total){
      i = min(which(ind_rough >= k))
      if (k %in% ind_rough) {
        if (k+1 <= n_total) {j = k + 1} else {j = 1}
        if (i+1 <= n_pop) {i_ = i + 1} else {i_ = 1}
        p = sqrt(n_neur[i_]*sum(n_neur))
        Z[k, n + 1] = Z[k,n] + delta*(-nu[i]*Z[k, n] + c_rate[i]*rate_function(x = Z[j,n], const = K[i_])) + c_rate[i]*sqrt(delta)*rnorm(1, mean = 0, sd = 1)*sqrt(rate_function(x = Z[j,n], const = K[i_])/p)
      } else { 
        Z[k, n + 1] = Z[k,n] + delta*(-nu[i]*Z[k, n] + Z[k + 1,n])
      }
    }
  }
  return(Z)
}



build_plot <- function(Z, time, ind_rough){
  # Build plot for s.a. of Hawkes process
  dimZ <- dim(Z)
  ylim_max <- 10 #max(Z)
  ylim_min <- -40 #min(Z)
  plot(time, Z[1,], type = "l", xlab = "", ylab = "", ylim = c(ylim_min, ylim_max), col = "grey")
  for (i in 2:dimZ[1]){
    if (i %in% ind_rough) {col_i = "black"} else {col_i = "grey"}
    lines(time, Z[i,], col = col_i)
  }
}

