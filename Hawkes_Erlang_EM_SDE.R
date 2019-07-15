rate_function <- function(x, const=1){
  # here for simplicity we use the same family of functions
  if (x<log(20)){
    value <- const*exp(x)
  } else {
    value <- 40*const/(1+400*exp(-2*x))
  }
  return(value)
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
        p = n_neur[i_] #/sum(n_neur)
        Z[k, n + 1] = Z[k,n] + delta*(-nu[i]*Z[k, n] + c_rate[i]*rate_function(x = Z[j,n], const = K[i_])) + c_rate[i]*sqrt(delta)*rnorm(1, mean = 0, sd = 1)*sqrt(rate_function(x = Z[j,n], const = K[i_])/p)
      } else { 
        Z[k, n + 1] = Z[k,n] + delta*(-nu[i]*Z[k, n] + Z[k + 1,n])
      }
    }
  }
  return(Z)
}

build_plot <- function(Z, ind_rough){
  # Build plot for s.a. of Hawkes process
  dimZ <- dim(Z)
  ylim_max <- 10 #max(Z)
  ylim_min <- -40 #min(Z)
  plot(Z[1,], type = "l", xlab = "", ylab = "", ylim = c(ylim_min, ylim_max), col = "grey")
  for (i in 2:dimZ[1]){
    if (i %in% ind_rough) {col_i = "black"} else {col_i = "grey"}
    lines(Z[i,], col = col_i)
  }
}

