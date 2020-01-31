library(futile.logger)
library(latex2exp)
setwd("University/neuroscience/Hawkes_with_Erlang_kernels/logs") # comment/uncomment/change the directory if necessary
date <- Sys.Date()          # Create a subdirectory with the current date
wd <- getwd()               # Save the name of the working directory 

# parameters of logging
path_to_logs <- file.path(wd,date)
dir.create(path_to_logs)
file <- paste(path_to_logs, "logfile", sep = "/")
# flog.appender(appender.file(file))
flog.appender(appender.tee(file)) # write both to console and to file 
flog.threshold(DEBUG)        # By default set threshold to INFO (because I can)
flog.debug("Debugging is on!") 

nb_pop = 2 # number of populations
nb_neur = c(10, 10) # number of neurons in population, vector of integers of length n_pop
eta_vec = c(3,2) # c(3,2) # number of memory variables, vector of integers of length n_pop
nu_vec = c(0.9,0.9) # auxilliary constants
c_vec = c(-1,1) # rates of population
K = c(10, 1) # constants for the rate functions
delta_gen = 1e-1
N_gen = round(1e2/delta_gen)

  
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
#------------------------------------------
#    Point process 
#------------------------------------------
# test = simul_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec+1, X_init, nb_points = 1e6, bound_method = "exact", stop_time = 30)
stop_time = 100
eta_vec = c(3,2) # c(3,2) # number of memory variables, vector of integers of length n_pop
nu_vec = c(1,1) # auxilliary constants
for (neur in c(10,50,100)){
  # name = paste(path_to_logs,"/X_thinning_",stop_time,"_neur",neur, ".csv", sep="")
  start_time <- Sys.time()
  test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = 0.001, stop_time, bound_method = "polyroot")
  end_time <- Sys.time()
  print(paste("Neur=", neur, " \n Time elapsed:",end_time-start_time,"s"))
  # write.table(cbind(test$X,test$spike_trains), file = name, sep = ",", append = FALSE, col.names = FALSE, row.names = FALSE)
}

###########################################################
# Plot the differences between different guys 
############################################################
test <- list()
X_split <- list()
stop_time = 300
neur <- c(10, 50, 100, 500)
delta_gen = 0.01
# First we generate some trajectories
for (j in length(neur)){
  test[[j]] = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur[j], 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = 0.01, stop_time, bound_method = "polyroot", name)
  X_split[[j]] = hawkes_splitting(N = as.integer(stop_time/delta_gen), delta = delta_gen, nb_pop, nb_neur=rep(neur[j], 2), eta_vec, nu_vec, c_vec, K)
}
# Save the trajectories
for (j in 1:length(neur)){
  name = paste(path_to_logs,"/X_thinning_",stop_time,"_neur",neur[j], ".csv", sep="")
  write.table(cbind(test[[j]]$X,test[[j]]$spike_trains), file = name, sep = ",", append = FALSE, col.names = FALSE, row.names = FALSE)
  name = paste(path_to_logs,"/X_splitting_",stop_time,delta_gen,"_neur",neur[j], ".csv", sep="")
  write.table(X_split[[j]], file = name, sep = ",", append = FALSE, col.names = FALSE, row.names = FALSE)
}

x_lim = c(15, 270)
# pdf("second_population_comparison.pdf", width = 10, height = 7) 
par(mfrow=c(length(neur),1), mai = c(0.3, 0.3, 0.1, 0.1))
for (i in 1:length(neur)){
  time_elapsed= seq(0, stop_time, delta_gen)
  if (i!=length(neur)){
    plot(test[[i]]$spike_train, test[[i]]$X[,5], type = "l", lwd = 2, xlim = x_lim, ylim = c(0,2), xlab = "", ylab = "", xaxt='n')
  } else { plot(test[[i]]$spike_train, test[[i]]$X[,5], type = "l", lwd = 2, xlim = x_lim, ylim = c(0,2), xlab = "", ylab = "")}
  # lines(test[[i]]$spike_train, test[[i]]$X[,1])
  lines(time_elapsed[1:length(X_split[[i]][,5])], X_split[[i]][,5], lty = 2, lwd = 2)
  legend(x_lim[2]-25,1.8, bg = "white", text.font = 20, legend = TeX(sprintf("$N_2=%d$",neur[i])))
}
# dev.off()
  
pdf("first_population_comparison.pdf", width = 10, height = 7)
par(mfrow=c(length(neur),1), mai = c(0.3, 0.3, 0.1, 0.1))
for (i in 1:length(neur)){
  time_elapsed = seq(0, stop_time, delta_gen)
  if (i!=length(neur)){
    plot(test[[i]]$spike_train, test[[i]]$X[,1], type = "l", xlim = x_lim, ylim = c(-5.5,-1.5), xlab = "", ylab = "",lwd = 2, xaxt='n')
  } else { plot(test[[i]]$spike_train, test[[i]]$X[,1], type = "l", xlim = x_lim, ylim = c(-5.5,-1.5), xlab = "",lwd = 2, ylab = "")}
  # lines(test[[i]]$spike_train, test[[i]]$X[,1])
  lines(time_elapsed[1:(length(time_elapsed)-1)], X_split[[i]][,1], lty = 2, lwd = 2)
  legend(x_lim[2]-25,-2, bg = "white", text.font = 20, legend = TeX(sprintf("$N_1=%d$",neur[i])))
}
dev.off()

# plot(test$spike_train, test$intensity[1,], type='l', ylim = c(0,40), col='blue')
# lines(test$spike_train, test$intensity[2,],col='black')
n = length(test$spike_train)
j = min(which(test$spike>100))
plot(test$spike_train[j:n], test$X[j:n,1], type = 'n', xlab = "t", ylab =TeX('$X_t$'), ylim = c(-10,10))
# for (i in c(2:(eta_vec[1]+1), (eta_vec[1]+3):(sum(eta_vec)+2)))
# {
#   lines(test$spike_train, test$X[,i], col='grey')
# }

lines(test$spike_train[j:n], test$X[j:n,1])
lines(test$spike_train[j:n], test$X[j:n,(eta_vec[1]+2)])

plot(test$X[1:stop_time,1], test$X[1:stop_time,5], type = "l")

plot(density(test$X[,1]), xlim = c(-10,10), ylim = c(0,0.7), ylab = TeX('Density of X'), type = "l", main = "t = 1000")
lines(density(test$X[,5]))
for (i in c(2:4, 6:7))
{
  lines(density(test$X[,i]), col='grey')
}

# Plot intensities

plot(test$spike_train, test$intensity[,1], type = "n", ylim = c(0,1), ylab = TeX('Density of X'), main = "")
lines(test$spike_train, test$intensity[,1])
lines(test$spike_train, test$intensity[,2], col = "blue")
#------------------------------------
# Plot densities
#------------------------------------
pprocess <- read.csv(file="X_thinning_1e+05_neur10.csv")
pprocess <- Thinning
n <- length(pprocess$X0)
# pdf("HawkesPDMP_10n_1e5.pdf") 
pdf("density_10n_1e5.pdf") 
plot(density(pprocess$X0), type = "l", xlim = c(-5,4), ylim = c(0,1), xlab = "x", ylab = TeX('$\\pi_x(x)$'), main = "Empirical density N = 20", lwd = 2)
lines(density(pprocess$X0.1), col = "grey", lwd = 2)
lines(density(pprocess$X0.2), col = "grey", lwd = 2)
lines(density(pprocess$X0.3), col = "grey", lwd = 2)
lines(density(pprocess$X0.4), col = "black", lwd = 2)
lines(density(pprocess$X0.5), col = "grey", lwd = 2)
lines(density(pprocess$X0.6), col = "grey", lwd = 2)
abline(v=mean(pprocess$X0), col = "red", lwd = 1)
abline(v = mean(pprocess$X0.4), col = "red", lwd = 1)
dev.off()

#------------------------------------
# Code for splitting
#-------------------------------------
stop_time = 1e5
delta_gen = 0.05
Z_split = hawkes_splitting(N = as.integer(stop_time/delta_gen), delta = delta_gen, nb_pop, nb_neur, eta_vec, nu_vec, c_vec, K)
time_elapsed = seq(0, stop_time, delta_gen)
write.csv(Z_split, "zsplit_1e5_n10.csv", row.names = F)
# build_plot(Z_split, time_elapsed, c(1,5))
# lines(time_elapsed, Z_split[1,], col = "black")
plot(time_elapsed, Z_split[,1], type = 'n', xlab = "t", ylab =TeX('$X_t$'), ylim = c(-40,10))
for (i in c(2:4, 6:7))
{
  lines(time_elapsed, Z_split[,i], col='grey')
}
lines(time_elapsed[100/delta_gen:stop_time/delta_gen], Z_split[100/delta_gen:N_gen,1], col = "red")
lines(time_elapsed[100/delta_gen:length(time_elapsed)], Z_split[100/delta_gen:N_gen,eta_vec[1]+2], col = "red")


pdf("HawkesSplit_10n_1e5.pdf") 
plot(density(Z_split[,1]), type = "l", xlim = c(-5,4), ylim = c(0,1), xlab = "x", ylab = TeX('$\\pi_x(x)$'), main = "", lwd = 2)
lines(density(Z_split[,1]), col = "black", lwd = 2, lty = 2)
lines(density(Z_split[,5]), col = "black", lwd = 2, lty = 2)
for (j in c(2:4, 6:7)){
  lines(density(Z_split[,j]), col = "grey", lwd = 2, lty = 2)
}
abline(v = mean(Z_split[,1]), col = "red", lty = 2)
abline(v = mean(Z_split[,5]), col = "red", lty = 2)
dev.off()

#-------------------------------------
# To plot intensities (EM)
#-------------------------------------
Z = hawkes_approximation(N = N_gen, delta = delta_gen, n_pop = nb_pop, n_neur = nb_neur, eta = eta_vec, nu = nu_vec, c_rate = c_vec, K = K)
lines(time_elapsed, sapply(Z[5,], rate_function, const = 1), type = "l", lty = 2)
lines(time_elapsed, sapply(Z[1,], rate_function, const = 10), col = "blue", lty=2)


# Now I will build cumulative intensities 
plot(test$spike_train, test$intensity[1,], type='l', ylim = c(0,40), col='blue')
lines(test$spike_train, test$intensity[2,],col='black')
lines(time_elapsed, sapply(Z[5,], rate_function, const = 1), type = "l", lty = 2)
lines(time_elapsed, sapply(Z[1,], rate_function, const = 10), col = "blue", lty=2)
lines(time_elapsed, sapply(Z_split[5,], rate_function, const = 1), type = "l", lty = 3)
lines(time_elapsed, sapply(Z_split[1,], rate_function, const = 10), col = "blue", lty=3)
legend(40, 40, legend=c("Point process", "EM scheme", "Splitting"),
       col=c("black", "black", "black"), lty=1:3, cex=1)


#==========================================================================
##### Computation of big amount of stuffs and measuring the time witht clusters
###########################################################################
library(tictoc)
library(parallel)
nneur = c(10,25,50,75,100)
time_discrete <- matrix(nrow = 100, ncol = 5)
time_cont <- matrix(nrow = 100, ncol = 5)
stop_time = 100
K = 100

int_f <- function(x){ min(0.1+max(0,x), 10)}
intensity_function <- list(int_f, int_f)

for (neur in nneur){
cl <- makeCluster(4, outfile = file)
clusterCall(cl = cl, function() {
  require("futile.logger")
  require("tictoc")
  require("expm")
  })
clusterExport(cl = cl, varlist = c("nb_pop", "nb_neur", "nneur", "neur", "intensity_function", "c_vec", "nu_vec", "eta_vec", "stop_time", "simul_Hawkes_Erlang_Anna", "time_discrete", "K", "linear_ODE_matrix","time_evaluation", "X_bound", "X_bound_root", "intensity_bound", "intensity"))
ind = which(nneur==neur)
time_discrete[,ind] = parSapply(cl = cl, 1:K, function(k){
    ind = which(nneur==neur)
    flog.debug("Launching the stuff... ")
    print(paste("neur=", neur, "k = ", k))
    tic("Launching the stuff")
    flog.debug("=====================================")
    flog.debug(paste("For neur=", neur, " and discrete bound"))
    test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = 1e-5, stop_time, bound_method = "polyroot")
    time = toc()
    flog.debug(paste("Elapsed time: ", time$toc - time$tic, "s"))
    return(round(time$toc - time$tic, digits = 5))
})
stopCluster(cl)
}
for (neur in nneur){
  cl <- makeCluster(4, outfile = file)
  clusterCall(cl = cl, function() {
    require("futile.logger")
    require("tictoc")
    require("expm")
  })
  clusterExport(cl = cl, varlist = c("nb_pop", "nb_neur", "nneur", "neur", "intensity_function", "c_vec", "nu_vec", "eta_vec", "stop_time", "simul_Hawkes_Erlang_Anna", "time_cont", "K", "linear_ODE_matrix","time_evaluation", "X_bound", "X_bound_root", "intensity_bound", "intensity"))
  ind = which(nneur==neur)
  time_cont[,ind] = parSapply(cl = cl, 1:K, function(k){
    ind = which(nneur==neur)
    flog.debug("Launching the stuff... ")
    tic("Launching the stuff")
    flog.debug("=====================================")
    flog.debug(paste("For neur=", neur, " and discrete bound"))
    test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = 1e-5, stop_time, bound_method = "exact")
    time = toc()
    flog.debug(paste("Elapsed time: ", time$toc - time$tic, "s"))
    return(round(time$toc - time$tic, digits = 5))
  })
  stopCluster(cl)
}
time_diff <- numeric(100)

for (k in 1:100){
  flog.debug("------------------------------------")
  tic("Launching the stuff")
  # flog.debug(paste("For neur=", neur, " and continuous bound"))
  # test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = 0.0001, stop_time, bound_method = "exact")
  Z <- hawkes_splitting(N = as.integer(stop_time/delta_gen), delta = delta_gen, nb_pop, nb_neur=rep(100, 2), eta_vec, nu_vec, c_vec, K)
  time = toc()
  time_diff[k] = round(time$toc - time$tic, digits = 5)
  flog.debug(paste("Elapsed time: ", time_diff[k], "s"))
  flog.debug("=====================================")
}


time_discrete <- read.csv("discrete_time_100.csv")
time_cont <- read.csv("cont_time_100.csv")
par(mai = c(0.5, 0.5, 0.5, 0.1))
neurons = c(10,25,50,75,100)*2
TD_Means = apply(time_discrete, 2, mean)
TD_SDs = apply(time_discrete, 2, sd)
TC_Means = apply(time_cont, 2, mean)
TC_SDs = apply(time_cont, 2, sd)

pdf("Elapsed_time_with_diff.pdf")
plot(neurons, TD_Means, lwd = 2, main = TeX("Exponential intensity"), type = "b", pch = 0, ylim = c(0,400), xlab = "N", ylab = "s")
# lines(neurons, TD_Means - TD_SDs, type = "l", lty = 2)
# lines(neurons, TD_Means + TD_SDs, type = "l", lty = 2)
lines(neurons, TC_Means, lwd = 2, type = "b", pch = 2, col = "blue")
# lines(neurons, TC_Means - TC_SDs, type = "l", lty = 2, col = "blue")
# lines(neurons, TC_Means + TC_SDs, type = "l", lty = 2, col = "blue")
abline(h = mean(time_diff), lwd = 2, col = "red", lty = 2)
legend(20, 400, bg = "white", lwd = 2, pch = c(0, 2, NA), lty = c(1,1,2), col = c("black", "blue", "red"), text.font = 20, legend = c(TeX('$\\tilde{f}^{\\Delta}(x)$'),TeX('$\\tilde{f}(x)$'),TeX("\\tilde{X}")))
dev.off()


for (Delta in c(0.001, 0.05, 0.5)){
  for (stop_time in c(10,100)){
    for (neur in c(10,100)){
      flog.debug("Launching the stuff... ")
      tic("Launching the stuff")
      flog.debug("=====================================")
      flog.debug(paste("We launch for Percent=",Delta," stop_time=", stop_time, "neur = ", neur ))
      test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = Delta, stop_time, bound_method = "polyroot")
      time = toc()
      flog.debug(paste("Elapsed time: ", time$toc - time$tic, "s"))
      flog.debug("=====================================")
    }
  }
}
       
  
  
#==========================================================================
##### Try to plot different bounds
###########################################################################
  
eta_vec <- c(3,2)
test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = 0.001, stop_time, bound_method = "polyroot")
  
t <- test$spike_train
X2 <- test$X[,5:7]
plot(t, X, type = "l")

lambda_1 <- function(Xk, nu){
  eta = length(Xk)
  power = (1:eta) - 1
  Y = Xk/(nu^power)
  alltime_bound = max(Y, 0)
  return(alltime_bound)
}
lambda_2 <- function(Xk, nu, Tmax){
  X <- Xk
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

l11 <- sapply(1:length(t), function(i) {lambda_1(X2[i,], nu_vec[2])} )
l12 <- sapply(1:length(t), function(i) {lambda_2(X2[i,], nu_vec[2], 1)} )
l21 <- sapply(1:length(t), function(i) {lambda_1(test$X[i,1:4], nu_vec[1])} )
l22 <- sapply(1:length(t), function(i) {lambda_2(test$X[i,1:4], nu_vec[1], 1)} )
n <- length(t)
ind = which(t>10)[1]

pdf("intensities.pdf", width = 10, height = 7)
par(mfrow=c(2,1), mai = c(0.3, 0.4, 0.1, 0.1))
plot(t[ind:n], test$intensity[ind:n,1], type = "n", ylim = c(0,4), lwd = 2, xlab = "", xaxt = "n", ylab = "")
lines(t[ind:n], f1(l21)[ind:n], col = "grey", lwd = 2)
lines(t[ind:n], f1(l22)[ind:n], col = "black", lwd = 2)
lines(t[ind:n], test$intensity[ind:n,1], col = "red", lty = 2, lwd = 3)
plot(t[ind:n], test$intensity[ind:n,2], type = "n", ylim = c(0,20), lwd = 2, xlab = "",  ylab = "")
lines(t[ind:n], f2(l11)[ind:n], col = "grey", lwd = 2)
lines(t[ind:n], f2(l12)[ind:n], col = "black", lwd = 2)
lines(t[ind:n], test$intensity[ind:n,2], col = "red", lty = 2, lwd = 3)

# legend(80, 20, bg = "white", lwd = 2, col = c("black", "blue", "red"), lty = c(1,1,2), text.font = 20, legend = c(TeX('$\\tilde{f}^{\\Delta}(x)$'),TeX('$\\tilde{f}(x)$', TeX('$f(x)$'))))

dev.off()
