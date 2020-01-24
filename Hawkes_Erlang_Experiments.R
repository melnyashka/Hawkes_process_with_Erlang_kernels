source("Hawkes_Erlang_new.R")
source("Hawkes_Erlang_EM_SDE.R")

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
stop_time = 1e5
eta_vec = c(3,2) # c(3,2) # number of memory variables, vector of integers of length n_pop
nu_vec = c(0.8,0.8) # auxilliary constants
for (neur in c(10, 50, 100)){
  name = paste(path_to_logs,"/X_thinning_",stop_time,"_neur",neur, ".csv", sep="")
  test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(neur, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = rep(0, sum(eta_vec+1)), percent_small_ISI = 0.1, stop_time, bound_method = "polyroot", name)
  write.table(cbind(test$X,test$spike_trains), file = name, sep = ",", append = FALSE, col.names = FALSE, row.names = FALSE)
}

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
# pprocess <- read.csv("HawkesPDMP_10n_1e4.csv")
n <- length(pprocess$V1)
pdf("HawkesPDMP_10n_1e4.pdf") 
plot(density(pprocess$V1), type = "l", xlim = c(-5,4), ylim = c(0,1), xlab = "x", ylab = TeX('$\\pi_x(x)$'), main = "")
lines(density(pprocess$V2), col = "grey")
lines(density(pprocess$V3), col = "grey")
lines(density(pprocess$V4), col = "grey")
lines(density(pprocess$V5), col = "black")
lines(density(pprocess$V6), col = "grey")
lines(density(pprocess$V7), col = "grey")
dev.off()

#------------------------------------
# Code for splitting
#-------------------------------------
Z_split = hawkes_splitting(N = as.integer(stop_time/delta_gen), delta = delta_gen, nb_pop, nb_neur, eta_vec, nu_vec, c_vec, K)
time_elapsed = seq(0, stop_time, delta_gen)
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
plot(density(Z_split[,1]), type = "l", xlim = c(-5,4), ylim = c(0,1), xlab = "x", ylab = TeX('$\\pi_x(x)$'), main = "")
lines(density(Z_split[,5]), col = "black")
for (j in c(2:4, 6:7)){
  lines(density(Z_split[,j]), col = "grey")
}
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
       