source("Hawkes_Erlang_new.R")
source("Hawkes_Erlang_EM_SDE.R")

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
nb_neur = c(100, 100) # number of neurons in population, vector of integers of length n_pop
eta_vec = c(3,2) # c(3,2) # number of memory variables, vector of integers of length n_pop
nu_vec = c(1,1) # auxilliary constants
c_vec = c(-1,1) # rates of population
K = c(10, 1) # constants for the rate functions
delta_gen = 1e-1
N_gen = round(70/delta_gen)
  

#------------------------------------------
#    Point process 
#------------------------------------------
# test = simul_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec+1, X_init, nb_points = 1e6, bound_method = "exact", stop_time = 30)
test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec+1, X_init, stop_time = 1e4, nb_points = 1e6, bound_method = "exact")

X <- t(test$X)
# plot(test$spike_train, test$intensity[1,], type='l', ylim = c(0,40), col='blue')
# lines(test$spike_train, test$intensity[2,],col='black')

plot(test$spike_train, X[1,], type = 'n', ylim = c(-40,10))
for (i in c(2:4, 6:7))
{
  lines(test$spike_train, X[i,], col='grey')
}
lines(test$spike_train, X[1,])
lines(test$spike_train, X[5,])

plot(density(test$X[1,]), xlim = c(-10,10), ylim = c(0,0.7), ylab = TeX('Density of X'), type = "l", main = "N=2000")
lines(density(test$X[5,]))
for (i in c(2:4, 6:7))
{
  lines(density(test$X[i,]), col='grey')
}
#------------------------------------
# Code for splitting
#-------------------------------------
Z_split = hawkes_splitting(N = N_gen, delta = delta_gen, nb_pop, nb_neur, eta_vec, nu_vec, c_vec, K)
time_elapsed = c(1:N_gen)*delta_gen
build_plot(Z_split, time_elapsed, c(1,5))
lines(time_elapsed, Z_split[1,], col = "black")

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