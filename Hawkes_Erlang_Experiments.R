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
eta_vec = c(3,2) # number of memory variables, vector of integers of length n_pop
nu_vec = c(1,1) # auxilliary constants
c_vec = c(-1,1) # rates of population
K = c(10, 1) # constants for the rate functions
delta_gen = 1e-2
N_gen = round(70/delta_gen)
  
Z = hawkes_approximation(N = N_gen, delta = delta_gen, n_pop = nb_pop, n_neur = nb_neur, eta = eta_vec, nu = nu_vec, c_rate = c_vec, K = K)


Z = hawkes_splitting(N = N_gen, delta = delta_gen, nb_pop, nb_neur, eta_vec, nu_vec, c_vec, K)
time_elapsed = c(1:N_gen)*delta_gen
build_plot(Z, time_elapsed, c(1,5))
lines(time_elapsed, Z[1,], col = "black")
# To plot intensities

plot(time_elapsed, sapply(Z[5,], rate_function, const = 1), type = "l", lty = 1)
lines(time_elapsed, sapply(Z[1,], rate_function, const = 10), col = "blue", lty=1)

