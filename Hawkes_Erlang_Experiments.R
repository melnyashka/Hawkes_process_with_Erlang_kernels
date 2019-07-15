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

n_pop = 2 # number of populations
n_neur = c(100, 100) # number of neurons in population, vector of integers of length n_pop
eta = c(3,2) # number of memory variables, vector of integers of length n_pop
nu = c(1,1) # auxilliary constants
c_rate = c(-1,1) # rates of population
K = c(10, 1) # constants for the rate functions
delta_gen = 1e-4
N_gen = round(60/delta_gen)

Z = hawkes_approximation(N = N_gen, delta = delta_gen, n_pop = n_pop, n_neur = n_neur, eta = eta, nu = nu, c_rate = c_rate, K = K)
build_plot(Z, c(4,7))
