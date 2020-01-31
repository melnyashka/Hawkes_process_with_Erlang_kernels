library(parallel)
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

# Setting up the parameters
nb_pop = 2 # number of populations
nb_neur = c(50, 50) # number of neurons in population, vector of integers of length n_pop
eta_vec = c(3,2) # c(3,2) # number of memory variables, vector of integers of length n_pop
nu_vec = c(1,1) # auxilliary constants
c_vec = c(-1,1) # rates of population
delta_gen = 1e-1

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

### Making the cluster

cl <- makeCluster(4, outfile = file)
clusterCall(cl = cl, function() {
  require("futile.logger")
  require("expm")
})
clusterExport(cl = cl, varlist = c("nb_pop", "nb_neur", "intensity_function", "c_vec", "nu_vec", "eta_vec", "simul_Hawkes_Erlang_Anna",  "linear_ODE_matrix","time_evaluation", "X_bound", "X_bound_root", "intensity_bound", "intensity"))
HawkesProcess = parLapply(cl = cl, 1:4, function(k){
  flog.debug("=====================================")
  flog.debug(paste("For neur=", neur, " and discrete bound"))
  X_val <- c(-4.8896667, -5.5195462, -5.5856640, -4.9379179,  1.1969114,  0.7250284,  0.2692858) # to put ourselves in stationary regime
  test = simul_Hawkes_Erlang_Anna(nb_pop, nb_neuron = rep(50, 2), intensity_function, c_vec, nu_vec, eta_vec+1, X_init = X_val, percent_small_ISI = 1e-3, stop_time=20, bound_method = "polyroot")
  return(test$X)
})
stopCluster(cl)