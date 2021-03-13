

# ---------------------------------------
# examples of various methods of parameter estimation  -- this is for testing


# loadfunctions("adapt")

# remotes::install_github( "jae0/adapt" )

require(adapt)

require(cmdstanr )
require(posterior )
require(bayesplot)

# parameter estimation via STAN

options(mc.cores = parallel::detectCores())



# ---------------------------------------
# 1. continuous ODE form via STAN/MCMC

i0 = 25/300
Npop = 300   # total population size (at start of simulation)
times = 0:99
Npreds = length(times)
params = list(beta=0.6, gamma=0.1)
inits = c(1-i0, i0, 0) # S, I, R, proportions

sim = simulate_data( selection="sir", Npop=Npop, inits=inits, times=times, params=params, plotdata=TRUE )  # Npop=1 forces results to be in proportions

# "partial" record dimensions
Nobs = 25  # length of observations time series
subsample = sort( sample.int( length(times), Nobs, replace=FALSE ) )
sim = sim[subsample,]

# data to pass to STAN
stan_data = list(
  Npop = Npop,
  Nobs = Nobs,
  Npreds = Npreds,   # here, only total number of predictions for output
  Sobs = sim$Sobs,
  Iobs = sim$Iobs,
  Robs = sim$Robs,
  time = as.integer(sim$time),
  time_pred = c(0:max(sim$time, times) ),
  t0 = -0.01
)



stancode = stan_initialize( stan_code=sir_stan_model_code( selection="continuous" ) )
stancode$compile()


fit = stancode$sample(        
  data = stan_data,
  #pars = c("y0", "params", "Ipred"),
  init = function() list( params=runif(2), S0=runif(1) ),
  chains = 3,
  warmup = 5000,
  iter = 6000
)

M = stan_extract( as_draws_df( fit$draws() ) )

# Check median estimates of parameters and initial conditions:
apply(M$params, 2, median)  # params
apply(M$y0, 2, median)[1:2]

plot_model_fit( stan_data=stan_data, M=M )





# ---------------------------------------
# 2. discrete form using recursive latent model via STAN/MCMC

i0 = 25/300
Npop = 300   # total population size (at start of simulation)
times = 0:99
params = list(beta=0.75, gamma=0.1)
inits = c(1-i0, i0, 0) # S, I, R, proportions

sim = simulate_data( selection="sir", Npop=Npop, inits=inits, times=times, params=params, plotdata=TRUE,  nsample=25 )  # Npop=1 forces results to be in proportions

# data to pass to STAN
stan_data = list(
    Npop = Npop,
    Nobs = length(sim$time) ,
    Npreds = 5,  # here this is additional time slices
    Sobs = sim$Sobs,
    Iobs = sim$Iobs,
    Robs = sim$Robs
)



stancode = stan_initialize( stan_code=sir_stan_model_code( selection="discrete_basic" ) )
stancode$compile()

fit = stancode$sample(        
  data = stan_data,
  adapt_delta = 0.95, 
  max_treedepth=15,
  chains = 3,
  warmup = 10000,
  iter = 15000
)

M = stan_extract( as_draws_df( fit$draws() ) )

# Check median estimates of parameters and initial conditions:
median(M$BETA)  # params
median(M$GAMMA)
median(M$K)

plot_model_fit( stan_data=stan_data, M=M )



# ---------------------------------------
# 3. discrete form using recursive latent model via STAN/MCMC and time-variable encounter rates



i0 = 25/300
Npop = 300   # total population size (at start of simulation)
times = 0:99
Npreds = length(times)
params = list(beta=0.75, gamma=0.1)
inits = c(1-i0, i0, 0) # S, I, R, proportions

sim = simulate_data( selection="sir", Npop=Npop, inits=inits, times=times, params=params, plotdata=TRUE, nsample=25 )  # Npop=1 forces results to be in proportions


# data to pass to STAN
stan_data = list(
  Npop = Npop,
  Nobs =  length(sim$time) ,
  Npreds = 5,  # here additional time slices
  Sobs = sim$Sobs,
  Iobs = sim$Iobs,
  Robs = sim$Robs,
  BETA_prior = 0.5,
  GAMMA_prior = 0.5,
  ER_prior = 0.5
)


stancode = stan_initialize( stan_code=sir_stan_model_code( selection="discrete_variable_encounter_rate" ) )
stancode$compile()

fit = stancode$sample(        
  data = stan_data,
   adapt_delta = 0.9, 
   max_treedepth=14 ,
#  pars = c("y0", "params", "Ipred"),
# init = function() list( params=runif(2), S0=runif(1) ),
  chains = 3,
  warmup = 1000,
  iter = 1500
)


M = stan_extract( as_draws_df( fit$draws() ) )

# Check median estimates of parameters and initial conditions:
plot( apply(M$K, 2, median) ) # params

plot( apply(M$I, 2, median), type="l" )
points( Iobs~ time, sim, col="red")

plot_model_fit( stan_data=stan_data, M=M )



# ---------------------------------------
# 3.5. stochastic discrete form via STAN/MCMC "discrete_binomial_process"

# set.seed(42)

i0 = 25/300

Npop = 300
times = 0:99
params = list(prob_infection=0.001, time_shedding=7)
# // prob_infection = probability of an individual infecting another in 1 unit of time (1 day)
# // time_shedding = 1/ probability of transition to recovered state (~ duration is units of time) .. ie., simple geometric  .. 14 days
inits = c(1-i0, i0, 0) # S, I, R, proportions

sim = simulate_data( selection="sir_stochastic", Npop=Npop, inits=inits, times=times, params=params, plotdata=TRUE )
  # prob_infection =  transmission probability (per person per 'close contact', per day)
  # time_shedding  =  mean shedding duration for an individual

stan_data = list(
    Npop = Npop,
    Nobs = nrow(sim),
    Sobs = sim$Sobs,
    Iobs = sim$Iobs,
    Robs = sim$Robs
)

stancode = stan_initialize( stan_code=sir_stan_model_code( selection="discrete_binomial_process" ) )
stancode$compile()

fit = stancode$sample(        
  data = stan_data,
  chains = 1,
  warmup = 10000,
  iter = 15000
)

M = stan_extract( as_draws_df( fit$draws() ) )


# ---------------------------------------
# 4. stochastic discrete form via STAN/MCMC

# CREDIT: Arie Voorman, April 2017
# https://rstudio-pubs-static.s3.amazonaws.com/270496_e28d8aaa285042f2be0c24fc915a68b2.html

require(dplyr)
require(ggplot2)
require(tidyr)
require(Matrix)


options(mc.cores = parallel::detectCores())

# set.seed(42)

i0 = 25/300

Npop = 300
times = 0:99
params = list(prob_infection=0.001, time_shedding=7)
# // prob_infection = probability of an individual infecting another in 1 unit of time (1 day)
# // time_shedding = 1/ probability of transition to recovered state (~ duration is units of time) .. ie., simple geometric  .. 14 days
inits = c(1-i0, i0, 0) # S, I, R, proportions

sim = simulate_data( selection="sir_stochastic", Npop=Npop, inits=inits, times=times, params=params, plotdata=TRUE )
  # prob_infection =  transmission probability (per person per 'close contact', per day)
  # time_shedding  =  mean shedding duration for an individual

stan_data = list(
    Npop = Npop,
    Nobs = nrow(sim),
    Sobs = sim$Sobs,
    Iobs = sim$Iobs,
    Robs = sim$Robs
)


stancode = stan_initialize( stan_code=sir_stan_model_code( selection="discrete_voorman_2017_basic" ) )
stancode$compile()

fit = stancode$sample(        
  data = stan_data,
  chains = 1,
  warmup = 10000,
  iter = 15000
)

M = stan_extract( as_draws_df( fit$draws() ) )

stan_dens(fit)

print(fit, digits = 4)

plot_model_fit( M=M, stan_data=stan_data )





# ---------------------------------------
# 5. discrete stochastic form via STAN/MCMC with covariates


Npop = 300 # number of subjects
times = 0:200
params = list(
  prob_infection=c(family=0.1, community=0.0005),
  time_shedding=c(child=10, adult=3)
)
  # lambda_family = 0.1 # transmission probability (per person per 'close contact', per day)
  # lambda_community = 0.0005# transmission probability (per person per 'person in community', per day)
  # gamma_child = 10 # mean shedding duration for child (naive child)
  # gamma_adult = 3 # mean shedding duration for an adult (prior immunity)

I0 = c(25/300, 0) # number of c(children,adults) initially infected

inits = list(
  children = c(1-i0[1], i0[1], 0), # S, I, R, proportions for children
  adults   = c(1-i0[2], i0[2], 0)  # S, I, R, proportions for adults
)

sim = simulate_data( selection="sir_stochastic_family", Npop=Npop, inits=inits, times=times, params=params, plotdata=TRUE )
  # prob_infection =  transmission probability (per person per 'close contact', per day)
  # time_shedding  =  mean shedding duration for an individual


stan_data = list(
    Npop = Npop,
    Nobs = length(sim$time),
    Sobs = sim$Sobs,
    Iobs = sim$Iobs,
    Robs = sim$Robs,
    fmat = fmat,
    adult= adult,
    lambda_max = 1,
    dt = rep(1,t+1)
)


stancode = stan_initialize( stan_code=sir_stan_model_code( selection="discrete_voorman_2017_covariates" ) )
stancode$compile()

fit = stancode$sample(        
  data = stan_data,
  chains = 3,
  warmup = 1000,
  iter = 1500
)

M = stan_extract( as_draws_df( fit$draws() ) )

stan_dens(fit)

print(fit, digits = 4)



# ---------------------------------------
# 6. discrete stochastic form via STAN/MCMC with variable time points

sim = simulate_data( selection="sir_stochastic_family", Npop=Npop, inits=inits, times=times, params=params, plotdata=TRUE )
  # prob_infection =  transmission probability (per person per 'close contact', per day)
  # time_shedding  =  mean shedding duration for an individual


# subsample ..
time.pts = c(1,3, 6, 10)

sim = sim[ time.pts, ]

stan_data = list(
    Npop = Npop,
    Nobs = length(sim$time),
    Sobs = sim$Sobs,
    Iobs = sim$Iobs,
    Robs = sim$Robs,
    fmat = as.matrix(fmat),
    adult = adult,
    lambda_max = 1,
    dt = c(diff(time.pts),0)
)


stancode = stan_initialize( stan_code=sir_stan_model_code( selection="discrete_voorman_2017_covariates_irregular_time" ) )
stancode$compile()

fit = stancode$sample(        
  data = stan_data,
  chains = 3,
  warmup = 1000,
  iter = 1500
)

M = stan_extract( as_draws_df( fit$draws() ) )

stan_dens(fit)

print(fit,digits = 4)

stan_dens(fit)
