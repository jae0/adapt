
# NOTE: 17 May 2020, JSC
# NOTE: this script is now subsumed by "example_parameter_estimation_SIR_provinces_of_Canada.R"
# NOTE: it is here for reference to show how to oeprate upon a single areal unit


# ---------------------------------------
# examples of various methods of parameter estimation using NS data

# remotes::install_github( "jae0/adapt" )
require(adapt)
require(cmdstanr )
require(posterior )
require(bayesplot)


options(mc.cores = parallel::detectCores())

# loadfunctions("adapt")

stan_results = list( stan_inputs =  data_nova_scotia(
  # interpolate_missing_data=TRUE,  # linear interpolation of missing data as a preprocessing step or estimate via imputation inside stan
  Npop = 971395,  # total population
  Npreds = 15,   # number of days for forward projectins
  BNP = 4,       # AR(BNP) process for beta number ( BNP= no of days lag) and also no days to average for forward projections
  BETA_max = 10,       # max rate param for S -> I  # approx number of contacts per person per time (day) multiplied by the probability of disease transmission in a contact between a susceptible and an infectious subject;  ~ 1/( typical time in days between contacts)
  GAMMA_max = 0.1,    # max rate param for I -> R  # ~ 1/(typical time until removal = 14) = 0.07
  EPSILON_max = 0.1,  # max rate param for I -> M  # about 5% seem to die ..
  modelname = "default" # "discrete_autoregressive_structured_beta_mortality_hybrid"  # splitting recovered and mortalities
) )



time_relaxation = as.numeric(stan_results$stan_inputs$time_relaxation - stan_results$stan_inputs$time_start)
time_distancing = as.numeric(stan_results$stan_inputs$time_distancing - stan_results$stan_inputs$time_start)

stancode = stan_initialize( stan_code=sir_stan_model_code( selection=stan_results$stan_inputs$modelname ) )
stancode$compile()

fit = stancode$sample(        
  data=stan_results$stan_inputs, 
  iter_warmup = 7000,
  iter_sampling = 3000,
  seed = 123,
  chains = 3,
  parallel_chains = 3,  # The maximum number of MCMC chains to run in parallel.
  max_treedepth = 15,
  adapt_delta = 0.95,
  refresh = 500
 )

fit$save_object( file = fnfit )   #  save this way due to R-lazy loading
fit = readRDS( fnfit )


 
  stan_results$posteriors = stan_extract( as_draws_df( fit$draws() ) )
 
 
  province = "Nova Scotia"
  outdir = file.path( "~", "bio", "adapt", "inst", "doc", province )
  to.screen = TRUE
      # to.screen = FALSE

  plot_model_fit( selection="susceptible", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="infected", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="recovered", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="deaths", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="reproductive_number", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="reproductive_number_histograms", stan_results=stan_results, outdir=outdir, to.screen=to.screen )


# --- now some simplistic stochastic simulations using joint posterior distributions from current day estimates:

sim = simulate( stan_results, nsims=2000, nprojections=150 )
plot_model_fit( selection="forecasts", stan_results=stan_results, outdir=outdir, sim=sim )

 