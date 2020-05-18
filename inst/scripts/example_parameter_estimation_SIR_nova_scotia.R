
# NOTE: 17 May 2020, JSC
# NOTE: this script is now subsumed by "example_parameter_estimation_SIR_provinces_of_Canada.R"
# NOTE: it is here for reference to show how to oeprate upon a single areal unit


# ---------------------------------------
# examples of various methods of parameter estimation using NS data

# remotes::install_github( "jae0/adapt" )
require(adapt)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# loadfunctions("adapt")

nscovid = data_nova_scotia(
  # interpolate_missing_data=TRUE,  # linear interpolation of missing data as a preprocessing step or estimate via imputation inside stan
  Npop = 971395,  # total population
  Npreds = 15,   # number of days for forward projectins
  BNP = 4,       # AR(BNP) process for beta number ( BNP= no of days lag) and also no days to average for forward projections
  BETA_max = 10,       # max rate param for S -> I  # approx number of contacts per person per time (day) multiplied by the probability of disease transmission in a contact between a susceptible and an infectious subject;  ~ 1/( typical time in days between contacts)
  GAMMA_max = 0.1,    # max rate param for I -> R  # ~ 1/(typical time until removal = 14) = 0.07
  EPSILON_max = 0.1,  # max rate param for I -> M  # about 5% seem to die ..
  modelname = "default" # "discrete_autoregressive_structured_beta_mortality_hybrid"  # splitting recovered and mortalities
)

time_relaxation = as.numeric(nscovid$time_relaxation - nscovid$time_start)
time_distancing = as.numeric(nscovid$time_distancing - nscovid$time_start)


stancode_compiled = rstan::stan_model( model_code=sir_stan_model_code( selection=nscovid$modelname ) )  # compile the code


f = rstan::sampling( stancode_compiled, data=nscovid, chains=3, warmup=7000, iter=10000, control= list(adapt_delta = 0.95, max_treedepth=15 ))


if (0) {
  fn = file.path("~", "tmp", paste( nscovid$modelname, "rdata", sep="."))
  save( f, file=fn, compress=TRUE)
  load(fn)
}

  M = extract(f)

  province = "Nova Scotia"
  outdir = file.path( "~", "bio", "adapt", "inst", "doc", province )
  to.screen = TRUE
      # to.screen = FALSE

  plot_model_fit( selection="susceptible", stan_data=nscovid, M=M, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="infected", stan_data=nscovid, M=M, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="recovered", stan_data=nscovid, M=M, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="deaths", stan_data=nscovid, M=M, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="reproductive_number", stan_data=nscovid, M=M, outdir=outdir, to.screen=to.screen )
  plot_model_fit( selection="reproductive_number_histograms", stan_data=nscovid, M=M, outdir=outdir, to.screen=to.screen )


# --- now some simplistic stochastic simulations using joint posterior distributions from current day estimates:

sim = simulate( M, istart=nscovid$Nobs-1, nsims=2000, nprojections=150 )
plot_model_fit( selection="forecasts", stan_data=nscovid, M=M, outdir=outdir, sim=sim )


if (0) {

    plot(f)
    plot(f, pars="I")
    print(f)

    traceplot(f)
    e = rstan::extract(f, permuted = TRUE) # return a list of arrays
    m2 = as.array(f)
    traceplot(f, pars=c("GAMMA"))
    traceplot(f, pars=c("MSErrorI"))
    traceplot(f, pars=c("MSErrorR"))
    traceplot(f, pars="lp__")
    summary(f)$summary[,"K1"]
    est=colMeans(M)
    prob=apply(M,2,function(x) I(length(x[x>0.10])/length(x) > 0.8)*1)

}
