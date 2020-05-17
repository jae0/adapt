

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

require(SimInf)

nsims = nrow(M$BETA)
today = nscovid$Nobs
nprojections = 120
sim = array( NA, dim=c(nsims, 3, nprojections) )


if (nscovid$modelname=="discrete_autoregressive_with_observation_error_structured_beta_mortality") {
  u0=data.frame(S=M$S[,today], I=M$I[,today], R=M$R[,today] + M$M[,today], beta=M$BETA[,today-1], gamma=M$GAMMA[] )
} else {
  u0=data.frame(S=M$S[,today], I=M$I[,today], R=M$R[,today], beta=M$BETA[,today-1], gamma=M$GAMMA[] )
}


for (i in 1:nsims) {
  sim[i,,] = run( SIR(
    u0=u0[i,c("S","I","R")],
    tspan=1:nprojections,
    beta=u0$beta[i],
    gamma=u0$gamma[i] )
  )@U[]
}

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
