

# ---------------------------------------

require(adapt)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
require(SimInf)

# loadfunctions("adapt")
# remotes::install_github( "jae0/adapt" )


workdir = file.path( "~", "tmp" )
setwd( workdir )
fn = file.path( workdir, "Covid19Canada.rdata")


tasks = c("download.data", "model", "plot", "forecast")


to.screen = FALSE
# to.screen = TRUE


if ("download.data" %in% tasks) res = data_provinces_of_canada( selection="download", fn=fn )

can = data_provinces_of_canada(
  fn = fn,
  Npreds = 20,   # number of days for ode-based forward projections
  BNP = 1,       # beta dynamics is AR(BNP) ; also the number of days to average for forward ode-based projections (incubation time is ~ 5-7 days) .. higher than 1 can cause problems ... high var in reporting causes + and - corrs
  BETA_max = 1.0,     # max rate param for S -> I  # approx number of contacts per person per time (day) multiplied by the probability of disease transmission in a contact between a susceptible and an infectious subject;  ~ 1/( typical time in days between contacts)
  # BETA_max is very important: seldom does this value go > 1 for Covid-19 in Canada,
  # if BETA_max is set too large, and due to the long tails of a cauchy, convergence can be slow and error distributions become very wide when infected numbers -> 0, 1 is a good upper bound
  GAMMA_max = 0.1,    # max rate param for I -> R  # ~ 1/(typical time until removal = 14) = 0.07
  EPSILON_max = 0.1,  # max rate param for I -> M  # > recovery time; < rate ..
  modelname="default"
)

provinces = names(can)


if ("model" %in% tasks ) {
# compile code

  stancode_compiled = rstan::stan_model( model_code=sir_stan_model_code( selection="default" ) )  # compile the code

  for (au in  provinces) {
    print(au)
    fn_model = file.path( workdir, paste( au, can[[au]]$modelname, "rdata", sep=".") )
    outdir = file.path( "~", "bio", "adapt", "inst", "doc", au)
    control.stan = list(adapt_delta = 0.95, max_treedepth=15 )
    # some  au's  have longer and more complex dynamics (i.e. parameter space)and likely reporting issues ... requires additional stabilzation
      if ( au %in% c( "Quebec", "Ontario" ) ) {
      #  control.stan = list(adapt_delta = 0.975, max_treedepth=15 )
       #  can[[au]]$BNP = 3
      }

    f = rstan::sampling( stancode_compiled, data=can[[au]], chains=3, warmup=6000, iter=8000, control=control.stan  )
    save(f, file=fn_model, compress=TRUE)
  }
}


if ( "plot" %in% tasks ) {

  for (au in  provinces) {
    print(au)
    fn_model = file.path( workdir, paste( au, can[[au]]$modelname, "rdata", sep=".") )
    outdir = file.path( "~", "bio", "adapt", "inst", "doc", au)
    load(fn_model)
    M = extract(f)
    plot_model_fit( selection="susceptible", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="infected", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="infected_effective", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="recovered", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="deaths", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="reproductive_number", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="reproductive_number_histograms", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="effective_number", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
  }
}


if ("forecast" %in% tasks ) {
  for (au in  provinces) {
    print(au)
    fn_model = file.path( workdir, paste( au, can[[au]]$modelname, "rdata", sep=".") )
    outdir = file.path( "~", "bio", "adapt", "inst", "doc", au)
    load(fn_model)
    M = extract(f)
    # --- simplistic stochastic simulations using joint posterior distributions from "current" day estimates:, if BNP is provided, this uses the average in the period specified
    sim = simulate( M, istart=can[[au]]$Nobs, nsims=1400, nprojections=200  )
    plot_model_fit( selection="forecasts", stan_data=can[[au]], M=M, outdir=outdir, sim=sim, to.screen=to.screen )
  }
}



## Comparisons across provinces: normalize to unit population
fn.summary = file.path( workdir, "Covid19Canada_summary.rdata")

res = summary_adapt( "summary.create", can=can, fn=fn.summary )
res = summary_adapt( "plot_reproductive_number_histograms", can=can, fn=fn.summary, to.screen=TRUE )
res = summary_adapt( "plot_all", can=can, fn=fn.summary, to.screen=TRUE )

