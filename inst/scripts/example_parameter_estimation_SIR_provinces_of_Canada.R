

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
  Npreds = 30,   # number of days for ode-based forward projections
  BNP = 5,       # beta dynamics is AR(BNP) ; also the number of days to average for forward ode-based projections (incubation time is ~ 5-7 days) .. higher than 1 can cause problems ... high var in reporting causes + and - corrs
  BETA_max = 5,       # max rate param for S -> I  # approx number of contacts per person per time (day) multiplied by the probability of disease transmission in a contact between a susceptible and an infectious subject;  ~ 1/( typical time in days between contacts)
  GAMMA_max = 0.1,    # max rate param for I -> R  # ~ 1/(typical time until removal = 14) = 0.07
  EPSILON_max = 0.1,  # max rate param for I -> M  # > recovery time; < rate ..
  modelname="default"
)
# str(can)



  # au = "Nova Scotia"
provinces = c(setdiff( sort( names(can) ), c("Ontario", "Quebec") ), c("Ontario", "Quebec"))  # do Ontario and Quebec last as they are very slow
# provinces = setdiff(provinces, "Nova Scotia")

if ("model" %in% tasks ) {
# compile code

  stancode_compiled = rstan::stan_model( model_code=sir_stan_model_code( selection="default" ) )  # compile the code

  for (au in  provinces) {
    print(au)
    fn_model = file.path( workdir, paste( au, can[[au]]$modelname, "rdata", sep=".") )
    outdir = file.path( "~", "bio", "adapt", "inst", "doc", au)
    control.stan = list(adapt_delta = 0.95, max_treedepth=14 )
    # if ( au %in% c("Quebec", "Ontario" ) ) {
    #    # these aus seem to have longer and more complex dynamics (i.e. parameter space) ... requires additional stabilzation
    #    control.stan = list(adapt_delta = 0.95, max_treedepth=15 )
    #    can[[au]]$BNP = 7
    # }
    f = rstan::sampling( stancode_compiled, data=can[[au]], chains=3, warmup=5000, iter=6000, control=control.stan  )
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
    plot_model_fit( selection="recovered", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="deaths", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="reproductive_number", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="reproductive_number_histograms", stan_data=can[[au]], M=M, outdir=outdir, to.screen=to.screen )
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
    sim = simulate( M, istart=can[[au]]$Nobs, nsims=2000, nprojections=150  )
    plot_model_fit( selection="forecasts", stan_data=can[[au]], M=M, outdir=outdir, sim=sim, to.screen=to.screen )
  }
}


## Comparisons across provinces: normalize to unit population
fn.summary = file.path( workdir, "Covid19Canada_summary.rdata")

# res = summary_adapt( can=can, fn=fn.summary )

summary_adapt( "plot.all", can=can, fn=fn.summary, to.screen=TRUE )

summary_adapt( "plot.all", can=can, fn=fn.summary, to.screen=FALSE )

