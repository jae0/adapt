

# ---------------------------------------

require(adapt)
require(cmdstanr )
require(posterior )
require(bayesplot)

options(mc.cores = parallel::detectCores())
require(SimInf)

library(scales)


# loadfunctions("adapt")
# remotes::install_github( "jae0/adapt" )


workdir = file.path( "~", "tmp" )
setwd( workdir )
fn = file.path( workdir, "Covid19Canada.rdata")


tasks = c("download.data", "model", "plot", "forecast")



if ("download.data" %in% tasks) res = data_provinces_of_canada( selection="download", fn=fn )

can = data_provinces_of_canada(
  fn = fn,
  Npreds = 7,   # number of days for ode-based forward projections
  BNP = 1,       # beta dynamics is AR(BNP) ; also the number of days to average for forward ode-based projections (incubation time is ~ 5-7 days) .. higher than 1 can cause problems ... high var in reporting causes + and - corrs
  BETA_max = 1.2,     # max rate param for S -> I  # approx number of contacts per person per time (day) multiplied by the probability of disease transmission in a contact between a susceptible and an infectious subject;  ~ 1/( typical time in days between contacts)
  # BETA_max is very important: seldom does this value go > 1 for Covid-19 in Canada,
  GAMMA_max = 0.2,    # max rate param for I -> R  # ~ 1/(typical time until removal = 14) = 0.07
  EPSILON_max = 0.1,  # max rate param for I -> M  # > recovery time; < rate ..
  modelname="default"
)

provinces = names(can)
# provinces = setdiff( provinces, c("Nova Scotia", "Quebec", "Ontario", "Alberta", "BC" ) )
# provinces = c( "Nova Scotia", "New Brunswick", "Quebec", "Ontario", "BC")

# compile code

stancode = stan_initialize( stan_code=sir_stan_model_code( selection="default_asymptomatic" ) )
stancode$compile()

to.screen = FALSE
# to.screen = TRUE

for (au in  provinces) {
    # au = "Nova Scotia"

    print(au)
    stan_results = list( stan_inputs=can[[au]] )
    fn_model = file.path( workdir, paste( au, stan_results$stan_inputs$modelname, "rdata", sep=".") )
    outdir = file.path( "~", "bio", "adapt", "inst", "doc", au)

    if ("model" %in% tasks ) {
      fit = stancode$sample(
        data=stan_results$stan_inputs, 
        chains=3, 
        warmup=5000, 
        iter=6000, 
        adapt_delta = 0.9, 
        max_treedepth=12
      )        
  
      stan_results$posteriors = stan_extract( as_draws_df( fit$draws() ) )
  
      save(stan_results, file=fn_model, compress=TRUE)
    } else {
      load(fn_model)
    }

    if ( "plot" %in% tasks ) {
      plot_model_fit( selection="susceptible", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
      plot_model_fit( selection="infected", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
      plot_model_fit( selection="recovered", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
      plot_model_fit( selection="deaths", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
      plot_model_fit( selection="reproductive_number", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
      plot_model_fit( selection="reproductive_number_histograms", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
      plot_model_fit( selection="infected_affected", stan_results=stan_results, outdir=outdir, to.screen=to.screen )
    }

    if ("forecast" %in% tasks ) {
      # --- simplistic stochastic simulations using joint posterior distributions from "current" day estimates:, if BNP is provided, this uses the average in the period specified
      sim = simulate( stan_results=stan_results, nsims=1400, nprojections=200, nthreads=1  ) # predictions are conditioned on beta estimates from t-1 by default
      plot_model_fit( selection="forecasts", stan_results=stan_results, sim=sim, outdir=outdir,
        to.screen=to.screen )
    }

}


## Comparisons across provinces: normalize to unit population
fn.summary = file.path( workdir, "Covid19Canada_summary.rdata")
modelname = can[[1]]$modelname

res = summary_adapt( "summary.create", aus=provinces, fn=fn.summary, workdir=workdir, modelname=modelname )

res = summary_adapt( "plot_reproductive_number_histograms", aus=provinces, fn=fn.summary, workdir=workdir, modelname=modelname, to.screen=TRUE )

res = summary_adapt( "plot_all", aus=provinces, fn=fn.summary, workdir=workdir, modelname=modelname, to.screen=FALSE )

