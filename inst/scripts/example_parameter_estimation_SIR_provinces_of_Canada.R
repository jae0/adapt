

# ---------------------------------------
# loadfunctions("adapt")
# remotes::install_github( "jae0/adapt" )


require(adapt)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
require(SimInf)


workdir = file.path( "~", "tmp" )
setwd( workdir )
fn = file.path( workdir, "Covid19Canada.rdata")

res = data_provinces_of_canada( selection="download", fn=fn )

can = data_provinces_of_canada(
  fn = fn,
  Npreds = 14,   # number of days for forward projections
  BNP = 3,        # beta number of days to average for forward projections
  modelname="default"
)

str(can)

tasks = c("model", "plot", "forecast")
to.screen = FALSE
# to.screen = TRUE

# compile code

if ("model" %in% tasks ) stancode_compiled = rstan::stan_model( model_code=sir_stan_model_code( selection="default" ) )  # compile the code

for (province in names(can) ) {

  # province = names(can)[i]
  print(province)

  fn_model = file.path( workdir, paste( province, can[[province]]$modelname, "rdata", sep=".") )
  outdir = file.path( "~", "bio", "adapt", "inst", "doc", province)

  if ("model" %in% tasks ) {
    f = rstan::sampling( stancode_compiled, data=can[[province]], chains=3, warmup=8000, iter=10000, control = list(adapt_delta = 0.95, max_treedepth=14 ) )
    save(f, file=fn_model, compress=TRUE)
  }


  load(fn_model)
  M = extract(f)

  if ( "plot" %in% tasks ) {

    plot_model_fit( selection="susceptible", stan_data=can[[province]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="infected", stan_data=can[[province]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="recovered", stan_data=can[[province]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="deaths", stan_data=can[[province]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="reproductive_number", stan_data=can[[province]], M=M, outdir=outdir, to.screen=to.screen )
    plot_model_fit( selection="reproductive_number_histograms", stan_data=can[[province]], M=M, outdir=outdir, to.screen=to.screen )
  }

  if ("forecast" %in% tasks ) {
    # --- simplistic stochastic simulations using joint posterior distributions from current day estimates:
    nsims = nrow(M$BETA)
    today = can[[province]]$Nobs
    nprojections = 120
    sim = array( NA, dim=c(nsims, 3, nprojections) )
    u0=data.frame(S=M$S[,today], I=M$I[,today], R=M$R[,today] + M$M[,today], beta=M$BETA[,today-1], gamma=M$GAMMA[] )

    for (i in 1:nsims) {
      sim[i,,] = run( SIR(
        u0=u0[i,c("S","I","R")], tspan=1:nprojections,
        beta=u0$beta[i], gamma=u0$gamma[i] )
      )@U[]
    }

    plot_model_fit( selection="forecasts", stan_data=can[[province]], M=M, outdir=outdir, sim=sim, to.screen=to.screen )

  }

}







