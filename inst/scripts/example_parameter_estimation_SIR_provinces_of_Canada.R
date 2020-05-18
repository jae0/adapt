

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


tasks = c("download.data", "model", "plot", "forecast")


to.screen = FALSE
# to.screen = TRUE


if ("download.data" %in% tasks) res = data_provinces_of_canada( selection="download", fn=fn )

can = data_provinces_of_canada(
  fn = fn,
  Npreds = 14,   # number of days for forward projections
  BNP = 3,        # beta number of days to average for forward projections
  modelname="default"
)

str(can)


# compile code

if ("model" %in% tasks ) stancode_compiled = rstan::stan_model( model_code=sir_stan_model_code( selection="default" ) )  # compile the code

for (province in names(can) ) {

  # province = names(can)[i]
  print(province)

  fn_model = file.path( workdir, paste( province, can[[province]]$modelname, "rdata", sep=".") )
  outdir = file.path( "~", "bio", "adapt", "inst", "doc", province)

  if ("model" %in% tasks ) {
    f = rstan::sampling( stancode_compiled, data=can[[province]], chains=4, warmup=5000, iter=7000, control = list(adapt_delta = 0.95, max_treedepth=14 ) )
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



## Comparisons across provinces: normalize to unit population

res = list()

for (province in names(can) ) {

  # province = names(can)[i]
  print(province)
  Npop = can[[province]]$Npop

  fn_model = file.path( workdir, paste( province, can[[province]]$modelname, "rdata", sep=".") )
  outdir = file.path( "~", "bio", "adapt", "inst", "doc", province)

  load(fn_model)
  M = extract(f)

  nsims = nrow(M$BETA)
  today = can[[province]]$Nobs
  nprojections = 120
  sim = array( NA, dim=c(nsims, 3, nprojections) )

  )

  res[[province]]$time = can[[province]]$time
  res[[province]]$Npop = can[[province]]$Npop
  res[[province]]$Nobs = can[[province]]$Nobs
  res[[province]]$Npreds = can[[province]]$Npreds
  res[[province]]$Nts = res[[province]]$Nobs + res[[province]]$Npreds
  res[[province]]$timeall = 1:res[[province]]$Nts
  res[[province]]$S = data.frame( cbind(
    median = apply(M$S/res[[province]]$Npop, 2, median, na.rm=TRUE),
    low = apply(M$S/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
    high = apply(M$S/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
  ))
  res[[province]]$I = data.frame( cbind(
    median = apply(M$I/res[[province]]$Npop, 2, median, na.rm=TRUE),
    low = apply(M$I/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
    high = apply(M$I/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
  ))
  res[[province]]$R = data.frame( cbind(
    median = apply(M$R/res[[province]]$Npop, 2, median, na.rm=TRUE),
    low = apply(M$R/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
    high = apply(M$R/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
  ))
  res[[province]]$M = data.frame( cbind(
    median = apply(M$M/res[[province]]$Npop, 2, median, na.rm=TRUE),
    low = apply(M$M/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
    high = apply(M$M/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
  ))
  res[[province]]$GAMMA = data.frame( cbind(
    median = apply(t(M$GAMMA), 1 median, na.rm=TRUE),
    low = apply(t(M$GAMMA), 1, quantile, probs=c(0.025), na.rm=TRUE),
    high = apply(t(M$GAMMA), 1, quantile, probs=c(0.975), na.rm=TRUE)
  ))
  res[[province]]$EPSILON = data.frame( cbind(
    median = apply(t(M$EPSILON), 1, median, na.rm=TRUE),
    low = apply(t(M$EPSILON), 1, quantile, probs=c(0.025), na.rm=TRUE),
    high = apply(t(M$EPSILON), 1, quantile, probs=c(0.975), na.rm=TRUE)
  ))
  res[[province]]$K = data.frame( cbind(
    median = apply(M$K, 2, median, na.rm=TRUE),
    low = apply(M$K, 2, quantile, probs=c(0.025), na.rm=TRUE),
    high = apply(M$K, 2, quantile, probs=c(0.975), na.rm=TRUE)
  ))

}


for (province in names(can) ) {

  plot( 0,0, xlim=...)
  lines( res[[province]]$I$median ~ res[[province]]$timeall, col///  )

}