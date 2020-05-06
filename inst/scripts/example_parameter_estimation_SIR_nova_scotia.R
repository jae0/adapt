

# ---------------------------------------
# examples of various methods of parameter estimation using NS data


# convert to data.frame
# loadfunctions("adapt")

# remotes::install_github( "jae0/adapt" )
require(adapt)

require(rstan)
require(SimInf)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


stan_data = data_nova_scotia(
  output = "stan_data",
  Npop = 971395,  # total population
  Npreds = 30,   # number of days for forward projectins
  BNP = 3,        # beta number of days to average for forward projections
  GAMMA_prior = 1/28,  # approx scale of GAMMA ("effective recovery rate") ~ 1/ recovery time (about 21 to 28 days)
  BETA_prior = 0.9,    # approx scale of BETA  ("effective infection rate", 1 -> 100%)
  # modelname = "discrete_autoregressive_with_observation_error_structured_beta"
  # modelname = "discrete_autoregressive_with_observation_error_unstructured_beta"  # slow
  # modelname = "continuous" # ODE form  really slow ... and does not work that well .. huge error bars
  # modelname = "discrete_basic"   # poor fit ..  latent .. similar to continuous .. ie. model is too simple and params too static
  modelname = "discrete_autoregressive_without_observation_error"  # no obsrevation error, only process error
)


stancode_compiled = rstan::stan_model( model_code=sir_stan_model_code( selection=stan_data$modelname ) )  # compile the code

f = rstan::sampling( stancode_compiled, data=stan_data, chains=3, warmup=4000, iter=6500, control= list(adapt_delta = 0.95, max_treedepth=14 ))


if (0) {
  fn = file.path("~", "tmp", paste( stan_data$modelname, "rdata", sep="."))
  save( f, file=fn, compress=TRUE)
  load(fn)
}


M = extract(f)


plot_model_fit( stan_data=stan_data, M=M )

outdir = file.path( "~/bio/adapt/inst/doc/")
nx = stan_data$Nobs + stan_data$Npreds - 1


io = stan_data$Iobs
io[io < 0] = NA

png(filename = file.path(outdir, "fit_with_projections_infected.png"))
  xrange = c(0, nx)
  yrange = c(0, max(M$I[, 1:nx]))
  plot( io ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$I, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( io ~ stan_data$time, col="darkgray", cex=1.2 )
  abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
  legend( "topright", "", paste( "Current date: ", Sys.Date(), "   "), bty="n")
dev.off()


ro = stan_data$Robs
ro[ro < 0] = NA
png(filename = file.path(outdir, "fit_with_projections_recovered.png"))
  xrange = c(0, nx)
  yrange = c(0, max(M$R[, 1:nx]))
  plot( ro ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Recovered", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$R, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$R, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$R, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( ro ~ stan_data$time, col="darkgray", cex=1.2 )
  abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
  legend( "topright", "", paste( "Current date: ", Sys.Date(), "   "), bty="n")
dev.off()


so = stan_data$Sobs
so[so < 0] = NA
png(filename = file.path(outdir, "fit_with_projections_susceptible.png"))
  xrange = c(0, nx)
  yrange = c(min(M$S[, 1:nx]), max(M$S[, 1:nx]))
  plot( so ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Susceptible", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$S, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$S, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$S, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( so ~ stan_data$time, col="darkgray", cex=1.2 )
  abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
  legend( "topright", "", paste( "Current date: ", Sys.Date(), "   "), bty="n")
dev.off()



png(filename = file.path(outdir, "reproductive_number.png"))
  plot( apply(M$K[,1:nx], 2, median) ~ seq(1,nx), lwd =3, col="darkgray", ylim=c(0,10), ylab="Reproductive number", xlab="Days (day 1 is 2020-03-17)" )
  lines( apply(M$K, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed" )
  lines( apply(M$K, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed" )
  abline( h=1, col="red", lwd=3 )
  abline( v=stan_data$Nob, , col="grey", lty="dashed" )
  legend( "topright", "", paste( "Current date: ", Sys.Date(), "   "), bty="n")
dev.off()


png(filename = file.path(outdir, "reproductive_number_today.png"))
  hist( M$K[,stan_data$Nobs] , "fd", xlab="Current Reproductive number", ylab="Frequency", main="")  # today's K
  abline( v=1, col="red", lwd=3 )
  legend( "topright", "", paste( "Current date: ", Sys.Date(), "   "), bty="n")
dev.off()



# --- now some simplistic stochastic simulations using joint posterior distributions from current day estimates:

require(SimInf)

nsims = nrow(M$BETA)
today = stan_data$Nobs
nprojections = 120
sim = array( NA, dim=c(nsims, 3, nprojections) )

for (i in 1:nsims) {
  sim[i,,] = run( SIR(
    u0=data.frame(S=M$S[i,today], I=M$I[i,today], R=M$R[i,today]),
    tspan=1:nprojections,
    beta=M$BETA[i,today],
    gamma=M$GAMMA[i] )
  )@U[]
}


library(scales)
png(filename = file.path(outdir, "fit_with_projections_and_stochastic_simulations.png"))
  simxval = stan_data$Nobs + c(1:nprojections)
  xrange = c(0, max(nx, simxval) )
  yrange = c(0, max(M$I[, 1:nx]))
  yrange = c(0, 1250)
  plot( io ~ stan_data$time, xlim=xrange, ylim=yrange,  type="n", ylab="Infected", xlab="Days (day 1 is 2020-03-17)")
  for ( i in 1:min(nsims, 1500)) lines( sim[i,2,] ~ simxval, col=alpha("slategray", 0.1), ltyp="dashed" )
  lines( apply(M$I, 2, median)[1:nx] ~ seq(1,nx), lwd = 3, col="slateblue" )
  lines( apply(M$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
  points( io ~ stan_data$time, xlim=xrange, ylim=yrange, col="darkgray", cex=1.2 )
  legend( "topright", "", paste( "Current date: ", Sys.Date(), "   "), bty="n")

dev.off()




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




