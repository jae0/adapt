

# ---------------------------------------
# examples of various methods of parameter estimation using NS data

# remotes::install_github( "jae0/adapt" )
# loadfunctions("adapt")
require(adapt)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


nscovid = data_nova_scotia(
  # interpolate_missing_data=TRUE,  # linear interpolation of missing data as a preprocessing step or estimate via imputation inside stan
  Npop = 971395,  # total population
  Npreds = 30,   # number of days for forward projectins
  BNP = 3,        # beta number of days to average for forward projections
  BETA_prior = 0.8,    # approx scale of BETA  ("effective infection rate", 1 -> 100%)
  GAMMA_prior = 1/25,  # approx scale of GAMMA ("effective recovery rate") ~ 1/ recovery time (about 21 to 28 days)
  EPSILON_prior = (1/25)/20,  # 1/20 of GAMMA .. porportion of positive tested that die is ~5%, so 1/20 of GAMMA_prior
  # modelname = "continuous" # ODE form  really slow ... and does not work that well .. huge error bars
  # modelname = "discrete_basic"   # poor fit ..  latent .. similar to continuous .. ie. model is too simple and params too static
  # modelname = "discrete_binomial_autoregressive"
  # modelname = "discrete_autoregressive_structured_beta_mortality_hybrid"  # splitting recovered and mortalities
  # modelname = "discrete_autoregressive_structured_beta_mortality_poisson"
  modelname = "discrete_autoregressive_structured_beta_mortality"
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

plot_model_fit( stan_data=nscovid, M=M )

outdir = file.path( "~/bio/adapt/inst/doc/")
nx = nscovid$Nobs + nscovid$Npreds - 1


io = nscovid$Iobs
io[io < 0] = NA

png(filename = file.path(outdir, "fit_with_projections_infected.png"))
  xrange = c(0, nx)
  yrange = c(0, quantile(M$I[, 1:nx], probs=0.99) )
  plot( io ~ nscovid$time, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$I, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( io ~ nscovid$time, col="darkgray", cex=1.2 )
  abline( v=nscovid$time[nscovid$Nobs], col="grey", lty="dashed" )
  abline( v=nscovid$time[time_distancing], col="orange", lty="dotted" )
  abline( v=nscovid$time[time_relaxation], col="green", lty="dotted" )
  legend( "bottomleft", "", "       [-- Social distancing -->", bty="n" )
  legend( "bottom", "", "                                  [-- Parks open -->", bty="n" )
  legend( "topright", "", paste( "Current date: ", nscovid$timestamp ), bty="n")
dev.off()


ro = nscovid$Robs
ro[ro < 0] = NA
png(filename = file.path(outdir, "fit_with_projections_recovered.png"))
  xrange = c(0, nx)
  yrange = c(0, quantile(M$R[, 1:nx], probs=0.99) )
  plot( ro ~ nscovid$time, xlim=xrange, ylim=yrange, ylab="Recovered", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$R, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$R, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$R, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( ro ~ nscovid$time, col="darkgray", cex=1.2 )
  abline( v=nscovid$time[nscovid$Nobs], col="grey", lty="dashed" )
  abline( v=nscovid$time[time_distancing], col="orange", lty="dotted" )
  abline( v=nscovid$time[time_relaxation], col="green", lty="dotted" )
  legend( "topleft", "", "\n\n       [-- Social distancing -->", bty="n" )
  legend( "top", "", "\n\n                                  [-- Parks open -->", bty="n" )
  legend( "topright", "", paste( "Current date: ", nscovid$timestamp ), bty="n")
dev.off()


mo = nscovid$Mobs
mo[mo < 0] = NA
png(filename = file.path(outdir, "fit_with_projections_mortalities.png"))
  xrange = c(0, nx)
  yrange = c(0, quantile(M$M[, 1:nx], probs=0.99) )
  plot( mo ~ nscovid$time, xlim=xrange, ylim=yrange, ylab="Mortalities", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$M, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$M, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$M, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( mo ~ nscovid$time, col="darkgray", cex=1.2 )
  abline( v=nscovid$time[nscovid$Nobs], col="grey", lty="dashed" )
  abline( v=nscovid$time[time_distancing], col="orange", lty="dotted" )
  abline( v=nscovid$time[time_relaxation], col="green", lty="dotted" )
  legend( "topleft", "", "\n\n       [-- Social distancing -->", bty="n" )
  legend( "top", "", "\n\n                                  [-- Parks open -->", bty="n" )
  legend( "topright", "", paste( "Current date: ", nscovid$timestamp, " "), bty="n")
dev.off()


so = nscovid$Sobs
so[so < 0] = NA
png(filename = file.path(outdir, "fit_with_projections_susceptible.png"))
  xrange = c(0, nx)
  yrange = c(min(M$S[, 1:nx]), max(M$S[, 1:nx]))
  plot( so ~ nscovid$time, xlim=xrange, ylim=yrange, ylab="Susceptible", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$S, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$S, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$S, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( so ~ nscovid$time, col="darkgray", cex=1.2 )
  abline( v=nscovid$time[nscovid$Nobs], col="grey", lty="dashed" )
  abline( v=nscovid$time[time_distancing], col="orange", lty="dotted" )
  abline( v=nscovid$time[time_relaxation], col="green", lty="dotted" )
  legend( "bottomleft", "", "       [-- Social distancing -->\n\n", bty="n" )
  legend( "bottom", "", "                                  [-- Parks open -->\n\n", bty="n" )
  legend( "topright", "", paste( "Current date: ", nscovid$timestamp ), bty="n")
dev.off()



png(filename = file.path(outdir, "reproductive_number.png"))
  plot( apply(M$K[,1:nx], 2, median) ~ seq(1,nx), lwd =3, col="darkgray", ylim=c(0,10), ylab="Reproductive number", xlab="Days (day 1 is 2020-03-17)" )
  lines( apply(M$K, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed" )
  lines( apply(M$K, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed" )
  abline( h=1, col="red", lwd=3 )
  abline( v=nscovid$time[nscovid$Nobs], col="grey", lty="dashed" )
  abline( v=nscovid$time[time_distancing], col="orange", lty="dotted" )
  abline( v=nscovid$time[time_relaxation], col="green", lty="dotted" )
  legend( "topleft", "", "\n\n       [-- Social distancing -->", bty="n" )
  legend( "top", "", "\n\n                                  [-- Parks open -->", bty="n" )
  legend( "topright", "", paste( "Current date: ", nscovid$timestamp ), bty="n")
dev.off()


png(filename = file.path(outdir, "reproductive_number_today.png"))
  hist( M$K[,nscovid$Nobs] , "fd", xlab="Current Reproductive number", ylab="Frequency", main="", xlim=c(0, max(1.2, M$K[,nscovid$Nobs])) )  # today's K
  abline( v=1, col="red", lwd=3 )
  legend( "topright", "", paste( "Current date: ", nscovid$timestamp, " "), bty="n")
dev.off()




# --- now some simplistic stochastic simulations using joint posterior distributions from current day estimates:

require(SimInf)

nsims = nrow(M$BETA)
today = nscovid$Nobs
nprojections = 120
sim = array( NA, dim=c(nsims, 3, nprojections) )


if (nscovid$modelname=="discrete_autoregressive_with_observation_error_structured_beta_mortality") {
  u0=data.frame(S=M$S[,today], I=M$I[,today], R=M$R[,today] + M$M[,today], beta=M$BETA[,today], gamma=M$GAMMA[] )
} else {
  u0=data.frame(S=M$S[,today], I=M$I[,today], R=M$R[,today], beta=M$BETA[,today], gamma=M$GAMMA[] )
}


for (i in 1:nsims) {
  sim[i,,] = run( SIR(
    u0=u0[i,c("S","I","R")],
    tspan=1:nprojections,
    beta=u0$beta[i],
    gamma=u0$gamma[i] )
  )@U[]
}



# library(scales)
png(filename = file.path(outdir, "fit_with_projections_and_stochastic_simulations.png"))
  simxval = nscovid$Nobs + c(1:nprojections)
  xrange = c(0, max(nx, simxval) )
  yrange = c(0, max(M$I[, 1:nx]))
  yrange = c(0, 800)
  plot( io ~ nscovid$time, xlim=xrange, ylim=yrange,  type="n", ylab="Infected", xlab="Days (day 1 is 2020-03-17)")
  for ( i in 1:min(nsims, 1500)) lines( sim[i,2,] ~ simxval, col=alpha("slategray", 0.1), ltyp="dashed" )
  lines( apply(M$I, 2, median)[1:nx] ~ seq(1,nx), lwd = 3, col="slateblue" )
  lines( apply(M$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( io ~ nscovid$time, xlim=xrange, ylim=yrange, col="darkgray", cex=1.2 )
  abline( v=nscovid$time[nscovid$Nobs], col="grey", lty="dashed" )
  abline( v=nscovid$time[time_distancing], col="orange", lty="dotted" )
  abline( v=nscovid$time[time_relaxation], col="green", lty="dotted" )
  legend( "topleft", "", "\n\n   [-- Social distancing -->", bty="n" )
  legend( "topleft", "", "\n\n\n                               [-- Parks open -->", bty="n" )
  legend( "top", "", paste( "Current date: ", nscovid$timestamp ), bty="n")

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
