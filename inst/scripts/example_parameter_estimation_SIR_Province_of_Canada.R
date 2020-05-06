

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



storage_data_file = file.path( "~", "tmp", "Covid19Canada.rdata")

modelname = "discrete_simple_autoregressive" # sim to discrete_variable_encounter_rate_autoregressive but reduced params .. operate directly on BETA (vs IR*ER)

stan_data = data_provinces_of_canada( output="stan_data", fn=storage_data_file,  Npreds=150, BNP=5, GAMMA_prior=0.3,  modelname=modelname )



# BNP = no of days of BETA estimates to use for forward projection

str(stan_data)
plot(stan_data$Iobs ~ stan_data$time )


stancode_compiled = rstan::stan_model( model_code=sir_stan_model_code( selection=stan_data$modelname ) )  # compile the code


f = rstan::sampling( stancode_compiled, data=stan_data, chains=3, warmup=6000, iter=8000, control = list(adapt_delta = 0.9, max_treedepth=12 ) )


# save( f, file=file.path("/home", "jae", "tmp", paste( stan_data$modelname, "rdata", sep=".")), compress=TRUE)

M = extract(f)


outdir = file.path( "~/bio/adapt/inst/doc/")

nx = stan_data$Nobs + trunc( stan_data$Npreds *0.2 )

png(filename = file.path(outdir, "fit_with_projections.png"))
  xrange = c(0, nx)
  yrange = c(0, max(M$I[, 1:nx]))
  plot( stan_data$Iobs ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$I[,1:nx], 2, median) ~ seq(1,nx), lwd =4, col="blue" )
  lines( apply(M$I[,1:nx], 2, quantile, probs=0.025) ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 3 )
  lines( apply(M$I[,1:nx], 2, quantile, probs=0.975) ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 3 )
  points( stan_data$Iobs ~ stan_data$time, col="black", pch=20 )
  abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
dev.off()


plot_model_fit( stan_data=stan_data, M=M )

png(filename = file.path(outdir, "reproductive_number.png"))
  plot( apply(M$K[,1:nx], 2, median) ~ seq(1,nx), lwd =3, col="darkgray", ylim=c(0,10), ylab="Reproductive number", xlab="Days (day 1 is 2020-03-17)" )
  lines( apply(M$K[,1:nx], 2, quantile, probs=0.025) ~ seq(1,nx), col="darkorange", lty="dashed" )
  lines( apply(M$K[,1:nx], 2, quantile, probs=0.975) ~ seq(1,nx), col="darkorange", lty="dashed" )
  abline( h=1, col="red", lwd=3 )
  abline( v=stan_data$Nob, , col="grey", lty="dashed" )
dev.off()


png(filename = file.path(outdir, "reproductive_number_today.png"))
  hist( M$K[,stan_data$Nobs] , "fd", xlab="Reproductive number", ylab="Frequency", main="Posterior distribution for today")  # today's K
  abline( v=1, col="red", lwd=3 )
dev.off()



# --- now some simplistic stochastic simulations using posterior distributions from current day estimates:

require(SimInf)

nsims = nrow(M$BETA)
today = stan_data$Nobs
nprojections = 30
sim = array( NA, dim=c(nsims, 3, nprojections) )

for (i in 1:nsims) {
  res = run( SIR(
    u0=data.frame(S=M$S[i,today], I=M$I[i,today], R=M$R[i,today]),
    tspan=1:nprojections,
    beta=M$BETA[i,today],
    gamma=M$GAMMA[i] )
  )
  sim[i,,] = res@U[]
}


library(scales)
png(filename = file.path(outdir, "fit_with_projections_and_stochastic_simulations.png"))
  xrange = c(0, nx)
  yrange = c(0, max(M$I[, 1:nx]))
  plot( stan_data$Iobs ~ stan_data$time, xlim=xrange, ylim=yrange, col="darkblue", type="n", ylab="Infected", xlab="Days (day 1 is 2020-03-17)")
  for ( i in 1:min(nsims, 2000)) lines( sim[i,2,] ~ simxval, col=alpha("slategray", 0.1), ltyp="dashed" )
  lines( apply(M$I[,1:nx], 2, median) ~ seq(1,nx), lwd = 5, col="blue" )
  lines( apply(M$I[,1:nx], 2, quantile, probs=0.025) ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 3 )
  lines( apply(M$I[,1:nx], 2, quantile, probs=0.975) ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 3 )
  abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
  simxval = stan_data$Nobs + c(1:nprojections)
  points( stan_data$Iobs ~ stan_data$time, xlim=xrange, ylim=yrange, col="darkblue")

dev.off()




if (0) {
    plot(f)
    plot(f, pars="I")
    print(f)

    traceplot(f)
    e = rstan::extract(f, permuted = TRUE) # return a list of arrays
    m2 = as.array(f)
    traceplot(f, pars=c("params", "y0"))
    traceplot(f, pars="lp__")
    summary(f)$summary[,"Rhat"]
    est=colMeans(M)
    prob=apply(M,2,function(x) I(length(x[x>0.10])/length(x) > 0.8)*1)
}




