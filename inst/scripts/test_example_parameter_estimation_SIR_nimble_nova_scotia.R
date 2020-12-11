
# NOTE: this attempts to implement a BUGS/Nimble version of a simplified form of the STAN-based model
# NOTE: It fails to converge and will require more tweaking .. just here for posterity

# remotes::install_github( "jae0/adapt" )
require(adapt)
require(nimble)
require(SimInf)

# loadfunctions("adapt")

nscovid = data_nova_scotia(
  # interpolate_missing_data=TRUE,  # linear interpolation of missing data as a preprocessing step or estimate via imputation inside stan
  Npop = 971395,  # total population
  Npreds = 30,   # number of days for forward projectins
  BNP = 3,        # beta number of days to average for forward projections
  BETA_max = 10,       # max rate param for S -> I  # approx number of contacts per person per time (day) multiplied by the probability of disease transmission in a contact between a susceptible and an infectious subject;  ~ 1/( typical time in days between contacts)
  GAMMA_max = 0.1,    # max rate param for I -> R  # ~ 1/(typical time until removal = 14) = 0.07
  EPSILON_max = 0.1,  # max rate param for I -> M  # about 5% seem to die ..
  modelname = "discrete_autoregressive_structured_beta_mortality_hybrid_nimble"  # splitting recovered and mortalities
)

time_relaxation = as.numeric(nscovid$time_relaxation - nscovid$time_start)
time_distancing = as.numeric(nscovid$time_distancing - nscovid$time_start)


if (0) {
  fn = file.path("~", "tmp", paste( nscovid$modelname, "rdata", sep="."))
  save( f, file=fn, compress=TRUE)
  load(fn)
}

#nimdata = nscovid[c("Iobs", "Sobs", "Robs", "Mobs")]
nscovid$Iobs[ which(nscovid$Iobs < 0) ] = NA
nscovid$Sobs[ which(nscovid$Sobs < 0) ] = NA
nscovid$Robs[ which(nscovid$Robs < 0) ] = NA
nscovid$Mobs[ which(nscovid$Mobs < 0) ] = NA

nscovid$dSI = -diff( nscovid$Sobs )
nscovid$dIR =  diff( nscovid$Robs )
nscovid$dIM =  diff( nscovid$Mobs )

vconstants = nscovid[
    c("Npop", "Nobs", "Npreds", "BNP", "BETA_prior", "GAMMA_prior", "EPSILON_prior" )]

    # ,
    # "Iobs", "Sobs", "Robs", "Mobs") ]

vinits = list(
#  dSImu=rep(1, nscovid$Nobs-1),
#  dIRmu=rep(1, nscovid$Nobs-1),
 # dIMmu=rep(1, nscovid$Nobs-1),
  Smu=runif(nscovid$Nobs, 0.99, 1),
  Imu=runif(nscovid$Nobs, 0, 0.01),
  Rmu=runif(nscovid$Nobs, 0, 0.01),
  Mmu=runif(nscovid$Nobs, 0, 0.01),
  Ssd=runif(1, 0, 0.25),
  Isd=runif(1, 0, 0.25),
  Rsd=runif(1, 0, 0.25),
  Msd=runif(1, 0, 0.25),
  BETA=runif(nscovid$Nobs-1, 0, 1) ,
  GAMMA=runif(1, 0, 0.5),
  EPSILON=runif(1, 0, 0.5) ,
  ar1=runif(1, -1,1),
  ar1sd=runif(1, 0, 0.1),
  ar1k=runif(1, -1, 1)
)

# vdata = list( Y=as.matrix(cbind( nscovid$dSI, nscovid$dIR, nscovid$dIM )) )
vdata = list(  Y = cbind(
  nscovid[["Sobs"]] / nscovid[["Npop"]] ,
  nscovid[["Iobs"]] / nscovid[["Npop"]] ,
  nscovid[["Robs"]] / nscovid[["Npop"]] ,
  nscovid[["Mobs"]] / nscovid[["Npop"]]
))


v <- nimbleCode({

  # basic sir in latent state space form
  ar1 ~  T( dnorm( mean=0, sd=1), -1, 1)
  ar1k ~ T( dnorm( mean=0, sd=1), -1, 1)
  ar1sd ~ T( dt( mu=0, tau=1, df=1), 1e-9, 1)  # student-t with df=1 is approximately a cauchy distribution

  for (i in 1:4) {
    Lsd[i] ~ T( dt( mu=0, tau=1, df=1), 1e-9, 0.25)  # student-t with df=1 is approximately a cauchy distribution
  }

  GAMMA ~ T( dnorm( GAMMA_prior, sd=GAMMA_prior ), 1e-9, 1) # recovery of I ... always < 1
  EPSILON ~ T( dnorm( EPSILON_prior, sd=EPSILON_prior ), 1e-9, 1) # recovery of I ... always < 1

  # pr_ir <- 1/GAMMA;
  # pr_ir <- 1.0 - exp(-GAMMA)
  # pr_im <- 1.0 - exp(-EPSILON)

  BETA[1]  ~ T( dnorm( BETA_prior, sd=BETA_prior ), 1e-9, 1) # recovery of I ... always < 1
  for (i in 1:(Nobs-2)){
    # BETA[i+1] ~ T( dnorm( mean=BETA_prior, sd=1 ), 0, )
    BETA[i+1] ~ T( dnorm(  ar1k + ar1*BETA[i], sd=ar1sd ), 1e-9, 1)
  }
  BETAproj ~ dnorm( mean(BETA[(Nobs-1-BNP):(Nobs-1)] ), sd=sd(BETA[(Nobs-1-BNP):(Nobs-1)] ) )

  for (i in 1:4) {
    Lmu[1,i] ~ T( dnorm( 1, sd=Lsd[i] ), 0, 1)
  }

  for (i in 1:(Nobs-1)){
    # pr_si[i] <- 1.0 - exp(-BETA[i+1] * Lmu[i,1] / Npop ); # approximation
    # dSLmu[i,2] ~ dbin( size=Lmu[i,1], prob=pr_si[i]); # prob of being infected
    # dILmu[i,3] ~ dbin( size=Lmu[i,2], prob=pr_ir); # prob of being infected
    # dILmu[i,4] ~ dbin( size=Lmu[i,2], prob=pr_im); # prob of being infected
    Lmu[i+1,1] ~ T( dnorm( Lmu[i,1] - Lmu[i,1] * Lmu[i,2]* BETA[i], sd=Lsd[1] ), 0, 1 )
    Lmu[i+1,2] ~ T( dnorm( Lmu[i,2] + Lmu[i,1] * Lmu[i,2]* BETA[i] - Lmu[i,2] * GAMMA  -Lmu[i,2] * EPSILON, sd=Lsd[2]), 0, 1 )
    Lmu[i+1,3] ~ T( dnorm( Lmu[i,3] + Lmu[i,2] * GAMMA, sd=Lsd[3] ), 0, 1 )
    Lmu[i+1,4] ~ T( dnorm( Lmu[i,4] + Lmu[i,2] * EPSILON, sd=Lsd[4] ), 0, 1 )
  }

  for(i in 1:(Nobs)){
    for (j in 1:4) {
      Y[i,j] ~ dnorm( Lmu[i,j], sd=Lsd[j]  )
    }
  }

})


vm <- nimbleModel( code=v, constants=vconstants )
vm$setData( vdata  )
vm$setInits( vinits )

mc <- buildMCMC(vm, monitors = c('BETA','GAMMA', 'EPSILON','Lmu', "Lsd", "ar1", "ar1k", "ar1sd"))
vc <- compileNimble(vm, mc, showCompilerOutput = TRUE )
vc$mc$run(1000, nburnin = 100, thin=10 )  # quick test

vc$mc$mvSaved["Lmu"]
vc$mc$mvSaved["pr_si"]
vc$mc$mvSaved["ar1"]
vc$mc$mvSaved["ar1sd"]
vc$mc$mvSaved["ar1k"]


vc$mc$run(500000, nburnin = 10000, thin=100 )

vss <- as.matrix(vc$mc$mvSamples)


if (0) {
  vm$getNodeNames()
  vm$plotGraph()
  vm$getDependencies(c("GAMMA", "BETA"))
  vm$getLogProb("BETA")
  # vm$calculate( vm$getDependencies(c("Imu")))

  # check that values makes sense
  vn = "GAMMA";   vm$simulate(vn); vm[[vn]]
  vn = "EPSILON";   vm$simulate(vn); vm[[vn]]
  vn = "ar1";   vm$simulate(vn); vm[[vn]]
  vn = "ar1k";   vm$simulate(vn); vm[[vn]]
  vn = "ar1sd";   vm$simulate(vn); vm[[vn]]

  vn = "pr_ir";   vm$simulate(vn); vm[[vn]]
  vn = "pr_im";   vm$simulate(vn); vm[[vn]]
  vn = "S0";   vm$simulate(vn); vm[[vn]]
  vn = "I0";   vm$simulate(vn); vm[[vn]]
  vn = "R0";   vm$simulate(vn); vm[[vn]]
  vn = "M0";   vm$simulate(vn); vm[[vn]]

  # simulate to check values
  for (i in 1:nscovid$Nobs) {
     vn = "BETA";   vm$simulate(vn); vm[[vn]]
    vn = "pr_si";   vm$simulate(vn); vm[[vn]]

    vn = "dSImu";   vm$simulate(vn); vm[[vn]]
    vn = "dIRmu";   vm$simulate(vn); vm[[vn]]
    vn = "dIMmu";   vm$simulate(vn); vm[[vn]]

    vn = "Smu";   vm$simulate(vn); vm[[vn]]
    vn = "Imu";   vm$simulate(vn); vm[[vn]]
    vn = "Rmu";   vm$simulate(vn); vm[[vn]]
    vn = "Mmu";   vm$simulate(vn); vm[[vn]]
  }

}

vss <- as.matrix(vc$mc$mvSamples)

## examine samples array, determine it hasn't converged yet...
if (0) {
  library(coda)
  vcoda <- as.mcmc(vss)
  vn = "BETA"
  vn = "GAMMA"
  plot(autocorr(vcoda[,vn]))
  quantile(vss[,vn], c(.025,.975))
  HPDinterval(as.mcmc(vss[,vn]), prob = .95)
  plot(vss[,vn], main = "Trace and Density")
  plot(density(vss[,'BETA']))

  vss <- as.matrix(vss)
  raftery.diag(vss)
  effectiveSize(vss)

  vc$mc$run(500000,  nburnin = 10000, thin=100 , reset=FALSE) ## continue that same run of the MCMC for more iterations:

vss <- as.matrix(vc$mc$mvSamples)

  vs <- nimbleMCMC(code=v,
    constants = vconstants,
    data = vdata,
    inits = vinits,
    nchains = 2, niter = 1000,
    summary = TRUE, WAIC = FALSE,
    monitors = c('BETA','Imu','Rmu', "Smu", "Mmu")
  )

  plot(vs$summary$all.chains[ grepl("Imu", rownames( vs$summary$all.chains)),"Mean"])

  # vs$WAIC

  vmConf <- configureMCMC(vm, print = TRUE)
  vmMCMC <- buildMCMC(vmConf)
  vmc <- compileNimble(vmMCMC, project=v )

  vss <- runMCMC(vmc, niter = 1000)

  vmc <- compileNimble(vm )
  vn="BETA"; vmc[[vn]]

}

o =colMeans(vss)
plot( o[grepl("Smu", names(o))] )
plot( o[grepl("Imu", names(o))] )
plot( o[grepl("BETA", names(o))] )


vn = "BETA"

plot(vss[ , vn], type = "l", xlab = "iteration", ylab = expression(beta))
acf(vss[, vn])


# extract simulated data & plot to see results
df <- tibble(t = seq_along(vm$Imu), Smu=vm$Smu)
df %>% ggplot(aes(t,Smu)) + geom_point()




outdir = file.path( "~/bio/adapt/inst/doc/")
nx = nscovid$Nobs + nscovid$Npreds - 1


io = nscovid$Iobs
io[io < 0] = NA

png(filename = file.path(outdir, "fit_with_projections_infected.png"))
  xrange = c(0, nx)
  yrange = c(0, max(M$I[, 1:nx]))
  plot( io ~ nscovid$time, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Days (day 1 is 2020-03-17)", type="n" )
  lines( apply(M$I, 2, median)[1:nx] ~ seq(1,nx), lwd =3, col="slateblue" )
  lines( apply(M$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  lines( apply(M$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
  points( io ~ nscovid$time, col="darkgray", cex=1.2 )
  abline( v=nscovid$time[nscovid$Nobs], col="grey", lty="dashed" )
  abline( v=nscovid$time[time_distancing], col="orange", lty="dotted" )
  abline( v=nscovid$time[time_relaxation], col="green", lty="dotted" )
  legend( "topleft", "", "\n       [-- Social distancing -->", bty="n" )
  legend( "top", "", "\n                                  [-- Parks open -->", bty="n" )
  legend( "topright", "", paste( "Current date: ", nscovid$timestamp ), bty="n")
dev.off()


ro = nscovid$Robs
ro[ro < 0] = NA
png(filename = file.path(outdir, "fit_with_projections_recovered.png"))
  xrange = c(0, nx)
  yrange = c(0, max(M$R[, 1:nx]))
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
  yrange = c(0, max(M$M[, 1:nx]))
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





