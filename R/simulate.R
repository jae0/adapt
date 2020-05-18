
simulate = function( M, istart=istart, nsims=1, nprojections=10, model="stochastic.simulation.sir" ) {

  require(SimInf)

  if (model=="stochastic.simulation.sir") {
    nsims = min( nsims, nrow(M$S) )
    iss = sample.int( nrow(M$S), nsims )
    sim = array( NA, dim=c(nsims, 3, nprojections) )
    u0=data.frame(
      S=M$S[iss,istart],
      I=M$I[iss,istart],
      R=M$R[iss,istart] + M$M[iss,istart],
      beta=M$BETA[iss,istart],
      gamma=M$GAMMA[iss]
    )
    for (i in 1:nsims) {
      sim[i,,] = run( SIR(
        u0=u0[i,c("S","I","R")], tspan=1:nprojections,
        beta=u0$beta[i], gamma=u0$gamma[i] )
      )@U[]
    }
    return( sim )
  }
}


