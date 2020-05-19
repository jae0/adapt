
simulate = function( M, istart=istart, BNP=0, nsims=1, nprojections=10, model="stochastic.simulation.sir" ) {

  require(SimInf)

  if (model=="stochastic.simulation.sir") {
    nsims = min( nsims, nrow(M$S) )
    iss = sample.int( nrow(M$S), nsims )
    sim = array( NA, dim=c(nsims, 3, nprojections) )
    if ( BNP > 0 ) {
      ibnp = c( (istart-BNP):(istart) )
      u0=data.frame(
        S= M$S[iss, istart],
        I= M$I[iss, istart],
        R= M$R[iss, istart] + M$M[iss, istart],
        beta= rowMeans(M$BETA[iss, ibnp]),
        gamma=M$GAMMA[iss]
      )
    } else {
      ibnp = istart
      u0=data.frame(
        S= M$S[iss, istart],
        I= M$I[iss, istart],
        R= M$R[iss, istart] + M$M[iss, istart],
        beta= M$BETA[iss, ibnp],
        gamma=M$GAMMA[iss]
      )
    }

    for (i in 1:nsims) {
      sim[i,,] = run( SIR(
        u0=u0[i,c("S","I","R")], tspan=1:nprojections,
        beta=u0$beta[i], gamma=u0$gamma[i] )
      )@U[]
    }
    return( sim )
  }
}


