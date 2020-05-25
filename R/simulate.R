
simulate = function( M, istart=istart, nsims=1, nprojections=10, nthreads=1, model="stochastic.simulation.sir" ) {

  require(SimInf)

  if (model=="stochastic.simulation.sir") {

    ST = c(
      "S -> BETA*S*I/(S+I+R+M) -> I" ,
      "I -> GAMMA*I -> R",
      "I -> EPSILON*I -> M"
    )

    SC = c( "S", "I", "R", "M" )

    nsims = min( nsims, nrow(M$S) )
    iss = sample.int( nrow(M$S), nsims )
    sim = array( NA, dim=c(nsims, length(SC), nprojections) )

    u0=data.frame(
      S= M$S[iss, istart],
      R= M$R[iss, istart],
      M= M$M[iss, istart],
      BETA= M$BETA[iss, istart-1],  # note, BETA is conditioned on previous time step. .
      GAMMA=M$GAMMA[iss],
      EPSILON=M$EPSILON[iss]
    )

    if (exists("Q", M) ) {
      u0$I = trunc( M$I[iss, istart] * M$Q[iss, istart-1] )
    } else {
      u0$I = M$I[iss, istart]
    }

    for (i in 1:nsims) {

      sim[i,,]  = run( model=mparse(
        transitions = ST,
        compartments = SC,
        gdata = c( BETA=u0$BETA[i], GAMMA=u0$GAMMA[i], EPSILON=u0$EPSILON[i] ),
        u0 = u0[i, SC],
        tspan = 1:nprojections
        ), threads=nthreads
      )@U[]

    }
    return( sim )
  }

  if (0) {
    # events :

    statevars = c( "S", "I", "R", "M", "Q" )
    event_types = c( "external_source", "external_sink",  "internal_transfer")

    E = matrix( c(
        1, 1, 1,
        0, 1, 1,
        0, 1, 1,
        0, 0, 0,
        0, 1, 0
      ),
      ncol=3, nrow=5, byrow=TRUE,
      dimnames=list( statevars, event_types )
    )

    N = matrix( c( 4,3,2,0,0 ), ncol=1, dimnames=list(statevars, "quarantined") ) # move x forward, i.e, S->Q

    quarantine = data.frame(
      event="intTrans",
      time=rep(21:52, each=50),
      node=1:100,
      dest=0,
      n=0,
      proportion=0.5,
      select=3,  # select from 3rd column of E
      shift=1  # move to compartment specified by first column of N
    )



    transitions = c(
      "S -> BETA*S*I/(S+I+R+M+Q) -> I" ,
      "I -> GAMMA*I -> R",
      "I -> EPSILON*I -> M"
    )

    params = c( BETA=0.3, GAMMA=0.1, EPSILON=0.01)

    initial_conditions = data.frame( S=rep(99,n), I=rep(1,n), R=rep(0,n), M=rep(0,n) )

    SIRMQ = mparse(
      transitions= transitions,
      compartments = statevars,  # susceptible, infected, recovered, mortality, quarantined
      gdata = params,
      u0 = initial_conditions,
      events = events, E=E, N=N,
      tspan = 1:100
    )

    set.seed(1234)
    res = run( model=SIRMQ, threads=1 )
    plot(res)


    o = trajectories(res)


      ###

      # source: https://rpubs.com/bbolker/SIRgillespie

      # SiRext

      # Ben Bolker
      # 6 October 2015

      # I wanted to run some continuous-time, stochastic simulations of the basic SIR model. What’s here is stuff that I and others have probably written hundreds of versions of over the years, but I thought it might be interesting to others as a worked example. Here I’m using the Gillespie algorithm, the simplest brute-force method for discrete-state, continuous-time stochastic simulation. The Gillespie algorithm is implemented in many other places, e.g. the GillespieSSA package for R, and in Darren Wilkinson’s smfsb (“Stochastic Modeling for Systems Biology”) package, along with other more sophisticated/efficient variations such as tau-leap algorithms, but this is the quick and dirty version.



      library("plyr")     ## for rdply()
      library("reshape2") ## for melt()
      library("emdbook")  ## for lambertW()
      library("ggplot2"); theme_set(theme_bw())

      # Functions for computing the event rates and the transitions to be executed when the events occur:

      ratefun <- function(x,p) {
          with(as.list(c(x,p)),{
                  c(inf=beta*S*I/N,  ## scale inf by pop size
                    recover=gamma*I)
              })
      }

      transfun <- function(x,w) {
          switch(w,
                x + c(-1,1),   ## infection: S-1, I+1
                x + c(0,-1))   ## removal: I-1
      }

      # A wrapper function to run the simulation with specified parameters/starting values and return either the ending state or a matrix of event times and types:

      run <- function(p=c(beta=2,gamma=1,N=100),
                      I0=1,
                      itmax=1e5,
                      ret=c("final","all")) {
          ret <- match.arg(ret)
          if (ret=="all") {
              rmat <- matrix(NA,nrow=itmax,ncol=2,
                            dimnames=list(NULL,c("t","trans")))
          }
          x <- c(S=unname(p["N"])-I0,I=I0)
          it <- 1
          t <- 0
          trans <- c(0,0)
          while (x["I"]>0 & it<itmax) {
              r <- ratefun(x,p)
              t <- t+rexp(1,rate=sum(r))
              w <- sample(length(r),size=1,prob=r)
              x <- transfun(x,w)
              if (ret=="all") rmat[it,] <- c(t,w)
              it <- it+1
          }
          if (ret=="all") return(rmat[!is.na(rmat[,1]),])
          return(c(S=unname(x["S"]),t=t,it=it))
      }

      # In this particular case the epidemic dies out early, after 2 infections and 3 recoveries (we started with a single infected individual …)

      set.seed(101)
      ex0 <- run(ret="all")
      plot(trans~t,data=ex0)

      ## can use .progress="text" if running interactively
      ex1 <- rdply(1e3,run())
      ex2 <- rdply(1e3,run(p=c(beta=1.1,gamma=1,N=100)))

      # Some handy functions:

      ## convenience function: we only want to keep final
      ## fraction unaffected, time to extinction ...

      mm <- function(x) {
          melt(x[c("S","t")],id.var=character(0))
      }

      ## analytic computation of expected final size
      ## from ?lambertW

      finalsize <- function(R0) {
          1+1/R0*lambertW(-R0*exp(-R0))
      }

      # Results with R0=2R0=2 (default) and R0=1.1R0=1.1:

      (g0 <- ggplot(mm(ex1),aes(x=value))+geom_histogram()+
          facet_wrap(~variable,scale="free"))

      finalsize(2)  ## check final size (should add to plot)

      ## [1] 0.7968121

      g0 %+% mm(ex2)




      ----

      # SIR via GillespieSSA

      library(GillespieSSA)
      library(reshape2)
      # set.seed(42)

      a <- c("beta*S*I","gamma*I")
      nu <- matrix(c(-1,+1,0,0,-1,+1),nrow=3,ncol=2,byrow=FALSE)

      parms <- c(beta=0.1/1000,gamma=0.05)
      x0 <- c(S=999,I=1,R=0)
      tf <- 200
      sir_out <- ssa(x0,a,nu,parms,tf=tf,simName="SIR")
      while( sir_out$stats$nSteps==1 ){
          sir_out <- ssa(x0,a,nu,parms,tf=tf,simName="SIR")
      }

      head(sir_out$data)
      plot(sir_out$data[,c(1,3)])



      ----

      # Input parameters ####################

      # int; total population

      N = 500
      T = 100.0 # float; maximum elapsed time
      t = 0.0 # float; start time
      V = 100.0 # float; spatial parameter
      alpha = 1.1 # float; rate of infection after contact
      beta = 0.5 # float; rate of cure
      n_I = 1 # int; initial infected population
      # Compute susceptible population, set recovered to zero
      n_S = N - n_I
      n_R = 0
      # Initialize results list
      simtot = 1000
      SIR = matrix(NA, ncol=4, nrow=simtot)
      SIR[1,] = c(t, n_S, n_I, n_R)
      # Main loop
      for (i in 2:simtot) {
          if (t >= T) break
        if (n_I == 0) break
        w1 = alpha * n_S * n_I / V
        w2 = beta * n_I
        W = w1 + w2
        dt = -log(runif(1)) / W
        t = t + dt
        if (runif(1) < (w1 / W) ) {
            n_S = n_S - 1
          n_I = n_I + 1
        } else {
            n_I = n_I - 1
          n_R = n_R + 1
        }
        SIR[i,] = c(t, n_S, n_I, n_R)
      }
      SIR
      plot(SIR[,2] ~ SIR[,1])

  }
}


