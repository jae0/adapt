

# param estimation via STAN

BYM of params
BYM of status

# forecasts mech with SimInf
-- nf days forecasts
for each au {
  -- sample from param dists and simulate
}

BYM of forecasts

# remotes::install_github("stewid/SimInf")

require(SimInf)

n = 1000

SIRM = mparse(
  transitions= c(
    "S -> BETA*S*I/(S+I+R+M) -> I" ,
    "I -> GAMMA*I -> R",
    "I -> EPSILON*I -> M"
  ),
  compartments = c( "S", "I", "R", "M" ),
  gdata = c( BETA=0.3, GAMMA=0.1, EPSILON=0.01),
  u0 = data.frame( S=rep(99,n), I=rep(1,n), R=rep(0,n), M=rep(0,n) ),
  tspan = 1:100
)

set.seed(1234)
res = run( model=SIRM, threads=1 )
plot(res)



quarantine = data.frame(
   event="intTrans",
   time=rep(21:52, each=50),
   node=1:100,
   dest=0,
   n=0,
   proportion=0.5,
   select=3,
   shift=1
)

# events :

statevars = c( "S", "I", "R", "M", "Q" )
event_types = c( "external_source", "external_sink",  "internal_transfer")

events = matrix( c(
    1, 1, 1,
    0, 1, 1,
    0, 0, 0,
    0, 1, 1,
    0, 1, 0
  ),
  ncol=3, nrow=5, byrow=TRUE,
  dimnames=list( statevars, event_types )
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
  tspan = 1:100
)

set.seed(1234)
res = run( model=SIRMQ, threads=1 )
plot(res)


o = trajectories(res)

