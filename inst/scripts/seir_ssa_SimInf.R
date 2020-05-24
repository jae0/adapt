

# param estimation via STAN

BYM of params
BYM of status

# forecasts mech with SimInf
-- nf days forecasts
for each au {
  -- sample from param dists and simulate
}

BYM of forecasts

remotes::install_github("stewid/SimInf")

  ## Create a SEIR model object.
  ## SimInf ...
  require(SimInf)
     model <- SEIR(u0 = data.frame(S = 99, E = 0, I = 1, R = 0),
                   tspan = 1:500,
                   beta = 0.16,
                   epsilon = 0.25,
                   gamma = 0.077)

     ## Run the SEIR model and plot the result.
     #set.seed(3)
     result <- run(model)
     plot(result)


  SIR(u0, tspan, events = NULL, beta = NULL, gamma = NULL)

     S -- beta S I / N --> I

                             I -- gamma I --> R
   require(SimInf)

  ## Create an SIR model object.
     model <- SIR(u0 = data.frame(S = 99, I = 1, R = 0),
                  tspan = 1:100,
                  beta = 0.16,
                  gamma = 0.077)

     ## Run the SIR model and plot the result.
     set.seed(22)
     result <- run(model)
     plot(result)
