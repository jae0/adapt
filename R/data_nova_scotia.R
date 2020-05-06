data_nova_scotia = function( output="stan_data", Npop=971395, Npreds=5, ... ) {

  # install.packages("googlesheets4")
  library(googlesheets4)

  gsdata = read_sheet( "https://docs.google.com/spreadsheets/d/1tgf2H9gDmRnGDGeQE-fC9IPhrmNxP8-JC7Nnnob_vuY/edit?usp=sharing" )

  gsdata = gsdata[ is.finite(gsdata$InfectedCurrently), ]

  gsdata$Iobs = gsdata$InfectedCurrently
  gsdata$Robs = gsdata$Recoveries + gsdata$Deaths
  gsdata$Sobs = Npop - gsdata$Robs - gsdata$Iobs
  gsdata$Dobs = gsdata$Deaths

  gsdata$dayno = lubridate::date( gsdata$Date)
  gsdata$dayno = gsdata$dayno - min(gsdata$dayno) + 1

  if (output=="raw_data") return (gsdata)

  daily = expand.grid( dayno=1:max(gsdata$dayno ))
  daily[, c("Iobs", "Sobs", "Robs", "Dobs")] = NA
  i = match( gsdata$dayno, daily$dayno )
  daily[i, c("Sobs", "Iobs", "Robs", "Dobs")] = gsdata[,c("Sobs", "Iobs",  "Robs", "Dobs")]

  # these are cummulative counts .. linear approximation where missing
  j = which( !is.finite(daily$Robs) )
  if (length(j) > 0) {
    o = approx( x=daily$dayno , y=daily$Robs, xout = daily$dayno, method="linear")
    daily$Robs[j] = trunc(o$y[j])
  }

  j = which( !is.finite(daily$Sobs) )
  if (length(j) > 0) {
    o = approx( x=daily$dayno , y=daily$Sobs, xout = daily$dayno, method="linear")
    daily$Sobs[j] = trunc(o$y[j])
  }

  j = which( !is.finite(daily$Iobs) )
  if (length(j) > 0) daily$Iobs[j] = Npop - daily$Sobs[j] - daily$Robs[j]

  # final check
  j = which( !is.finite(daily$Sobs) ); if (length(j) > 0) daily$Sobs[j] = -1
  j = which( !is.finite(daily$Iobs) ); if (length(j) > 0) daily$Iobs[j] = -1
  j = which( !is.finite(daily$Robs) ); if (length(j) > 0) daily$Robs[j] = -1

  if (output=="daily_data") return (daily)


  stan_data = list(
    Npop = Npop,
    Nobs = nrow( daily ),
    Npreds = Npreds,   # here, only total number of predictions for output
    Sobs = daily$Sobs,
    Iobs = daily$Iobs,
    Robs = daily$Robs,
    time = as.integer(daily$dayno),
    time_pred = as.integer( c(daily$dayno, max(daily$dayno)+c(1:Npreds)) ) ,
    t0 = -0.01
  )
  stan_data = c( stan_data, list(...) )

  if (!exists("modelname", stan_data)) stan_data$modelname="discrete_autoregressive_with_observation_error"

  # add a few more flags for discrete_variable_encounter_rate and
  if (!exists("BETA_prior", stan_data)) stan_data$BETA_prior = 0.5
  if (!exists("GAMMA_prior", stan_data)) stan_data$GAMMA_prior = 1/28
  if (!exists("me_prior", stan_data)) stan_data$me_prior = 0.05  # % of Infected that are asymptomatic
  if (!exists("BNP", stan_data)) stan_data$BNP = 3 # % of Infected that are asymptomatic

  if (output=="stan_data") return (stan_data)

}