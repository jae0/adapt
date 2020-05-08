data_nova_scotia = function( output="stan_data", Npop=971395, Npreds=5, interpolate_missing_data=FALSE, ... ) {

  # install.packages("googlesheets4")
  library(googlesheets4)

  gsdata = read_sheet( "https://docs.google.com/spreadsheets/d/1tgf2H9gDmRnGDGeQE-fC9IPhrmNxP8-JC7Nnnob_vuY/edit?usp=sharing" )

  gsdata = gsdata[ is.finite(gsdata$InfectedCurrently), ]

  gsdata$Iobs = gsdata$InfectedCurrently
  gsdata$Robs = gsdata$Recoveries + gsdata$Deaths  # note: recovered = deaths+recoveries
  gsdata$Sobs = Npop - gsdata$Robs - gsdata$Iobs
  gsdata$Mobs = gsdata$Deaths   # mortalities .. this is redundant here but in some models, recoveries is split apart from mortalities and so useful

  gsdata$dayno = lubridate::date( gsdata$Date)
  gsdata$dayno = gsdata$dayno - min(gsdata$dayno) + 1

  if (output=="raw_data") return (gsdata)

  daily = expand.grid( dayno=1:max(gsdata$dayno ))
  daily[, c("Iobs", "Sobs", "Robs", "Mobs")] = NA
  i = match( gsdata$dayno, daily$dayno )
  daily[i, c("Sobs", "Iobs", "Robs", "Mobs")] = gsdata[,c("Sobs", "Iobs",  "Robs", "Mobs")]


  if (interpolate_missing_data) {
  # these are cummulative counts .. linear approximation where missing

    j = which( !is.finite(daily$Robs) )
    if (length(j) > 0) {
      o = approx( x=daily$dayno , y=daily$Robs, xout = daily$dayno, method="linear")
      daily$Robs[j] = trunc(o$y[j])
    }

    j = which( !is.finite(daily$Mobs) )
    if (length(j) > 0) {
      o = approx( x=daily$dayno , y=daily$Mobs, xout = daily$dayno, method="linear")
      daily$Mobs[j] = trunc(o$y[j])
    }

    j = which( !is.finite(daily$Sobs) )
    if (length(j) > 0) {
      o = approx( x=daily$dayno , y=daily$Sobs, xout = daily$dayno, method="linear")
      daily$Sobs[j] = trunc(o$y[j])
    }

    j = which( !is.finite(daily$Iobs) )
    if (length(j) > 0) daily$Iobs[j] = Npop - daily$Sobs[j] - daily$Robs[j] - daily$Mobs[j]
  }

  # final check .. set missing values as -1
  j = which( !is.finite(daily$Sobs) ); if (length(j) > 0) daily$Sobs[j] = -1
  j = which( !is.finite(daily$Iobs) ); if (length(j) > 0) daily$Iobs[j] = -1
  j = which( !is.finite(daily$Robs) ); if (length(j) > 0) daily$Robs[j] = -1
  j = which( !is.finite(daily$Mobs) ); if (length(j) > 0) daily$Mobs[j] = -1

  if (output=="daily_data") return (daily)


  stan_data = list(
    Npop = Npop,
    Nobs = nrow( daily ),
    Npreds = Npreds,   # here, only total number of predictions for output
    Sobs = daily$Sobs,
    Iobs = daily$Iobs,
    Robs = daily$Robs,
    Mobs = daily$Mobs,
    time = as.integer(daily$dayno)
  )


  stan_data = c( stan_data, list(...) )

  if (!exists("modelname", stan_data)) stan_data$modelname="discrete_autoregressive_with_observation_error_structured_beta_mortality"
  stan_data$timestamp = max( lubridate::date( gsdata$Date) )

  if ( stan_data$modelname %in% c("continuous") ) {
    # used by ODE-based methods for rk4 integration
    if (!exists("time_pred", stan_data)) stan_data$time_pred = as.integer( c(daily$dayno, max(daily$dayno)+c(1:Npreds)) )
    if (!exists("t0", stan_data)) stan_data$t0 = -0.01
  }

  # add a few more flags for discrete_variable_encounter_rate and
  if (!exists("BETA_prior", stan_data)) stan_data$BETA_prior = 0.5
  if (!exists("GAMMA_prior", stan_data)) stan_data$GAMMA_prior = 1/28
  if (!exists("EPSILON_prior", stan_data)) stan_data$me_prior = 0.05  # % of positively Infected that die
  if (!exists("BNP", stan_data)) stan_data$BNP = 3 # number of days to use for beta averaging for projections

  if (output=="stan_data") return (stan_data)

}