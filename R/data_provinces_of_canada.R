data_provinces_of_canada = function( output="stan_data", fn=NULL, Npreds=5, ... ) {

  # this accesses data from : https://github.com/ishaberry/Covid19Canada
  # Citation:  Berry I, Soucy J-PR, Tuite A, Fisman D. Open access epidemiologic data and an interactive dashboard to monitor the COVID-19 outbreak in Canada. CMAJ. 2020 Apr 14;192(15):E420. doi: https://doi.org/10.1503/cmaj.75262

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19Canada.rdata")  # default to current work directory

  if (output =="download") {
    library(curl)
    cases = read.csv(curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/cases.csv"), header=TRUE)
    mortality = read.csv(curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/mortality.csv"), header=TRUE)
    recovered = read.csv(curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/recovered_cumulative.csv"), header=TRUE)
    res = list(timestamp=Sys.Date(), cases=cases, mortality=mortality, recovered=recovered)
    save( res, file=fn, compress=TRUE )
  }

  if (!file.exists(fn)) data_provinces_of_canada( output="download", fn=fn )
  load(fn)

  if (!exists("res")) data_provinces_of_canada( output="download", fn=fn )
  if (res$timestamp != Sys.Date() )  data_provinces_of_canada( output="download", fn=fn )

  infected = as.data.frame.table( tapply( 1:nrow(res$cases), INDEX=list( date=res$cases$date_report, province=res$cases$province ), length ) )
  colnames(infected) = c("date", "province", "infected")
  infected$date = lubridate::dmy( infected$date )


  death = as.data.frame.table( tapply( 1:nrow(res$mortality), INDEX=list( date=res$mortality$date_death_report, province=res$mortality$province ), length ) )
  colnames(death) = c("date", "province", "death")
  death$date = lubridate::dmy( death$date )

  recovered = res$recovered
  colnames(recovered) = c("date", "province", "recovered")  # cummulative recovered
  recovered$date = lubridate::dmy( recovered$date )


  daterange = range( c(death$date, recovered$date, infected$date), na.rm=TRUE )
  daynos = 1:(diff(daterange)+1)


  au = sort( unique( cases$province))
  au = setdiff(au,  "Repatriated")  # not clear what to do with these .. droppoing for now

  infected = infected[ infected$province %in% au , ]
  death = death[ death$province %in% au , ]
  recovered = recovered[ recovered$province %in% au , ]

  gsdata = array( NA, dim=c( length(daynos), length(au), 4 ) )  # SIRD (D=death)

  infected$dayno = infected$date - daterange[1] + 1
  recovered$dayno = recovered$date - daterange[1] + 1
  death$dayno = death$date - daterange[1] + 1

  gsdata[ cbind( match(infected$dayno, daynos ), match(infected$province, au), rep(2, nrow(infected) ) ) ] = infected$infected

  gsdata[ cbind( match(recovered$dayno, daynos ), match(recovered$province, au), rep(3, nrow(recovered) ) ) ] = recovered$recovered

  gsdata[ cbind( match(death$dayno, daynos ), match(death$province, au), rep(4, nrow(death) ) ) ] = death$death


  if (output=="raw_data") return (gsdata)


  j = which( !is.finite(gsdata[,,2]) )
  if (length(j) > 0) {
    gsdata[j,2] = 0
  }

  j = which( !is.finite(gsdata$Robs) )
  if (length(j) > 0) {
    o = approx( x=gsdata$dayno , y=gsdata$Robs, xout = gsdata$dayno, method="linear")
    gsdata$Robs[j] = trunc(o$y[j])
  }

  j = which( !is.finite(gsdata$Sobs) )
  if (length(j) > 0) {
    gsdata$Sobs[j] = Npop - gsdata$Iobs[j] - gsdata$Robs[j]
  }

  gsdata$managementmeasures = "none"
  gsdata$managementmeasures[i] = gsdata$managementmeasures

  if (output=="daily_data") return (gsdata)


  stan_data = list(
    Npop = Npop,
    Nobs = nrow( gsdata ),
    Npreds = Npreds,   # here, only total number of predictions for output
    Sobs = gsdata$Sobs,
    Iobs = gsdata$Iobs,
    Robs = gsdata$Robs,
    time = as.integer(gsdata$dayno),
    time_pred = as.integer( c(gsdata$dayno, max(gsdata$dayno)+c(1:Npreds)) ) ,
    t0 = -0.01
  )
  stan_data = c( stan_data, list(...) )

  if (!exists("modelname", stan_data)) stan_data$modelname="discrete_basic"

  # add a few more flags for discrete_variable_encounter_rate and
  if (!exists("BETA_prior", stan_data)) stan_data$BETA_prior = 0.5
  if (!exists("GAMMA_prior", stan_data)) stan_data$GAMMA_prior = 0.5
  if (!exists("ER_prior", stan_data)) stan_data$ER_prior = 0.5
  if (!exists("me_prior", stan_data)) stan_data$me_prior = 0.1  # % of Infected that are asymptomatic
  if (!exists("BNP", stan_data)) stan_data$BNP = 5  # % of Infected that are asymptomatic

  if (output=="stan_data") return (stan_data)

}