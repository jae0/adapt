data_provinces_of_canada = function( selection="default", fn=NULL, Npreds=5, ... ) {

  # this accesses data from : https://github.com/ishaberry/Covid19Canada
  # Citation:  Berry I, Soucy J-PR, Tuite A, Fisman D. Open access epidemiologic data and an interactive dashboard to monitor the COVID-19 outbreak in Canada. CMAJ. 2020 Apr 14;192(15):E420. doi: https://doi.org/10.1503/cmaj.75262

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19Canada.rdata")  # default to current work directory

  outdir = dirname(fn)
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )

  fn_pop = file.path( dirname(fn), "POP_Canada.rdata")

  if (selection =="download_pop") {
    pop = read.csv(curl(
    "https://www150.statcan.gc.ca/t1/tbl1/en/dtl!downloadDbLoadingData.action?pid=1710000901&latestN=0&startDate=20200101&endDate=20201001&csvLocale=en&selectedMembers=%5B%5B2%2C3%2C4%2C5%2C6%2C7%2C8%2C9%2C10%2C11%2C12%2C14%2C15%5D%5D&checkedLevels=0D1" ), header=TRUE)
    pop = pop[, c( "REF_DATE", "GEO", "VALUE")]
    names(pop) =c("date", "province", "Npop")
    pop$au = c(
      "Canada",
      "NL",
      "PEI",
      "Nova Scotia",
      "New Brunswick",
      "Quebec",
      "Ontario",
      "Manitoba",
      "Saskatchewan",
      "Alberta",
      "BC",
      "Yukon",
      "NWT",
      "Nunavut"
    )
    save( pop, file=fn_pop, compress=TRUE )  # default to current work directory
  }
  load( fn_pop )

  res = NULL
  if (selection =="download") {
    library(curl)
    cases = read.csv(curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/cases.csv"), header=TRUE)
    mortality = read.csv(curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/mortality.csv"), header=TRUE)
    recovered = read.csv(curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/recovered_cumulative.csv"), header=TRUE)

    res = list(timestamp=Sys.Date(), cases=cases, mortality=mortality, recovered=recovered)
    save( res, file=fn, compress=TRUE )
    return(res)
  }

  load(fn)

  if (!file.exists(fn)) res = data_provinces_of_canada( selection="download", fn=fn )
  if (!exists("res")) res = data_provinces_of_canada( selection="download", fn=fn )
  if (res$timestamp != Sys.Date() )  res = data_provinces_of_canada( selection="download", fn=fn )

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

  au = sort( unique( res$cases$province))
  au = setdiff(au,  "Repatriated")  # not clear what to do with these .. droppoing for now

  infected = infected[ infected$province %in% au , ]
  death = death[ death$province %in% au , ]
  recovered = recovered[ recovered$province %in% au , ]

  daily = array( NA, dim=c( length(daynos), length(au), 4 ) )  # SIRD (D=death)

  infected$dayno = infected$date - daterange[1] + 1
  recovered$dayno = recovered$date - daterange[1] + 1
  death$dayno = death$date - daterange[1] + 1

  daily[ cbind( match(infected$dayno, daynos ), match(infected$province, au), rep(2, nrow(infected) ) ) ] = infected$infected

  daily[ cbind( match(recovered$dayno, daynos ), match(recovered$province, au), rep(3, nrow(recovered) ) ) ] = recovered$recovered

  daily[ cbind( match(death$dayno, daynos ), match(death$province, au), rep(4, nrow(death) ) ) ] = death$death

  Npop = pop$Npop
  names(Npop) = pop$au

  if (selection=="daily_data") return (daily)

  Ntimes = dim(daily)[1]

  # fill start of TS with 0 until first case is observed
  for (k in 1:length(au)) {
  for (j in 1:Ntimes) {
    if (any(is.finite( daily[j,k,2:4] ))) break()
    daily[j,k,2:4] = 0
  }}

  totalRecoveries = daily[ ,,3]
  newRecoveries = totalRecoveries[] * 0
  totalRecoveries[ which( !is.finite(totalRecoveries))] = 0
  newRecoveries[2:Ntimes,] = totalRecoveries[ 2:Ntimes,] - totalRecoveries[ 1:(Ntimes-1),]
  newRecoveries[ which(newRecoveries < 0)] = 0  ## there is a typo in Alberta 19-05-2020       Alberta                 5854
  newRecoveries[ which( !is.finite(newRecoveries))] = 0

  # cummulative sums for Mortalities
  newDeaths = daily[ ,, 4]
  newDeaths[ which(!is.finite(newDeaths))] = 0
  newDeaths[ which(newDeaths < 0)] = 0  ## in case of data entry errors
  daily[ , , 4] = 0
  for (k in 1:length(au)) {
    daily[,k,4] = cumsum(newDeaths[,k])
  }

  # compute infecteds
  newInfecteds = daily[ ,, 2]
  newInfecteds[ which(!is.finite(newInfecteds))] = 0
  newInfecteds[ which(newInfecteds < 0)] = 0  ## in case of data entry errors
  daily[ , , 2] = 0
  for (k in 1:length(au)) {
  for (j in 1:(Ntimes-1)) {
    daily[j+1,k,2] = daily[j,k,2] + newInfecteds[j,k] - newRecoveries[j,k] - newDeaths[j,k]
  }}

  # compute Susceptibles
  # fill with initial pop
  daily[1,,1] = Npop[au]  # start of data
  for (k in 1:length(au)) {
  for (j in 1:(Ntimes-1) ) {
    daily[j+1,k,1] =  daily[j,k,1] - sum(daily[j,k,2:4], na.rm=TRUE )
  }}

  daily [ !is.finite(daily)] = -1

    # default is to return this:
  data_province = list()
  for ( i in 1:length(au) ) {
    prov = au[i]
    stan_data = NULL
    stan_data  = list(
      Npop = Npop[[ prov ]],
      Nobs = Ntimes,
      Npreds = Npreds,
      Sobs = daily[,i,1],
      Iobs = daily[,i,2],
      Robs = daily[,i,3],
      Mobs = daily[,i,4],
      daterange = daterange,
      daynos = daynos,
      au = prov,
      statevars = c("susceptible", "infected", "recovered", "dead"),  # names of last dim of "daily"
      time = daynos,
      time_start = daterange[1],
      timestamp = max( lubridate::as_date( daterange ) )
    )


    stan_data = c( stan_data, list(...) )
    if (!exists("modelname", stan_data)) stan_data$modelname="default"
    # add a few more flags for discrete_variable_encounter_rate and
    if (!exists("BETA_max", stan_data)) stan_data$BETA_max = 0.5
    if (!exists("GAMMA_max", stan_data)) stan_data$GAMMA_max = 0.5
    if (!exists("EPSILON_max", stan_data)) stan_data$EPSILON_max = 0.5
    if (!exists("BNP", stan_data)) stan_data$BNP = 5  # % of Infected that are asymptomatic

    if ( stan_data$modelname %in% c("continuous") ) {
      # used by ODE-based methods for rk4 integration
      if (!exists("time_pred", stan_data)) stan_data$time_pred = as.integer( c(daily$dayno, max(daily$dayno)+c(1:Npreds)) )
      if (!exists("t0", stan_data)) stan_data$t0 = -0.01
    }

    data_province[[prov]] = stan_data
  }
  return( data_province )
}