#' @title data_health_regions_of_canada
#' @description This is a placeholder for a description.
#' @param selection default is \code{"default"}
#' @param fn default is \code{NULL}
#' @param Npreds default is \code{5}
#' @param ... other arguments passed to methods
#' @return  This is a placeholder for what it returns.
#' @author Jae Choi, \email{choi.jae.seok@gmail.com}
#' @export
data_health_regions_of_canada = function( selection="default", fn=NULL, Npreds=5, ... ) {

  # this accesses data from : https://github.com/ishaberry/Covid19Canada
# Citation  Berry I, Soucy J-PR, Tuite A, Fisman D. Open access epidemiologic data and an interactive dashboard to monitor the COVID-19 outbreak in Canada. CMAJ. 2020 Apr 14;192(15):E420. doi: https://doi.org/10.1503/cmaj.75262

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19CanadaHR.rdata")  # default to current work directory

  outdir = dirname(fn)
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )


  fn_pop = file.path( dirname(fn), "POP_Canada.rdata")

  if (selection =="download_pop") {
    pop = utils::read.csv(curl::curl(
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
  if (!file.exists(fn_pop)) {
    pop = data_health_regions_of_canada( selection="download_pop", fn=fn )
  } else {
    load( fn_pop )
  }


  fn_pop_by_health_region = file.path( dirname(fn), "POP_Canada_health_regions.rdata")
  if (selection =="download_pop_health_regions") {
    pop_hr = utils::read.csv(curl::curl(
    "https://www150.statcan.gc.ca/t1/tbl1/en/dtl!downloadDbLoadingData.action?pid=1710013401&latestN=0&startDate=20150101&endDate=20190101&csvLocale=en&selectedMembers=%5B%5B%5D%2C%5B%5D%2C%5B%5D%5D&checkedLevels=0D1%2C0D2%2C0D3%2C0D4%2C1D1%2C2D1" ), header=TRUE)
    pop_hr = pop_hr[, c( "REF_DATE", "GEO", "VALUE")]
    names(pop_hr) =c("date", "AU", "Npop")
    o = which( sapply( strsplit( pop_hr$AU, ","), length) == 2 )
    pop_hr = pop_hr[o,]
    au = strsplit( pop_hr$AU, ", ")
    pop_hr$province = sapply( au, function(x) x[[2]])
    pop_hr$health_region = sapply( au, function(x) x[[1]])
    pop_hr$AU = paste( pop_hr$province, pop_hr$health_region, sep=" __ ")
    save( pop_hr, file=fn_pop_by_health_region, compress=TRUE )  # default to current work directory
  }
  if (!file.exists(fn_pop_by_health_region)) {
    pop_hr = data_health_regions_of_canada( selection="download_pop_health_regions", fn=fn )
  } else {
    load(fn_pop_by_health_region )
  }
  pop_hr = pop_hr[ pop_hr$date== 2019, ]

  res = NULL
  if (selection =="download") {
    #library(curl)
    cases = utils::read.csv(curl::curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/cases.csv"), header=TRUE)
    mortality = utils::read.csv(curl::curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/mortality.csv"), header=TRUE)
    recovered = utils::read.csv(curl::curl("https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/recovered_cumulative.csv"), header=TRUE)

    res = list(timestamp=Sys.Date(), cases=cases, mortality=mortality, recovered=recovered)
    save( res, file=fn, compress=TRUE )
    return(res)
  }

  if (!file.exists(fn)) {
    res = data_health_regions_of_canada( selection="download", fn=fn )
  } else {
    load (fn )
  }

  if (res$timestamp != Sys.Date() )  res = data_health_regions_of_canada( selection="download", fn=fn )


  daterange = range( c(
    lubridate::dmy( res$cases$date_report ),
    lubridate::dmy( res$mortality$date_death_report ),
    lubridate::dmy( res$recovered$date_recovered )),
    na.rm=TRUE
  )

  tus = 1:(diff(daterange)+1)
  aus = as.character( sort( unique( paste( res$cases$province, res$case$health_region, sep=" __ " ) )) )
  provinces = setdiff( sort( unique( res$cases$province )), "Repatriated")


  # newly reported infections
  au = paste( res$cases$province, res$case$health_region, sep=" __ " )
  infected = as.data.frame.table( tapply( 1:nrow(res$cases), INDEX=list( date=res$cases$date_report, au=au ), length ), stringsAsFactors=FALSE  )
  colnames(infected) = c("date", "au", "infected")
  infected$date = lubridate::dmy( infected$date )
  infected = infected[ infected$au %in% aus , ]
  infected$tu = infected$date - daterange[1] + 1
  infected$tu_index = match( infected$tu, tus )
  infected$au_index = match( infected$au, aus )

  # new mortalities
  au = paste( res$mortality$province, res$mortality$health_region, sep=" __ " )
  death = as.data.frame.table( tapply( 1:nrow(res$mortality), INDEX=list( date=res$mortality$date_death_report, au=au ), length ), stringsAsFactors=FALSE  )
  colnames(death) = c("date", "au", "death")
  death$date = lubridate::dmy( death$date )
  death = death[ death$au %in% aus , ]
  death$tu = death$date - daterange[1] + 1
  death$tu_index = match(death$tu, tus )
  death$au_index = match(death$au, aus )


  # newly reported daily as a matrix
  I_daily = array( NA, dim=c( length(tus), length(aus)  ) )
  I_daily[ cbind( infected$tu_index, infected$au_index) ] = infected$infected

  # newly reported daily as a matrix
  D_daily = array( NA, dim=c( length(tus), length(aus)  ) )
  D_daily[ cbind( death$tu_index, death$au_index ) ] = death$death

  # fill start of TS with 0 until first case is observed
  for (k in 1:length(aus)) {
  for (j in 1:length(tus)) {
    if (any( is.finite( c( I_daily[j,k], D_daily[j,k]) ))) break()
    I_daily[j,k] = 0
    D_daily[j,k] = 0
  }}

  I_daily[ which(!is.finite(I_daily))] = 0
  I_daily[ which(I_daily < 0)] = 0  ## in case of data entry errors

  D_daily[ which(!is.finite(D_daily))] = 0
  D_daily[ which(D_daily < 0)] = 0  ## in case of data entry errors

  # convert daily new deaths to cummulative sums for Mortalities
  D_cumsum = D_daily[] * 0
  for (k in 1:length(aus)) {
    D_cumsum[,k] = cumsum(D_daily[,k])
  }


  # RECOVERED .... recovered data only exist by province but we want it by health region:
  # we can estimate by health region based upon the fraction of incidence in a health region to the total incidence in a province
  # use infected as the au X time template and fill with INF_LOCAL and INF_PROV


  # newly reported daily by province as a matrix
  infected_province = as.data.frame.table( tapply( 1:nrow(res$cases), INDEX=list( date=res$cases$date_report, province=res$cases$province ), length ), stringsAsFactors=FALSE  )
  colnames(infected_province) = c("date", "province", "infected_province")
  infected_province$date = lubridate::dmy( infected_province$date )
  infected_province = infected_province[ infected_province$province %in% provinces , ]
  infected_province$tu = infected_province$date - daterange[1] + 1
  infected_province$tu_index = match( infected_province$tu, tus )
  infected_province$province_index = match( infected_province$province, provinces )

  # new mortalities
  death_province = as.data.frame.table( tapply( 1:nrow(res$mortality), INDEX=list( date=res$mortality$date_death_report, province=res$mortality$province ), length ), stringsAsFactors=FALSE )
  colnames(death_province) = c("date", "province", "death_province")
  death_province$date = lubridate::dmy( death_province$date )
  death_province = death_province[ death_province$province %in% provinces , ]
  death_province$tu = death_province$date - daterange[1] + 1
  death_province$tu_index = match(death_province$tu, tus )
  death_province$province_index = match(death_province$province, provinces )


  recovered_province = res$recovered
  colnames(recovered_province) = c("date", "province", "recovered")  # cummulative recovered_province
  recovered_province$date = lubridate::dmy( recovered_province$date )
  recovered_province = recovered_province[ recovered_province$province %in% provinces , ]
  recovered_province$tu = recovered_province$date - daterange[1] + 1
  recovered_province$tu_index = match(recovered_province$tu, tus )
  recovered_province$province_index = match(recovered_province$province, provinces )


  I_daily_province = array( NA, dim=c( length(tus), length(provinces)  ) )
  I_daily_province[ cbind( infected_province$tu_index, infected_province$province_index) ] = infected_province$infected

  # newly reported daily as a matrix
  D_daily_province = array( NA, dim=c( length(tus), length(provinces)  ) )
  D_daily_province[ cbind( death_province$tu_index, death_province$province_index ) ] = death_province$death

  R_cumsum_province = array( NA, dim=c( length(tus), length(provinces)  ) )
  R_cumsum_province[ cbind( recovered_province$tu_index, recovered_province$province_index ) ] = recovered_province$recovered

  for (k in 1:length(provinces)) {
  for (j in 1:length(tus)) {
    if (any( is.finite( c( I_daily_province[j,k] )))) break()
    I_daily_province[j,k] = 0
    R_cumsum_province[j,k] = 0
    D_daily_province[j,k] = 0
  }}


  Ntimes = length(tus)
  R_cumsum_province[ which( !is.finite(R_cumsum_province))] = 0
  R_daily_province = R_cumsum_province[] * 0
  R_daily_province[2:Ntimes,] = R_cumsum_province[ 2:Ntimes,] - R_cumsum_province[ 1:(Ntimes-1),]
  R_daily_province[ which(R_daily_province < 0)] = 0  ## there is a typo in Alberta 19-05-2020       Alberta                 5854
  R_daily_province[ which( !is.finite(R_daily_province))] = 0

  # convert daily new deaths to cummulative sums for Mortalities
  D_daily_province[ which(!is.finite(D_daily_province))] = 0
  D_daily_province[ which(D_daily_province < 0)] = 0  ## in case of data entry errors
  D_cumsum_province = D_daily_province[] * 0
  for (k in 1:length(provinces)) {
    D_cumsum_province[,k] = cumsum(D_daily_province[,k])
  }

  I_daily_province[ which(!is.finite(I_daily_province))] = 0
  I_daily_province[ which(I_daily_province < 0)] = 0  ## in case of data entry errors
  I_current_province = I_daily_province[] * 0
  for (k in 1:length(provinces)) {
  for (j in 1:(length(tus)-1)) {
    I_current_province[j+1,k] = I_current_province[j,k] + I_daily_province[j,k]  - D_daily_province[j,k]  - R_daily_province[j,k]
  }}


  # cummulative infections by hr and province
  # compute infecteds
  I_cumsum = I_daily[] * 0
  for (k in 1:length(aus)) {
  for (j in 1:(length(tus)-1)) {
    I_cumsum[j+1,k] = I_cumsum[j,k] + I_daily[j,k]
  }}

  rownames(I_cumsum) = tus
  colnames(I_cumsum) = aus
  infected_hr_cumm = as.data.frame.table( I_cumsum, stringsAsFactors=FALSE )
  colnames(infected_hr_cumm) = c("tu", "au", "infected_hr")
  infected_hr_cumm$province = sapply( strsplit( infected_hr_cumm$au, " __ "), function(x) x[[1]])
  infected_hr_cumm$tu = as.numeric( infected_hr_cumm$tu )

  # cummulative infections by hr and province

  I_cumsum_province = I_daily_province[] * 0
  for (k in 1:length(provinces)) {
  for (j in 1:(length(tus)-1)) {
    I_cumsum_province[j+1,k] = I_cumsum_province[j,k] + I_daily_province[j,k]
  }}

  rownames(I_cumsum_province) = tus
  colnames(I_cumsum_province) = provinces
  infected_province_cumm = as.data.frame.table( I_cumsum_province, stringsAsFactors=FALSE )
  colnames(infected_province_cumm) = c("tu", "province", "infected_province")
  infected_province_cumm$province = sapply( strsplit( infected_province_cumm$province, " __ "), function(x) x[[1]])
  infected_hr_cumm$tu = as.numeric( infected_hr_cumm$tu )  # cummulative infections by hr and province


  # cummulative infections by hr and province
  infected_province_cumm$tu = as.numeric( infected_province_cumm$tu )  # cummulative infections by hr and province
  # infected_province_cumm$tu = infected_province_cumm$tu - 7  # offset day ~ recovery
  infected_hr_cumm = merge( infected_hr_cumm, infected_province_cumm, by=c( "tu", "province" ), all.x=TRUE, all.y=FALSE )
  infected_hr_cumm =  infected_hr_cumm[, c("tu", "au", "province", "infected_hr", "infected_province") ]
  infected_hr_cumm$fraction_infected = infected_hr_cumm$infected_hr / infected_hr_cumm$infected_province  # fraction relative to province on a given day

  recovered = res$recovered
  recovered$date = lubridate::dmy( recovered$date_recovered )
  recovered$tu =  recovered$date - daterange[1] + 1
  recovered$date = as.character( recovered$date )

  recovered = merge( infected_hr_cumm, recovered, by=c("tu", "province"), all.x=TRUE, all.y=FALSE )
  recovered$fraction_infected[ which( !is.finite( recovered$fraction_infected)) ] = 0

  recovered$recovered_estimated = round( recovered$fraction_infected * recovered$cumulative_recovered )
  recovered$au_index = match( recovered$au, aus )
  recovered$tu_index = match( recovered$tu, tus )

  R_cumsum = array( 0, dim=c( length(tus), length(aus)  ) )
  R_cumsum[ cbind( recovered$tu_index, recovered$au_index ) ] = recovered$recovered

  R_daily = R_cumsum[] * 0
  R_daily[2:length(tus),] = R_cumsum[ 2:length(tus),] - R_cumsum[ 1:(length(tus)-1),]
  R_daily[ which(R_daily < 0)] = 0  ## there is a typo in Alberta 19-05-2020       Alberta                 5854
  R_daily[ which( !is.finite(R_daily))] = 0

  I_current = I_daily[] * 0
  for (k in 1:length(aus)) {
  for (j in 1:(length(tus)-1)) {
    I_current[j+1,k] = max(0, I_current[j,k] + I_daily[j,k]  - D_daily[j,k] - R_daily[j,k] )
  }}

  # compute Susceptibles
  # fill with initial pop


# note: Stats Can and Covid19 name regions and provinces differently ... using lookup defined in: au_lookup_canada.R
  Npop = pop_hr$Npop
  names(Npop) = pop_hr$AU

  S_cumsum = array( NA, dim=c( length(tus), length(aus)  ) )
  S_cumsum[1,] = Npop[aus]  # start of data
  for (k in 1:length(aus)) {
  for (j in 1:(length(tus)-1) ) {
    S_cumsum[j+1,k] = S_cumsum[j,k] - sum( I_current[j,k] + R_cumsum[j,k] + D_cumsum[j,k], na.rm=TRUE )
  }}

    # default is to return this:
  data_au = list()
  for ( i in 1:length(aus) ) {
    au = aus[i]
    au_statscan = au_lookup_canada( covid19=aus[i] )
    if (is.na(au_statscan)) next()  # no match
    stan_data = NULL
    stan_data  = list(
      Npop = Npop[[ au_statscan ]],
      Nobs = length(tus),
      Npreds = Npreds,
      Sobs = floor( S_cumsum[,i] ),
      Iobs = floor( I_current[,i] ),
      Robs = floor( R_cumsum[,i] ),
      Mobs = D_cumsum[,i],
      daterange = daterange,
      tu = tus,
      au = au,
      statevars = c("susceptible", "infected", "recovered", "dead"),  # names of last dim of "daily"
      time = tus,
      time_start = daterange[1],
      timestamp = max( lubridate::as_date( daterange ) )
    )


    #stan_data = c( stan_data, list(...) )
    if (!exists("modelname", stan_data)) stan_data$modelname="default"
    # add a few more flags for discrete_variable_encounter_rate and
    if (!exists("BETA_max", stan_data)) stan_data$BETA_max = 0.5
    if (!exists("GAMMA_max", stan_data)) stan_data$GAMMA_max = 0.5
    if (!exists("EPSILON_max", stan_data)) stan_data$EPSILON_max = 0.5
    if (!exists("BNP", stan_data)) stan_data$BNP = 5  # % of Infected that are asymptomatic

    if ( stan_data$modelname %in% c("continuous") ) {
      # used by ODE-based methods for rk4 integration
      if (!exists("time_pred", stan_data)) stan_data$time_pred = as.integer( c(tus, max(tus)+c(1:Npreds)) )
      if (!exists("t0", stan_data)) stan_data$t0 = -0.01
    }

    data_au[[au]] = stan_data
  }
  return( data_au )
}
