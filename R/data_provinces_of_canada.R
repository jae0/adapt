#' @title data_provinces_of_canada
#' @description This is a placeholder for a description.
#' @param selection default is \code{"default"}
#' @param fn default is \code{NULL}
#' @param Npreds default is \code{5}
#' @param ... other arguments passed to methods
#' @return  This is a placeholder for what it returns.
#' @author Jae Choi, \email{choi.jae.seok@gmail.com}
#' @export
data_provinces_of_canada = function( selection="default", fn=NULL, Npreds=5, ... ) {

  # this accesses data from : https://github.com/ishaberry/Covid19Canada
  # Citation:  Berry I, Soucy J-PR, Tuite A, Fisman D. Open access epidemiologic data and an interactive dashboard to monitor the COVID-19 outbreak in Canada. CMAJ. 2020 Apr 14;192(15):E420. doi: https://doi.org/10.1503/cmaj.75262

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19Canada.rdata")  # default to current work directory

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
    res = data_provinces_of_canada( selection="download", fn=fn )
  } else {
    load (fn )
  }
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
  tus = 1:(diff(daterange)+1)

  aus = sort( unique( res$cases$province))
  aus = setdiff(aus,  "Repatriated")  # not clear what to do with these .. droppoing for now
  aus = c(setdiff( aus, c("Ontario", "Quebec") ), c("Ontario", "Quebec"))  # do Ontario and Quebec last as they are slow

  infected = infected[ infected$province %in% aus , ]
  death = death[ death$province %in% aus , ]
  recovered = recovered[ recovered$province %in% aus , ]

  infected$dayno = infected$date - daterange[1] + 1
  recovered$dayno = recovered$date - daterange[1] + 1
  death$dayno = death$date - daterange[1] + 1

  Ntimes = length(tus)


  I_daily = array( 0, dim=c( length(tus), length(aus)  ) )  # SIRD (D=death)
  I_daily[ cbind( match(infected$dayno, tus ), match(infected$province, aus) ) ] = infected$infected
  I_daily[ which(!is.finite(I_daily))] = 0
  I_daily[ which(I_daily < 0)] = 0  ## in case of data entry errors

  D_daily = array( 0, dim=c( length(tus), length(aus)  ) )  # SIRD (D=death)
  D_daily[ cbind( match(death$dayno, tus ), match(death$province, aus) ) ] = death$death
  D_daily[ which(!is.finite(D_daily))] = 0
  D_daily[ which(D_daily < 0)] = 0  ## in case of data entry errors

  R_daily = array( 0, dim=c( length(tus), length(aus)  ) )  # SIRD (D=death)
  R_cumsum = array( 0, dim=c( length(tus), length(aus)  ) )  # SIRD (D=death)
  R_cumsum[ cbind( match(recovered$dayno, tus ), match(recovered$province, aus) ) ] = recovered$recovered
  R_cumsum[ which( !is.finite(R_cumsum))] = 0


  Npop = pop$Npop
  names(Npop) = pop$au

  R_daily[2:Ntimes,] = R_cumsum[ 2:Ntimes,] - R_cumsum[ 1:(Ntimes-1),]
  R_daily[ which(R_daily < 0)] = 0  ## there is a typo in Alberta 19-05-2020       Alberta                 5854
  R_daily[ which( !is.finite(R_daily))] = 0

  # cummulative sums for Mortalities
  D_cumsum = D_daily[] * 0
  for (k in 1:length(aus)) D_cumsum[,k] = cumsum( D_daily[,k])

  # compute Susceptibles
  # fill with initial pop
  S_cumsum = array( 0, dim=c( length(tus), length(aus)  ) )  # SIRD (D=death)
  S_cumsum[1,] = Npop[aus]  # start of data
  for (k in 1:length(aus)) {
  for (j in 1:(Ntimes-1) ) {
    S_cumsum[j+1,k] =  S_cumsum[j,k] - I_daily[j,k]
  }}
  S_cumsum[ which(!is.finite(S_cumsum))] = 0
  S_cumsum[ which(S_cumsum < 0)] = 0  ## in case of data entry errors


  I_active = I_daily[] * 0
  for (k in 1:length(aus)) {
  for (j in 1:(length(tus)-1)) {
    I_active[j+1,k] = max(0, I_active[j,k] + I_daily[j,k]  - D_daily[j,k] - R_daily[j,k] )
  }}
  I_active[ which(!is.finite(I_active))] = 0
  I_active[ which(I_active < 0)] = 0  ## in case of data entry errors

    # default is to return this:
  data_province = list()
  for ( i in 1:length(aus) ) {
    prov = aus[i]
    stan_data = NULL
    stan_data  = list(
      Npop = Npop[[ prov ]],
      Nobs = Ntimes,
      Npreds = Npreds,
      Sobs = S_cumsum[,i],
      Iobs = I_active[,i],
      Robs = R_cumsum[,i],
      Mobs = D_cumsum[,i],
      daterange = daterange,
      tus = tus,
      au = prov,
      statevars = c("susceptible", "infected", "recovered", "dead"),  # names of last dim of "daily"
      time = tus,
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
