#' @title simulate_data
#' @description This is a placeholder for a description.
#' @param selection default is \code{"sir"}
#' @param Npop default is \code{1}
#' @param inits default is \code{c(0.99, 0.01, 0)}
#' @param times default is \code{0:99}
#' @param params default is \code{list(beta=0.6, gamma=0.1)}
#' @param plotdata default is \code{TRUE}
#' @param sample_init default is \code{NULL}
#' @param nsample default is \code{NULL}
#' @return  This is a placeholder for what it returns.
#' @importFrom magrittr "%>%"
#' @author Jae Choi, \email{choi.jae.seok@gmail.com}
#' @export
simulate_data = function( selection="sir", Npop=1, inits=c(0.99, 0.01, 0), times=0:99, params=list(beta=0.6, gamma=0.1), plotdata=TRUE, sample_init=NULL, nsample=NULL ) {

subject <- status <- time <- start_shed <- NA
. = NULL #this is to prevent build check from warning about valid dplyr use 
  if (selection=="sir") {
    # CREDIT to:
    # https://jrmihalj.github.io/estimating-transmission-by-fitting-mechanistic-models-in-Stan/

    # inits = initial proportions in each SIR category (S0, I0, R0)
    # params = c( transmission and pathogen-induced death rates)

    # library(deSolve)

    # ODE function
    delta_sir = function(t, y, params) {
      with( as.list(c(params, y) ), {
        dS = - beta * y[1] * y[2]
        dI = beta * y[1] * y[2] - gamma * y[2]
        dR = gamma * y[2]
        res = c(dS,dI,dR)
        list(res)
      })
    }

    # Run the integration:
    res = as.data.frame( deSolve::ode(inits, times, delta_sir, params, method="ode45") )
    colnames(res) = c("time", "Sobs", "Iobs", "Robs")
    res[, c("Iobs", "Sobs", "Robs")] = Npop * res[, c("Iobs", "Sobs", "Robs")]
    res[, "Sobs" ] = as.integer( res[, "Sobs" ])
    res[, "Iobs" ] = as.integer( res[, "Iobs" ])
    res[, "Robs" ] = as.integer( res[, "Robs" ])


    if (!is.null(nsample)) {
      if (is.null( sample_init)) sample_init = 1:max( floor(length(times)/2), nsample )

      tokeep = sort( sample( sample_init, nsample, replace=FALSE ) )
      ss = min( tokeep ) : max( tokeep )
      todrop = setdiff( ss, tokeep )
      res$todrop = FALSE
      res$todrop[todrop] = TRUE
      res$dS = NA
      res$dI = NA
      res$dR = NA
      res$dS[2:nrow(res)] = diff(res$S )
      res$dI[2:nrow(res)] = diff(res$I )
      res$dR[2:nrow(res)] = diff(res$R )

      sim = res[ss,]

      accum = c(0,0,0)
      for ( i in 1:nrow(sim) ) {
        if ( ! sim$todrop[i] ) {
          sim$Sobs[i] = max(0, sim$Sobs[i] + accum[1])
          sim$Iobs[i] = max(0, sim$Iobs[i] + accum[2])
          sim$Robs[i] = max(0, sim$Robs[i] + accum[3])
          accum = c(0,0,0)
        } else {
          sim$Sobs[i] = sim$Sobs[i-1]
          sim$Iobs[i] = sim$Iobs[i-1]
          sim$Robs[i] = sim$Robs[i-1]
          accum = accum + c( sim$dS[i-1], sim$dI[i-1], sim$dR[i-1] )
        }

      }

      sim$Sobs[ which(!is.finite(sim$Sobs)) ] = -1 # reset NAs to Inf as stan does not take NAs
      sim$Iobs[ which(!is.finite(sim$Iobs)) ] = -1 # reset NAs to Inf as stan does not take NAs
      sim$Robs[ which(!is.finite(sim$Robs)) ] = -1 # reset NAs to Inf as stan does not take NAs

      # missing data
      # add "missing counts" to adjacent (subsequent) time period to mimic irregular counting periods

    }

    if (plotdata) {

      graphics::plot(NA,NA, xlim=range(times), ylim=range( 0, max(1, Npop) ), xlab = "Time", ylab="Population")

      graphics::lines(res$Sobs ~ res$time, col="black")
      graphics::lines(res$Iobs ~ res$time, col="red")
      graphics::lines(res$Robs ~ res$time, col="green")

      if (!is.null(nsample)) {
        graphics::points(sim$Sobs ~ sim$time, col="black")
        graphics::points(sim$Iobs ~ sim$time, col="red")
        graphics::points(sim$Robs ~ sim$time, col="green")
      }
      graphics::legend(x = 0.3*max(res$time), y = 0.8*Npop, legend = c("Susceptible", "Infected", "Recovered"),
            col = c("black", "red", "green"), lty = c(1, 1, 1), lwd=c(3,3,3), bty="n")

    }

    if (!exists("sim")) sim= res

    return(sim)
  }


  # ----------


  if (selection=="sir_stochastic") {

    # CREDIT to Arie Voorman, April 2017
    # adapted by Jae Choi, April 2020
    # https://rstudio-pubs-static.s3.amazonaws.com/270496_e28d8aaa285042f2be0c24fc915a68b2.html

    # Npop = 300 # number of subjects
    # I0 = 25 # number of children initially infected
    # prob_infection = 0.001 # transmission probability (per person per 'close contact', per day)
    # time_shedding = 7 # mean I->R duration for an individual

    # require(dplyr)
    # require(ggplot2)
    # require(tidyr)

    prob_infection = params$prob_infection
    time_shedding = params$time_shedding

    #data frame to store results
    data = expand.grid(subject = 1:Npop, time=times) %>% dplyr::tbl_df()
    wdata = data %>%
      dplyr::mutate(status = NA_character_) %>%
      tidyr::spread(subject,status) %>%
      stats::setNames(names(.) %>%
      make.names
    )

    # initial state
    I0 = floor( Npop * inits[2] )

    infected = c(rep(1,I0),rep(0,Npop-I0))
    susceptible = 1-infected
    recovered = rep(0,Npop)

    wdata[1,-1][ which( infected == 1)] = 'I'
    wdata[1,-1][ which( (1-infected) ==1)] = 'S'
    wdata[1,-1][ which( recovered == 1) ] = 'R'

    #run through the simulation:
    for (t in 1:max(times)){
      wdata[t+1,-1] = wdata[t,-1]

      i = which( wdata[t+1,-1] == "I" )
      s = which( wdata[t+1,-1] == "S" )

      pr.infection = rep(0, Npop )
      pr.infection[s] = 1 - (1 - prob_infection) ^ length(i)

      pr.recovery =  1/time_shedding

      s_to_i = intersect( s, which( stats::rbinom(Npop, 1, pr.infection) == 1  ) )
      if (length(s_to_i) > 0)  wdata[t+1,-1][ s_to_i ] = "I"

      i_to_r = intersect( i, which( stats::rbinom(Npop, 1, pr.recovery) == 1 ) )
      if (length(i_to_r) > 0)  wdata[t+1,-1][ i_to_r ] = "R"
    }

    data = tidyr::gather(wdata, subject, status, -time) %>% dplyr::tbl_df()

    data = data %>%
      dplyr::group_by(subject) %>%
      dplyr::mutate(start_shed = as.numeric(min(time[status == 'I']))) %>%
      dplyr::ungroup %>%
    dplyr::mutate(subject = stats::reorder(subject,-start_shed), status=factor(status, levels = c('S','I','R'),ordered = T))

    sim = data %>% dplyr::group_by(time) %>%
      plyr::summarise(Sobs = sum(status =='S') %>% as.integer,
                Iobs = sum(status =='I') %>% as.integer,
                Robs = sum(status =='R') %>% as.integer)

    sim = as.data.frame(sim)

    if (plotdata) {
      #plot the SIR status for each individual over time:
      ggplot2::ggplot( data, ggplot2::aes(x=time,y=as.numeric(subject),fill = status)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual('Status',values = c('Sobs'='green3','Iobs'='indianred','Robs'='blue3')) +
        ggplot2::ylab('Subject') + ggplot2::xlab('Time')
    }

    return(sim)
  }



  if (selection=="sir_stochastic_family") {

    # CREDIT to Arie Voorman, April 2017
    # adapted by Jae Choi, April 2020
    # https://rstudio-pubs-static.s3.amazonaws.com/270496_e28d8aaa285042f2be0c24fc915a68b2.html

    # Npop = 300 # number of subjects
    # I0 = 25 # number of children initially infected

    # require(dplyr)
    # require(ggplot2)
    # require(tidyr)
    # 
    # require(Matrix)

    prob_infection = params$prob_infection
    time_shedding = params$time_shedding
    # lambda_family = 0.1 # transmission probability (per person per 'close contact', per day)
    # lambda_community = 0.0005# transmission probability (per person per 'person in community', per day)

    # gamma_child = 10 # mean recovery duration for child (naive child)
    # gamma_adult = 3 # mean recovery duration for an adult (prior immunity)

    fmat = Matrix::bdiag( replicate( Npop/4, matrix(1,4,4), simplify = FALSE ) ) #matrix indicating family membership
    diag(fmat) = 0

    adult = rep(c(1,1,0,0), Npop/4) # two adults and two children per HH

    data = expand.grid(subject = 1:Npop, time=times) %>% dplyr::tbl_df()

    wdata = data %>%
      dplyr::mutate(status = NA_character_) %>%
      tidyr::spread(subject,status) %>%
      stats::setNames(names(.) %>%
      make.names
    )

  # initial state
    #MM - I0 was i0 - which was never defined.
    inits = list(
      c(1-I0, I0, 0), # S, I, R, proportions for children 
      c(1, 0, 0)      # S, I, R, proportions for adults
    )

    I0 = floor( Npop * inits[[1]][2] )

    # keep track of immune and naive individuals
    infected = c( rep(c(0,0,0,1), I0 ), rep(0, Npop-I0*4 )) %>% as.numeric
    susceptible = 1-infected
    recovered = rep(0,Npop)

    wdata[1,-1][ which( infected == 1)] = 'I'
    wdata[1,-1][ which( (1-infected) ==1)] = 'S'
    wdata[1,-1][ which( recovered == 1) ] = 'R'


    #run through the simulation:
    for(t in 1:max(times)){

      #probability of infection = 1- probability of not getting infected, from either within HH, or within community
      #set to zero for those that aren't susceptible

      wdata[t+1,-1] = wdata[t,-1]

      i = which( wdata[t+1,-1] == "I" )
      s = which( wdata[t+1,-1] == "S" )

      pr.infection = rep(0, Npop )
      pr.infection[s] = 1 - ((1-prob_infection$family)^(fmat%*%status$infected)) * ((1-prob_infection$community)^(sum(status$infected) - fmat%*%status$infected) )

      pr.recovery =  1/time_shedding

      s_to_i = intersect( s, which( stats::rbinom(Npop, 1, pr.infection) == 1  ) )
      if (length(s_to_i) > 0)  wdata[t+1,-1][ s_to_i ] = "I"

      i_to_r = intersect( i, which( stats::rbinom(Npop, 1, pr.recovery) == 1 ) )
      if (length(i_to_r) > 0)  wdata[t+1,-1][ i_to_r ] = "R"
    }


    sim = list(
      Sobs = as.numeric(wdata[,-1] == 'S') %>% matrix(t+1,Npop),
      Iobs = as.numeric(wdata[,-1] == 'I') %>% matrix(t+1,Npop),
      Robs = as.numeric(wdata[,-1] == 'R') %>% matrix(t+1,Npop),
      fmat = as.matrix(fmat),
      adult = adult
    )

    if (plotdata) {
      out = tidyr::gather(wdata,subject,status, -time) %>% dplyr::tbl_df()

      out = out %>% dplyr::group_by(subject) %>% dplyr::mutate(start_shed = as.numeric(min(time[status == 'I']))) %>%
        dplyr::ungroup %>%
        dplyr::mutate(subject = stats::reorder(subject,-start_shed), status = factor(status, levels = c('S','I','R'),ordered = T))

      out = as.data.frame(out)

      #plot the SIR status for each individual over time:
      ggplot2::ggplot(out, ggplot2::aes(x=time,y=as.numeric(subject),fill = status)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual('Status',values = c('Sobs'='green3','Iobs'='indianred','Robs'='blue3')) +
        ggplot2::ylab('Subject') + ggplot2::xlab('Time')
    }

    return(sim)
  }

}


