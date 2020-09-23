#' @title plot_model_fit
#' @description This is a placeholder for a description.
#' @param selection default is sim
#' @param stan_results default is \code{NULL}
#' @param sim default is \code{NULL}
#' @param nsimlines default is \code{2000}
#' @param outdir default is \code{""}
#' @param nx default is \code{NULL}
#' @param to.screen default is \code{FALSE}
#' @return  This is a placeholder for what it returns.
#' @author Jae Choi, \email{choi.jae.seok@gmail.com}
#' @export
plot_model_fit = function( selection="default", stan_results=NULL,
  sim=NULL, nsimlines=2000, outdir="", nx=NULL, to.screen =FALSE ) {
  
  # create variable that will be used
sample_prop <- ts <- mod_time <- mod_median <- mod_high <- mod_low <-NA
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )

  stan_data = stan_results$stan_inputs
  posteriors = rstan::extract( stan_results$stan_samples )  # posteriors = mcmc posteriors from STAN

  # nx = stan_data$Nobs + trunc( stan_data$Npreds *0.2 )
  if (is.null (nx)) nx = stan_data$Nobs + stan_data$Npreds - 1
  ndat = stan_data$Nobs
  xrange = c(0, nx)

  if (selection=="default") {
    # data to plot
    df_sample = data.frame( cbind(
      sample_prop = stan_data$Iobs / stan_data$Npop,
      ts=stan_data$time
    ))

    if (!is.null(posteriors)) {
      # Model predictions across the sampling time period.
      df_fit = data.frame( cbind(
        mod_median = apply(posteriors$I/stan_data$Npop, 2, stats::median, na.rm=TRUE),
        mod_low = apply(posteriors$I/stan_data$Npop, 2, stats::quantile, probs=c(0.025), na.rm=TRUE),
        mod_high = apply(posteriors$I/stan_data$Npop, 2, stats::quantile, probs=c(0.975), na.rm=TRUE),
        mod_time = c( stan_data$time, max(stan_data$time) + c(1:stan_data$Npreds) )
      ))
      # Median and 95% Credible Interval
      ggplot2::ggplot(df_sample, ggplot2::aes(x=ts, y=sample_prop)) +
        ggplot2::geom_point(col="black", shape = 19, size = 1.5) +
        ggplot2::geom_line(data = df_fit, ggplot2::aes(x=mod_time, y=mod_median), color = "red") +
        ggplot2::geom_line(data = df_fit, ggplot2::aes(x=mod_time, y=mod_high), color = "red", linetype=3) +
        ggplot2::geom_line(data = df_fit, ggplot2::aes(x=mod_time, y=mod_low), color = "red", linetype=3) +
        ggplot2::labs(x = "Time (days)", y = "Proportion Infected") +
        ggplot2::scale_x_continuous(limits=c(0, dim(posteriors$I)[2]), breaks=c(0,25,50, 75, 100)) +
        ggplot2::scale_y_continuous(limits=c(0, max(df_fit$mod_high) ), breaks=c(0,.2, .4, .6, .8, 1)) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.line.x = ggplot2::element_line(color="black"),
              axis.line.y = ggplot2::element_line(color="black")
      )
    } else {
      # Median and 95% Credible Interval
      ggplot2::ggplot(df_sample, ggplot2::aes(x=ts, y=sample_prop)) +
        ggplot2::geom_point(col="black", shape = 19, size = 1.5) +
        ggplot2::labs(x = "Time (days)", y = "Proportion Infected") +
        ggplot2::scale_x_continuous(limits=c(0, 50), breaks=c(0,25,50)) +
        ggplot2::scale_y_continuous(limits=c(0,1), breaks=c(0,.5,1)) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.line.x = ggplot2::element_line(color="black"),
              axis.line.y = ggplot2::element_line(color="black")
      )
    }
  }


  if (selection=="infected_affected") {

    ap = (posteriors$R + posteriors$I+ posteriors$M)[,1:ndat]
    ip = posteriors$I[,1:ndat]

    col_palette = scales::alpha( rep(c( "slateblue", "darkgreen", "cyan", "magenta", "gray", "darkorange" ), length.out=nx), 0.02 )
    cols = matrix( col_palette, ncol=nx, nrow=nrow(ip), byrow=TRUE )
    xrange = c(0, log10( max( ap, na.rm=TRUE) ) )
    yrange = c(0, log10( max( ip[,1:ndat], na.rm=TRUE) ) )

    xticks = seq( xrange[1], xrange[2], length.out=7)
    yticks = seq( yrange[1], yrange[2], length.out=7)

    xvals = round( 10^xticks )
    yvals = round( 10^yticks )

    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "infected_affected.png"))
    } else {
      grDevices::dev.new()
    }
    graphics::plot( 0 , 0, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Total affected", type="n", axes=FALSE )
      graphics::points( log10(ip) ~ log10(ap), col=cols, cex=0.8, pch=19 )
      graphics::lines( log10(apply( ip, 2, median)) ~ log10( apply( ap, 2,stats::median) ), col="orange", lty="solid", lwd=2 )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
      graphics::axis( 1, at=xticks, labels=xvals )
      graphics::axis( 2, at=yticks, labels=yvals )

    if (!to.screen) grDevices::dev.off()
  }


  if (selection=="susceptible") {
    so = stan_data$Sobs
    so[so < 0] = NA
    sp = apply(posteriors$S, 2, stats::median)[1:nx]
    xrange = range(1:nx)
    yrange = range( c(so, sp[1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2] )
    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "fit_with_projections_susceptible.png"))
    } else {
      grDevices::dev.new()
    }
      cols = scales::alpha( "magenta", 0.005 )
      graphics::plot( so ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Susceptible", xlab="Days", type="n" )
      graphics::lines( sp ~ seq(1,nx), lwd = 3, col="slateblue" )

      nr = nrow(posteriors$S)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        graphics::points( posteriors$S [i,1:nx] ~ c(1:nx) , col=cols, cex=0.8, pch=19 )
      }
      graphics::lines( apply(posteriors$S, 2, stats::quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::lines( apply(posteriors$S, 2, stats::quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::points( so ~ stan_data$time, col="darkgray", cex=1.2 )
      graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) grDevices::dev.off()
  }


  if (selection=="infected_effective") {
    io = stan_data$Iobs
    io[io < 0] = NA

    if (! exists("Q", posteriors)) return("Nothing to do")  # nothing to do

    ip = posteriors$I * 0
    ip[ , c(1:stan_data$Nobs-1)] = posteriors$I[ , c(1:stan_data$Nobs-1)] * posteriors$Q[, c(1:stan_data$Nobs-1)]
    ip[ , c(stan_data$Nobs: nx)] = posteriors$I[ , c(stan_data$Nobs: nx)] * posteriors$Q[, stan_data$Nobs-1]

    ipm = apply(ip, 2, stats::median)[1:nx]
    ipl = apply(ip, 2, stats::quantile, probs=0.025)[1:nx]
    ipu = apply(ip, 2, stats::quantile, probs=0.975)[1:nx]

    yrange = range( c(io, ipm[1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "fit_with_projections_infected_effective.png"))
    } else {
      grDevices::dev.new()
    }

      graphics::plot( io ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Effective nmber of infected", xlab="Days", type="n" )
      graphics::lines( ipm ~ seq(1,nx), lwd =3, col="slateblue" )
      graphics::lines( ipu ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::lines( ipl ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::points( io ~ stan_data$time, col="darkgray", cex=1.2 )
      graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
   if (!to.screen) grDevices::dev.off()
  }


  if (selection=="infected") {
    io = stan_data$Iobs
    io[io < 0] = NA

    ip = apply(posteriors$I, 2, stats::median)[1:nx]

    yrange = range( c(io, posteriors$I[,1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "fit_with_projections_infected.png"))
    } else {
      grDevices::dev.new()
    }
      cols = scales::alpha( "magenta", 0.005 )
      graphics::plot( io ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Days", type="n" )

      nr = nrow(posteriors$I)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        graphics::points( posteriors$I[i,1:nx] ~ c(1:nx) , col=cols, cex=1, pch=19 )
      }
      graphics::lines( ip ~ seq(1,nx), lwd =3, col="lightblue" )
      graphics::lines( apply(posteriors$I, 2, stats::quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::lines( apply(posteriors$I, 2, stats::quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::points( io ~ stan_data$time, col="slateblue", cex=1.2 )
      graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
   if (!to.screen) grDevices::dev.off()
  }

  if (selection=="recovered") {
    ro = stan_data$Robs
    ro[ro < 0] = NA
    rp = apply(posteriors$R, 2, stats::median)[1:nx]

    yrange = range( c(ro, rp[1:ndat]),  na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])

    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "fit_with_projections_recovered.png"))
    } else {
      grDevices::dev.new()
    }
      cols = scales::alpha( "magenta", 0.002 )

      graphics::plot( ro ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Recovered", xlab="Days", type="n" )

      nr = nrow(posteriors$R)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        graphics::points( posteriors$R[i,1:nx] ~ c(1:nx) , col=cols, cex=0.8, pch=19 )
      }
      graphics::lines( rp ~ seq(1,nx), lwd =3, col="slateblue" )
      graphics::lines( apply(posteriors$R, 2, stats::quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::lines( apply(posteriors$R, 2, stats::quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::points( ro ~ stan_data$time, col="darkgray", cex=1.2 )
      graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) grDevices::dev.off()
  }


  if (selection=="deaths") {
    mo = stan_data$Mobs
    mo[mo < 0] = NA
    mp = apply(posteriors$M, 2, stats::median)[1:nx]
    yrange = range( c(mo, mp[1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "fit_with_projections_mortalities.png"))
    } else {
      grDevices::dev.new()
    }
      cols = scales::alpha( "magenta", 0.002 )

      graphics::plot( mo ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Mortalities", xlab="Days", type="n" )

      nr = nrow(posteriors$M)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        graphics::points( posteriors$M[i,1:nx] ~ c(1:nx) , col=cols, cex=0.8, pch=19 )
      }
      graphics::lines( mp ~ seq(1,nx), lwd =3, col="slateblue" )
      graphics::lines( apply(posteriors$M, 2, stats::quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::lines( apply(posteriors$M, 2, stats::quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::points( mo ~ stan_data$time, col="darkgray", cex=1.2 )
      graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) grDevices::dev.off()
  }


  if (selection=="effective_number") {
    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "effective_number.png"))
    } else {
      grDevices::dev.new()
    }

    nxr = stan_data$Nobs-1
    dr = c(1:nxr)
    ipm = apply(posteriors$Q[,dr], 2, stats::median)
    ipl = apply(posteriors$Q[,dr], 2, stats::quantile, probs=0.025)
    ipu = apply(posteriors$Q[,dr], 2, stats::quantile, probs=0.975)

    yrange = range( posteriors$Q )
    # yrange = c(0.5, 1.5)  #range( c(ipm) ) #, ipl, ipu) )

    graphics::plot( ipm ~ dr, type="l", lwd =3, col="slateblue", ylim=yrange, ylab="Effective number (proportion of infected)", xlab="Days" )
    graphics::lines( ipl ~ dr, col="darkorange", lty="dashed" )
    graphics::lines( ipu ~ dr, col="darkorange", lty="dashed" )
    graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed", lwd=2 )
    graphics::abline( h=1, col="grey", lty="dashed", lwd=2 )
    graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) grDevices::dev.off()
  }



  if (selection=="reproductive_number") {
    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "reproductive_number.png"))
    } else {
      grDevices::dev.new()
    }
      cols = scales::alpha( "magenta", 0.002 )

      nr = nrow(posteriors$K)
      ny = min( nr, 1200)

      nxr = stan_data$Nobs + stan_data$Npreds - 4  # in case data are not upto-date
      rp = apply(posteriors$K[,1:nxr], 2, stats::median)
      yrange = range(c(0, rp, posteriors$K[,1:nxr], 3))

      graphics::plot( rp ~ seq(1,nxr), type="l", lwd =3, col="slateblue", ylim=yrange, ylab="Reproductive number", xlab="Days" )

      for (i in sample.int(nr, ny)) {
        graphics::points( posteriors$K[i,1:nxr] ~ c(1:nxr) , col=cols, cex=0.8, pch=19 )
      }

      graphics::lines( apply(posteriors$K, 2, stats::quantile, probs=0.025)[1:nxr] ~ seq(1,nxr), col="darkorange", lty="dashed" )
      graphics::lines( apply(posteriors$K, 2, stats::quantile, probs=0.975)[1:nxr] ~ seq(1,nxr), col="darkorange", lty="dashed" )
      graphics::abline( h=1, col="red", lwd=2 )
      graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) grDevices::dev.off()
  }

  if (selection=="reproductive_number_histograms") {
    brks = 30
    days0 = graphics::hist( posteriors$K[,stan_data$Nobs-1] , breaks=brks, plot=FALSE)  # "today" is actually today-1 as K has no data to condition for estimation
    days1 = graphics::hist( posteriors$K[,stan_data$Nobs-2] , breaks=brks, plot=FALSE )  # today's K
    days7 = graphics::hist( posteriors$K[,stan_data$Nobs-8] , breaks=brks, plot=FALSE )  # today's K
    yrange = range( c(days0$density, days1$density, days7$density))
    yrange[2] = yrange[2] * 1.15
    xrange = range( c(0, days0$mides, days1$mids, days7$mids, 1.3))
    xrange[2] = xrange[2] * 1.15

    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "reproductive_number_today.png"))
    } else {
      grDevices::dev.new()
    }
      graphics::plot(  days0$density ~ days0$mids, col="green", lwd=3, xlab="Reproductive number", ylab="Probability density", main="", xlim=xrange, ylim=yrange, type ="l")
      graphics::lines( days1$density ~ days1$mids, col="slateblue", lwd=3, lty="dotted")
      graphics::lines( days7$density ~ days7$mids, col="darkorange", lwd=3, lty="dashed")
      graphics::abline( v=1, col="red", lwd=3 )
      graphics::abline( h=0, col="gray", lwd=1 )
      graphics::legend( "topright", "", paste( "Current date: ", stan_data$timestamp, " "), bty="n")
      graphics::legend( "topleft", legend=c("Current", "Yesterday", "7 days ago"), lty=c("solid", "dotted", "dashed"), col=c("green", "slateblue", "darkorange"), lwd=c(3,3,3), bty="n")
      graphics::title( main= paste( stan_data$au, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) grDevices::dev.off()

  }

  if (selection=="forecasts") {
    if (!to.screen) {
      grDevices::png(filename = file.path(outdir, "fit_with_projections_and_stochastic_simulations.png"))
    } else {
      grDevices::dev.new()
    }
      simxval = stan_data$Nobs + c(1:dim(sim)[3])
      xrange = c(0, max(nx, simxval) )
      io = stan_data$Iobs
      io[io < 0] = NA
      ip = apply(posteriors$I, 2, stats::median)[1:nx]
      yrange = range( c(io, ip), na.rm=TRUE)
      yrange[2] = yrange[2] * 2
      graphics::plot( io ~ stan_data$time, xlim=xrange, ylim=yrange,  type="n", ylab="Infected", xlab="Days")
      nsimlines = min( nsimlines, dim(sim)[1])
      for ( i in 1:min( dim(sim)[1], nsimlines)) graphics::lines( sim[i,2,] ~ simxval, col=scales::alpha("slategray", 0.1), lty="solid" )
      graphics::lines( ip ~ seq(1,nx), lwd = 3, col="slateblue" )
      graphics::lines( apply(posteriors$I, 2, stats::quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::lines( apply(posteriors$I, 2, stats::quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      graphics::points( io ~ stan_data$time, xlim=xrange, ylim=yrange, col="darkgray", cex=1.2 )
      graphics::abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      graphics::title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) grDevices::dev.off()

  }


}

