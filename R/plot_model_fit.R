
plot_model_fit = function( selection="default", stan_results=NULL,
  sim=NULL, nsimlines=2000, outdir="", nx=NULL, to.screen =FALSE ) {

  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )

  stan_data = stan_results$stan_inputs
  posteriors = rstan::extract( stan_results$stan_samples )  # posteriors = mcmc posteriors from STAN

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
        mod_median = apply(posteriors$I/stan_data$Npop, 2, median, na.rm=TRUE),
        mod_low = apply(posteriors$I/stan_data$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        mod_high = apply(posteriors$I/stan_data$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE),
        mod_time = c( stan_data$time, max(stan_data$time) + c(1:stan_data$Npreds) )
      ))
      # Median and 95% Credible Interval
      ggplot(df_sample, aes(x=ts, y=sample_prop)) +
        geom_point(col="black", shape = 19, size = 1.5) +
        geom_line(data = df_fit, aes(x=mod_time, y=mod_median), color = "red") +
        geom_line(data = df_fit, aes(x=mod_time, y=mod_high), color = "red", linetype=3) +
        geom_line(data = df_fit, aes(x=mod_time, y=mod_low), color = "red", linetype=3) +
        labs(x = "Time (days)", y = "Proportion Infected") +
        scale_x_continuous(limits=c(0, dim(posteriors$I)[2]), breaks=c(0,25,50, 75, 100)) +
        scale_y_continuous(limits=c(0, max(df_fit$mod_high) ), breaks=c(0,.2, .4, .6, .8, 1)) +
        theme_classic() +
        theme(axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black")
      )
    } else {
      # Median and 95% Credible Interval
      ggplot(df_sample, aes(x=ts, y=sample_prop)) +
        geom_point(col="black", shape = 19, size = 1.5) +
        labs(x = "Time (days)", y = "Proportion Infected") +
        scale_x_continuous(limits=c(0, 50), breaks=c(0,25,50)) +
        scale_y_continuous(limits=c(0,1), breaks=c(0,.5,1)) +
        theme_classic() +
        theme(axis.line.x = element_line(color="black"),
              axis.line.y = element_line(color="black")
      )
    }
  }


  if (selection=="infected_affected") {

    ap = (posteriors$R + posteriors$I+ posteriors$M)[,1:ndat]
    ip = posteriors$I[,1:ndat]

    col_palette = alpha( rep(c( "slateblue", "darkgreen", "cyan", "magenta", "gray", "darkorange" ), length.out=nx), 0.02 )
    cols = matrix( col_palette, ncol=nx, nrow=nrow(ip), byrow=TRUE )
    xrange = c(0, log10( max( ap, na.rm=TRUE) ) )
    yrange = c(0, log10( max( ip[,1:ndat], na.rm=TRUE) ) )

    xticks = seq( xrange[1], xrange[2], length.out=7)
    yticks = seq( yrange[1], yrange[2], length.out=7)

    xvals = round( 10^xticks )
    yvals = round( 10^yticks )

    if (!to.screen) {
      png(filename = file.path(outdir, "infected_affected.png"))
    } else {
      dev.new()
    }
      plot( 0 , 0, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Total affected", type="n", axes=FALSE )
      points( log10(ip) ~ log10(ap), col=cols, cex=0.8, pch=19 )
      lines( log10(apply( ip, 2, median)) ~ log10( apply( ap, 2, median) ), col="orange", lty="solid", lwd=2 )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
      axis( 1, at=xticks, labels=xvals )
      axis( 2, at=yticks, labels=yvals )

    if (!to.screen) dev.off()
  }


  if (selection=="susceptible") {
    so = stan_data$Sobs
    so[so < 0] = NA
    sp = apply(posteriors$S, 2, median)[1:nx]
    xrange = range(1:nx)
    yrange = range( c(so, sp[1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2] )
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_susceptible.png"))
    } else {
      dev.new()
    }
      cols = alpha( "magenta", 0.005 )
      plot( so ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Susceptible", xlab="Days", type="n" )
      lines( sp ~ seq(1,nx), lwd = 3, col="slateblue" )

      nr = nrow(posteriors$S)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        points( posteriors$S [i,1:nx] ~ c(1:nx) , col=cols, cex=0.8, pch=19 )
      }
      lines( apply(posteriors$S, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(posteriors$S, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( so ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }


  if (selection=="infected_effective") {
    io = stan_data$Iobs
    io[io < 0] = NA

    if (! exists("Q", posteriors)) return("Nothing to do")  # nothing to do

    ip = posteriors$I * 0
    ip[ , c(1:stan_data$Nobs-1)] = posteriors$I[ , c(1:stan_data$Nobs-1)] * posteriors$Q[, c(1:stan_data$Nobs-1)]
    ip[ , c(stan_data$Nobs: nx)] = posteriors$I[ , c(stan_data$Nobs: nx)] * posteriors$Q[, stan_data$Nobs-1]

    ipm = apply(ip, 2, median)[1:nx]
    ipl = apply(ip, 2, quantile, probs=0.025)[1:nx]
    ipu = apply(ip, 2, quantile, probs=0.975)[1:nx]

    yrange = range( c(io, ipm[1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_infected_effective.png"))
    } else {
      dev.new()
    }

      plot( io ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Effective nmber of infected", xlab="Days", type="n" )
      lines( ipm ~ seq(1,nx), lwd =3, col="slateblue" )
      lines( ipu ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( ipl ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( io ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
   if (!to.screen) dev.off()
  }


  if (selection=="infected") {
    io = stan_data$Iobs
    io[io < 0] = NA

    ip = apply(posteriors$I, 2, median)[1:nx]

    yrange = range( c(io, posteriors$I[,1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_infected.png"))
    } else {
      dev.new()
    }
      cols = alpha( "magenta", 0.005 )
      plot( io ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Days", type="n" )

      nr = nrow(posteriors$I)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        points( posteriors$I[i,1:nx] ~ c(1:nx) , col=cols, cex=1, pch=19 )
      }
      lines( ip ~ seq(1,nx), lwd =3, col="lightblue" )
      lines( apply(posteriors$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(posteriors$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( io ~ stan_data$time, col="slateblue", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
   if (!to.screen) dev.off()
  }

  if (selection=="recovered") {
    ro = stan_data$Robs
    ro[ro < 0] = NA
    rp = apply(posteriors$R, 2, median)[1:nx]

    yrange = range( c(ro, rp[1:ndat]),  na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])

    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_recovered.png"))
    } else {
      dev.new()
    }
      cols = alpha( "magenta", 0.002 )

      plot( ro ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Recovered", xlab="Days", type="n" )

      nr = nrow(posteriors$R)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        points( posteriors$R[i,1:nx] ~ c(1:nx) , col=cols, cex=0.8, pch=19 )
      }
      lines( rp ~ seq(1,nx), lwd =3, col="slateblue" )
      lines( apply(posteriors$R, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(posteriors$R, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( ro ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }


  if (selection=="deaths") {
    mo = stan_data$Mobs
    mo[mo < 0] = NA
    mp = apply(posteriors$M, 2, median)[1:nx]
    yrange = range( c(mo, mp[1:ndat]), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_mortalities.png"))
    } else {
      dev.new()
    }
      cols = alpha( "magenta", 0.002 )

      plot( mo ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Mortalities", xlab="Days", type="n" )

      nr = nrow(posteriors$M)
      ny = min( nr, 1000)

      for (i in sample.int(nr, ny)) {
        points( posteriors$M[i,1:nx] ~ c(1:nx) , col=cols, cex=0.8, pch=19 )
      }
      lines( mp ~ seq(1,nx), lwd =3, col="slateblue" )
      lines( apply(posteriors$M, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(posteriors$M, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( mo ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }


  if (selection=="effective_number") {
    if (!to.screen) {
      png(filename = file.path(outdir, "effective_number.png"))
    } else {
      dev.new()
    }

    nxr = stan_data$Nobs-1
    dr = c(1:nxr)
    ipm = apply(posteriors$Q[,dr], 2, median)
    ipl = apply(posteriors$Q[,dr], 2, quantile, probs=0.025)
    ipu = apply(posteriors$Q[,dr], 2, quantile, probs=0.975)

    yrange = range( posteriors$Q )
    # yrange = c(0.5, 1.5)  #range( c(ipm) ) #, ipl, ipu) )

    plot( ipm ~ dr, type="l", lwd =3, col="slateblue", ylim=yrange, ylab="Effective number (proportion of infected)", xlab="Days" )
    lines( ipl ~ dr, col="darkorange", lty="dashed" )
    lines( ipu ~ dr, col="darkorange", lty="dashed" )
    abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed", lwd=2 )
    abline( h=1, col="grey", lty="dashed", lwd=2 )
    title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }



  if (selection=="reproductive_number") {
    if (!to.screen) {
      png(filename = file.path(outdir, "reproductive_number.png"))
    } else {
      dev.new()
    }
      cols = alpha( "magenta", 0.002 )

      nr = nrow(posteriors$K)
      ny = min( nr, 1200)

      nxr = stan_data$Nobs + stan_data$Npreds - 4  # in case data are not upto-date
      rp = apply(posteriors$K[,1:nxr], 2, median)
      yrange = range(c(0, rp, posteriors$K[,1:nxr], 3))

      plot( rp ~ seq(1,nxr), type="l", lwd =3, col="slateblue", ylim=yrange, ylab="Reproductive number", xlab="Days" )

      for (i in sample.int(nr, ny)) {
        points( posteriors$K[i,1:nxr] ~ c(1:nxr) , col=cols, cex=0.8, pch=19 )
      }

      lines( apply(posteriors$K, 2, quantile, probs=0.025)[1:nxr] ~ seq(1,nxr), col="darkorange", lty="dashed" )
      lines( apply(posteriors$K, 2, quantile, probs=0.975)[1:nxr] ~ seq(1,nxr), col="darkorange", lty="dashed" )
      abline( h=1, col="red", lwd=2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }

  if (selection=="reproductive_number_histograms") {
    brks = 30
    days0 = hist( posteriors$K[,stan_data$Nobs-1] , breaks=brks, plot=FALSE)  # "today" is actually today-1 as K has no data to condition for estimation
    days1 = hist( posteriors$K[,stan_data$Nobs-2] , breaks=brks, plot=FALSE )  # today's K
    days7 = hist( posteriors$K[,stan_data$Nobs-8] , breaks=brks, plot=FALSE )  # today's K
    yrange = range( c(days0$density, days1$density, days7$density))
    yrange[2] = yrange[2] * 1.15
    xrange = range( c(0, days0$mides, days1$mids, days7$mids, 1.3))
    xrange[2] = xrange[2] * 1.15

    if (!to.screen) {
      png(filename = file.path(outdir, "reproductive_number_today.png"))
    } else {
      dev.new()
    }
      plot(  days0$density ~ days0$mids, col="green", lwd=3, xlab="Reproductive number", ylab="Probability density", main="", xlim=xrange, ylim=yrange, type ="l")
      lines( days1$density ~ days1$mids, col="slateblue", lwd=3, lty="dotted")
      lines( days7$density ~ days7$mids, col="darkorange", lwd=3, lty="dashed")
      abline( v=1, col="red", lwd=3 )
      abline( h=0, col="gray", lwd=1 )
      legend( "topright", "", paste( "Current date: ", stan_data$timestamp, " "), bty="n")
      legend( "topleft", legend=c("Current", "Yesterday", "7 days ago"), lty=c("solid", "dotted", "dashed"), col=c("green", "slateblue", "darkorange"), lwd=c(3,3,3), bty="n")
      title( main= paste( stan_data$au, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()

  }

  if (selection=="forecasts") {
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_and_stochastic_simulations.png"))
    } else {
      dev.new()
    }
      simxval = stan_data$Nobs + c(1:dim(sim)[3])
      xrange = c(0, max(nx, simxval) )
      io = stan_data$Iobs
      io[io < 0] = NA
      ip = apply(posteriors$I, 2, median)[1:nx]
      yrange = range( c(io, ip), na.rm=TRUE)
      yrange[2] = yrange[2] * 2
      plot( io ~ stan_data$time, xlim=xrange, ylim=yrange,  type="n", ylab="Infected", xlab="Days")
      nsimlines = min( nsimlines, dim(sim)[1])
      for ( i in 1:min( dim(sim)[1], nsimlines)) lines( sim[i,2,] ~ simxval, col=alpha("slategray", 0.1), lty="solid" )
      lines( ip ~ seq(1,nx), lwd = 3, col="slateblue" )
      lines( apply(posteriors$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(posteriors$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( io ~ stan_data$time, xlim=xrange, ylim=yrange, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      title( main= paste( stan_data$au, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()

  }


}





