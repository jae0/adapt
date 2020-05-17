
plot_model_fit = function( selection="default", stan_data, M,
  sim=NULL,
  outdir="",
  nx=stan_data$Nobs + stan_data$Npreds - 1,
  is.nova.scotia=FALSE,
  to.screen =FALSE
  ) {

  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )

  # M = mcmc posteriors from STAN
  # nx = stan_data$Nobs + trunc( stan_data$Npreds *0.2 )
  xrange = c(0, nx)

  if (selection=="default") {
    # data to plot
    df_sample = data.frame( cbind(
      sample_prop = stan_data$Iobs / stan_data$Npop,
      ts=stan_data$time
    ))

    if (!is.null(M)) {
      # Model predictions across the sampling time period.
      df_fit = data.frame( cbind(
        mod_median = apply(M$I/stan_data$Npop, 2, median, na.rm=TRUE),
        mod_low = apply(M$I/stan_data$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        mod_high = apply(M$I/stan_data$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE),
        mod_time = c( stan_data$time, max(stan_data$time) + c(1:stan_data$Npreds) )
      ))
      # Median and 95% Credible Interval
      ggplot(df_sample, aes(x=ts, y=sample_prop)) +
        geom_point(col="black", shape = 19, size = 1.5) +
        geom_line(data = df_fit, aes(x=mod_time, y=mod_median), color = "red") +
        geom_line(data = df_fit, aes(x=mod_time, y=mod_high), color = "red", linetype=3) +
        geom_line(data = df_fit, aes(x=mod_time, y=mod_low), color = "red", linetype=3) +
        labs(x = "Time (days)", y = "Proportion Infected") +
        scale_x_continuous(limits=c(0, dim(M$I)[2]), breaks=c(0,25,50, 75, 100)) +
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


  if (selection=="susceptible") {
    so = stan_data$Sobs
    so[so < 0] = NA
    sp = apply(M$S, 2, median)[1:nx]
    yrange = range( c(so, sp), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_susceptible.png"))
    } else {
      dev.new()
    }
      plot( so ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Susceptible", xlab="Days", type="n" )
      lines( sp ~ seq(1,nx), lwd =3, col="slateblue" )
      lines( apply(M$S, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(M$S, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( so ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      if (is.nova.scotia) {
        abline( v=stan_data$time[time_distancing], col="orange", lty="dotted" )
        abline( v=stan_data$time[time_relaxation], col="green", lty="dotted" )
        legend( "bottomleft", "", "      [-- Social distancing -->\n\n", bty="n" )
        legend( "bottom", "", "                                 [-- Parks open -->\n\n", bty="n" )
      }
      title( main= paste( stan_data$plotlabel, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }

  if (selection=="infected") {
    io = stan_data$Iobs
    io[io < 0] = NA
    ip = apply(M$I, 2, median)[1:nx]
    yrange = range( c(io, ip), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_infected.png"))
    } else {
      dev.new()
    }
      plot( io ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Infected", xlab="Days", type="n" )
      lines( ip ~ seq(1,nx), lwd =3, col="slateblue" )
      lines( apply(M$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(M$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( io ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      if (is.nova.scotia) {
        abline( v=stan_data$time[time_distancing], col="orange", lty="dotted" )
        abline( v=stan_data$time[time_relaxation], col="green", lty="dotted" )
        legend( "bottomleft", "", "      [-- Social distancing -->", bty="n" )
        legend( "bottom", "", "                                 [-- Parks open -->", bty="n" )
      }
      title( main= paste( stan_data$plotlabel, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
   if (!to.screen) dev.off()
  }

  if (selection=="recovered") {
    ro = stan_data$Robs
    ro[ro < 0] = NA
    rp = apply(M$R, 2, median)[1:nx]
    yrange = range( c(ro, rp),  na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_recovered.png"))
    } else {
      dev.new()
    }
      plot( ro ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Recovered", xlab="Days", type="n" )
      lines( rp ~ seq(1,nx), lwd =3, col="slateblue" )
      lines( apply(M$R, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(M$R, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( ro ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      if (is.nova.scotia) {
        abline( v=stan_data$time[time_distancing], col="orange", lty="dotted" )
        abline( v=stan_data$time[time_relaxation], col="green", lty="dotted" )
        legend( "left", "", "\n\n      [-- Social distancing -->", bty="n" )
        legend( "center", "", "\n\n                                [-- Parks open -->", bty="n" )
      }
      title( main= paste( stan_data$plotlabel, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }


  if (selection=="deaths") {
    mo = stan_data$Mobs
    mo[mo < 0] = NA
    mp = apply(M$M, 2, median)[1:nx]
    yrange = range( c(mo, mp), na.rm=TRUE)
    yrange = c(yrange[1], yrange[2])
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_mortalities.png"))
    } else {
      dev.new()
    }
      plot( mo ~ stan_data$time, xlim=xrange, ylim=yrange, ylab="Mortalities", xlab="Days", type="n" )
      lines( mp ~ seq(1,nx), lwd =3, col="slateblue" )
      lines( apply(M$M, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(M$M, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( mo ~ stan_data$time, col="darkgray", cex=1.2 )
      abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
      if (is.nova.scotia) {
        abline( v=stan_data$time[time_distancing], col="orange", lty="dotted" )
        abline( v=stan_data$time[time_relaxation], col="green", lty="dotted" )
        legend( "topleft", "", "\n\n      [-- Social distancing -->", bty="n" )
        legend( "top", "", "\n\n                                 [-- Parks open -->", bty="n" )
      }
      title( main= paste( stan_data$plotlabel, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }


  if (selection=="reproductive_number") {
    if (!to.screen) {
      png(filename = file.path(outdir, "reproductive_number.png"))
    } else {
      dev.new()
    }
      rp = apply(M$K[,1:nx], 2, median)
      yrange = range(rp)
      plot( rp ~ seq(1,nx), type="l", lwd =3, col="darkgray", ylim=yrange, ylab="Reproductive number", xlab="Days" )
      lines( apply(M$K, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed" )
      lines( apply(M$K, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed" )
      abline( h=1, col="red", lwd=3 )
      abline( v=stan_data$time[stan_data$Nobs-1], col="grey", lty="dashed" )
      if (is.nova.scotia) {
        abline( v=stan_data$time[time_distancing], col="orange", lty="dotted" )
        abline( v=stan_data$time[time_relaxation], col="green", lty="dotted" )
        legend( "topleft", "", "\n\n      [-- Social distancing -->", bty="n" )
        legend( "top", "", "\n\n                                 [-- Parks open -->", bty="n" )
      }
      title( main= paste( stan_data$plotlabel, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()
  }

  if (selection=="reproductive_number_histograms") {
    brks = 30
    days0 = hist( M$K[,stan_data$Nobs-1] , breaks=brks, plot=FALSE)  # today's K
    days1 = hist( M$K[,stan_data$Nobs-2] , breaks=brks, plot=FALSE )  # today's K
    days7 = hist( M$K[,stan_data$Nobs-8] , breaks=brks, plot=FALSE )  # today's K
    yrange = range( c(days0$density, days1$density, days7$density))
    xrange = range( c(0, days0$mides, days1$mids, days7$mids, 1.3))
    if (!to.screen) {
      png(filename = file.path(outdir, "reproductive_number_today.png"))
    } else {
      dev.new()
    }

      plot(  days0$density ~ days0$mids, col="green", lwd=3, xlab="Current Reproductive number", ylab="Probability density", main="", xlim=xrange, ylim=yrange, type ="l")
      lines( days1$density ~ days1$mids, col="slateblue", lwd=3, lty="dotted")
      lines( days7$density ~ days7$mids, col="darkorange", lwd=3, lty="dashed")
      abline( v=1, col="red", lwd=3 )
      legend( "topright", "", paste( "Current date: ", stan_data$timestamp, " "), bty="n")
      legend( "right", legend=c("Today", "Yesterday", "7 days ago"), lty=c("solid", "dotted", "dashed"), col=c("green", "slateblue", "darkorange"), lwd=c(3,3,3), bty="n")
      title( main= paste( stan_data$plotlabel, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()

  }

  if (selection=="forecasts") {
    if (!to.screen) {
      png(filename = file.path(outdir, "fit_with_projections_and_stochastic_simulations.png"))
    } else {
      dev.new()
    }
      simxval = stan_data$Nobs + c(1:nprojections)
      xrange = c(0, max(nx, simxval) )
      io = stan_data$Iobs
      io[io < 0] = NA
      ip = apply(M$I, 2, median)[1:nx]
      yrange = range( c(io, ip), na.rm=TRUE)
      plot( io ~ stan_data$time, xlim=xrange, ylim=yrange,  type="n", ylab="Infected", xlab="Days")
      for ( i in 1:min(nsims, 2000)) lines( sim[i,2,] ~ simxval, col=alpha("slategray", 0.1), ltyp="dashed" )
      lines( ip ~ seq(1,nx), lwd = 3, col="slateblue" )
      lines( apply(M$I, 2, quantile, probs=0.025)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      lines( apply(M$I, 2, quantile, probs=0.975)[1:nx] ~ seq(1,nx), col="darkorange", lty="dashed", lwd = 2 )
      points( io ~ stan_data$time, xlim=xrange, ylim=yrange, col="darkgray", cex=1.2 )
      if (is.nova.scotia) {
        abline( v=stan_data$time[stan_data$Nobs], col="grey", lty="dashed" )
        abline( v=stan_data$time[time_distancing], col="orange", lty="dotted" )
        abline( v=stan_data$time[time_relaxation], col="green", lty="dotted" )
        legend( "topleft", "", "\n\n   [-- Social distancing -->", bty="n" )
        legend( "topleft", "", "\n\n\n                              [-- Parks open -->", bty="n" )
       }
      title( main= paste( stan_data$plotlabel, "Start: ", stan_data$time_start, "  Current date: ", stan_data$timestamp ) )
    if (!to.screen) dev.off()

  }


}





