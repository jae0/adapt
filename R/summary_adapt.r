
summary_adapt = function( selection="summary.load", aus, fn=NULL, workdir=getwd(), modelname="default", brks=30, to.screen=TRUE ) {

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19Canada_summary.rdata")  # default to current work directory
  outdir = dirname(fn)
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )

  res = NULL

  if (selection=="summary.create" ) {

    res = list()

    for (au in aus ) {

      print(au)

      fn_model = file.path( workdir, paste( au, modelname, "rdata", sep=".") )
      outdir = file.path( "~", "bio", "adapt", "inst", "doc", au)

      load(fn_model)

      posteriors = rstan::extract( stan_results$stan_samples )  # posteriors = mcmc posteriors from STAN

      res[[au]]$timestamp = stan_results$stan_inputs$timestamp
      res[[au]]$time = stan_results$stan_inputs$time
      res[[au]]$Npop = stan_results$stan_inputs$Npop
      res[[au]]$Nobs = stan_results$stan_inputs$Nobs
      res[[au]]$Npreds = stan_results$stan_inputs$Npreds
      res[[au]]$Nts = res[[au]]$Nobs + res[[au]]$Npreds
      res[[au]]$timeall = 1:res[[au]]$Nts

      res[[au]]$S = data.frame( cbind(
        median = apply(posteriors$S/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(posteriors$S/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(posteriors$S/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$I = data.frame( cbind(
        median = apply(posteriors$I/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(posteriors$I/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(posteriors$I/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$R = data.frame( cbind(
        median = apply(posteriors$R/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(posteriors$R/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(posteriors$R/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$M = data.frame( cbind(
        median = apply(posteriors$M/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(posteriors$M/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(posteriors$M/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$GAMMA = data.frame( cbind(
        median = apply(t(posteriors$GAMMA), 1, median, na.rm=TRUE),
        low = apply(t(posteriors$GAMMA), 1, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(t(posteriors$GAMMA), 1, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$EPSILON = data.frame( cbind(
        median = apply(t(posteriors$EPSILON), 1, median, na.rm=TRUE),
        low = apply(t(posteriors$EPSILON), 1, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(t(posteriors$EPSILON), 1, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$K = data.frame( cbind(
        median = apply(posteriors$K, 2, median, na.rm=TRUE),
        low = apply(posteriors$K, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(posteriors$K, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$histogram_K0 = hist( posteriors$K[,res[[au]]$Nobs] , breaks=brks, plot=FALSE)
      res[[au]]$histogram_K1 = hist( posteriors$K[,res[[au]]$Nobs-1] , breaks=brks, plot=FALSE)
      res[[au]]$histogram_K7 = hist( posteriors$K[,res[[au]]$Nobs-7] , breaks=brks, plot=FALSE)
    }

    save ( res, file=fn, compress=TRUE )
  }


  if (selection=="summary.load" ) {
    if (file.exists(fn)) {
      load(fn)
    } else {
      res = summary_adapt( "summary.create", aus=aus, fn=fn, workdir=workdir, modelname=modelname )
    }
    return(res)
  }

  if (is.null(res)) res = summary_adapt( "summary.load", aus=aus, fn=fn, workdir=workdir, modelname=modelname )

  if ( grepl("plot", selection))  {

    colours = 1:length(aus)
    ltypes = 1:length(aus)
    pchs = 1:length(aus)

    nx = 0
    for (i in 1:length(aus) ) {
      nx = max( nx, length(res[[aus[i]]]$time) + 3 )
    }

    if (selection %in% c( "plot_all", "plot_susceptible") ) {
      if (!to.screen) {
        png(filename = file.path(outdir, "susceptible.png"))
      } else {
        dev.new()
      }
        # using cubic folded root to visualize
        frp = 1/2  # folded root power
        xvals = seq(0, nx, by=20)
        yrange = folded_root(c(0.96, 1), frp)
        yticks = seq( yrange[1], yrange[2], length.out=8)
        yvals = round( folded_root( yticks, frp, inverse=TRUE) *100 , 2)  # convert to %
        plot( 0,0, type="n", xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
        for (i in 1:length(aus) ) {
          au = aus[i]
          lines( folded_root(res[[au]]$S$median, frp) ~ res[[au]]$timeall, col=alpha(colours[i], 0.9), lty=ltypes[i]   )
        }
        axis( 1, at=xvals )
        axis( 2, at=yticks, labels=yvals )
        legend( "bottomleft", legend=aus, col=colours, lty=ltypes, bty="n" )
        title( xlab="Time (days)", ylab="Percent susceptible (Folded root=1/2)" )
      if (!to.screen) dev.off()
    }

    if (selection %in% c( "plot_all", "plot_infected") ) {
      if (!to.screen) {
        png(filename = file.path(outdir, "infected.png"))
      } else {
        dev.new()
      }
        # using cubic folded root to visualize
        frp = 1/2  # folded root power
        xvals = seq(0, nx, by=20)
        yrange = folded_root(c(0, 0.003), frp)
        yticks = seq( yrange[1], yrange[2], length.out=8)
        yvals = round( folded_root( yticks, frp, inverse=TRUE) *100 , 2)  # convert to %
        plot( 0,0, type="n", xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
        for (i in 1:length(aus) ) {
          au = aus[i]
          lines( folded_root(res[[au]]$I$median, frp) ~ res[[au]]$timeall, col=alpha(colours[i], 0.9), lty=ltypes[i]   )
        }
        axis( 1, at=xvals )
        axis( 2, at=yticks, labels=yvals )
        legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
        title( xlab="Time (days)", ylab="Percent infected (Folded root=1/2)" )
      if (!to.screen) dev.off()
    }


    if (selection %in% c( "plot_all", "plot_recovered") ) {
      if (!to.screen) {
        png(filename = file.path(outdir, "recovered.png"))
      } else {
        dev.new()
      }
        # using cubic folded root to visualize
        frp = 1  # folded root power
        xvals = seq(0, nx, by=20)
        yrange = folded_root(c(0, 0.0012), frp)
        yticks = seq( yrange[1], yrange[2], length.out=8)
        yvals = round( folded_root( yticks, frp, inverse=TRUE) *100 , 2)  # convert to %
        plot( 0,0, type="n", xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
        for (i in 1:length(aus) ) {
          au = aus[i]
          lines( folded_root(res[[au]]$R$median, frp) ~ res[[au]]$timeall, col=alpha(colours[i], 0.9), lty=ltypes[i]   )
        }
        axis( 1, at=xvals )
        axis( 2, at=yticks, labels=yvals )
        legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
        title( xlab="Time (days)", ylab="Percent recovered (Folded root=1)" )
      if (!to.screen) dev.off()
    }


    if (selection %in% c( "plot_all", "plot_mortalities") ) {

      if (!to.screen) {
        png(filename = file.path(outdir, "mortalities.png"))
      } else {
        dev.new()
      }
        # using cubic folded root to visualize
        frp = 1/2  # folded root power
        xvals = seq(0, nx, by=20)
        yrange = folded_root(c(0, 0.00045), frp)
        yticks = seq( yrange[1], yrange[2], length.out=8)
        yvals = round( folded_root( yticks, frp, inverse=TRUE) *100 , 2)  # convert to %
        plot( 0,0, type="n", xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
        for (i in 1:length(aus) ) {
          au = aus[i]
          lines( folded_root(res[[au]]$M$median, frp) ~ res[[au]]$timeall, col=alpha(colours[i], 0.9), lty=ltypes[i]   )
        }
        axis( 1, at=xvals )
        axis( 2, at=yticks, labels=yvals )
        legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
        title( xlab="Time (days)", ylab="Percent mortality (Folded root=1/2)" )
        abline(h =1, col="lightgray")
      if (!to.screen) dev.off()
    }


    if (selection %in% c( "plot_all", "plot_EPSILON_GAMMA") ) {

      if (!to.screen) {
        png(filename = file.path(outdir, "EPSILON_GAMMA.png"))
      } else {
        dev.new()
      }
        xrange = NULL
        yrange = NULL
        for ( i in 1:length(aus)) {
          au = aus[i]
          xrg0 = range( res[[au]]$GAMMA$median, na.rm=TRUE )
          xrange = range( c(xrange, xrg0 )) #, yrg1, yrg7 ) )
          yrg0 = range( res[[au]]$EPSILON$median, na.rm=TRUE )
          yrange = range( c(yrange, yrg0 )) #, yrg1, yrg7 ) )
        }
        yrange[1] =yrange[1]*0.8
        yrange[2] =yrange[2]*1.2
        xrange[1] =xrange[1]*0.8
        xrange[2] =xrange[2]*1.2

        xvals = round( seq(xrange[1], xrange[2], length.out=5 ), 3)
        yvals = round( seq(yrange[1], yrange[2], length.out=5 ), 3)
        plot( 0,0, type="n",  xlab="", ylab="", ylim=yrange, xlim=xrange, axes=FALSE)
        for (i in 1:length(aus) ) {
          au = aus[i]
          points( res[[au]]$EPSILON$median ~ res[[au]]$GAMMA$median, pch=pchs[i], col=alpha(colours[i], 0.95), cex=1.5  )
        }
        axis( 1, at=xvals )
        axis( 2, at=yvals )
        legend( "topright", legend=aus, col=colours, pch=pchs, bty="n", cex=1.25 )
        title( ylab="Mortality rate constant (EPSILON)", xlab="Recovery rate constant (GAMMA)" )
      if (!to.screen) dev.off()
    }


    if (selection %in% c( "plot_all", "plot_reproductive_number") ) {

      if (!to.screen) {
        png(filename = file.path(outdir, "reproductive_number.png"))
      } else {
        dev.new()
      }
        # using cubic folded root to visualize
        frp = 1/3  # folded root power
        KMAX = 100

        yrange = 0
        for ( i in 1:length(aus)) {
          au = aus[i]
          yrg0 = range( res[[au]]$K$median, na.rm=TRUE )
          yrange = range( c(yrange, yrg0 )) #, yrg1, yrg7 ) )
        }
        yrange[2] = yrange[2] * 1.1
        yrange = folded_root(yrange/KMAX, frp)

        xvals = seq(0, nx, by=20)

        yticks = seq( yrange[1], yrange[2], length.out=8)
        yvals = round( folded_root( yticks, frp, inverse=TRUE)* KMAX )

        plot( 0,0, type="n", xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
        for (i in 1:length(aus) ) {
          au = aus[i]
          lines( folded_root(res[[au]]$K$median/KMAX, frp) ~ res[[au]]$timeall[-1], col=alpha(colours[i], 0.9), lty=ltypes[i], lwd=2   )
       }
        axis( 1, at=xvals )
        axis( 2, at=yticks, labels=yvals )
        legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
        title( xlab="Time (days)", ylab="Reproductive number (Folded root=1/3)" )
        abline( h =folded_root(1/KMAX, frp), col="lightgray")
      if (!to.screen) dev.off()
    }

    if (selection %in% c( "plot_all", "plot_reproductive_number_histograms") ) {
      brks = 30
      yrange = 0
      for ( i in 1:length(aus)) {
        au = aus[i]
        yrg0 = range( res[[au]]$histogram_K0$density, na.rm=TRUE )
        # yrg1 = range( res[[au]]$histogram_K1$density, na.rm=TRUE )
        # yrg7 = range( res[[au]]$histogram_K7$density, na.rm=TRUE )
        yrange = range( c(yrange, yrg0 )) #, yrg1, yrg7 ) )
      }
      yrange[2] = yrange[2] * 1.1
      yrange[2] = min(10, yrange[2])
      xrange = range( c(0, res[[au]]$histogram_K0$mids) )
      xrange[2] = xrange[2] * 1.25
      xvals = round( seq(xrange[1], xrange[2], length.out=8 ), 3)
      yvals = round( seq(yrange[1], yrange[2], length.out=8 ), 3)

      if (!to.screen) {
        png(filename = file.path(outdir, "reproductive_number_histograms.png"))
      } else {
        dev.new()
      }
        plot(  0,0, type="n", xlab="Reproductive number", ylab="Probability density", main="", xlim=xrange, ylim=yrange, axes=FALSE )
        for (i in 1:length(aus)) {
          au = aus[i]
          lines( res[[au]]$histogram_K0$density ~ res[[au]]$histogram_K0$mids,
            col=alpha(colours[i], 0.95), lty=ltypes[i], cex=1.5)
        }
        abline( v=1, col="red" )
        abline( h=0, col="gray", lwd=1 )
        axis( 1, at=xvals )
        axis( 2, at=yvals )
        legend( "topright", legend=aus, col=colours, lty=ltypes, bty="n", cex=1.25, lwd=2 )
        title( main= paste( "  Current date: ", res[[1]]$timestamp ) )
      if (!to.screen) dev.off()

    }

  }


  return (res)

}

