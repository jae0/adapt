
summary_adapt = function( option="summary.load", can, fn=NULL, to.screen=TRUE ) {

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19Canada_summary.rdata")  # default to current work directory
  outdir = dirname(fn)
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )

  res = NULL

  if (option=="summary.create" ) {

    res = list()

    for (au in names(can) ) {

      # au = names(can)[i]
      print(au)
      Npop = can[[au]]$Npop

      fn_model = file.path( workdir, paste( au, can[[au]]$modelname, "rdata", sep=".") )
      outdir = file.path( "~", "bio", "adapt", "inst", "doc", au)

      load(fn_model)
      M = extract(f)

      res[[au]]$time = can[[au]]$time
      res[[au]]$Npop = can[[au]]$Npop
      res[[au]]$Nobs = can[[au]]$Nobs
      res[[au]]$Npreds = can[[au]]$Npreds
      res[[au]]$Nts = res[[au]]$Nobs + res[[au]]$Npreds
      res[[au]]$timeall = 1:res[[au]]$Nts
      res[[au]]$S = data.frame( cbind(
        median = apply(M$S/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$S/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$S/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$I = data.frame( cbind(
        median = apply(M$I/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$I/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$I/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$R = data.frame( cbind(
        median = apply(M$R/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$R/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$R/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$M = data.frame( cbind(
        median = apply(M$M/res[[au]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$M/res[[au]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$M/res[[au]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$GAMMA = data.frame( cbind(
        median = apply(t(M$GAMMA), 1, median, na.rm=TRUE),
        low = apply(t(M$GAMMA), 1, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(t(M$GAMMA), 1, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$EPSILON = data.frame( cbind(
        median = apply(t(M$EPSILON), 1, median, na.rm=TRUE),
        low = apply(t(M$EPSILON), 1, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(t(M$EPSILON), 1, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[au]]$K = data.frame( cbind(
        median = apply(M$K, 2, median, na.rm=TRUE),
        low = apply(M$K, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$K, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
    }

    save ( res, file=fn, compress=TRUE )
  }


  if (option=="summary.load" ) {
    if (file.exists(fn)) {
      load(fn)
    } else {
      res = summary_adapt( "summary.create", can=can, fn=fn )
    }
    return(res)
  }

  if (is.null(res)) res = summary_adapt( "summary.load", can=can, fn=fn )

  if (option %in% c( "plot", "plot.infected") ) {

    aus = names(can)
    colours = 1:length(aus)
    ltypes = 1:length(aus)
    pchs = 1:length(aus)

    nx = min( can[[1]]$Nobs, length(res[[1]]$time) )


    if (!to.screen) {
      png(filename = file.path(outdir, "infected.png"))
    } else {
      dev.new()
    }
      # using cubic folded root to visualize
      frp = 1/3  # folded root power
      xvals = seq(0, nx, by=10)
      yrange = folded_root(c(0, 0.0135), frp)
      yticks = seq( yrange[1], yrange[2], length.out=8)
      yvals = round( folded_root( yticks, frp, inverse=TRUE) *100 , 2)  # convert to %
      plot( 0,0, xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
      for (i in 1:length(aus) ) {
        au = aus[i]
        lines( folded_root(res[[au]]$I$median[1:nx], frp) ~ res[[au]]$timeall[1:nx], col=alpha(colours[i], 0.8), lty=ltypes[i], lwd=3   )
      }
      axis( 1, at=xvals )
      axis( 2, at=yticks, labels=yvals )
      legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
      title( xlab="Time (days)", ylab="Percent infected (Folded root=1/3)" )
    if (!to.screen) dev.off()


    if (!to.screen) {
      png(filename = file.path(outdir, "recovered.png"))
    } else {
      dev.new()
    }
      # using cubic folded root to visualize
      frp = 1  # folded root power
      xvals = seq(0, nx, by=10)
      yrange = folded_root(c(0, 0.0012), frp)
      yticks = seq( yrange[1], yrange[2], length.out=8)
      yvals = round( folded_root( yticks, frp, inverse=TRUE) *100 , 2)  # convert to %
      plot( 0,0, xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
      for (i in 1:length(aus) ) {
        au = aus[i]
        lines( folded_root(res[[au]]$R$median[1:nx], frp) ~ res[[au]]$timeall[1:nx], col=alpha(colours[i], 0.8), lty=ltypes[i], lwd=3   )
      }
      axis( 1, at=xvals )
      axis( 2, at=yticks, labels=yvals )
      legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
      title( xlab="Time (days)", ylab="Percent recovered (Folded root=1)" )
    if (!to.screen) dev.off()


    if (!to.screen) {
      png(filename = file.path(outdir, "mortalities.png"))
    } else {
      dev.new()
    }
      # using cubic folded root to visualize
      frp = 1/2  # folded root power
      xvals = seq(0, nx, by=10)
      yrange = folded_root(c(0, 0.00045), frp)
      yticks = seq( yrange[1], yrange[2], length.out=8)
      yvals = round( folded_root( yticks, frp, inverse=TRUE) *100 , 2)  # convert to %
      plot( 0,0, xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
      for (i in 1:length(aus) ) {
        au = aus[i]
        lines( folded_root(res[[au]]$M$median[1:nx], frp) ~ res[[au]]$timeall[1:nx], col=alpha(colours[i], 0.8), lty=ltypes[i], lwd=3   )
      }
      axis( 1, at=xvals )
      axis( 2, at=yticks, labels=yvals )
      legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
      title( xlab="Time (days)", ylab="Percent mortality (Folded root=1/2)" )
      abline(h =1, col="lightgray", lwd=3)
    if (!to.screen) dev.off()



    if (!to.screen) {
      png(filename = file.path(outdir, "reproductive_number.png"))
    } else {
      dev.new()
    }
      xvals = seq(0, nx, by=10)
      yrange = c(0, 10)
      yvals = seq(yrange[1], yrange[2], by =2)
      plot( 0,0, xlab="", ylab="", ylim=yrange, xlim=c(0, nx), axes=FALSE)
      for (i in 1:length(aus) ) {
        au = aus[i]
        lines( res[[au]]$K$median[1:nx] ~ res[[au]]$timeall[1:nx], col=alpha(colours[i], 0.6), lty=ltypes[i], lwd=2   )
      }
      axis( 1, at=xvals )
      axis( 2, at=yvals )
      legend( "topleft", legend=aus, col=colours, lty=ltypes, bty="n" )
      title( xlab="Time (days)", ylab="Reproductive number" )
      abline(h =1, col="lightgray", lwd=3)
    if (!to.screen) dev.off()



    if (!to.screen) {
      png(filename = file.path(outdir, "EPSILON_GAMMA.png"))
    } else {
      dev.new()
    }
      yrange = c(0,0.0057)
      xrange = c(0, 0.114)
      xvals = round( seq(xrange[1], xrange[2], length.out=8 ), 3)
      yvals = round( seq(yrange[1], yrange[2], length.out=8 ), 3)
      plot( 0,0,  xlab="", ylab="", ylim=yrange, xlim=xrange, axes=FALSE)
      for (i in 1:length(aus) ) {
        au = aus[i]
        points( res[[au]]$EPSILON$median ~ res[[au]]$GAMMA$median, pch=pchs[i], col=alpha(colours[i], 0.75), cex=1.5  )
      }
      axis( 1, at=xvals )
      axis( 2, at=yvals )
      legend( "topright", legend=aus, col=colours, pch=pchs, bty="n", cex=1.25 )
      title( ylab="EPSILON (Mortality rate constant)", xlab="GAMMA (Recovery rate constant)" )
    if (!to.screen) dev.off()

  }


  return (res)

}

