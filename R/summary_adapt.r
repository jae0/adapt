
summary_adapt = function( option="summary.load", can, fn=NULL ) {

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19Canada_summary.rdata")  # default to current work directory
  outdir = dirname(fn)
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )


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
    return(res)
  }


  if (option=="summary.load" ) {
    res = NULL
    if (file.exists(fn)) {
      load(fn)
    } else {
      res = summary_adapt( "summary.create", can=can, fn=fn )
    }
    return(res)
  }



  if (option %in% c( "plot", "plot.XXX") ) {

    res = summary_adapt( "summary.load", can=can, fn=fn )

    for (au in names(can) ) {

      plot( 0,0)
      lines( res[[province]]$I$median ~ res[[province]]$timeall  )
    }

  }

  return (res)

}

