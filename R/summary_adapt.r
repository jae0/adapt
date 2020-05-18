
summary_adapt = function( option="summary", can, fn=NULL ) {

  if (is.null(fn) ) fn = file.path(getwd(), "Covid19Canada_summary.rdata")  # default to current work directory
  outdir = dirname(fn)
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings=FALSE, recursive=TRUE )


  if (option=="summary" ) {

    res = list()

    for (province in names(can) ) {

      # province = names(can)[i]
      print(province)
      Npop = can[[province]]$Npop

      fn_model = file.path( workdir, paste( province, can[[province]]$modelname, "rdata", sep=".") )
      outdir = file.path( "~", "bio", "adapt", "inst", "doc", province)

      load(fn_model)
      M = extract(f)

      res[[province]]$time = can[[province]]$time
      res[[province]]$Npop = can[[province]]$Npop
      res[[province]]$Nobs = can[[province]]$Nobs
      res[[province]]$Npreds = can[[province]]$Npreds
      res[[province]]$Nts = res[[province]]$Nobs + res[[province]]$Npreds
      res[[province]]$timeall = 1:res[[province]]$Nts
      res[[province]]$S = data.frame( cbind(
        median = apply(M$S/res[[province]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$S/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$S/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[province]]$I = data.frame( cbind(
        median = apply(M$I/res[[province]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$I/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$I/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[province]]$R = data.frame( cbind(
        median = apply(M$R/res[[province]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$R/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$R/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[province]]$M = data.frame( cbind(
        median = apply(M$M/res[[province]]$Npop, 2, median, na.rm=TRUE),
        low = apply(M$M/res[[province]]$Npop, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$M/res[[province]]$Npop, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[province]]$GAMMA = data.frame( cbind(
        median = apply(t(M$GAMMA), 1, median, na.rm=TRUE),
        low = apply(t(M$GAMMA), 1, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(t(M$GAMMA), 1, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[province]]$EPSILON = data.frame( cbind(
        median = apply(t(M$EPSILON), 1, median, na.rm=TRUE),
        low = apply(t(M$EPSILON), 1, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(t(M$EPSILON), 1, quantile, probs=c(0.975), na.rm=TRUE)
      ))
      res[[province]]$K = data.frame( cbind(
        median = apply(M$K, 2, median, na.rm=TRUE),
        low = apply(M$K, 2, quantile, probs=c(0.025), na.rm=TRUE),
        high = apply(M$K, 2, quantile, probs=c(0.975), na.rm=TRUE)
      ))
    }

    save ( res, file=fn, compress=TRUE )
  }


  load( fn)

  if (option %in% c( "plot", "plot.XXX") ) {

    for (province in names(can) ) {

      plot( 0,0)
      lines( res[[province]]$I$median ~ res[[province]]$timeall  )
    }

  }

  return (res)

}

