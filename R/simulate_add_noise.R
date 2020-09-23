#' @title simulate_add_noise
#' @description This is a placeholder for a description.
#' @param sim default is sim
#' @param nobs default is \code{30}
#' @return  This is a placeholder for what it returns.
#' @author Jae Choi, \email{choi.jae.seok@gmail.com}
#' @export
simulate_add_noise = function( sim, nobs=30 ) {

  subsample = sort( sample.int( length(1:50), nobs, replace=FALSE ) )

  # missing data
  missing = which(!is.finite(sim$S[subsample] ))
  trange = range(sim$time[subsample])
  out = data.frame(time=seq(trange[1], trange[2], by=1), S=NA, I=NA, R=NA )
  with_data = match( sim$time[subsample], out$time)
  out[with_data, "S" ] = as.integer(sim$S[subsample])
  out[with_data, "I" ] = as.integer(sim$I[subsample])
  out[with_data, "R" ] = as.integer(sim$R[subsample])

  out$S[ which(!is.finite(out$S)) ] = -1 # reset NAs to Inf as stan does not take NAs
  out$I[ which(!is.finite(out$I)) ] = -1 # reset NAs to Inf as stan does not take NAs
  out$R[ which(!is.finite(out$R)) ] = -1 # reset NAs to Inf as stan does not take NAs

  return(out)

}
