plot_model_fit = function( stan_data, M=NULL, selection="default"  ) {
  # M = mcmc posteriors from STAN

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

}





