
---

# Current COVID-19 status for Nova Scotia:

Here we show a few elements of COVID-19 status in Nova Scotia based on data compiled from publicly available information (https://novascotia.ca/coronavirus/data/). This has been kindly compiled by J. Strang and N. St-Jacques and stored at: https://docs.google.com/spreadsheets/d/1tgf2H9gDmRnGDGeQE-fC9IPhrmNxP8-JC7Nnnob_vuY/edit#gid=1323236978

The following figures are generated from https://github.com/jae0/adapt/blob/master/inst/scripts/example_parameter_estimation_SIR_nova_scotia.R.

The use of smaller areal units in a spatially explicit mode would of course be more informative, however, this is all that is currently available.

The number of cases are modelled as a latent variable in a recursive SIR model and then projected forward using parameter estimates from the last few days of observations, in this case, three days.

I will keep updating this as much as possible.

---


# Infected number of people with simple projections

![](./inst/doc/fit_with_projections_infected.png)

The number of infected people as a function of time (days; day 1 is 2020-03-17) in circles. Vertical line represents "today". The blue line shown is the model fit to a modified SIR model with 95% Credible Intervals in orange. Simple deterministic (mean-field) forecasts from the recursive model are shown.

---

# Recovered number of people with simple projections

![](./inst/doc/fit_with_projections_recovered.png)

The number of recovered (and dead) people as a function of time (days; day 1 is 2020-03-17) in circles. Vertical line represents "today". The blue line shown is the model fit to a modified SIR model with 95% Credible Intervals in orange. Simple deterministic (mean-field) forecasts from the recursive model are shown.



---

# Reproductive number

![](./inst/doc/reproductive_number.png)

How the reproductive number has been changing over the course of the epidemic. Currently it is below the critical value of 1. Above this number, an epidemic is more likely.

![](./inst/doc/reproductive_number_today.png)

The current estimate of the reproductive number (posterior distribution). There is still a non-zero probablity that the Reproductive number is larger than 1 (red line), but the majority of the distribution is currently below 1.


---
# Forecast with stochastic simulations

![](./inst/doc/fit_with_projections_and_stochastic_simulations.png)

Finally, individual trajectories of stochastic simulations are shown. These are based upon the joint posterior distributions of the parameter estimates for the most "current day", obtained from the above analysis. These trajectories represent possible futures, accounting for small number stochasticity (unlike the mean-field ODE-based "forward" predictions), assuming the current "situation" remains constant/consistent (i.e., control measures and population behaviours encapsualted in the joint-posterior distributions of the model parameters) . The majority of the trajectories suggest a continued decline in the number of infected people. Nonetheless, there still exists some scenarios of rapid disease outbreak.



---

# What is "adapt"?

*adapt* (Areal Disease Analysis and Predicion Tools) is a set of routines to analyze publicly available disease epidemic data such as COVID-19 that you can customize for your town, province or state or country. This data tends to be rather crude counts of cases and recovered people and deaths. Your area of interest probably has these announcements and the information is likely captured by concerned citizens. To make sense of this information, beyond the daily ups and downs, you need to model it. To help, *adapt* attempts to assimilate disease spread data for COVID-19 and other similar diseases and estimate model parameters of classical epidemiological models using Bayesian methods (MCMC via STAN: https://mc-stan.org/) and then simulate/forcast disease progression using stochastic simulation modelling approaches (via the Gillespie method, implemented beautifully by SimInf: https://github.com/stewid/SimInf). Areal unit-based approaches using CAR/BYM models (via INLA) are also (soon) being developed to help understand and model spatial patterns.

A basic model appropriate for such crude data is the SIR (Susceptible-Infected-Recovered) compartmental model. This is a well understood model that has its share of limitations but still sufficient to get a crude sense of what is going on. Look it up if you want to know details. The programs here are used to fit a variation of this model to the crude data, as best we can. It relies upon some advanced statistical and mathematical engines developed by opensourced, cutting edge software (R, STAN, SimInf, and many others, see below), and approaches the problem as a "latent, state-space" model, often encountered in fisheries-modelling problems. Use of "adapt", however, requires only some minimal understanding of programming, mostly R (https://cran.r-project.org/).

If you just want to get a sense of what things are like for your area of interest, you need to change only the data (see below). Right now, it is accessing a Google sheet that stores the required information:

https://docs.google.com/spreadsheets/d/1tgf2H9gDmRnGDGeQE-fC9IPhrmNxP8-JC7Nnnob_vuY/edit#gid=1323236978

This data was manually compiled and upated by J. Strang and N. St-Jacques. You can use this as a template. Alternatively, you can estimate the numbers required (see example below).

Please note: No guarantees are being made here. There are always errors in models, programs that implement such models and in the data itself. However, this is a functional way of helping make sense of information such that we can engage in more informed discussions with your community on next steps in these trying times.

Good luck,

Jae

---

# Installation

To install you need to install R, and then bootstrap from github directly:

```
  remotes::install_github( "jae0/adapt" )
```

and also the Rpackages, "rstan" and "SimInf". They will pull in their own dependencies.

Ultimately, you just need to create a data list with the information required: number of infected people on a daily basis ("InfectedCurrently" in the spreadsheet) as well as the cummulative number of "Recoveries" and "Deaths" and dates. You will also need the total population size of your area of interest. Look inside the function (https://github.com/jae0/adapt/blob/master/R/data_nova_scotia.R) to see how it is done here. Use the Nova Scotia example as a template. Thereafter, you can probably run the short code in https://github.com/jae0/adapt/blob/master/inst/scripts/example_parameter_estimation_SIR_nova_scotia.R with minimal modification.


Here is an example of the data structure that is expected in "stan_data":


```

R> str(stan_data)
List of 14
 $ Npop       : num 971395
 $ Nobs       : int 51
 $ Npreds     : num 30
 $ Sobs       : num [1:51] 971394 971392 971390 971390 -1 ...
 $ Iobs       : num [1:51] 1 3 5 5 -1 28 41 51 68 73 ...
 $ Robs       : num [1:51] 0 0 0 0 -1 0 0 0 0 0 ...
 $ time       : int [1:51] 1 2 3 4 5 6 7 8 9 10 ...
 $ time_pred  : int [1:81] 1 2 3 4 5 6 7 8 9 10 ...
 $ t0         : num -0.01
 $ BNP        : num 3
 $ GAMMA_prior: num 0.0357
 $ BETA_prior : num 0.9
 $ modelname  : chr "discrete_autoregressive_without_observation_error"

```
