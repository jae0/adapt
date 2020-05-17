
# Current COVID-19 status by province in Canada:

You can see the current status of COVID-19 disease progression by province, where data are available, by choosing from the links below.

- ![Alberta](./inst/doc/Alberta/README.md)
- ![British Columbia](./inst/doc/BC/README.md)
- ![Manitoba](./inst/doc/Manitoba/README.md)
- ![New Brunswick](./inst/doc/New%20Brunswick/README.md)
- ![Newfoundland and Labrador](./inst/doc/NL/README.md)
- ![Nova Scotia](./inst/doc/Nova%20Scotia/README.md)
- ![Ontario](./inst/doc/Ontario/README.md)
- ![Quebec](./inst/doc/Quebec/README.md)
- ![Yukon](./inst/doc/Yukon/README.md)
- ![Northwest Territories](./inst/doc/NWT/README.md)


More information about the models and data can be found on the ![main page](https://github.com/jae0/adapt/blob/master/README.md)

Please note that these results are generated from an automated process. There might be problems due to unforeseen issues. I will keep tweaking and updating this as much as possible.


---


# The data

The data are obtained directly from: https://github.com/ishaberry/Covid19Canada. This is a data repository maintained by:

Berry I, Soucy J-PR, Tuite A, Fisman D. Open access epidemiologic data and an interactive dashboard to monitor the COVID-19 outbreak in Canada. CMAJ. 2020 Apr 14;192(15):E420. doi: https://doi.org/10.1503/cmaj.75262

and assimilated in https://github.com/jae0/adapt/blob/master/R/data_provinces_of_canada.R .


---

# The model

The number of cases are modelled as a latent state-space variable in a variant of the compartmental SIR model. The variation is that Mortalities are separated from Recovered people and the infection rate parameter is modelled as an autoregressive AR(K) process. For the purposes of these results, a K=3 day lag is used, a balance between computational time and stabilization of the estimates.

Mean-field projections from this recursive model are presented with 95% posterior credible intervals, based upon an average and standard deviation of the infection rate parameter over the last K days (i.e., 3 days).

A stochastic simulation (Master Equation-based) from these parameters are also presented.



---

# What is "adapt"?

*adapt* (Areal Disease Analysis and Predicion Tools) is a set of routines to analyze publicly available disease epidemic data such as COVID-19 that you can customize for your town, province or state or country. This data tends to be rather crude counts of cases and recovered people and deaths. Your area of interest probably has these announcements and the information is likely captured by concerned citizens. To make sense of this information, beyond the daily ups and downs, you need to model it. To help, *adapt* attempts to assimilate disease spread data for COVID-19 and other similar diseases and estimate model parameters of classical epidemiological models using Bayesian methods (MCMC via STAN: https://mc-stan.org/) and then simulate/forcast disease progression using stochastic simulation modelling approaches (via the Gillespie method, implemented beautifully by SimInf: https://github.com/stewid/SimInf). Areal unit-based approaches using CAR/BYM models (via INLA) are also (soon) being developed to help understand and model spatial patterns.

A basic model appropriate for such crude data is the SIR (Susceptible-Infected-Recovered) compartmental model. This is a well understood model that has its share of limitations but still sufficient to get a crude sense of what is going on. Look it up if you want to know details. The programs here are used to fit a variation of this model to the crude data, as best we can. It relies upon some advanced statistical and mathematical engines developed by opensourced, cutting edge research groups (R, STAN, SimInf, and many others, see below), and approaches the problem as a "latent, state-space" problem, often encountered also in fisheries assessment problems. Use of "adapt", however, requires only some minimal understanding of programming, mostly R (https://cran.r-project.org/).

If you just want to get a sense of what things are like for your area of interest, you need to change the input data (see below for examples). For Nova Scotia's status, it is accessing a Google sheet that stores the required information:

https://docs.google.com/spreadsheets/d/1tgf2H9gDmRnGDGeQE-fC9IPhrmNxP8-JC7Nnnob_vuY/edit#gid=1323236978

This data was/is compiled and updated by Jennifer Strang and Nathalie Saint-Jacques. You can use this as a template.

Alternatively, you can directly estimate the numbers required and manually create the data structures (see example in and assimilated in https://github.com/jae0/adapt/blob/master/R/data_provinces_of_canada.R ).

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

Ultimately, you just need to create a data list with the information required: number of infected people on a daily basis ("InfectedCurrently" in the spreadsheet) as well as the cummulative number of "Recoveries" and "Deaths" on a daily basis. You will also need the total population size of your area of interest. Look inside the function (https://github.com/jae0/adapt/blob/master/R/data_nova_scotia.R) to see how it is done here. Use the Nova Scotia example as a template. Thereafter, you can probably run the short code in https://github.com/jae0/adapt/blob/master/inst/scripts/example_parameter_estimation_SIR_nova_scotia.R with minimal modification.


Here is an example of the data structure that is expected in "stan_data":


```

R> str(stan_data)
List of 14
 $ Npop       : num 971395
 $ Nobs       : int 51
 $ Npreds     : num 30
 $ time       : int [1:51] 1 2 3 4 5 6 7 8 9 10 ...
 $ Sobs       : num [1:51] 971394 971392 971390 971390 -1 ...
 $ Iobs       : num [1:51] 1 3 5 5 -1 28 41 51 68 73 ...
 $ Robs       : num [1:51] 0 0 0 0 -1 0 0 0 0 0 ...
 $ Mobs       : num [1:51] 0 0 0 0 -1 0 0 0 0 0 ...
 $ BNP        : num 3
 $ modelname  : chr "discrete_autoregressive_without_observation_error"
 $ plotlabel  : chr "Nova Scotia"

```
