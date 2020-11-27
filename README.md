
---

# Updates have restarted due to the current outbreak in Nova Scotia.

---

# Current COVID-19 status by province in Canada, based on an *adapt* analysis:

You can see *adapt*'s assessement of the current status of COVID-19 by choosing from the links below.

- ![Alberta](./inst/doc/Alberta/README.md)
- ![British Columbia](./inst/doc/BC/README.md)
- ![Manitoba](./inst/doc/Manitoba/README.md)
- ![New Brunswick](./inst/doc/New%20Brunswick/README.md)
- ![Newfoundland and Labrador](./inst/doc/NL/README.md)
- ![Northwest Territories](./inst/doc/NWT/README.md)
- ![Nova Scotia](./inst/doc/Nova%20Scotia/README.md)
- ![Ontario](./inst/doc/Ontario/README.md)
- ![Prince Edward Island](./inst/doc/PEI/README.md)
- ![Saskatchewan](./inst/doc/Saskatchewan/README.md)
- ![Quebec](./inst/doc/Quebec/README.md)
- ![Yukon](./inst/doc/Yukon/README.md)


Some advice: These assessments would be more informative for smaller areal units as they would make management decisions more relevant. This is because dynamics in a city or a rural area is likely very different with spatial patterns related to how connected each area is to neighbouring areas and population behaviour in terms of density, mobility, distancing, service support, etc. So if you have access to information for your city/town, you might want to run this for your city/town to obtain the most relevant information.

---


# Data

Currently, the data are obtained from: https://github.com/ishaberry/Covid19Canada. This is a data repository maintained by:

Berry I, Soucy J-PR, Tuite A, Fisman D. Open access epidemiologic data and an interactive dashboard to monitor the COVID-19 outbreak in Canada. CMAJ. 2020 Apr 14;192(15):E420. doi: https://doi.org/10.1503/cmaj.75262

and assimilated in https://github.com/jae0/adapt/blob/master/R/data_provinces_of_canada.R .


---

# Analytical model

The number of cases are modelled as a latent state-space variable in a variant of the compartmental SIR model. The former is a technique we often use in fishery status and estimation problems. We deviate from the standard SIR model in that Mortalities are separated from Recovered people and the infection rate parameter is modelled as an autoregressive AR(K=1) process. A K=1 day lag is used as there exists temporal autocorelation in disease progression which if ignored can introduce biased parameter estimates; this represents a functional balance between computational time and stabilization of the estimates. No claims are made that these are "true" rate parameters, even though latent in formulation, as we do not have a good understanding of asymptomatic cases.

Nonetheless, they represent an "effective" (that is, a pragmatic and consistent) estimate of the observed state of disease progression where all dynamics are absorbed by the time-varying "effective" infection rate. We can assess this "effective" infection rate relative to the other "effective" rate processes (death and recovery) by examining the manner in which the "effective" Reproductive number evolves over time. If this is greater than 1, then it means it is spreading exponentially; if below 1, then it is decreasing exponentially; and if it is near 1, this means the disease is steady where each new infection is balanced by a recovery or death. Keeping this value low by wearing masks, social distancing, handwashing, etc. is critical and indeed can be reasonably used as an indicator of the success of such measures and eventually identification of hot and cold-spots of infection spread that require additional measures (to do next).

"Simple" projections (i.e., deterministic, "mean-field" predictions) from this recursive model are presented with 95% posterior credible intervals. A stochastic simulation from these parameters are also presented; these are based on a master equation formulation with a Gillespie approximation.

The model options and main results are created by running:

  https://github.com/jae0/adapt/blob/master/inst/scripts/example_parameter_estimation_SIR_provinces_of_Canada.R.




---

# What is "adapt"?

*adapt* (Areal Disease Analysis and Predicion Tools) is a set of routines to analyze publicly available disease epidemic data such as COVID-19 that you can customize for your town, province or state or country. This data tends to be rather crude counts of cases and recovered people and deaths. Your area of interest probably has these announcements and the information is likely captured by concerned citizens. To make sense of this information, beyond the daily ups and downs, you need to model it. To help, *adapt* attempts to assimilate disease spread data for COVID-19 and other similar diseases and estimate model parameters of classical epidemiological models using Bayesian methods (MCMC via STAN: https://mc-stan.org/) and then simulate/forecast disease progression using stochastic simulation modelling approaches (via the Gillespie method, implemented beautifully by SimInf: https://github.com/stewid/SimInf). Areal unit-based approaches using CAR/BYM models (via INLA) are also used to help understand and model spatial patterns (todo).

A basic model appropriate for such crude data is the SIR (Susceptible-Infected-Recovered) compartmental model. This is a well understood model that has its share of limitations but still sufficient to get a crude sense of what is going on. Look it up if you want to know details. The programs here are used to fit a variation of this model to the crude data, as best we can. It diverges from standard usage in that the default model's infection rate parameter is treated as an autoregressive AR(K=1) process, rather than as a static, disease-specific constant. It, therefore, absorbs the temporal variability in disease dynamics that a simple SIR cannot model. This "effective" infection rate and associated "effective" reproductive number can help us understand disease progression. More formally, *adapt*  approaches the estimation of these and associated parameters as a "latent, state-space" problem, also often encountered also in fisheries assessment problems, where due to the size of the oceans, reality can not generally be directly observed.

Use of *adapt*, requires only some minimal understanding of programming, mostly R (https://cran.r-project.org/). If you just want to get a sense of what things are like for your area of interest, you need to change the input data (see below for examples). For Nova Scotia's status, it is accessing a Google sheet that stores the required information:

https://docs.google.com/spreadsheets/d/1tgf2H9gDmRnGDGeQE-fC9IPhrmNxP8-JC7Nnnob_vuY/edit#gid=1323236978

This was kindly compiled and updated by Jennifer Strang and Nathalie Saint-Jacques, though no longer maintained. You can use this as a template.

Alternatively, you can directly estimate the numbers required and manually create the data structures (see example in and assimilated in https://github.com/jae0/adapt/blob/master/R/data_provinces_of_canada.R ).

Please note: No guarantees are being made here. There are always errors in models, programs that implement such models and in the data itself. However, this is a functional way of helping make sense of information such that we can engage in more informed discussions with your community on next steps in these trying times.

Finally, computations require a lot of time. Currently, due to the length of the timeseries and the complexity of dynamics in some provinces, it can take upwards of 10 hours to finish some estimates (on my 4 year old laptop). Shorter term models runs updating only the recent dynamics are likely better operationally and so you might want to truncate the timeseries. I will try to get this done soon.

Please send a note if you are using it for your area of interest.

Best wishes,
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
 $ Npop       : num 971395  # total population size
 $ Nobs       : int 51      # number of observed time slices (days)
 $ Npreds     : num 30      # number of prediction time slices (days)
 $ time       : int [1:51] 1 2 3 4 5 6 7 8 9 10 ... # time index
 $ Sobs       : num [1:51] 971394 971392 971390 971390 -1 ... # total susceptible population size
 $ Iobs       : num [1:51] 1 3 5 5 -1 28 41 51 68 73 ...      # infected population size
 $ Robs       : num [1:51] 0 0 0 0 -1 0 0 0 0 0 ...           # total recovered population size
 $ Mobs       : num [1:51] 0 0 0 0 -1 0 0 0 0 0 ...           # total mortalities population size
 $ BNP        : num 1                                         # number of lags in BETA (infection rate); AR(K=BNP)

```
