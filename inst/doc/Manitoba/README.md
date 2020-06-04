
---

# Current COVID-19 status for Manitoba:

Here we show a few elements of COVID-19 status by province based on data compiled from publicly available information.

The figures are generated from https://github.com/jae0/adapt/blob/master/inst/scripts/example_parameter_estimation_SIR_provinces_of_Canada.R. More information about the models and data can be found on the ![main page](https://github.com/jae0/adapt/blob/master/README.md)

Please note that these results are generated from an automated process. There might be problems due to unforeseen issues. I will keep tweaking and updating this as much as possible.

---

# Infected number of people with simple projections

![](./fit_with_projections_infected.png)

The number of infected people as a function of time (days) in circles. Vertical line represents "today". The blue line shown is the model fit to a modified SIR model with 95% Credible Intervals in orange  and posterior distributions in cyan. Simple deterministic (mean-field) forecasts from the recursive model are shown.

---

# Recovered number of people with simple projections

![](./fit_with_projections_recovered.png)

The number of recovered people as a function of time (days) in circles. Vertical line represents "today". The blue line shown is the model median fit to a modified SIR model with posterior distributions in cyan. Simple deterministic (mean-field) forecasts from the recursive model are shown.

---

# Number of deaths with simple projections

![](./fit_with_projections_mortalities.png)

The number of deaths as a function of time (days) in circles. Vertical line represents "today". The blue line shown is the model fit to a modified SIR model with 95% Credible Intervals in orange  and posterior distributions in cyan. Simple deterministic (mean-field) forecasts from the recursive model are shown.

---

# Infected vs total affected population

![](./infected_affected.png)

The number of infected people as a function of the total number of people affected by Covid-19 (infected, recovered, mortality), on log10 - log10 scales. This is commonly used to pinpoint the moment of the flattening of the infection curve (i.e., when it deviates from an exponential increase). To help identify days, they are coloured in alternating sequences.  These are posterior estimates derived from the model fit. The median value is shown as a solid line.

---

# Reproductive number

![](./reproductive_number.png)

How the reproductive number has been changing over the course of the epidemic. If this value is below the critical value of 1, then disease spread is being controlled. If it is above 1, an epidemic is more likely. The blue line shown is the model fit to a modified SIR model with 95% Credible Intervals in orange and posterior distributions in cyan.

![](./reproductive_number_today.png)

The current and recent estimates of the reproductive number (posterior distribution) in relation to the critical value of  1 (red line)!

---
# Forecast with stochastic simulations

![](./fit_with_projections_and_stochastic_simulations.png)

Here, individual trajectories of stochastic simulations are shown. These are based upon the joint posterior distributions of the parameter estimates for the most "current day", obtained from the above analysis. Similar to the simple projections, "current" is defined as the mean over the last BNP=3 days. These trajectories represent possible futures, accounting for small number stochasticity (unlike the mean-field ODE-based "simple" predictions). This essentially  amounts to assuming that the current "situation" remains constant/consistent (i.e., control measures and population behaviours encapsualted in the joint-posterior distributions of the model parameters).

