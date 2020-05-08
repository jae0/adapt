sir_stan_model_code = function( selection="discrete_autoregressive_with_observation_error"  ) {

  if ( selection=="continuous" ) {

  # CREDIT to (copied from with minor changes):
  # https://jrmihalj.github.io/estimating-transmission-by-fitting-mechanistic-models-in-Stan/
  # priors are marginally informative

    return(
"

functions {
  real[] SI(
    real time,
    real[] Y,
    real[] params,
    real[] x_r,
    int[] x_i ) {
      real dYdt[3];
      dYdt[1] = - params[1] * Y[1] * Y[2];
      dYdt[2] = params[1] * Y[1] * Y[2] - params[2] * Y[2];
      dYdt[3] = params[2] * Y[2];
      return dYdt;
  }
}

data {
  int<lower = 1> Nobs;   // Number of days sampled
  int<lower = 1> Npop;       // Number of hosts sampled at each time point.
  int<lower = 1> Npreds;  // This is to generate predictions
  real t0;                // Initial time point
  int  Iobs[Nobs];          // Data
  real time[Nobs];         // Time points of data
  real time_pred[Nobs + Npreds];   // Time points for predictions
}

transformed data {
  real x_r[0];    // containers for ode integration poorly documented in STAN ...
  int  x_i[0];    // containers for ode integration poorly documented in STAN ...
  int Ntimeall;
  Ntimeall = Nobs + Npreds;
}

parameters {
  real<lower = 0> params[2];      // Model parameters = 2 in SIR
  real<lower = 0, upper = 1> S0;  // Initial fraction of hosts susceptible
}

transformed parameters{
  real <lower =0, upper = 1> Imu[Nobs, 3]; // Output from the ODE solver SI = 3
  real <lower =0, upper = 1> y0[3];           // Initial conditions for both S and I ; dim=3
  y0[1] = S0;     // S
  y0[2] = 1 - S0; // I
  y0[3] = 0;      // R
  Imu = integrate_ode_rk45(SI, y0, t0, time, params, x_r, x_i);
  for( j in 1:Nobs) Imu[j,2] = Imu[j,2]  + 1e-4;  // positive offset to balance overshooting in solver
}

model {
  params ~ cauchy(0, 0.5) ;
  S0 ~ beta( 0.5, 0.5 );
  Iobs ~ binomial( Npop, Imu[, 2]); //Imu[,2] are the fractions infected from the ODE solver
}

generated quantities {
  int<lower = 0, upper = Npop> I[Ntimeall]; // latent I
  real<lower= 0, upper = 1>  preds[Ntimeall, 3];  // 3  derivatives
  real<lower=0> K;
  // sample from  mean process (proportions to counts)
  // Generate predicted data over the whole time series:
  preds = integrate_ode_rk45(SI, y0, t0, time_pred, params, x_r, x_i);
  for( j in 1:Ntimeall) preds[j,2] = preds[j,2]  + 1e-4;  // positive offset to balance overshooting in solver
  I = binomial_rng( Npop, preds[,2] );
  K = params[1]/params[2]; // the 'contact number, the fraction of S in contact with I == basic reproductive number
}

"
    )  # end return
  } # end selection


  if ( selection=="discrete_basic" ) {
    return(

"
data {
 //declare variables
 int<lower=0> Npop;  // Npop total
 int<lower=0> Nobs;  //number of time slices
 int<lower=0> Npreds;  //additional number of time slices for prediction
 int Sobs[Nobs]; // observed S
 int Iobs[Nobs]; // observed I
 int Robs[Nobs]; // observed R
}


transformed data {
  int Ntimeall;
  real<lower = 0, upper =1> Sprop[Nobs]; // observed S
  real<lower = 0, upper =1> Iprop[Nobs]; // observed I
  real<lower = 0, upper =1> Rprop[Nobs]; // observed R

  Ntimeall = Nobs + Npreds;

  // convert to proportions
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0) {
      Sprop[i] = (Sobs[i] * 1.0)/ (Npop * 1.0)  ; // fiddling to get extra precision
    }
    if (Iobs[i] >= 0) {
      Iprop[i] = (Iobs[i]* 1.0)/ ( Npop * 1.0) ;
    }
    if (Robs[i] >= 0) {
      Rprop[i] = (Robs[i]* 1.0)/ (Npop* 1.0) ;
    }
  }

}

parameters {
  real<lower=1e-9, upper=1> GAMMA;
  real<lower=0> BETA;
  real<lower = 1e-9, upper =0.5>  Ssd;
  real<lower = 1e-9, upper =0.5>  Isd ;
  real<lower = 1e-9, upper =0.5>  Rsd;
}

transformed parameters{
  real<lower = 0, upper =1> Smu[Nobs]; // mean process S
  real<lower = 0, upper =1> Imu[Nobs]; // mean process I
  real<lower = 0, upper =1> Rmu[Nobs]; // mean process R

  // recursive model eq for discrete form of SIR; set intial conditions
  Smu[1] = Sprop[1] ;
  Imu[1] = Iprop[1] ;
  Rmu[1] = Rprop[1] ;

  for (i in 1:(Nobs-1) ) {
    Smu[i+1] = Smu[i] - BETA * Smu[i] * Imu[i];
    Imu[i+1] = Imu[i] + BETA * Smu[i] * Imu[i] - GAMMA * Imu[i];
    Rmu[i+1] = Rmu[i] + GAMMA * Imu[i];
  }
 }


model {

  Ssd ~ beta( 0.5, 0.5 );
  Isd ~ beta( 0.5, 0.5 );
  Rsd ~ beta( 0.5, 0.5 );

  // diffuse .. sim to a long tailed t-distrib with no variance
  BETA  ~ cauchy(0.0, 1.0);
  GAMMA ~ beta(0.5, 0.5);

  // data likelihoods, if *obs ==-1, then data was missing
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0) {  // this is to handle missing values
      Sprop[i] ~ normal( Smu[i], Ssd );
    }
    if (Iobs[i] >= 0) {
      Iprop[i] ~ normal( Imu[i], Isd );
    }
    if (Robs[i] >= 0) {
      Rprop[i] ~ normal( Rmu[i], Rsd );
    }
  }
 }


generated quantities {
  real<lower=0> K;
  int<lower = 0, upper =Npop> S[Ntimeall]; // latent S
  int<lower = 0, upper =Npop> I[Ntimeall]; // latent I
  int<lower = 0, upper =Npop> R[Ntimeall]; // latent R
  // sample from  mean process (proportions to counts)

  real<lower = 0, upper =1> Sp[Npreds+1]; // mean process S (predictive)
  real<lower = 0, upper =1> Ip[Npreds+1]; // mean process I
  real<lower = 0, upper =1> Rp[Npreds+1]; // mean process R

  // copy
  Sp[1] = Smu[Nobs] ;
  Ip[1] = Imu[Nobs] ;
  Rp[1] = Rmu[Nobs] ;
  for (i in 1:Npreds ) {
    Sp[i+1] = Sp[i] - BETA * Sp[i] * Ip[i];
    Ip[i+1] = Ip[i] + BETA * Sp[i] * Ip[i] - GAMMA * Ip[i];
    Rp[i+1] = Rp[i] + GAMMA * Ip[i];
  }

  S[(Nobs+1):(Nobs+Npreds)] = binomial_rng( Npop, Sp[2:(Npreds+1)] );
  I[(Nobs+1):(Nobs+Npreds)] = binomial_rng( Npop, Ip[2:(Npreds+1)] );
  R[(Nobs+1):(Nobs+Npreds)] = binomial_rng( Npop, Rp[2:(Npreds+1)] );

  S[1:Nobs] = binomial_rng( Npop, Smu );
  I[1:Nobs] = binomial_rng( Npop, Imu );
  R[1:Nobs] = binomial_rng( Npop, Rmu );

  K = BETA/GAMMA; // the 'contact number, the fraction of S in contact with I == basic reproductive number as this is normalized to 1
}

"

    )  # end return
  } # end selection


# ------------------




  if ( selection=="discrete_autoregressive_with_observation_error_structured_beta" ) {
    return(

"
data {
  //declare variables
  int<lower=0> Npop;  // Npop total
  int<lower=0> Nobs;  //number of time slices
  int<lower=0> Npreds;  //additional number of time slices for prediction
  int<lower=0> BNP; // the last no days to use for BETA to project forward
  real<lower=0> BETA_prior;
  real<lower=0> GAMMA_prior;
  int Sobs[Nobs]; // observed S
  int Iobs[Nobs]; // observed I
  int Robs[Nobs]; // observed R
}

transformed data {
  int Ntimeall;
  real<lower = 0, upper =1> Sprop[Nobs]; // observed S in proportion of total pop
  real<lower = 0, upper =1> Iprop[Nobs]; // observed I
  real<lower = 0, upper =1> Rprop[Nobs]; // observed R

  Ntimeall = Nobs + Npreds;

  // * 1.0 is fiddling to comvert int to real
  // checking for > 0 is to check for missing values == -1
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 && Iobs[i] >= 0) {
      Sprop[i] = (Sobs[i] * 1.0 )/ (Npop * 1.0)  ; // observation error .. some portion of infected is not captured
    } else {
      Sprop[i]=0.0; //dummy value
    }
    if (Sobs[i] >= 0 && Iobs[i] >= 0 && Robs[i] >= 0) {
      Iprop[i] = (Iobs[i]* 1.0 )/ ( Npop * 1.0) ;
    } else {
      Iprop[i]=0.0; //dummy value
    }
    if (Robs[i] >= 0) {
      Rprop[i] = (Robs[i]* 1.0 )/ (Npop* 1.0) ;
    } else {
      Rprop[i]=0.0; //dummy value
    }
  }
}

parameters {
  real<lower=0.0, upper =1> GAMMA;     // recovery rate .. proportion of infected recovering
  real<lower = 0, upper=1> BETA[Ntimeall-1];  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate
  real<lower=-0.2, upper =0.2>  MSErrorI;  // fractional mis-specification error due to latent, asymptompatic cases, reporting irregularities
  real<lower=-0.2, upper =0.2>  MSErrorR;  // fractional mis-specification error due to variable disease progression/secondary infections
  real<lower = 1e-9, upper =0.2>  Ssd;  // these are fractional .. i.e CV's
  real<lower = 1e-9, upper =0.2>  Isd ;
  real<lower = 1e-9, upper =0.2>  Rsd;
  real<lower = -1, upper =1> ar1;
  real<lower = 0 > ar1sd;
  real ar1k;
  real<lower=0, upper =1> sir0[3];
}

transformed parameters{
  real<lower = 0, upper =1> Smu[Ntimeall]; // mean process S
  real<lower = 0, upper =1> Imu[Ntimeall]; // mean process I
  real<lower = 0, upper =1> Rmu[Ntimeall]; // mean process R

  // process model: SIR ODE in Euler difference form
  // recursive model eq for discrete form of SIR;
  //set intial conditions
  Smu[1] = sir0[1] ; //latent S
  Imu[1] = sir0[2] ; //latent I
  Rmu[1] = sir0[3] ; //latent R
   for (i in 1:(Ntimeall-1) ) {
    real dS = BETA[i] * Smu[i] * Imu[i];
    real dI = GAMMA * Imu[i];
    Smu[i+1] = Smu[i] - dS ;
    Imu[i+1] = Imu[i] + dS - dI;
    Rmu[i+1] = Rmu[i] + dI;
  }
}

model {

  // non informative hyperpriors
  Ssd ~ cauchy(0, 0.5);
  Isd ~ cauchy(0, 0.5);
  Rsd ~ cauchy(0, 0.5);

  sir0[1] ~ normal(Sprop[1], 0.2) ;
  sir0[2] ~ normal(Iprop[1], 0.2) ;
  sir0[3] ~ normal(Rprop[1], 0.2) ;

  ar1 ~ normal( 0, 1.0 ); // autoregression
  ar1sd ~ normal( 0, 1.0 );
  ar1k ~ normal( 0, 1.0 );

  // .. MSErrorI  is the mis-specification dur to asymptomatic cases
  MSErrorI ~ normal( 0, 0.1 );  // proportion of I that are asymtomatic
  MSErrorR ~ normal( 0, 0.1 );  // proportion of R that are misclassified

  GAMMA ~ normal( GAMMA_prior, 1.0 );  // recovery of I ... always < 1
  BETA[1] ~ normal( BETA_prior, 1.0 );  // # 10% CV
  for (i in 1:(Nobs-1)) {
    BETA[i+1] ~ normal( ar1k + ar1 * BETA[i], ar1sd );
  }
  for ( i in (Nobs+1):(Ntimeall-1) ) {
    BETA[i] ~ normal( mean( BETA[(Nobs-1-BNP):(Nobs-1)] ), sd( BETA[(Nobs-1-BNP):(Nobs-1)] ) );
  }

  // data likelihoods, if *obs ==-1, then data was missing  . same conditions as in transformed parameters
  // observation model:
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 && Iobs[i] >= 0) {  // to handle missing values in SI
      (Sprop[i] + Iprop[i]*MSErrorI) ~ normal( Smu[i] , Ssd );
    }
    if (Sobs[i] >= 0 && Iobs[i] >= 0 && Robs[i] >= 0) {
      (Iprop[i] - Iprop[i]*MSErrorI - Rprop[i]*MSErrorR) ~ normal( Imu[i], Isd );
    }
    if (Robs[i] >= 0 ) {
      (Rprop[i] + Rprop[i]*MSErrorR)  ~ normal( Rmu[i], Rsd );  // assume no observation error / mis-specification error
    }
  }
}

generated quantities {
  real<lower=0> K[Ntimeall-1];
  int<lower = 0, upper =Npop> S[Ntimeall]; // latent S
  int<lower = 0, upper =Npop> I[Ntimeall]; // latent I
  int<lower = 0, upper =Npop> R[Ntimeall]; // latent R

  S = binomial_rng( Npop, Smu );
  I = binomial_rng( Npop, Imu );
  R = binomial_rng( Npop, Rmu );

  // sample from  mean process (proportions to counts)
  for (i in 1:(Ntimeall-1) ) {
    K[i] = BETA[i] / GAMMA; // the contact number = fraction of S in contact with I
  }

}
"
    )  # end return
  } # end selection


# ------------------



  if ( selection=="discrete_autoregressive_with_observation_error_unstructured_beta" ) {
    return(

"

data {
  //declare variables
  int<lower=0> Npop;  // Npop total
  int<lower=0> Nobs;  //number of time slices
  int<lower=0> Npreds;  //additional number of time slices for prediction
  int<lower=0> BNP; // the last no days to use for BETA to project forward
  real<lower=0> BETA_prior;
  real<lower=0> GAMMA_prior;
  int Sobs[Nobs]; // observed S
  int Iobs[Nobs]; // observed I
  int Robs[Nobs]; // observed R
}

transformed data {
  int Ntimeall;
  real<lower = 0, upper =1> Sprop[Nobs]; // observed S in proportion of total pop
  real<lower = 0, upper =1> Iprop[Nobs]; // observed I
  real<lower = 0, upper =1> Rprop[Nobs]; // observed R
  Ntimeall = Nobs + Npreds;
  // * 1.0 is fiddling to comvert int to real
  // checking for > 0 is to check for missing values == -1
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 && Iobs[i] >= 0) {
      Sprop[i] = (Sobs[i] * 1.0 )/ (Npop * 1.0)  ; // observation error .. some portion of infected is not captured
    } else {
      Sprop[i]=0.0; //dummy value
    }
    if (Sobs[i] >= 0 && Iobs[i] >= 0 && Robs[i] >= 0) {
      Iprop[i] = (Iobs[i]* 1.0 )/ ( Npop * 1.0) ;
    } else {
      Iprop[i]=0.0; //dummy value
    }
    if (Robs[i] >= 0) {
      Rprop[i] = (Robs[i]* 1.0 )/ (Npop* 1.0) ;
    } else {
      Rprop[i]=0.0; //dummy value
    }
  }
}

parameters {
  real<lower=0.0, upper =1> GAMMA;     // recovery rate .. proportion of infected recovering
  real<lower = 0, upper=1> BETA[Ntimeall-1];  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate
  real<lower=-0.2, upper =0.2>  MSErrorI;  // fractional mis-specification error due to latent, asymptompatic cases, reporting irregularities
  real<lower=-0.2, upper =0.2>  MSErrorR;  // fractional mis-specification error due to variable disease progression/secondary infections
  real<lower = 1e-9, upper =0.25>  Ssd;  // these are fractional .. i.e CV's
  real<lower = 1e-9, upper =0.25>  Isd ;
  real<lower = 1e-9, upper =0.25>  Rsd;
  real<lower = -1, upper =1> ar1;
  real<lower = 0, upper =0.25> ar1sd;
  real<lower = -1, upper =1> ar1k;
  real<lower=0, upper =1> sir0[3];
}

transformed parameters{
  real<lower = 0, upper =1> Smu[Ntimeall]; // mean process S
  real<lower = 0, upper =1> Imu[Ntimeall]; // mean process I
  real<lower = 0, upper =1> Rmu[Ntimeall]; // mean process R


  // recursive model eq for discrete form of SIR;
  //set intial conditions
  Smu[1] = sir0[1] ; //latent S
  Imu[1] = sir0[2] ; //latent I
  Rmu[1] = sir0[3] ; //latent R

  // process model: SIR ODE in Euler difference form
   for (i in 1:(Ntimeall-1) ) {
    real dS = BETA[i] * Smu[i] * Imu[i];
    real dI = GAMMA * Imu[i];
    Smu[i+1] = Smu[i] - dS ;
    Imu[i+1] = Imu[i] + dS - dI;
    Rmu[i+1] = Rmu[i] + dI;
  }

}


model {

  // non informative hyperpriors
  Ssd ~ cauchy(0, 0.5);
  Isd ~ cauchy(0, 0.5);
  Rsd ~ cauchy(0, 0.5);
  sir0 ~ beta(0.5, 0.5) ;

  // .. MSErrorI  is the mis-specification dur to asymptomatic cases
  MSErrorI ~ normal( 0, 0.1 );  // proportion of I that are asymtomatic
  MSErrorR ~ normal( 0, 0.1 );  // proportion of R that are misclassified

  GAMMA ~ normal( GAMMA_prior, GAMMA_prior/3.0 );  // recovery of I ... always < 1

  BETA[1:(Nobs-1)] ~ normal( BETA_prior, BETA_prior/3.0 );
  for ( i in (Nobs):(Ntimeall-1) ) {
    BETA[i] ~ normal( mean( BETA[(Nobs-1-BNP):(Nobs-1)] ), sd( BETA[(Nobs-1-BNP):(Nobs-1)] ) );
  }

  // data likelihoods, if *obs ==-1, then data was missing  . same conditions as in transformed parameters
  // observation model:
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 && Iobs[i] >= 0) {  // to handle missing values in SI
      (Sprop[i] + Iprop[i]*MSErrorI) ~ normal( Smu[i] , Ssd );
    }
    if (Sobs[i] >= 0 && Iobs[i] >= 0 && Robs[i] >= 0) {
      (Iprop[i] - Iprop[i]*MSErrorI - Rprop[i]*MSErrorR) ~ normal( Imu[i], Isd );
    }
    if (Robs[i] >= 0 ) {
      (Rprop[i] + Rprop[i]*MSErrorR)  ~ normal( Rmu[i], Rsd );  // assume no observation error / mis-specification error
    }
  }
}

generated quantities {
  real<lower=0> K[Ntimeall-1];
  int<lower = 0, upper =Npop> S[Ntimeall]; // latent S
  int<lower = 0, upper =Npop> I[Ntimeall]; // latent I
  int<lower = 0, upper =Npop> R[Ntimeall]; // latent R

  S = binomial_rng( Npop, Smu );
  I = binomial_rng( Npop, Imu );
  R = binomial_rng( Npop, Rmu );

  // sample from  mean process (proportions to counts)
  for (i in 1:(Ntimeall-1) ) {
    K[i] = BETA[i] / GAMMA; // the contact number = fraction of S in contact with I
  }
}

"
    )  # end return
  } # end selection


# ------------------

# ------------------



  if ( selection=="discrete_autoregressive_without_observation_error" ) {
    return(
"
data {
  //declare variables
  int<lower=0> Npop;  // Npop total
  int<lower=0> Nobs;  //number of time slices
  int<lower=0> Npreds;  //additional number of time slices for prediction
  int<lower=0> BNP; // the last no days to use for BETA to project forward
  real<lower=0> BETA_prior;
  real<lower=0> GAMMA_prior;
  int Sobs[Nobs]; // observed S
  int Iobs[Nobs]; // observed I
  int Robs[Nobs]; // observed R
}

transformed data {
  int Ntimeall;
  real<lower = 0, upper =1> Sprop[Nobs]; // observed S in proportion of total pop
  real<lower = 0, upper =1> Iprop[Nobs]; // observed I
  real<lower = 0, upper =1> Rprop[Nobs]; // observed R

  Ntimeall = Nobs + Npreds;

  // * 1.0 is fiddling to comvert int to real
  // checking for > 0 is to check for missing values == -1
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 ) {
      Sprop[i] = (Sobs[i] * 1.0 )/ (Npop * 1.0)  ; // observation error .. some portion of infected is not captured
    } else {
      Sprop[i]=0.0; //dummy value
    }
    if ( Iobs[i] >= 0  ) {
      Iprop[i] = (Iobs[i]* 1.0 )/ ( Npop * 1.0) ;
    } else {
      Iprop[i]=0.0; //dummy value
    }
    if (Robs[i] >= 0) {
      Rprop[i] = (Robs[i]* 1.0 )/ (Npop* 1.0) ;
    } else {
      Rprop[i]=0.0; //dummy value
    }
  }
}

parameters {
  real<lower=0.0, upper =1> GAMMA;     // recovery rate .. proportion of infected recovering
  real<lower = 0, upper=1> BETA[Ntimeall-1];  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate
  real<lower = 1e-9, upper =0.2>  Ssd;  // these are fractional .. i.e CV's
  real<lower = 1e-9, upper =0.2>  Isd ;
  real<lower = 1e-9, upper =0.2>  Rsd;
  real<lower = -1, upper =1> ar1;
  real<lower = 0 > ar1sd;
  real ar1k;
  real<lower=0, upper =1> sir0[3];
}

transformed parameters{
  real<lower = 0, upper =1> Smu[Ntimeall]; // mean process S
  real<lower = 0, upper =1> Imu[Ntimeall]; // mean process I
  real<lower = 0, upper =1> Rmu[Ntimeall]; // mean process R

  // process model: SIR ODE in Euler difference form
  // recursive model eq for discrete form of SIR;
  //set intial conditions
  Smu[1] = sir0[1] ; //latent S
  Imu[1] = sir0[2] ; //latent I
  Rmu[1] = sir0[3] ; //latent R
   for (i in 1:(Ntimeall-1) ) {
    real dS = BETA[i] * Smu[i] * Imu[i];
    real dI = GAMMA * Imu[i];
    Smu[i+1] = Smu[i] - dS ;
    Imu[i+1] = Imu[i] + dS - dI;
    Rmu[i+1] = Rmu[i] + dI;
  }
}

model {

  // non informative hyperpriors
  Ssd ~ cauchy(0, 0.5);
  Isd ~ cauchy(0, 0.5);
  Rsd ~ cauchy(0, 0.5);

  sir0[1] ~ normal(Sprop[1], 0.2) ;
  sir0[2] ~ normal(Iprop[1], 0.2) ;
  sir0[3] ~ normal(Rprop[1], 0.2) ;

  ar1 ~ normal( 0, 1.0 ); // autoregression
  ar1sd ~ normal( 0, 1.0 );
  ar1k ~ normal( 0, 1.0 );

  GAMMA ~ normal( GAMMA_prior, 1.0 );  // recovery of I ... always < 1
  BETA[1] ~ normal( BETA_prior, 1.0 );  // # 10% CV
  for (i in 1:(Nobs-1)) {
    BETA[i+1] ~ normal( ar1k + ar1 * BETA[i], ar1sd );
  }
  for ( i in (Nobs+1):(Ntimeall-1) ) {
    BETA[i] ~ normal( mean( BETA[(Nobs-1-BNP):(Nobs-1)] ), sd( BETA[(Nobs-1-BNP):(Nobs-1)] ) );
  }

  // data likelihoods, if *obs ==-1, then data was missing  . same conditions as in transformed parameters
  // observation model:
  for (i in 1:Nobs) {
    if ( Sobs[i] >= 0 ) {  // to handle missing values in SI
      (Sprop[i] ) ~ normal( Smu[i] , Ssd );
    }
    if ( Iobs[i] >= 0 ) {
      (Iprop[i] ) ~ normal( Imu[i], Isd );
    }
    if (Robs[i] >= 0 ) {
      (Rprop[i] ) ~ normal( Rmu[i], Rsd );  // assume no observation error / mis-specification error
    }
  }
}

generated quantities {
  real<lower=0> K[Ntimeall-1];
  int<lower = 0, upper =Npop> S[Ntimeall]; // latent S
  int<lower = 0, upper =Npop> I[Ntimeall]; // latent I
  int<lower = 0, upper =Npop> R[Ntimeall]; // latent R

  S = binomial_rng( Npop, Smu );
  I = binomial_rng( Npop, Imu );
  R = binomial_rng( Npop, Rmu );

  // sample from  mean process (proportions to counts)
  for (i in 1:(Ntimeall-1) ) {
    K[i] = BETA[i] / GAMMA; // the contact number = fraction of S in contact with I
  }

}
"

    )  # end return
  } # end selection

# -----------------------

  if ( selection=="discrete_autoregressive_with_observation_error_structured_beta_mortality" ) {
    return(

"
data {
  //declare variables
  int<lower=0> Npop;  // Npop total
  int<lower=0> Nobs;  //number of time slices
  int<lower=0> Npreds;  //additional number of time slices for prediction
  int<lower=0> BNP; // the last no days to use for BETA to project forward
  real<lower=0> BETA_prior;
  real<lower=0> GAMMA_prior;
  real<lower=0> EPSILON_prior;
  int Sobs[Nobs]; // observed S
  int Iobs[Nobs]; // observed I
  int Robs[Nobs]; // observed Recovered (including deaths)
  int Mobs[Nobs]; // observed mortality (Deaths only)
}

transformed data {
  int Ntimeall;
  real<lower = 0, upper =1> Sprop[Nobs]; // observed S in proportion of total pop
  real<lower = 0, upper =1> Iprop[Nobs]; // observed I
  real<lower = 0, upper =1> Rprop[Nobs]; // observed R excluding deaths  ..  chaning meaning of R here (vs Robs)
  real<lower = 0, upper =1> Mprop[Nobs]; // observed mortalities

  Ntimeall = Nobs + Npreds;

  // * 1.0 is fiddling to comvert int to real
  // checking for > 0 is to check for missing values == -1
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 ) {
      Sprop[i] = (Sobs[i] * 1.0 )/ (Npop * 1.0)  ; // observation error .. some portion of infected is not captured
    } else {
      Sprop[i]=0.0; //dummy value
    }
    if ( Iobs[i] >= 0) {
      Iprop[i] = (Iobs[i]* 1.0 )/ ( Npop * 1.0) ;
    } else {
      Iprop[i]=0.0; //dummy value
    }
    if (Mobs[i] >= 0) {
      Mprop[i] = (Mobs[i]* 1.0 )/ (Npop* 1.0) ;  // deaths
    } else {
      Mprop[i]=0.0; //dummy value
    }
    if (Robs[i] >= 0 && Mobs[i] >= 0) {
      Rprop[i] = ( (Robs[i] - Mobs[i])*1.0)/ (Npop* 1.0) ;  // recoveries only (with no deaths)
    } else {
      Rprop[i]=0.0; //dummy value
    }

  }
}

parameters {
  real<lower=0.0, upper =1> GAMMA;     // recovery rate .. proportion of infected recovering
  real<lower=0.0, upper =1> EPSILON;   // death rate .. proportion of infected dying
  real<lower=0.0, upper  =1> BETA[Ntimeall-1];  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate
  real<lower=0, upper =1>  MSErrorI;  // fractional mis-specification error due to latent, asymptompatic cases, reporting irregularities
  real<lower = 1e-9, upper =0.2>  Ssd;  // these are fractional .. i.e CV's
  real<lower = 1e-9, upper =0.2>  Isd;
  real<lower = 1e-9, upper =0.2>  Rsd;
  real<lower = 1e-9, upper =0.2>  Msd;
  real<lower = -1, upper =1> ar1;
  real<lower = 1e-9, upper =1 > ar1sd;
  real ar1k;
  real<lower=0, upper =1> sir0[4];
}

transformed parameters{
  real<lower = 0, upper =1> Smu[Ntimeall]; // mean process S
  real<lower = 0, upper =1> Imu[Ntimeall]; // mean process I
  real<lower = 0, upper =1> Rmu[Ntimeall]; // mean process Recoveries only (no deaths)
  real<lower = 0, upper =1> Mmu[Ntimeall]; // mean process Mortalities

  // process model: SIR ODE in Euler difference form
  // recursive model eq for discrete form of SIR;
  //set intial conditions
  Smu[1] = sir0[1] ; //latent S
  Imu[1] = sir0[2] ; //latent I
  Rmu[1] = sir0[3] ; //latent R
  Mmu[1] = sir0[4] ; //latent Mortalities

  for (i in 1:(Ntimeall-1) ) {
    real dSI = BETA[i] * Smu[i] * Imu[i];
    real dIR = GAMMA * Imu[i];
    real dIM = EPSILON * Imu[i] ;
    Smu[i+1] = Smu[i] - dSI ;
    Imu[i+1] = Imu[i] + dSI - dIR - dIM;
    Rmu[i+1] = Rmu[i] + dIR;
    Mmu[i+1] = Mmu[i] + dIM;
  }
}

model {

  // non informative hyperpriors
  Ssd ~ cauchy(0, 0.5);
  Isd ~ cauchy(0, 0.5);
  Rsd ~ cauchy(0, 0.5);
  Msd ~ cauchy(0, 0.5);

  sir0[1] ~ cauchy(0, 0.5) ;
  sir0[2] ~ cauchy(0, 0.5) ;
  sir0[3] ~ cauchy(0, 0.5) ;
  sir0[4] ~ cauchy(0, 0.5) ;

  ar1 ~ cauchy(0, 0.5); // autoregression
  ar1sd ~ cauchy(0, 0.5);
  ar1k ~ cauchy(0, 0.5);

  // .. MSErrorI  is the mis-specification dur to asymptomatic cases
  MSErrorI ~ cauchy(0, 0.5);  // proportion of I that are asymtomatic

  GAMMA ~ normal( GAMMA_prior, GAMMA_prior );  // recovery of I ... always < 1
  EPSILON ~ normal( EPSILON_prior, EPSILON_prior );  // recovery of I ... always < 1
  BETA[1] ~ normal( BETA_prior, BETA_prior);  // # 10% CV
  for (i in 1:(Nobs-1)) {
    BETA[i+1] ~ normal( ar1k + ar1 * BETA[i], ar1sd );
  }
  for ( i in (Nobs+1):(Ntimeall-1) ) {
    BETA[i] ~ normal( mean( BETA[(Nobs-1-BNP):(Nobs-1)] ), sd( BETA[(Nobs-1-BNP):(Nobs-1)] ) );
  }

  // data likelihoods, if *obs ==-1, then data was missing  . same conditions as in transformed parameters
  // observation model:
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 ) {  // to handle missing values in SI
      (Sprop[i] - Iprop[i]*MSErrorI) ~ normal( Smu[i] , Ssd );
    }
    if (Iobs[i] >= 0 ) {
      (Iprop[i] + Iprop[i]*MSErrorI ) ~ normal( Imu[i], Isd );
    }
    if (Robs[i] >= 0 ) {
      Rprop[i]  ~ normal( Rmu[i], Rsd );  // assume no observation error / mis-specification error
    }
    if (Mobs[i] >= 0 ) {
      Mprop[i]  ~ normal( Mmu[i], Msd );  // assume no observation error / mis-specification error
    }
  }
}

generated quantities {
  real<lower=0> K[Ntimeall-1];
  int<lower = 0, upper =Npop> S[Ntimeall]; // latent S
  int<lower = 0, upper =Npop> I[Ntimeall]; // latent I
  int<lower = 0, upper =Npop> R[Ntimeall]; // latent R (no mortality)
  int<lower = 0, upper =Npop> M[Ntimeall]; // latent M (mortality)

  S = binomial_rng( Npop, Smu );
  I = binomial_rng( Npop, Imu );
  R = binomial_rng( Npop, Rmu );
  M = binomial_rng( Npop, Mmu );

  // sample from  mean process (proportions to counts)
  for (i in 1:(Ntimeall-1) ) {
    K[i] = BETA[i] / GAMMA; // the contact number = fraction of S in contact with I
  }

}
"
    )  # end return
  } # end selection


# ------------------


# ------------------


  if ( selection=="discrete_voorman_2017_basic" ) {

  # CREDIT to (copied from with minor changes):
  # Arie Voorman, April 2017
  # https://rstudio-pubs-static.s3.amazonaws.com/270496_e28d8aaa285042f2be0c24fc915a68b2.html
  # note this is using "discrete transition probability" for each individual rather than "continuous rates" for a population

    return(
    "
data {
  //declare variables
  int<lower=0> Npop; //number of individuals
  int<lower=0> Nobs; //number of time points
  int<lower = 0, upper =Npop> Sobs[Nobs]; //SIR
  int<lower = 0, upper =Npop> Iobs[Nobs];
  int<lower = 0, upper =Npop> Robs[Nobs];
}
transformed data {
//Calculate transitions, based on SIR status
  int<lower=0, upper = Npop> z_si[Nobs-1];
  int<lower=0, upper = Npop> z_ir[Nobs-1];
  for(i in 1:(Nobs-1)){
    z_si[i] = Sobs[i] - Sobs[i+1];
    z_ir[i] = Robs[i+1] - Robs[i];
  }
}
parameters {
  real<lower=1e-9> gamma;  // probability of transition to recovered state = 1/() duration is γ units of time) .. ie., simple geometric
  real lambda;  // probability of an individual infecting another in 1 unit of time
}

model {
  lambda ~ beta(0.5, 0.5);
  for(i in 1:(Nobs-1)){
    if(Iobs[i] > 0){ //only define z_si when there are infections - otherwise distribution is degenerate and STAN has trouble
      z_si[i] ~ binomial(Sobs[i], 1-(1-lambda)^Iobs[i]); // prob of being infected
    }
    z_ir[i] ~ binomial(Iobs[i], 1/gamma);
  }
}

generated quantities {
   int I[Nobs];
   int<lower=0> r_0; // simulate a single infected individual in the population, and count secondary infections (i.e. R_0).
  r_0 = 0;
  while (1) {
    r_0 = r_0 + binomial_rng(Npop-r_0,lambda);
    if (bernoulli_rng(1/gamma)) break;
  }
  for (i in 1:Nobs) {
    I[i] = binomial( Npop, )
  }
}
"
    )  # end return
  } # end selection



  if ( selection=="discrete_voorman_2017_covariates" ) {

  # CREDIT to (copied from with minor changes):

  # Arie Voorman, April 2017
  # https://rstudio-pubs-static.s3.amazonaws.com/270496_e28d8aaa285042f2be0c24fc915a68b2.html
  # note this is using "discrete transition probability" rather than "continuous rates"

    return(
"
data {
  int<lower=0> Npop;
  int<lower=0> Nobs;
  int<lower = 0,upper=1> Sobs[Nobs,Npop];
  int<lower = 0,upper=1> Iobs[Nobs,Npop];
  int<lower = 0,upper=1> Robs[Nobs,Npop];
  matrix[Npop,Npop] fmat;
  int<lower=0,upper=1> adult[Npop];
  real<lower=0> lambda_max;
}
transformed data {
  int<lower=0, upper = 1> z_si[Nobs-1,Npop];
  int<lower=0, upper = 1> z_ir[Nobs-1,Npop];
  matrix[Nobs,Npop] I_mat;
  matrix[Nobs,Npop] fam_I;
  matrix[Nobs,Npop] com_I;
  I_mat = to_matrix(Iobs);
  fam_I = I_mat*fmat;

  for(t in 1:(Nobs-1)){
    for(i in 1:Npop){
      z_si[t,i] = Sobs[t,i] - Sobs[t+1,i];
      z_ir[t,i] = Robs[t+1,i] - Robs[t,i];
      com_I[t,i] = sum(I_mat[t]) - I_mat[t]*fmat[,i];
    }
  }
}
parameters {
  real<lower=1> gamma_adult;
  real<lower=1> gamma_child;
  real<lower=0, upper = lambda_max> lambda_family;
  real<lower=0, upper = lambda_max> lambda_community;
}
model {
real prob_infection;

for(i in 1:Npop){
  for(t in 1:(Nobs-1)){
    if( sum(I_mat[t]) > 0 && Sobs[t,i]){
      prob_infection = 1- (  pow( 1-lambda_family, fam_I[t,i])*pow(1-lambda_community,com_I[t,i]));
      z_si[t,i] ~ bernoulli(prob_infection);
    }
    if(Iobs[t,i]){
      if(adult[i]) z_ir[t,i] ~ bernoulli(1/gamma_adult);
      if(1-adult[i]) z_ir[t,i] ~ bernoulli(1/gamma_child);
    }
   }
}
}
"
    )  # end return
  } # end selection




  if ( selection=="discrete_voorman_2017_covariates_irregular_time" ) {

  # CREDIT to (copied from with minor changes):

  # Arie Voorman, April 2017
  # https://rstudio-pubs-static.s3.amazonaws.com/270496_e28d8aaa285042f2be0c24fc915a68b2.html
  # note this is using "discrete transition probability" rather than "continuous rates"


    return(
    "
data {
  int<lower=0> Npop;
  int<lower=0> Nobs;
  int<lower = 0,upper=1> Sobs[Nobs,Npop];
  int<lower = 0,upper=1> Iobs[Nobs,Npop];
  int<lower = 0,upper=1> Robs[Nobs,Npop];
  matrix[Npop,Npop] fmat;
  int<lower=0,upper=1> adult[Npop];
  real<lower=0> lambda_max;
  real<lower=0> dt[Nobs];
}
transformed data {
  int<lower=0, upper = 1> z_si[Nobs-1,Npop];
  int<lower=0, upper = 1> z_ir[Nobs-1,Npop];
  matrix[Nobs,Npop] I_mat;
  matrix[Nobs,Npop] fam_I;
  matrix[Nobs,Npop] com_I;
  I_mat = to_matrix(Iobs);
  fam_I = I_mat*fmat;

  for(t in 1:(Nobs-1)){
    for(i in 1:Npop){
      z_si[t,i] = Sobs[t,i] - Sobs[t+1,i];
      z_ir[t,i] = Robs[t+1,i] - Robs[t,i];
      com_I[t,i] = sum(I_mat[t]) - I_mat[t]*fmat[,i];
    }
  }
}
parameters {
  real<lower=1> gamma_adult;
  real<lower=1> gamma_child;
  real<lower=0, upper = lambda_max> lambda_family;
  real<lower=0, upper = lambda_max> lambda_community;
}
model {
real prob_infection;

for(i in 1:Npop){
  for(t in 1:(Nobs-1)){
    if( sum(I_mat[t]) > 0 && Sobs[t,i]){
      prob_infection = 1- (  pow( 1-lambda_family*dt[t], fam_I[t,i])*pow(1-lambda_community*dt[t],com_I[t,i]));
      z_si[t,i] ~ bernoulli(prob_infection);
    }
    if(Iobs[t,i]){
      if(adult[i]) z_ir[t,i] ~ bernoulli( 1-(1-(1/gamma_adult))^dt[t] );
      if(1-adult[i]) z_ir[t,i] ~ bernoulli( 1- (1-(1/gamma_child))^dt[t] );
    }
   }
}
}
"
    )  # end return
  } # end selection



}
