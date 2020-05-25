sir_stan_model_code = function( selection="default"  ) {


  if ( selection %in% c("default", "discrete_autoregressive_structured_beta_mortality_hybrid") ) {

    ## this tried to add the binomial data costraints but STAN really does not permit integers as rando variables and
    ## so the probabilities are computed only for post-processing

    return(

"
data {
  //declare variables
  int<lower=0> Npop;  // Npop total
  int<lower=0> Nobs;  //number of time slices
  int<lower=0> Npreds;  //additional number of time slices for prediction
  int<lower=0> BNP; // the last no days to use for BETA to project forward
  real<lower=0> BETA_max; // this value is important, seldom does this value go > 1 for Covid-19 in Canada, if too large then convergence is slow and error distributions become very wide when infected numbers ->0
  real<lower=0> GAMMA_max; //
  real<lower=0> EPSILON_max; //
  int Sobs[Nobs]; // observed S
  int Iobs[Nobs]; // observed I
  int Robs[Nobs]; // observed Recovered (including deaths  .. because this generally tends to be how it is reported)
  int Mobs[Nobs]; // observed mortality (Deaths only)
}

transformed data {
  int Ntimeall;
  int Nobs_1;
  int BNP1;
  real<lower = 0.0, upper =1.0> Sprop[Nobs]; // observed S in proportion of total pop
  real<lower = 0.0, upper =1.0> Iprop[Nobs]; // observed I
  real<lower = 0.0, upper =1.0> Rprop[Nobs]; // observed R excluding deaths  ..  chaning meaning of R here (vs Robs)
  real<lower = 0.0, upper =1.0> Mprop[Nobs]; // observed mortalities

  Nobs_1 = Nobs - 1;
  Ntimeall = Nobs + Npreds;
  BNP1 = BNP+1;

  // * 1.0 is fiddling to convert int to real
  // checking for > 0 is to check for missing values == -1
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 ) {
      Sprop[i] = fmin(1.0, fmax( 0.0, (Sobs[i] * 1.0) / (Npop * 1.0) )) ; // observation error .. some portion of infected is not captured
    } else {
      Sprop[i]=0.0; //dummy value
    }
    if ( Iobs[i] >= 0) {
      Iprop[i] = fmin(1.0, fmax( 0.0, (Iobs[i]* 1.0 )/ ( Npop * 1.0) )) ;
    } else {
      Iprop[i]=0.0; //dummy value
    }
    if (Mobs[i] >= 0) {
      Mprop[i] = fmin(1.0, fmax( 0.0, (Mobs[i]* 1.0 )/ (Npop* 1.0) )) ;  // deaths
    } else {
      Mprop[i]=0.0; //dummy value
    }
    if (Robs[i] >= 0 && Mobs[i] >= 0) {
      Rprop[i] = fmin(1.0, fmax( 0.0, ( (Robs[i] - Mobs[i])*1.0)/ (Npop* 1.0) ));  // recoveries only (with no deaths)
    } else {
      Rprop[i]=0.0; //dummy value
    }
  }
}

parameters {
  real<lower=1.0e-9, upper =GAMMA_max> GAMMA;     // recovery rate .. proportion of infected recovering
  real<lower=1.0e-9, upper =EPSILON_max> EPSILON;   // death rate .. proportion of infected dying
  real<lower=0.0, upper =BETA_max> BETA[Nobs_1];  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate
  real<lower = -1.0, upper =1.0> BETAar[BNP];
  real<lower = -1.0, upper =1.0> BETAark;  // BETA of AR(0) >=0 is sensible
  real<lower = 1.0e-9, upper =0.1 > BETAsd;
  real<lower = 1.0e-9, upper =0.2>  Ssd;  // these are fractional .. i.e CV's
  real<lower = 1.0e-9, upper =0.2>  Isd;
  real<lower = 1.0e-9, upper =0.2>  Rsd;
  real<lower = 1.0e-9, upper =0.2>  Msd;
  real<lower = 0.0, upper =1.0> Smu[Nobs]; // mean process S
  real<lower = 0.0, upper =1.0> Imu[Nobs]; // mean process I
  real<lower = 0.0, upper =1.0> Rmu[Nobs]; // mean process Recoveries only (no deaths)
  real<lower = 0.0, upper =1.0> Mmu[Nobs]; // mean process Mortalities
  real<lower = 0.0, upper =2.0> Q[Nobs_1]; // unobserved
}

transformed parameters{
}

model {
  // non informative hyperpriors (process error)
  Ssd ~ cauchy(0.0, 0.1);
  Isd ~ cauchy(0.0, 0.1);
  Rsd ~ cauchy(0.0, 0.1);
  Msd ~ cauchy(0.0, 0.1);

  GAMMA ~ normal(0, 0.1/4.0);;  // recovery of I ... always < 1, shrinks towards 0
  EPSILON ~ normal(0, 0.1/4.0);;  // recovery of I ... always < 1, shrinks towards 0

  // AR(k=BNP) model for BETA
  BETAar ~ normal( 0.0, 0.25 ); // autoregression (AR(k=BNP)) ..  shrink to 0
  BETAark ~ cauchy( 0.0, 0.1 ); //, shrinks towards 0

  BETAsd ~ cauchy( 0.0, 0.1 ); // , shrinks towards 0
  BETA[1:BNP] ~ normal( 0.0, 1.0/4.0 );  //  centered on 0, shrink towards 0

  Q ~ normal( 1.0, 0.1 );  //  centered on 0, shrink towards 0

  //set intial conditions
  Smu[1] ~ normal(Sprop[1], Ssd) ;
  Imu[1] ~ normal(Iprop[1], Isd) ;
  Rmu[1] ~ normal(Rprop[1], Rsd) ;
  Mmu[1] ~ normal(Mprop[1], Msd) ;

  // process model
  for ( i in 1:Nobs_1 ) {
    if ( i >= BNP1 ) {
      real BETAmu = BETAark;
      for ( j in 1:BNP) {
        BETAmu += BETAar[j] * step( Imu[i] ) * BETA[i-j];
      }
      BETA[i] ~ normal( BETAmu, BETAsd );
    }
    Smu[i+1] ~ normal( Smu[i] - BETA[i] *  Smu[i] * Imu[i] * Q[i] , Ssd)  ;
    Imu[i+1] ~ normal( Imu[i] + BETA[i] *  Smu[i] * Imu[i] * Q[i] - GAMMA * Imu[i]* Q[i] - EPSILON * Imu[i]* Q[i] , Isd);
    Rmu[i+1] ~ normal( Rmu[i] + GAMMA * Imu[i]* Q[i] , Rsd ) ;
    Mmu[i+1] ~ normal( Mmu[i] + EPSILON * Imu[i]* Q[i] , Msd) ;
  }

  // data likelihoods, if *obs ==-1, then data was missing  . same conditions as in transformed parameters
  // observation model with binomial observation error: slow .. swithcing to normal

  for (i in 1:Nobs) {
    if (Sobs[i] >= 0  ) {  // to handle missing values in SI
      Sprop[i] ~ normal( Smu[i] , Ssd );
      // Sobs[i] ~ binomial( Npop, Smu[i] );  // slow
    }
    if (Iobs[i] >= 0 ) {
      Iprop[i] ~ normal( Imu[i], Isd );
      // Iobs[i] ~ binomial( Npop, Imu[i] );
    }
    if (Robs[i] >= 0 ) {
      Rprop[i]  ~ normal( Rmu[i], Rsd );
      // Robs[i] ~ binomial( Npop, Rmu[i] );
    }
    if (Mobs[i] >= 0 ) {
      Mprop[i]  ~ normal( Mmu[i], Msd );
      // Mobs[i] ~ binomial( Npop, Mmu[i] );
    }
  }
}


generated quantities {
  real<lower=0> K[Ntimeall-1];
  int<lower = 0, upper =Npop> S[Ntimeall]; // latent S
  int<lower = 0, upper =Npop> I[Ntimeall]; // latent I
  int<lower = 0, upper =Npop> R[Ntimeall]; // latent R (no mortality)
  int<lower = 0, upper =Npop> M[Ntimeall]; // latent M (mortality)
  real<lower = 0.0, upper =1.0> Spp[Npreds+1]; // mean process S
  real<lower = 0.0, upper =1.0> Ipp[Npreds+1]; // mean process I
  real<lower = 0.0, upper =1.0> Rpp[Npreds+1]; // mean process Recoveries only (no deaths)
  real<lower = 0.0, upper =1.0> Mpp[Npreds+1]; // mean process Mortalities

  // predicted observations
  for (i in 1:Nobs) {
    S[i] = binomial_rng( Npop, Smu[i] );
    I[i] = binomial_rng( Npop, Imu[i] );
    R[i] = binomial_rng( Npop, Rmu[i] );
    M[i] = binomial_rng( Npop, Mmu[i] );
  }

  // initial conditions
  Spp[1] = Smu[Nobs];
  Ipp[1] = Imu[Nobs];
  Rpp[1] = Rmu[Nobs];
  Mpp[1] = Mmu[Nobs];

  for ( i in 1:Npreds ) {
    real dsi = BETA[Nobs_1] * Spp[i] * Ipp[i] * Q[Nobs_1] ;
    real dir = GAMMA * Ipp[i]* Q[Nobs_1] ;
    real dim = EPSILON * Ipp[i] * Q[Nobs_1] ;
    Spp[i+1] = fmax(0, fmin( 1, Spp[i] - dsi ) )  ;
    Ipp[i+1] = fmax(0, fmin( 1, Ipp[i] + dsi - dir - dim  ));
    Rpp[i+1] = fmax(0, fmin( 1, Rpp[i] + dir )) ;
    Mpp[i+1] = fmax(0, fmin( 1, Mpp[i] + dim )) ;
  }

  // predicted observations
  for ( i in 1:Npreds ) {
    S[Nobs+i] = binomial_rng( Npop, Spp[i+1] );
    I[Nobs+i] = binomial_rng( Npop, Ipp[i+1] );
    R[Nobs+i] = binomial_rng( Npop, Rpp[i+1] );
    M[Nobs+i] = binomial_rng( Npop, Mpp[i+1] );
  }

  // sample from  mean process (proportions to counts)
  for (i in 1:Nobs_1 ) {
    K[i] = BETA[i] / (GAMMA+EPSILON); // the contact number = fraction of S in contact with I
  }
  for (i in Nobs:(Ntimeall-1) ) {
    K[i] = BETA[Nobs_1] / (GAMMA+EPSILON); // the contact number = fraction of S in contact with I
  }

}
"
    )  # end return
  } # end selection


# ------------------


# ------------------



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
  real<lower = 0.0> params[2];      // Model parameters = 2 in SIR
  real<lower = 0.0, upper = 1> S0;  // Initial fraction of hosts susceptible
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
  real<lower = 0.0, upper =1.0> Sprop[Nobs]; // observed S
  real<lower = 0.0, upper =1.0> Iprop[Nobs]; // observed I
  real<lower = 0.0, upper =1.0> Rprop[Nobs]; // observed R

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
  real<lower=1.0e-9, upper=1> GAMMA;
  real<lower=0> BETA;
  real<lower = 1.0e-9, upper =0.5>  Ssd;
  real<lower = 1.0e-9, upper =0.5>  Isd ;
  real<lower = 1.0e-9, upper =0.5>  Rsd;
}

transformed parameters{
  real<lower = 0.0, upper =1.0> Smu[Nobs]; // mean process S
  real<lower = 0.0, upper =1.0> Imu[Nobs]; // mean process I
  real<lower = 0.0, upper =1.0> Rmu[Nobs]; // mean process R

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

  real<lower = 0.0, upper =1.0> Sp[Npreds+1]; // mean process S (predictive)
  real<lower = 0.0, upper =1.0> Ip[Npreds+1]; // mean process I
  real<lower = 0.0, upper =1.0> Rp[Npreds+1]; // mean process R

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



# -----------------------

  if ( selection=="discrete_autoregressive_structured_beta_mortality" ) {
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
  real<lower = 0.0, upper =1.0> Sprop[Nobs]; // observed S in proportion of total pop
  real<lower = 0.0, upper =1.0> Iprop[Nobs]; // observed I
  real<lower = 0.0, upper =1.0> Rprop[Nobs]; // observed R excluding deaths  ..  chaning meaning of R here (vs Robs)
  real<lower = 0.0, upper =1.0> Mprop[Nobs]; // observed mortalities

  Ntimeall = Nobs + Npreds;

  // * 1.0 is fiddling to comvert int to real
  // checking for > 0 is to check for missing values == -1
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0 ) {
      Sprop[i] = (Sobs[i] * 1.0 )/ (Npop * 1.0)  ;
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
  real<lower=0.0, upper =0.1> GAMMA;     // recovery rate .. proportion of infected recovering
  real<lower=0.0, upper =0.01> EPSILON;   // death rate .. proportion of infected dying
  real<lower=0.0, upper  =1> BETA[Ntimeall-1];  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate
  real<lower = 1.0e-9, upper =0.2>  Ssd;  // these are fractional .. i.e CV's
  real<lower = 1.0e-9, upper =0.2>  Isd;
  real<lower = 1.0e-9, upper =0.2>  Rsd;
  real<lower = 1.0e-9, upper =0.2>  Msd;
  real<lower = -1.0, upper =1.0> ar1;
  real<lower = 1.0e-9, upper =1 > BETAsd;
  real ar1k;
  real <lower =0, upper =1 > latent0[4];
}

transformed parameters{
  real<lower = 0.0, upper =1.0> Smu[Ntimeall]; // mean process S
  real<lower = 0.0, upper =1.0> Imu[Ntimeall]; // mean process I
  real<lower = 0.0, upper =1.0> Rmu[Ntimeall]; // mean process Recoveries only (no deaths)
  real<lower = 0.0, upper =1.0> Mmu[Ntimeall]; // mean process Mortalities

  // process model: SIR ODE in Euler difference form
  // recursive model eq for discrete form of SIR;
  // some fraction of recovered die (rather than directly from infected),
  // this is due to large a latency between infection and death (30 days+),
  // using Recovered as it is closer to the timescale of the deaths

  Smu[1] = latent0[1];
  Imu[1] = latent0[2];
  Rmu[1] = latent0[3];
  Mmu[1] = latent0[4];

  for (i in 1:(Ntimeall-1) ) {
    real dSI = BETA[i] * Smu[i] * Imu[i];
    real dIR = GAMMA * Imu[i];
    real dIM = EPSILON * Imu[i];
    Smu[i+1] = Smu[i] - dSI ;
    Imu[i+1] = Imu[i] + dSI - dIR - dIM;
    Rmu[i+1] = Rmu[i] + dIR ;
    Mmu[i+1] = Mmu[i] + dIM ;
  }
}

model {

  // non informative hyperpriors
  Ssd ~ cauchy(0, 0.5);
  Isd ~ cauchy(0, 0.5);
  Rsd ~ cauchy(0, 0.5);
  Msd ~ cauchy(0, 0.5);

  //set intial conditions
  latent0[1] ~ normal(Sprop[1], Ssd) ;
  latent0[2] ~ normal(Iprop[1], Isd) ;
  latent0[3] ~ normal(Rprop[1], Rsd) ;
  latent0[4] ~ normal(Mprop[1], Msd) ;

  ar1 ~ normal(0, 1); // autoregression
  BETAsd ~ cauchy(0, 0.5);
  ar1k ~ cauchy(0, 0.5);

  BETA[1] ~ normal( BETA_prior, BETA_prior );  // # 10% CV
  for (i in 1:(Nobs-2)) {
    BETA[i+1] ~ normal( ar1k + ar1 * BETA[i], BETAsd );
  }
  for ( i in (Nobs):(Ntimeall-1) ) {
    BETA[i] ~ normal( mean( BETA[(Nobs-1-BNP):(Nobs-1)] ), sd( BETA[(Nobs-1-BNP):(Nobs-1)] ) );
  }

  GAMMA ~ normal( GAMMA_prior, GAMMA_prior );  // recovery of I ... always < 1
  EPSILON ~ normal( EPSILON_prior, EPSILON_prior );  // recovery of I ... always < 1

  // data likelihoods, if *obs ==-1, then data was missing  . same conditions as in transformed parameters
  // observation model:
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0  ) {  // to handle missing values in SI
      Sprop[i] ~ normal( Smu[i] , Ssd );
    }
    if (Iobs[i] >= 0 ) {
      Iprop[i]  ~ normal( Imu[i], Isd );
    }
    if (Robs[i] >= 0 ) {
      Rprop[i]  ~ normal( Rmu[i], Rsd );
    }
    if (Mobs[i] >= 0 ) {
      Mprop[i]  ~ normal( Mmu[i], Msd );
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


# -----------------------


# -----------------------

  if ( selection=="discrete_autoregressive_structured_beta_mortality_testing" ) {

    ## this tried to add the binomial data costraints but STAN really does not permit integers as rando variables and
    ## so the probabilities are computed only for post-processing

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
  real<lower = 0.0, upper =1.0> Sprop[Nobs]; // observed S in proportion of total pop
  real<lower = 0.0, upper =1.0> Iprop[Nobs]; // observed I
  real<lower = 0.0, upper =1.0> Rprop[Nobs]; // observed R excluding deaths  ..  chaning meaning of R here (vs Robs)
  real<lower = 0.0, upper =1.0> Mprop[Nobs]; // observed mortalities

  int<lower=0, upper = Npop> z_si[Nobs-1];
  int<lower=0, upper = Npop> z_ir[Nobs-1];
  int<lower=0, upper = Npop> z_im[Nobs-1];

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


  for ( i in 1:(Nobs-1)){
    if (Sobs[i] >=0 && Sobs[i+1] >= 0) {
     z_si[i] = Sobs[i] - Sobs[i+1]; // infected dynamics (and Susceptibles)
    } else {
      z_si[i] = 0;  // dummy value
    }
    if (Robs[i] >=0 && Robs[i+1] >=0 && Mobs[i] >=0 && Mobs[i+1] >=0) {
      z_ir[i] = Robs[i+1] - Robs[i] - Mobs[i+1] + Mobs[i];  // recoveries only (excluding deaths)
    } else {
      z_ir[i] = 0;
    }
    if (Mobs[i] >=0 && Mobs[i+1] >=0) {
      z_im[i] = Mobs[i+1] - Mobs[i]; // death dynamics
    } else {
      z_im[i] = 0;
    }
  }


}

parameters {
  real<lower=0.0, upper =0.5> GAMMA;     // recovery rate .. proportion of infected recovering
  real<lower=0.0, upper =0.1> EPSILON;   // death rate .. proportion of infected dying
  real<lower=0.0, upper  =1> BETA[Ntimeall-1];  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate
  real<lower = 1.0e-9, upper =0.2>  Ssd;  // these are fractional .. i.e CV's
  real<lower = 1.0e-9, upper =0.2>  Isd;
  real<lower = 1.0e-9, upper =0.2>  Rsd;
  real<lower = 1.0e-9, upper =0.2>  Msd;

  real<lower = -1.0, upper =1.0> ar1;
  real<lower = 1.0e-9, upper =0.2 > BETAsd;
  real ar1k;
  real<lower = 0.0, upper =1.0> latent0[4];
}

transformed parameters{
  real<lower = 0.0, upper =1.0> Smu[Ntimeall]; // mean process S
  real<lower = 0.0, upper =1.0> Imu[Ntimeall]; // mean process I
  real<lower = 0.0, upper =1.0> Rmu[Ntimeall]; // mean process Recoveries only (no deaths)
  real<lower = 0.0, upper =1.0> Mmu[Ntimeall]; // mean process Mortalities

  real<lower = 0.0, upper =0.5> dSI[Ntimeall];
  real<lower = 0.0, upper =0.5> dIR[Ntimeall];
  real<lower = 0.0, upper =0.5> dIM[Ntimeall];

  real<lower=0.0, upper =1> pr_si[Ntimeall-1];
  real<lower=0.0, upper =1> pr_ir;
  real<lower=0.0, upper =1> pr_im;

  //set intial conditions
  Smu[1] = latent0[1] ;
  Imu[1] = latent0[2] ;
  Rmu[1] = latent0[3] ;
  Mmu[1] = latent0[4] ;

  for (i in 1:(Ntimeall-1) ) {
    dSI[i] = BETA[i] * Smu[i] * Imu[i] ;
    dIR[i] = GAMMA * Imu[i]  ;
    dIM[i] = EPSILON  * Imu[i] ;
    Smu[i+1] = Smu[i] - dSI[i]   ;
    Imu[i+1] = Imu[i] + dSI[i] - dIR[i] - dIM[i] ;
    Rmu[i+1] = Rmu[i] + dIR[i] ;
    Mmu[i+1] = Mmu[i] + dIM[i] ;
  }

  // pr_ir[i] = 1/GAMMA;
  pr_ir = 1.0 - exp(-GAMMA);
  pr_im = 1.0 - exp(-EPSILON);

  // process model: SIR ODE in Euler difference form
  // recursive model eq for discrete form of SIR;
  // some fraction of recovered die (rather than directly from infected),
  // this is due to large a latency between infection and death (30 days+),
  // using Recovered as it is closer to the timescale of the deaths
  for (i in 1:(Nobs-1) ) {
    // pr_si[i] = 1-(1-BETA[i])^Imu[i]*Npop;  // per capita probability
    pr_si[i] = 1.0 - exp( -BETA[i] * Imu[i] ); // approximation
  }
  for ( i in (Nobs):(Ntimeall-1) ) {
    pr_si[i] = mean( pr_si[(Nobs-1-BNP):(Nobs-1)] ) ;
  }

}

model {

  // non informative hyperpriors
  Ssd ~ cauchy(0, 0.5);
  Isd ~ cauchy(0, 0.5);
  Rsd ~ cauchy(0, 0.5);
  Msd ~ cauchy(0, 0.5);

  ar1 ~ normal(0, 1); // autoregression
  BETAsd ~ cauchy(0, 0.5);
  ar1k ~ cauchy(0, 0.5);

  BETA[1] ~ normal( BETA_prior, BETA_prior );  // # 10% CV
  for (i in 1:(Nobs-2)) {
    BETA[i+1] ~ normal( ar1k + ar1 * BETA[i], BETAsd );
  }
  for ( i in (Nobs):(Ntimeall-1) ) {
    BETA[i] ~ normal( mean( BETA[(Nobs-1-BNP):(Nobs-1)] ), sd( BETA[(Nobs-1-BNP):(Nobs-1)] ) );
  }

  GAMMA ~ normal( GAMMA_prior, GAMMA_prior );  // recovery of I ... always < 1
  EPSILON ~ normal( EPSILON_prior, EPSILON_prior );  // recovery of I ... always < 1

  //set intial conditions
  latent0[1] ~ normal(Sprop[1], Ssd) ;
  latent0[2] ~ normal(Iprop[1], Isd) ;
  latent0[3] ~ normal(Rprop[1], Rsd) ;
  latent0[4] ~ normal(Mprop[1], Msd) ;

  // data likelihoods, if *obs ==-1, then data was missing  . same conditions as in transformed parameters
  // observation model:
  for (i in 1:Nobs) {
    if (Sobs[i] >= 0  ) {  // to handle missing values in SI
      Sprop[i] ~ normal( Smu[i] , Ssd );
    }
    if (Iobs[i] >= 0 ) {
      Iprop[i]  ~ normal( Imu[i], Isd );
    }
    if (Robs[i] >= 0 ) {
      Rprop[i]  ~ normal( Rmu[i], Rsd );
    }
    if (Mobs[i] >= 0 ) {
      Mprop[i]  ~ normal( Mmu[i], Msd );
    }
  }


  // additional likelihood constraints on direct incremental differences
  // .. not using binomial as STAN does not operate on integer random variables
  // use poisson instead on increments
    for (i in 1:(Nobs-1)){
      if (Sobs[i] >=0 && Sobs[i+1] >=0) {
        z_si[i] ~ binomial( Npop, pr_si[i]*Smu[i] );  // prob of being infected
        // z_si[i] ~ poisson(  Npop*Smu[i]* pr_si[i] ) ;
        // z_si[i] ~ poisson( dSI[i]*Npop );
      }
      if (Robs[i] >=0 && Robs[i+1] >=0 && Mobs[i] >=0 && Mobs[i+1] >=0) {
       // z_ir[i] ~ binomial( Npop, pr_ir*Imu[i] );
        // z_ir[i] ~ poisson( Imu[i]*Npop*pr_ir );
        // z_ir[i] ~ poisson( dIR[i]*Npop );
     }
     if (Mobs[i] >=0 && Mobs[i+1] >=0) {
       z_im[i] ~ binomial( Npop, pr_im*Mmu[i] );
        // z_im[i] ~ poisson( Imu[i]*Npop*pr_im) ;
        //z_im[i] ~ poisson( dIM[i]*Npop) ;
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


  if ( selection=="discrete_binomial_autoregressive_nonlatent_basic" ) {

    return(
    "
data {
  //declare variables
  int<lower=0> Npop;  // Npop total
  int<lower=0> Nobs;  //number of time slices
  int<lower=0> Npreds;  //additional number of time slices for prediction
  real<lower=0> BETA_prior;
  real<lower=0> GAMMA_prior;
  int Sobs[Nobs]; // observed S
  int Iobs[Nobs]; // observed I
  int Robs[Nobs]; // observed Recovered (including deaths)
}

transformed data {
//Calculate transitions, based on SIR status
  int Ntimeall;
  int<lower=0, upper = Npop> z_si[Nobs-1];
  int<lower=0, upper = Npop> z_ir[Nobs-1];

  Ntimeall = Nobs + Npreds;

  for(i in 1:(Nobs-1)){
    z_si[i] = Sobs[i] - Sobs[i+1]; // infected dynamics (and Susceptibles)
    z_im[i] = Mobs[i+1] - Mobs[i]; // death dynamics
  }

}

parameters {
  real<lower=1.0e-9, upper =0.1> GAMMA;  // probability of transition to recovered state = 1/( duration is γ units of time) .. ie., simple geometric
  real <lower=0.0, upper  =1> BETA;  // == beta in SIR , here we do *not* separate out the Encounter Rate from the infection rate  // probability of an individual infecting another in 1 unit of time
}

transformed parameters{
  real<lower=0.0, upper =1> pr_si[Ntimeall-1];
  real<lower=0.0, upper =1> pr_ir;

  // pr_ir[i] = 1/GAMMA;
  pr_ir = 1.0 - exp(-GAMMA);

  for (i in 1:(Nobs-1) ) {
    if ( Iobs[i] > 0){
      // pr_si[i] = 1-(1-BETA)^Iobs[i];  // per capita probability
      pr_si[i] = 1.0 - exp(-BETA * (Iobs[i] *1.0) / ( Npop * 1.0) ); // approximation
    } else {
      pr_si[i] = 0.0;
    }
  }
}

model {
  BETA  ~ normal( BETA_prior, BETA_prior );
  GAMMA ~ normal( GAMMA_prior, GAMMA_prior );  // recovery of I ... always < 1
  // likelihoods
  for (i in 1:(Nobs-1)){
    if(Iobs[i] > 0){ //only define z_si when there are infections - otherwise distribution is degenerate and STAN has trouble
      z_si[i] ~ binomial( Sobs[i], pr_si[i] ); // prob of being infected
    }
    z_ir[i] ~ binomial( Iobs[i], pr_ir );
  }
}

generated quantities {
}
"
    )  # end return
  } # end selection



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
  real<lower=1.0e-9> gamma;  // probability of transition to recovered state = 1/() duration is γ units of time) .. ie., simple geometric
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
//   for (i in 1:Nobs) {
//    I[i] = binomial( Npop, )
//  }
}
"
    )  # end return
  } # end selection


  # --------------------------


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
