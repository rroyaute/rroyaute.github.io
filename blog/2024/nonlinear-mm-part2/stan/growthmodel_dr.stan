data{
     vector[1260] l;
    array[1260] int t;
     vector[1260] Dose;
    array[1260] int ID;
}
parameters{
     vector<lower=0>[210] a0;
     vector<lower=0>[210] b0;
     vector<lower=0>[210] c0;
     real a1;
     real b1;
     real c1;
     real<lower=0> mu_a0;
     real<lower=0> mu_b0;
     real<lower=0> mu_c0;
     real<lower=0> sigma_a0;
     real<lower=0> sigma_b0;
     real<lower=0> sigma_c0;
     real<lower=0> sigma_obs;
}
model{
     vector[1260] mu;
     vector[1260] l0;
     vector[1260] linf;
     vector[1260] alpha;
    sigma_obs ~ exponential( 1 );
    sigma_c0 ~ exponential( 1 );
    sigma_b0 ~ exponential( 1 );
    sigma_a0 ~ exponential( 1 );
    mu_c0 ~ normal( 0.03 , 0.006 );
    mu_b0 ~ normal( 1 , 0.03 );
    mu_a0 ~ normal( 0.15 , 0.03 );
    c1 ~ normal( 0 , 0.05 );
    b1 ~ normal( 0 , 0.05 );
    a1 ~ normal( 0 , 0.05 );
    c0 ~ normal( mu_c0 , sigma_c0 );
    b0 ~ normal( mu_b0 , sigma_b0 );
    a0 ~ normal( mu_a0 , sigma_a0 );
    
    for ( i in 1:1260 ) {
        alpha[i] = c0[ID[i]] + c1 * Dose[i];
    }
    for ( i in 1:1260 ) {
        linf[i] = b0[ID[i]] + b1 * Dose[i];
    }
    for ( i in 1:1260 ) {
        l0[i] = a0[ID[i]] + a1 * Dose[i];
    }
    for ( i in 1:1260 ) {
        mu[i] = linf[i] * exp(log(l0[i]/linf[i]) * exp(-alpha[i] * t[i]));
    }
    l ~ normal( mu , sigma_obs );
}
