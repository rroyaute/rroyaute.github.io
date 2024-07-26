data{
     vector[264] dose_cat;
     vector[264] dose;
     vector[264] outcome;
     vector[264] time;
     vector[264] dose_adj;
    array[264] int id;
}
parameters{
     vector[24] a0;
     vector[24] b0;
     real mu_a;
     real mu_b;
     real<lower=0> sigma_a;
     real<lower=0> sigma_b;
     real a1;
     real b1;
     real sigma;
}
model{
     vector[264] mu;
     vector[264] alpha;
     vector[264] beta;
    sigma ~ cauchy( 0 , 1 );
    b1 ~ normal( -0.3 , 1 );
    a1 ~ normal( 0.3 , 1 );
    sigma_b ~ cauchy( 0 , 1 );
    sigma_a ~ cauchy( 0 , 1 );
    mu_b ~ normal( 0.5 , 1 );
    mu_a ~ normal( 2 , 1 );
    b0 ~ normal( mu_b , sigma_b );
    a0 ~ normal( mu_a , sigma_a );
    for ( i in 1:264 ) {
        beta[i] = b0[id[i]] + b1 * dose_adj[i];
    }
    for ( i in 1:264 ) {
        alpha[i] = a0[id[i]] + a1 * dose_adj[i];
    }
    for ( i in 1:264 ) {
        mu[i] = exp(alpha[i]) * log(time[i]) - exp(beta[i]) * time[i];
    }
    outcome ~ normal( mu , sigma );
}
