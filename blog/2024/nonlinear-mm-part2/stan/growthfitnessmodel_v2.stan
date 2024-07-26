data{
     vector[180] l;
    array[180] int t;
     vector[180] w_rel;
    array[180] int ID;
}
parameters{
     vector<lower=0>[30] alpha;
     vector<lower=0>[30] linf;
     vector<lower=0>[30] l0;
     real <lower=0> mu_l0;
     real <lower=0> mu_linf;
     real <lower=0> mu_alpha;
     vector<lower=0>[3] sigma_i;
     real<lower=0> sigma_l;
     corr_matrix[3] Rho;
     real b0;
     real b1;
     real b2;
     real b3;
     real g11;
     real g22;
     real g33;
     real g12;
     real g13;
     real g23;
     real<lower=0> sigma_w;
}
model{
     vector[180] mu_l;
     vector[180] mu_w;
    sigma_w ~ exponential( 1 );
    g23 ~ normal( 0 , 1 );
    g13 ~ normal( 0 , 1 );
    g12 ~ normal( 0 , 1 );
    g33 ~ normal( 0 , 1 );
    g22 ~ normal( 0 , 1 );
    g11 ~ normal( 0 , 1 );
    b3 ~ normal( 0 , 1 );
    b2 ~ normal( 0 , 1 );
    b1 ~ normal( 0 , 1 );
    b0 ~ normal( 1 , 0.5 );
    for ( i in 1:180 ) {
        mu_w[i] = b0 + b1 * l0[ID[i]] + b2 * linf[ID[i]] + b3 * alpha[ID[i]] +
        g11 * l0[ID[i]]^2 + g22 * linf[ID[i]]^2 + g33 * alpha[ID[i]]^2 + 
        g12 * (l0[ID[i]] * linf[ID[i]]) + 
        g13 * (l0[ID[i]] * alpha[ID[i]]) + 
        g23 * (linf[ID[i]] * alpha[ID[i]]);
    }
    w_rel ~ normal( mu_w , sigma_w );
    Rho ~ lkj_corr( 2 );
    sigma_l ~ exponential( 1 );
    sigma_i ~ exponential( 1 );
    mu_alpha ~ normal( log(0.03) , 0.1 );
    mu_linf ~ normal( log(1) , 0.001 );
    mu_l0 ~ normal( log(0.15) , 0.1 );
    {
    array[30] vector[3] YY;
    vector[3] MU;
    MU = [ mu_l0 , mu_linf , mu_alpha ]';
    for ( j in 1:30 ) YY[j] = [ l0[j] , linf[j] , alpha[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho , sigma_i) );
    }
    for ( i in 1:180 ) {
        mu_l[i] = linf[ID[i]] * 
        exp(log(l0[ID[i]] / linf[ID[i]]) * 
        exp(-alpha[ID[i]] * t[i]));
    }
    l ~ normal( mu_l , sigma_l );
}
