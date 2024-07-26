//
//
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] x; // trait x
  vector[N] y; // trait y
  vector[N] w; // fitness w
  vector[N] dose; // contaminant doses
}

parameters {
  real mu_x;
  real<lower=0> sigma_x;
  real mu_y;
  real<lower=0> sigma_y;
  real mu_w;
  real<lower=0> sigma_w;
  //real b0;
  real b1;
  real b2;
  real g11;
  real g22;
  real g12;
  real d0_w;
  real d1_w;
  real d0_x;
  real d1_x;
  real d0_y;
  real d1_y;
  real d0_x2;
  real d1_x2;
  real d0_y2;
  real d1_y2;
  real d0_xy;
  real d1_xy;
}

model {
  x ~ normal(mu_x, sigma_x);
  y ~ normal(mu_y, sigma_y);
  w ~ normal(mu_w, sigma_w);
  
  // Lande & Arnold equation linking fitness to x & y
  for (i in 1:N) {
  w[i] = mu_w + 
    b1 * x[i] + b2 * y[i] +
    g11 * x[i]^2 + g22 * y[i]^2 +
    g12 * (x[i]*y[i]);
    }
    // Dose-response on mean fitness
    //b0 = d0_w + d1_w * dose;
    
    // Dose-response on selection gradients
    // b1 = d0_x + d1_x * dose;
    // b2 = d0_y + d1_y * dose;
    // g11 = d0_x2 + d1_x2 * dose;
    // g22 = d0_y2 + d1_y2 * dose;
    // g12 = d0_xy + d1_xy * dose;


  
}

