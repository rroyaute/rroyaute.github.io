// Stan program for estimating selection gradients on two independent traits (x, y)
// The input data are vectors of length 'N'.
data {
  int<lower=0> N; // number of observations
  vector[N] x; // trait x
  vector[N] y; // trait y
  vector[N] w; // fitness w
}

parameters {
  real mu_x; // mean of trait x
  real<lower=0> sigma_x; // sd of trait x
  real mu_y;  // mean of trait y
  real<lower=0> sigma_y;  // sd of trait y
  real mu_w;  // mean fitness
  real<lower=0> sigma_w;  // sd fitness
  real b1; // linear selection on x
  real b2; // linear selection on y
  real g11; // quad. selection on x
  real g22; // quad. selection on x
  real g12; // correlational selection on x & y
}

model {
  vector[N] w_pred; // fitness stored as vector
  x ~ normal(mu_x, sigma_x); // x follows a normal distribution of mean mu_x and sd sd_x
  y ~ normal(mu_y, sigma_y);
  //w ~ normal(mu_w, sigma_w);
  
  // Lande & Arnold equation linking fitness to x & y
  
  w_pred = mu_w + 
    b1 * x + b2 * y +
    g11 * x^2 + g22 * y^2 +
    g12 * (x * y);
  
  w ~ normal(w_pred, sigma_w);  
}
