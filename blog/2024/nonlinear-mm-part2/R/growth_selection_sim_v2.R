library(tidyverse); library(patchwork); library(tidybayes)
library(brms); library(marginaleffects); library(viridis)
library(easystats); library(kableExtra)

# Setting theme
theme_set(theme_bw(16))

# Gompertz function ----
Gomp.fun = function(t, Linf, L0, alpha){
  logLhat = exp(-alpha * t) * (log(L0) - log(Linf)) + log(Linf) # predicted growth
  return(logLhat)}

Gomp.fun.2 = function(t, Linf, L0, alpha){
  Lhat = exp(exp(-alpha * t) * (log(L0) - log(Linf)) + log(Linf)) # predicted growth
  return(Lhat)}

# Data simulation ----
n_id = 10 # 10 individuals
t = seq(0,126, by = 24) # measure every h for 126 h 
Linf = 1370
L0 = 182
alpha = 0.028
sigma = 100 # random noise

set.seed(42)
df = crossing(t, 1:n_id) %>% 
  mutate(logLhat = Gomp.fun(t, Linf = Linf, L0 = L0, alpha = alpha)) %>% 
  mutate(Lhat = exp(logLhat)) %>% 
  mutate(L = rnorm(n(), Lhat, sigma))

df %>% 
  ggplot(aes(y = L, x = t)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = Gomp.fun.2, 
                args = list(Linf = Linf, L0 = L0, alpha = alpha),
                color = "red", size = 1) +
  ylim(0, 1600) +
  labs(x = "Hours since hatching", y = "Body-length (mum)")

# Log-scale model no vi ----
bf.l = bf(log(L) ~ 
            log(Linf*exp(log(L0/Linf) * exp(-alpha * t))), # Gompertz population curve
          L0 + Linf + alpha ~ 1,
          nl = T)
get_prior(bf.l, df)

priors = 
  # Intercept priors
  prior(normal(200, 40), nlpar = L0, class = b, lb = 0) +
  prior(normal(1500, 300), nlpar = Linf, class = b, lb = 0) +
  prior(normal(.03,.006), nlpar = alpha, class = b, lb = 0) +
  # Residual prior
  prior(exponential(1), class = sigma)

gomp.to.dose.prior = brm(bf.l, 
                         # family = gaussian(link = "exponential"),
                         df,
                         backend = "cmdstan",
                         prior = priors, 
                         # init = 0,
                         iter = 1000,
                         sample_prior = "only",
                         # file = "stan/gomp.to",
                         seed = 42, 
                         cores = 4,
                         threads = threading(3))
model_parameters(gomp.to.dose.prior, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.to.dose.prior))
plot(conditional_effects(gomp.to.dose.prior, 
                         ndraws = 100, spaghetti = T))

gomp.to.dose = brm(bf.l,
                   df,
                   backend = "cmdstan",
                   prior = priors, 
                   # init = 0,
                   warmup = 2000,
                   iter = 3000,
                   control = list(adapt_delta = .95,
                                  max_treedepth = 15),
                   sample_prior = "yes",
                   # file = "stan/gomp.to",
                   seed = 42, 
                   cores = 4,
                   threads = threading(3))
model_parameters(gomp.to.dose, effects = "all") %>%
  kable(digits = 2)
pp_check(gomp.to.dose, ndraws = 100)
plot(conditional_effects(gomp.to.dose), points = T)

# Log-scale model with vi ----

n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
sd = 10 # random noise
L0_mu = 182 # initial length (micrometers)
Linf_mu = 1370 # maximal length (micrometers)
alpha_mu = 0.028 # Growth rate (hour^-1)

rho = 0 # Suppose all parameters are independent

mu     = c(L0_mu, Linf_mu, alpha_mu)
sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, rho, rho,
                   rho, 1, rho,
                   rho, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = 1:n_id)

# Simulate individual growth
df = ID %>%
  tidyr::expand(nesting(ID, L0_i, Linf_i, alpha_i), 
                t = times) %>%
  mutate(Lhat = Linf_i * exp(log(L0_i/Linf_i)*exp(-alpha_i*t))) %>%
  mutate(L = rnorm(n(), Lhat, sd))

fig.title = df %>% 
  ggplot(aes(y = L, x = t, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .2, size = .5) +
  geom_function(fun = Gomp.fun.2, 
                args = list(Linf = Linf_mu, L0 = L0_mu, alpha = alpha_mu),
                color = "red", size = 1) +
  ylim(0, 1600) +
  labs(x = "Hours since hatching", 
       y = expression(paste("Body-length (", mu, "m)")))
fig.title

bf.l = bf(log(L) ~ 
            log(Linf*exp(log(L0/Linf) * exp(-alpha * t))), # Gompertz population curve
          L0 + Linf + alpha ~ 1 + (1|c|ID),
          nl = T)
get_prior(bf.l, df)

priors = 
  # Intercept priors
  prior(normal(200, 40), nlpar = L0, class = b, lb = 0) +
  prior(normal(1500, 300), nlpar = Linf, class = b, lb = 0) +
  prior(normal(.03,.006), nlpar = alpha, class = b, lb = 0) + 
  # Random effects priors (informative priors with 20 % CV)
  prior(exponential(1), nlpar = L0, class = sd, group = ID) +
  prior(exponential(1), nlpar = Linf, class = sd, group = ID) +
  prior(exponential(1), nlpar = alpha, class = sd, group = ID) +
  # Residual prior
  prior(exponential(1), class = sigma) +
  # Correlation prior
  prior(lkj(2), class = cor)


gomp.to.dose.prior = brm(bf.l, 
                         # family = gaussian(link = "exponential"),
                         df,
                         backend = "cmdstan",
                         prior = priors, 
                         # init = 0,
                         iter = 1000,
                         sample_prior = "only",
                         # file = "stan/gomp.to",
                         seed = 42, 
                         cores = 4,
                         threads = threading(3))
model_parameters(gomp.to.dose.prior, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.to.dose.prior, re_formula = NULL))
plot(conditional_effects(gomp.to.dose.prior, 
                         ndraws = 100, spaghetti = T))

gomp.to.dose = brm(bf.l,
                   df,
                   backend = "cmdstan",
                   prior = priors, 
                   init = 0,
                   warmup = 2000,
                   iter = 3000,
                   control = list(adapt_delta = .95,
                                  max_treedepth = 15),
                   sample_prior = "yes",
                   # file = "stan/gomp.to",
                   seed = 42, 
                   cores = 4,
                   threads = threading(3))
model_parameters(gomp.to.dose, effects = "all") %>%
  kable(digits = 2)
pp_check(gomp.to.dose, ndraws = 100)
plot(conditional_effects(gomp.to.dose, re_formula = NULL), points = T)

# Dose-response simulation ----
n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
sd = 10 # random noise
L0_mu = 182 # initial length (micrometers)
Linf_mu = 1370 # maximal length (micrometers)
alpha_mu = 0.028 # Growth rate (hour^-1)
Dose = c(0, .1, .3, .5, .9, 1.1, 1.2)
# a0_m = 1 # Mean trait in absence of contaminants
a1_L0 = -.25 * L0_mu # Change in mean trait per unit of contaminant concentration
a1_Linf = -.25 * Linf_mu # Change in mean trait per unit of contaminant concentration
a1_alpha = -.25 * alpha_mu # Change in mean trait per unit of contaminant concentration

rho = -.7 # Assume strong negative correlation between Linf and alpha

## Simulate Dose = 0 ----
mu = c(L0_mu, Linf_mu, alpha_mu)
sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 5 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID_d0 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = 1:n_id)

df = ID_d0 %>%
  tidyr::expand(nesting(ID, L0_i, Linf_i, alpha_i), 
                t = times) %>%
  mutate(Lhat = Linf_i * exp(log(L0_i/Linf_i)*exp(-alpha_i*t))) %>%
  mutate(L = rnorm(n(), Lhat, sd))

df %>% 
  ggplot(aes(y = L, x = t, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .2, linewidth = .5) +
  geom_function(fun = Gomp.fun.2, 
                args = list(Linf = Linf_mu,# + a1_L0 *Dose, 
                            L0 =  L0_mu, # + a1_Linf *Dose, 
                            alpha = alpha_mu), # + a1_alpha *Dose),
                color = "red", size = 1) +
  # ylim(0, 1600) +
  labs(x = "Hours since hatching", y = "Body-length (mum)")

## Simulate Dose = .1 ----
D = Dose[2]

mu = c(L0_mu + a1_L0 * D, 
       Linf_mu + a1_Linf * D, 
       alpha_mu + a1_alpha * D)

sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID_d1 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id + 1) : (n_id * 2))

## Simulate Dose = .3 ----
D = Dose[3]

mu = c(L0_mu + a1_L0 * D, 
       Linf_mu + a1_Linf * D, 
       alpha_mu + a1_alpha * D)

sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID_d2 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 2 + 1) : (n_id * 3))

## Simulate Dose = .5 ----
D = Dose[4]

mu = c(L0_mu + a1_L0 * D, 
       Linf_mu + a1_Linf * D, 
       alpha_mu + a1_alpha * D)

sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID_d3 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 3 + 1) : (n_id * 4))

## Simulate Dose = .9 ----
D = Dose[5]

mu = c(L0_mu + a1_L0 * D, 
       Linf_mu + a1_Linf * D, 
       alpha_mu + a1_alpha * D)

sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID_d4 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 4 + 1) : (n_id * 5))

## Simulate Dose = 1.1 ----
D = Dose[6]

mu = c(L0_mu + a1_L0 * D, 
       Linf_mu + a1_Linf * D, 
       alpha_mu + a1_alpha * D)

sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID_d5 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 5 + 1) : (n_id * 6))

## Simulate Dose = 1.2 ----
D = Dose[7]

mu = c(L0_mu + a1_L0 * D, 
       Linf_mu + a1_Linf * D, 
       alpha_mu + a1_alpha * D)

sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID_d6 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 6 + 1) : (n_id * 7))

## Combine ID data ----
ID = rbind(ID_d0, ID_d1, ID_d2, ID_d3, ID_d4, ID_d5, ID_d6)
ID$Dose = c(rep(Dose[1], n_id),
            rep(Dose[2], n_id),
            rep(Dose[3], n_id),
            rep(Dose[4], n_id),
            rep(Dose[5], n_id),
            rep(Dose[6], n_id),
            rep(Dose[7], n_id))
ID %>% 
  select(-ID) %>% 
  GGally::ggpairs() +
  theme_bw() 

## Simulate individual growth ----
df = ID %>%
  tidyr::expand(nesting(ID, L0_i, Linf_i, alpha_i, Dose), 
                t = times) %>%
  mutate(Lhat = Linf_i * exp(log(L0_i/Linf_i)*exp(-alpha_i*t))) %>%
  mutate(L = rnorm(n(), Lhat, sd))

df %>% 
  ggplot(aes(y = L, x = t, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .2, linewidth = .5) +
  geom_function(fun = Gomp.fun.2, 
                args = list(Linf = Linf_mu,# + a1_L0 *Dose, 
                            L0 =  L0_mu, # + a1_Linf *Dose, 
                            alpha = alpha_mu), # + a1_alpha *Dose),
                color = "red", size = 1) +
  facet_wrap(~Dose) +
  # ylim(0, 1600) +
  labs(x = "Hours since hatching", y = "Body-length (mum)")


# Log-scale model no vi ----
bf.l = bf(log(L) ~ 
            log(Linf*exp(log(L0/Linf) * exp(-alpha * t))), # Gompertz population curve
          nlf(L0 ~ a0 + doseL0 ),
          nlf(Linf ~ b0 + doseLinf),
          nlf(alpha ~ c0 + dosealpha),
          a0 + b0 + c0 ~ 0 + Intercept,
          lf(doseL0 + doseLinf + dosealpha ~ 0 + Dose, cmc = F),
          nl = T)
default_prior(bf.l, df)

priors = 
  # Intercept priors
  prior(normal(200, 40), nlpar = a0, class = b, coef = Intercept) + # , lb = 0
  prior(normal(1500, 300), nlpar = b0, class = b, coef = Intercept) +
  prior(normal(.03,.006), nlpar = c0, class = b, coef = Intercept) +
  # Regression priors
  prior(normal(0, .05), nlpar = doseL0, class = b, coef = Dose) +
  prior(normal(0, .05), nlpar = doseLinf, class = b, coef = Dose) +
  prior(normal(0, .05), nlpar = dosealpha, class = b, coef = Dose) +
  # Residual prior
  prior(exponential(1), class = sigma)

## Setting initial values ----
set_inits <- function(seed = 1) {
  set.seed(seed)
  list(
    b_a0 = rnorm(n = 1, 200, 40),
    b_b0 = rnorm(n = 1, 1500, 300),
    b_c0 = rnorm(n = 1, 0.03, 0.006),
    b_doseL0 = rnorm(n = 1, 20, 0.05),
    b_doseLinf = rnorm(n = 1, 0, 0.05),
    b_dosealpha = rnorm(n = 1, 0, 0.05),
    sigma = rexp(1, 1)
  )
}

list_of_inits <- list(
  set_inits(seed = 1),
  set_inits(seed = 2),
  set_inits(seed = 3),
  set_inits(seed = 4))

# what have we done?
str(list_of_inits)

## Fitting the model ----
gomp.to.dose.prior = brm(bf.l, 
                         # family = gaussian(link = "exponential"),
                         df,
                         backend = "cmdstan",
                         prior = priors, 
                         init = list_of_inits,
                         iter = 1000,
                         sample_prior = "only",
                         # file = "stan/gomp.to",
                         seed = 42, 
                         cores = 4,
                         threads = threading(3))
model_parameters(gomp.to.dose.prior, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.to.dose.prior))
plot(conditional_effects(gomp.to.dose.prior, 
                         ndraws = 100, spaghetti = T))

gomp.to.dose = brm(bf.l,
                   df,
                   backend = "cmdstan",
                   prior = priors, 
                   init = list_of_inits,
                   warmup = 2000,
                   iter = 3000,
                   control = list(adapt_delta = .95,
                                  max_treedepth = 15),
                   sample_prior = "yes",
                   # file = "stan/gomp.to",
                   seed = 42, 
                   cores = 4,
                   threads = threading(3))
model_parameters(gomp.to.dose, effects = "all") %>%
  kable(digits = 2)
pp_check(gomp.to.dose, ndraws = 100)
plot(conditional_effects(gomp.to.dose), points = T)




# Log-scale model no vi ----
bf.l = bf(log(L) ~ 
            log(Linf*exp(log(L0/Linf) * exp(-alpha * t))), # Gompertz population curve
          nlf(L0 ~ a0 + doseL0 ),
          nlf(Linf ~ b0 + doseLinf),
          nlf(alpha ~ c0 + dosealpha),
          a0 + b0 + c0 ~ 1 + (1 |c|ID),
          lf(doseL0 + doseLinf + dosealpha ~ 0 + Dose, cmc = F),
          nl = T)
default_prior(bf.l, df)

priors = 
  # Intercept priors
  prior(normal(200, 40), nlpar = a0, class = b, coef = Intercept) + # , lb = 0
  prior(normal(1500, 300), nlpar = b0, class = b, coef = Intercept) +
  prior(normal(.03,.006), nlpar = c0, class = b, coef = Intercept) +
  # Regression priors
  prior(normal(0, .05), nlpar = doseL0, class = b, coef = Dose) +
  prior(normal(0, .05), nlpar = doseLinf, class = b, coef = Dose) +
  prior(normal(0, .05), nlpar = dosealpha, class = b, coef = Dose) +
  # Random effects priors (informative priors with 20 % CV)
  prior(exponential(1), nlpar = a0, class = sd, group = ID) +
  prior(exponential(1), nlpar = b0, class = sd, group = ID) +
  prior(exponential(1), nlpar = c0, class = sd, group = ID) +
  # Residual prior
  prior(exponential(1), class = sigma) +
  # Correlation prior
  prior(lkj(2), class = cor)

## Setting initial values ----
set_inits <- function(seed = 1) {
  set.seed(seed)
  list(
    b_a0 = rnorm(n = 1, 200, 40),
    b_b0 = rnorm(n = 1, 1500, 300),
    b_c0 = rnorm(n = 1, 0.03, 0.006),
    b_doseL0 = rnorm(n = 1, 20, 0.05),
    b_doseLinf = rnorm(n = 1, 0, 0.05),
    b_dosealpha = rnorm(n = 1, 0, 0.05)
    # sd_a0 = rexp(1, 1),
    # sd_b0 = rexp(1, 1),
    # sd_c0 = rexp(1, 1),
    # sigma = rexp(1, 1)
  )
}

list_of_inits <- list(
  set_inits(seed = 1),
  set_inits(seed = 2),
  set_inits(seed = 3),
  set_inits(seed = 4))

# what have we done?
str(list_of_inits)

## Fitting the model ----
gomp.to.dose.prior = brm(bf.l, 
                         # family = gaussian(link = "exponential"),
                         df,
                         backend = "cmdstan",
                         prior = priors, 
                         init = list_of_inits,
                         iter = 1000,
                         sample_prior = "only",
                         # file = "stan/gomp.to",
                         seed = 42, 
                         cores = 4,
                         threads = threading(3))
gomp.to.dose.prior
plot(conditional_effects(gomp.to.dose.prior))
plot(conditional_effects(gomp.to.dose.prior, 
                         ndraws = 100, spaghetti = T))

gomp.to.dose = brm(bf.l,
                   df,
                   backend = "cmdstan",
                   prior = priors, 
                   init = list_of_inits,
                   warmup = 2000,
                   iter = 3000,
                   control = list(adapt_delta = .95,
                                  max_treedepth = 15),
                   sample_prior = "yes",
                   # file = "stan/gomp.to",
                   seed = 42, 
                   cores = 4,
                   threads = threading(3))
model_parameters(gomp.to.dose, effects = "all") %>%
  kable(digits = 2)
pp_check(gomp.to.dose, ndraws = 100)
plot(conditional_effects(gomp.to.dose), points = T)


