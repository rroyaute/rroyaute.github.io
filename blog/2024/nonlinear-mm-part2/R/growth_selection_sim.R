library(tidyverse); library(patchwork); library(tidybayes)
library(brms); library(marginaleffects); library(viridis)
library(easystats); library(kableExtra)

# Setting theme
theme_set(theme_bw(16))

# Gompertz function
Gomp.fun = function(t, Linf, L0, alpha){
  Lhat = exp(Linf) * exp(log(exp(L0)/exp(Linf))*exp(-exp(alpha)*t)) # predicted growth
  return(Lhat)}

# Simulation parameters
n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
sd = 10 # random noise
L0_mu = log(182) # initial length (micrometers)
Linf_mu = log(1370) # maximal length (micrometers)
alpha_mu = log(0.028) # Growth rate (hour^-1)
Dose = c(0, .1, .3, .5, .9, 1.1, 1.2)
# a0_m = 1 # Mean trait in absence of contaminants
a1_L0 = -.25 * L0_mu # Change in mean trait per unit of contaminant concentration
a1_Linf = -.25 * Linf_mu # Change in mean trait per unit of contaminant concentration
a1_alpha = -.25 * alpha_mu # Change in mean trait per unit of contaminant concentration

a0_s = 1 # Selection gradients in absence of contaminants
a1_s = .5 # Change in gradient per unit of contaminant concentration

rho = -.7 # Assume strong negative correlation between Linf and alpha

## Simulate Dose = 0 ----
mu = c(L0_mu, Linf_mu, alpha_mu)
sigmas = c(L0_mu*.01, Linf_mu*.01, alpha_mu*.01) # 1 % CV around the mean on log-scale
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
  mutate(Lhat = exp(Linf_i) * exp(log(exp(L0_i)/exp(Linf_i))*exp(-exp(alpha_i)*t))) %>%
  mutate(L = rnorm(n(), Lhat, sd))

df %>% 
  ggplot(aes(y = L, x = t, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .2, linewidth = .5) +
  geom_function(fun = Gomp.fun, 
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

sigmas = c(L0_mu*.01, Linf_mu*.01, alpha_mu*.01) # 10 % CV around the mean
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

sigmas = c(L0_mu*.01, Linf_mu*.01, alpha_mu*.01) # 10 % CV around the mean
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

sigmas = c(L0_mu*.01, Linf_mu*.01, alpha_mu*.01) # 10 % CV around the mean
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

sigmas = c(L0_mu*.01, Linf_mu*.01, alpha_mu*.01) # 10 % CV around the mean
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

sigmas = c(L0_mu*.01, Linf_mu*.01, alpha_mu*.01) # 10 % CV around the mean
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

sigmas = c(L0_mu*.01, Linf_mu*.01, alpha_mu*.01) # 10 % CV around the mean
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
  mutate(Lhat = exp(Linf_i) * exp(log(exp(L0_i)/exp(Linf_i))*exp(-exp(alpha_i)*t))) %>%
  mutate(L = rnorm(n(), Lhat, sd))

df %>% 
  ggplot(aes(y = L, x = t, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .2, linewidth = .5) +
  geom_function(fun = Gomp.fun, 
                args = list(Linf = Linf_mu,# + a1_L0 *Dose, 
                            L0 =  L0_mu, # + a1_Linf *Dose, 
                            alpha = alpha_mu), # + a1_alpha *Dose),
                color = "red", size = 1) +
  facet_wrap(~Dose) +
  # ylim(0, 1600) +
  labs(x = "Hours since hatching", y = "Body-length (mum)")

## Dose response growth modeling ----
df$l = with(df, L/exp(Linf_mu))
df$logl = log(df$l)

### No Vi model ----
bf.l = bf(exp(logl) ~ 
            exp(linf)*exp(log(exp(l0)/exp(linf)) * exp(-exp(alpha) * t)), # Gompertz population curve
          l0 + linf + alpha ~ 1 + Dose,
          # nlf(l0 ~ a0 + dosel0),
          # nlf(linf ~ b0 + doselinf),
          # nlf(alpha ~ c0 + dosealpha),
          # a0 + b0 + c0 ~ 1, # Intercepts
          # lf(dosel0 + doselinf + dosealpha ~ 0 + Dose, cmc = F), # Linear formulas for change in parameter with dose
          nl = T)
get_prior(bf.l, df)

priors = 
  # Intercept priors
  prior(normal(-1.89712, .5), class = b, coef = Intercept, nlpar = l0) +
  prior(normal(0, .5), class = b, coef = Intercept, nlpar = linf) +
  prior(normal(-3.575551, .5), class = b, coef = Intercept, nlpar = alpha) + 
  
  # Regression priors
  prior(normal(0, .05), class = b, coef = Dose, nlpar = l0) +
  prior(normal(0, .05), class = b, coef = Dose, nlpar = linf) +
  prior(normal(0, .05), class = b, coef = Dose, nlpar = alpha) +
  # prior(normal(0, .05), nlpar = doselinf) +
  # prior(normal(0, .05), nlpar = dosealpha) +
  
  # Residual prior
  prior(exponential(1), class = sigma)

gomp.to.dose.prior = brm(bf.l,
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
plot(conditional_effects(gomp.to.dose, 
                         ndraws = 100, spaghetti = T))
plot(conditional_effects(gomp.to.dose), points = T)

Doses <- data.frame(Dose = unique(df$Dose))
rownames(Doses) <- unique(df$Dose)
me_dose <- conditional_effects(
  gomp.to.dose, conditions = Doses,
  re_formula = NULL, method = "predict")
plot(me_dose, ncol = 5, points = TRUE)

### With Vi model ----
bf.l.vi = bf(l ~ linf*exp(log(l0/linf) * exp(-alpha * t)), # Gompertz population curve
             nlf(l0 ~ a0 + dosel0 ),
             nlf(linf ~ b0 + doselinf),
             nlf(alpha ~ c0 + dosealpha),
             a0 + b0 + c0 ~ 1 + (1 |c| ID), # Intercepts
             lf(dosel0 + doselinf + dosealpha ~ 0 + Dose, cmc = F), # Linear formulas for change in parameter with dose
             nl = T)
get_prior(bf.l.vi, df)

priors = 
  # Intercept priors
  prior(normal(.15, .03), nlpar = a0, class = b, lb = 0) +
  prior(normal(1, .03), nlpar = b0, class = b, lb = 0) +
  prior(normal(.03,.006), nlpar = c0, class = b, lb = 0) + 
  # Regression priors
  prior(normal(0, .05), nlpar = dosel0) +
  prior(normal(0, .05), nlpar = doselinf) +
  prior(normal(0, .05), nlpar = dosealpha) +
  # Random effects priors (informative priors with 20 % CV)
  prior(exponential(33), nlpar = a0, class = sd, group = ID) +
  prior(exponential(5), nlpar = b0, class = sd, group = ID) +
  prior(exponential(170), nlpar = c0, class = sd, group = ID) +
  
  # Residual prior
  prior(exponential(1), class = sigma) +
  # Correlation prior
  prior(lkj(2), class = cor)

gomp.to.dose.vi.prior = brm(bf.l.vi,
                         df,
                         backend = "cmdstan",
                         prior = priors, 
                         init = 0,
                         iter = 1000,
                         sample_prior = "only",
                         # file = "stan/gomp.to",
                         seed = 42, 
                         cores = 4,
                         threads = threading(3))
model_parameters(gomp.to.dose.vi.prior, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.to.dose.vi.prior, 
                         ndraws = 100, spaghetti = T))

gomp.to.dose.vi = brm(bf.l.vi,
                   df,
                   backend = "cmdstan",
                   prior = priors, 
                   init = 0,
                   # warmup = 4000,
                   # iter = 5000,
                   control = list(adapt_delta = .95,
                                  max_treedepth = 15),
                   sample_prior = "yes",
                   # file = "stan/gomp.to",
                   seed = 42, 
                   cores = 4,
                   threads = threading(3))
model_parameters(gomp.to.dose.vi, effects = "all") %>%
  kable(digits = 2)
pp_check(gomp.to.dose.vi, ndraws = 100)
plot(conditional_effects(gomp.to.dose.vi, 
                         ndraws = 100, spaghetti = T))


## Fitness modeling ----
bf.w = bf(w ~ b0 + b1 * l0 + b2 * linf + b3 * alpha + 
            g1 * l0^2 + g2 * linf^2 + g3 * alpha^2 + 
            g4 * l0:linf + g5 * l0:alpha + g6 * linf:alpha,
          b0 + b1 + b2 + b3 + g1 + g2 + g3 + g4 + g5 + g6 ~ a0_s + a1_s * Dose, 
          nl = T)