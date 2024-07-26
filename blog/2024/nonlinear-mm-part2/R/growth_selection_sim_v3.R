library(here); library(tidyverse); library(patchwork); library(tidybayes)
library(rethinking); library(marginaleffects); library(viridis); library(scales)
library(easystats); library(kableExtra); library(plot3D)

# Setting theme ----
theme_set(theme_bw(16))

# Function ----
Gomp.fun = function(t, Linf, L0, alpha){
  Lhat = Linf * exp(log(L0/Linf)*exp(-alpha*t)) # predicted growth
  return(Lhat)
}

# Store parameters for non-exposed populations ----
n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
sd = 10 # random noise
L0_mu = 182 # initial length (micrometers)
Linf_mu = 1370 # maximal length (micrometers)
alpha_mu = 0.028 # Growth rate (hour^-1)

rho = -.7 # Assume strong negative correlation between Linf and alpha

mu     = c(L0_mu, Linf_mu, alpha_mu)
sigmas = c(L0_mu*.1, Linf_mu*.1, alpha_mu*.1) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = 1:n_id)

# Simulate individual growth ----
df = ID %>%
  tidyr::expand(nesting(ID, L0_i, Linf_i, alpha_i), 
                t = times) %>%
  mutate(Lhat = Linf_i * exp(log(L0_i/Linf_i)*exp(-alpha_i*t))) %>%
  mutate(L = rnorm(n(), Lhat, sd))

fig.title = df %>% 
  ggplot(aes(y = L, x = t, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .2, size = .5) +
  geom_function(fun = Gomp.fun, 
                args = list(Linf = Linf_mu, L0 = L0_mu, alpha = alpha_mu),
                color = "red", size = 1) +
  ylim(0, 1600) +
  labs(x = "Hours since hatching", 
       y = expression(paste("Body-length (", mu, "m)")))
fig.title

# Ulam model ----
## LMM for test ----
lmm.1 = ulam(alist(
  L ~ normal(mu_L, sigma_obs),
  mu_L <- b0[ID] + b1[ID] * t,
  
  c(b0, b1)[ID] ~ multi_normal(c(mu_b0, mu_b1), 
                               Rho, 
                               sigma_i),
  Rho ~ lkjcorr(2),
  
  mu_b0 ~ normal(200, 40),
  mu_b1 ~ normal(10, 2),
  
  sigma_i ~ exponential(1),
  sigma_obs ~ exponential(1)), 
  
  data = df, 
  # start = inits, 
  # constraints = constraints,
  chains = 4, 
  cores = 4, 
  # threads = 3,
  seed = 42)
precis(lmm.1)
stancode(lmm.1)

## NLMM ----
# inits & constraints
set.seed(42)
inits = list(mu_L0 = rnorm(1, 200, 40), 
             mu_Linf = rnorm(1, 1500, 300),
             mu_alpha = rnorm(1, .03,.006))

constraints = list(L0="lower=0", 
                   Linf="lower=0",
                   alpha="lower=0",
                   mu_L0="lower=0", 
                   mu_Linf="lower=0",
                   mu_alpha="lower=0")

# Fit!

nlmm.1 = ulam(
  alist(
    L ~ normal(mu_L, sigma_obs),
    mu_L <- Linf[ID] * exp( log(L0[ID] / Linf[ID] )*exp( -alpha[ID] * t)),
    
    
    c(L0, Linf, alpha)[ID] ~ multi_normal(c(mu_L0, mu_Linf, mu_alpha), 
                                          Rho, 
                                          sigma_i),
    Rho ~ lkjcorr(2),
    
    mu_L0 ~ normal(200, 40),
    mu_Linf ~ normal(1500, 300),
    mu_alpha ~ normal(.03,.006),
    
    sigma_i ~ exponential(1),
    
    sigma_obs ~ exponential(1)), 
  data = df, 
  start = inits, 
  constraints = constraints,
  chains = 4, 
  cores = 4, 
  # threads = 3,
  seed = 42)
precis(nlmm.1)
stancode(nlmm.1)



# Ulam reparametrized model ----
df$l = with(df, L/Linf_mu)
set.seed(42)
inits = list(mu_l0 = rnorm(1, .15, .03), 
             mu_linf = rnorm(1, 1, .03),
             mu_alpha = rnorm(1, .03,.006))
# Positive values
inits = list(mu_l0 = rexp(1,1), 
             mu_linf = rexp(1,1),
             mu_alpha = rexp(1,1))


constraints = list(l0="lower=0", 
                   linf="lower=0",
                   alpha="lower=0",
                   mu_l0="lower=0", 
                   mu_linf="lower=0",
                   mu_alpha="lower=0")


nlmm.2 = ulam(
  alist(
    l ~ normal(mu_l, sigma_obs),
    mu_l <- linf[ID] * exp( log(l0[ID] / linf[ID] )*exp( -alpha[ID] * t)),
    
    
    c(l0, linf, alpha)[ID] ~ multi_normal(c(mu_l0, mu_linf, mu_alpha), 
                                          Rho, 
                                          sigma_i),
    Rho ~ lkjcorr(2),
    
    mu_l0 ~ normal(.15, .03),
    mu_linf ~ normal(1, .03),
    mu_alpha ~ normal(.03,.006),
    
    sigma_i ~ exponential(1),
    sigma_obs ~ exponential(1)), 
  data = df, 
  start = inits, 
  constraints = constraints,
  chains = 4, 
  cores = 4, 
  # threads = 3,
  seed = 42)
precis(nlmm.2, depth = 3)
stancode(nlmm.2)

 

# Ulam re-reparamtrized model ----
df$l = with(df, L/Linf_mu)
set.seed(42)
inits = list(mu_l0 = rnorm(1, log(.15), .1), 
             mu_linf = rnorm(1, log(1), .001),
             mu_alpha = rnorm(1, log(.03), .1))

# constraints = list(l0 = "lower=0", 
#                    linf = "lower=0",
#                    alpha = "lower=0")


nlmm.3 = ulam(
  alist(
    l ~ normal(mu_l, sigma_obs),
    mu_l <- exp(linf[ID]) * 
      exp( log( exp(l0[ID]) / exp(linf[ID]) ) * 
             exp( - exp(alpha[ID]) * t)),
    
    
    c(l0, linf, alpha)[ID] ~ multi_normal(c(mu_l0, mu_linf, mu_alpha), 
                                          Rho, 
                                          sigma_i),
    Rho ~ lkjcorr(2),
    
    mu_l0 ~ normal(log(.15), .1),
    mu_linf ~ normal(log(1), .001),
    mu_alpha ~ normal(log(.03), .1),
    
    sigma_i ~ exponential(1),
    sigma_obs ~ exponential(1)), 
  data = df, 
  start = inits, 
  # constraints = constraints, 
  control = list(adapt_delta = .999,
                 maxtreedepth = 15),
  chains = 4, 
  cores = 4, 
  # threads = 3,
  seed = 42)
precis(nlmm.3, depth = 3)
stancode(nlmm.3)


# Growth-fitness joint model ----
mu_w = 1000 # Mean fitness
sigma_w = 5 # measurement error
hist(rnorm(1e6, mu_w, sigma_w))
# mu_w = 1 # Mean relative fitness

b1 = -15 # Directional selection on l0
b2 = 15 # Directional selection on linf
b3 = 15 # Directional selection on alpha
g11 = 15 # quadratic selection on l0
g22 = -15 # quadratic selection on linf
g33 = -15 # quadratic selection on alpha
g12 = -15 # correlational selection on l0xlinf
g13 = -15 # correlational selection on l0xalpha
g23 = -15 # correlational selection on l0xalpha


# Plot x,y,w, z = 0
l0 = seq(.07, .2, length.out = 100)
linf = seq(.8, 1.2, length.out = 100)
alpha = seq(.005, .05, length.out = 100)
l0_sc = (l0 - mean(l0))/sd(l0)
linf_sc = (linf - mean(linf))/sd(linf)
alpha_sc = (alpha - mean(alpha))/sd(alpha)


fit_f = function(l0 = 0, linf = 0, alpha = 0){
  w = mu_w + 
    b1 * l0 + b2 * linf + b3 * alpha +
    g11 * l0^2 + g22 * linf^2 + g33 * alpha^2 +
    g12 * (l0 * linf) + 
    g13 * (l0 * alpha) + 
    g23 * (linf * alpha)
  return(w)
}
plot(l0_sc, fit_f(l0 = l0_sc), "l")
plot(linf_sc, fit_f(linf = linf_sc), "l")
plot(alpha_sc, fit_f(alpha = alpha_sc), "l")

plot(l0_sc*sd(l0) + mean(l0), fit_f(l0 = l0_sc, linf = 0, alpha = 0), "l")

z = outer(l0_sc, linf_sc, function(l0, linf, alpha) fit_f(l0_sc, linf_sc, alpha = 0))
persp(x, y, z,
      col = "lightblue", 
      theta = 300, 
      xlim = c(-2, 2), 
      ylim = c(-2, 2),
      zlim = c(900, 1010))

w = ID %>% 
  mutate(L0_i2 = L0_i^2, Linf_i2 = Linf_i^2, alpha_i2 = alpha_i^2, 
         L0Linf_i = L0_i * Linf_i, L0alpha_i = L0_i * alpha_i, 
         Linfalpha_i = Linf_i * alpha_i) %>% 
  mutate(L0_i_sd = as.numeric(scale(L0_i)),
         Linf_i_sd = as.numeric(scale(Linf_i)),
         alpha_i_sd = as.numeric(scale(alpha_i)),
         L0_i2_sd = as.numeric(scale(L0_i2)),
         Linf_i2_sd = as.numeric(scale(Linf_i2)),
         alpha_i2_sd = as.numeric(scale(alpha_i2)),
         L0Linf_i_sd = as.numeric(scale(L0Linf_i)),
         L0alpha_i_sd = as.numeric(scale(L0alpha_i)),
         Linfalpha_i_sd = as.numeric(scale(Linfalpha_i))) %>% 
  mutate(w_i = mu_w + 
           b1 * L0_i_sd + b2 * Linf_i_sd + b3 * alpha_i_sd   +
           g11 *L0_i2_sd + g22 * Linf_i2_sd + g33 * alpha_i2_sd +
           g12 * L0Linf_i_sd + 
           g13 * L0alpha_i_sd  + g23 * Linfalpha_i_sd) %>% 
  mutate(w_obs = rnorm(n(), w_i, sigma_w)) %>% 
  mutate(w_rel = w_i / mean(w_i))
hist(w$w_rel)
w %>% 
  ggplot(aes(x = L0_i_sd, y = w_rel)) +
  geom_point() +
  geom_smooth(method = "lm", 
              formula=y ~ poly(x, 2, raw=TRUE))

w %>% 
  ggplot(aes(x = Linf_i_sd, y = w_rel)) +
  geom_point() +
  geom_smooth(method = "lm", 
              formula=y ~ poly(x, 2, raw=TRUE))

w %>% 
  ggplot(aes(x = alpha_i_sd, y = w_rel)) +
  geom_point() +
  geom_smooth(method = "lm", 
              formula=y ~ poly(x, 2, raw=TRUE))

df = left_join(df, w)
df$l = with(df, L/Linf_mu)
dlist = list(l = df$l,
             w_rel = df$w_rel,
             t = df$t,
             ID = df$ID)
set.seed(42)
inits = list(mu_l = rexp(1, 1),
             mu_l0 = rexp(1, 1),
             mu_linf = rexp(1, 1),
             mu_alpha = rexp(1, 1),
             mu_w = rexp(1, 1),
             sigma_l = rexp(1, 1),
             sigma_l0 = rexp(1, 1),
             sigma_linf = rexp(1, 1),
             sigma_alpha = rexp(1, 1),
             sigma_w = rexp(1, 1))
constraints = list(mu_l = "lower=0",
                   mu_w = "lower=0")

form = alist(
  # Growth model
  l ~ normal(mu_l, sigma_l),
  mu_l <- linf[ID] * exp( log( l0[ID] / linf[ID] ) *
                           exp( - alpha[ID] * t ) ),
  c(l0, linf, alpha)[ID] ~ multi_normal(c(mu_l0, mu_linf, mu_alpha), 
                                        Rho, 
                                        sigma_i),
  # Priors
  # means
  mu_l0 ~ normal(.15, .1),
  mu_linf ~ normal(0, .001),
  mu_alpha ~ normal(.03, .1),
  # sd
  sigma_i ~ exponential(1),
  sigma_l ~ exponential(1),
  # correlations
  Rho ~ lkjcorr(2),
  
  # Transform l0, linf & alpha to standard deviation units
  # transdata> l0_sc <- (l0 - mu_l0)/sigma_i[1],
  # transdata> linf_sc <- (linf - mu_linf)/sigma_i[2],
  # transdata> alpha_sc <- (alpha - mu_alpha)/sigma_i[3],
  # Transform l0^2, linf^2 & alpha^2 to standard deviation units
  # transdata> l02_sc <- l0_sc^2,
  # transdata> linf2_sc <- linf_sc^2,
  # transdata> alpha2_sc <- alpha_sc^2,

  # Fitness error-in-variable
  w_rel ~ normal(mu_w, sigma_w),
  mu_w <- b0 +
    b1 * l0[ID] + b2 * linf[ID] + b3 * alpha[ID] +
    g11 * l0[ID]^2 + g22 * linf[ID]^2 + g33 * alpha[ID]^2 +
    g12 * (l0[ID] * linf[ID]) +
    g13 * (l0[ID] * alpha[ID]) +
    g23 * (linf[ID] * alpha[ID]),
  # mu_w <- b0 +
  #   b1 * l0_sc[ID] + b2 * linf_sc[ID] + b3 * alpha_sc[ID] +
  #   g11 * l02_sc[ID] + g22 * linf2_sc[ID] + g33 * alpha2_sc[ID] +
  #   g12 * (l0_sc[ID] * linf_sc[ID]) +
  #   g13 * (l0_sc[ID] * alpha_sc[ID]) +
  #   g23 * (linf_sc[ID] * alpha_sc[ID]),
  
  # Priors for fitness
  b0 ~ normal(1, .5), # mean relative fitness
  # selection gradients
  b1 ~ normal(0, 1),
  b2 ~ normal(0, 1),
  b3 ~ normal(0, 1),
  g11 ~ normal(0, 1),
  g22 ~ normal(0, 1),
  g33 ~ normal(0, 1),
  g12 ~ normal(0, 1),
  g13 ~ normal(0, 1),
  g23 ~ normal(0, 1),
  # sd
  sigma_w ~ exponential(1))

nlmm.5 = ulam(
  form, 
  data = dlist,
  start = inits,
  constraints = constraints,
  control = list(adapt_delta = .99,
                 max_treedepth = 15),
  chains = 4,
  cores = 4, 
  # threads = 3,
  seed = 42, sample = T)
precis(nlmm.5)
stancode(nlmm.5)
plot(nlmm.5)
plot(nlmm.5, 
     pars = c("b0","b1", "b2", "b3", 
              "g11", "g22", "g33", 
              "g12", "g13", "g23"))
plot(nlmm.5, 
     pars = c("sigma_i"))
# check with cmdstan
mod <- cmdstan_model("stan/growthfitnessmodel.stan")
mod$print()
fit_mcmc <- mod$sample(
  data = dlist,
  seed = 123,
  chains = 2,
  parallel_chains = 2)

mod.2 <- cmdstan_model("stan/growthfitnessmodel_v2.stan")
mod.2$print()
fit_mcmc <- mod.2$sample(
  data = dlist,
  seed = 123,
  chains = 2,
  parallel_chains = 2, 
  adapt_delta = .9, 
  max_treedepth = 15)
fit_mcmc$diagnostic_summary()

# Dose-response model ----
## Dose-response simulation ----
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

### Simulate Dose = 0 ----
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
### Simulate Dose = .1 ----
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

set.seed(43)
ID_d1 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id + 1) : (n_id * 2))

### Simulate Dose = .3 ----
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

set.seed(44)
ID_d2 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 2 + 1) : (n_id * 3))

### Simulate Dose = .5 ----
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

set.seed(45)
ID_d3 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 3 + 1) : (n_id * 4))

### Simulate Dose = .9 ----
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

set.seed(46)
ID_d4 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 4 + 1) : (n_id * 5))

### Simulate Dose = 1.1 ----
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

set.seed(47)
ID_d5 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 5 + 1) : (n_id * 6))

### Simulate Dose = 1.2 ----
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

set.seed(48)
ID_d6 = MASS::mvrnorm(n_id, mu, sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = (n_id * 6 + 1) : (n_id * 7))

### Combine ID data ----
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

### Simulate individual growth ----
df = ID %>%
  tidyr::expand(nesting(ID, L0_i, Linf_i, alpha_i, Dose), 
                t = times) %>%
  mutate(Lhat = Linf_i * exp(log(L0_i/Linf_i)*exp(-alpha_i*t))) %>%
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
  labs(x = "Hours since hatching", 
       expression(paste("Body-length (", mu, "m)")))



# Dose-response Growth model ----
df$l = with(df, L/Linf_mu)
dlist = list(l = df$l,
             Dose = df$Dose,
             t = df$t,
             ID = df$ID)
set.seed(42)
inits = list(mu_a0 = rexp(1, 1),
             mu_b0 = rexp(1, 1),
             mu_c0 = rexp(1, 1))
form = alist(
  l ~ normal(mu, sigma_obs),
  mu <- exp(linf) * exp( log( exp(l0) / exp(linf) ) *
                           exp( - exp(alpha) * t ) ),
  l0 <- a0[ID] + a1 * Dose,
  linf <- b0[ID] + b1 * Dose,
  alpha <- c0[ID] + c1 * Dose,
  a0[ID] ~ normal(mu_a0, sigma_a0),
  b0[ID] ~ normal(mu_b0, sigma_b0),
  c0[ID] ~ normal(mu_c0, sigma_c0),
  a1 ~ normal(0, .05),
  b1 ~ normal(0, .05),
  c1 ~ normal(0, .05),
  mu_a0 ~ normal(log(.15), .1),
  mu_b0 ~ normal(log(1), .001),
  mu_c0 ~ normal(log(.03), .1),
  sigma_a0 ~ exponential(1),
  sigma_b0 ~ exponential(1),
  sigma_c0 ~ exponential(1),
  sigma_obs ~ exponential(1))

nlmm.4 = ulam(
  form, 
  data = dlist,
  start = inits,
  # constraints = constraints,
  control = list(adapt_delta = .999, 
                 max_treedepth = 15),
  chains = 4,
  cores = 4,
  seed = 42)
precis(nlmm.4)
stancode(nlmm.4)

# Fit directly from stan because fuck it I'm done ----
# make Stan model. 
nlmm.3.stan <- cmdstanr::cmdstan_model(here("blog/2024/nonlinear-mm-part2/stan/growthmodel_dr.stan"), 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)
nlmm.3.stan.fit <- nlmm.3.stan$sample(
  data = dlist,
  chains = 1,
  init = list(inits),
  max_treedepth = 15, 
  adapt_delta = .999,
  seed = 42)
nlmm.3.stan.fit$diagnostic_summary()
pos_samples = nlmm.3.stan.fit$draws(inc_warmup = FALSE, format = "draws_df")
ypred_df = pos_samples %>% select(starts_with("ypred"))
ppcheck = bayesplot::ppc_dens_overlay(df$l, as.matrix(ypred_df))
plot(m1_p2)


# Fitness error-in-variable model ----
## Simulate fitness values ----
mu_w = 1000 # Mean fitness at dose = 0
a_w = -50 # Shift in mean fitness with dose

a0_s = -5 # Directional selection gradients in absence of contaminants
b0_s = -5 # quadratic selection gradients in absence of contaminants
c0_s = -5 # correlational selection gradients in absence of contaminants
a1_s = -2 # Change in gradient per unit of contaminant concentration
b1_s = -2 # Change in gradient per unit of contaminant concentration
c1_s = -2 # Change in gradient per unit of contaminant concentration

# Plot x,y,w, z = 0
x = seq(-2, 2, by = .1)
y = seq(-2, 2, by = .1)

fit_f = function(x, y, Dose){
  w = mu_w + a_w * Dose +
    (a0_s + a1_s * Dose) * (x + y) +
    (b0_s + b1_s * Dose) * (x^2 + y^2) +
    (c0_s + c1_s * Dose) * (x * y)
  return(w)
}
par(mfrow = c(2, 2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,2,2) + 0.1)
plot(x, fit_f(x,y, Dose = 0), "l", ylim = c(900, 1000))
plot(x, fit_f(x,y, Dose = .25), "l", ylim = c(900, 1000))
plot(x, fit_f(x,y, Dose = .75), "l", ylim = c(900, 1000))
plot(x, fit_f(x,y, Dose = 1), "l", ylim = c(900, 1000))

par(mfrow = c(2, 2),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,2,2) + 0.1)
z = outer(x, y, function(x,y) fit_f(x,y, Dose = 0))
persp(x, y, z, main = "Dose = 0",
      col = "lightblue", 
      theta = 300, 
      xlim = c(-2, 2), 
      ylim = c(-2, 2),
      zlim = c(900, 1010))

z = outer(x, y, function(x,y) fit_f(x,y, Dose = .5))
persp(x, y, z, main = "Dose = 0.5",
      col = "lightblue", 
      theta = 300, 
      xlim = c(-2, 2), 
      ylim = c(-2, 2), 
      zlim = c(900, 1010))

z = outer(x, y, function(x,y) fit_f(x,y, Dose = .75))
persp(x, y, z, main = "Dose = 0.5",
      col = "lightblue", 
      theta = 300, 
      xlim = c(-2, 2), 
      ylim = c(-2, 2), 
      zlim = c(900, 1010))

z = outer(x, y, function(x,y) fit_f(x,y, Dose = 1))
persp(x, y, z , main = "Dose = 1",
      col = "lightblue", 
      theta = 300, 
      xlim = c(-2, 2), 
      ylim = c(-2, 2), 
      zlim = c(900, 1010))

d = data.frame(x = x, y = y) %>% 
  expand_grid(Dose = Dose)

d$w = mu_w + a_w * Dose +
  (a0_s + a1_s * d$Dose) * (d$x + d$y) +
  (b0_s + b1_s * d$Dose) * (d$x^2 + d$y^2) +
  (c0_s + c1_s * d$Dose) * (d$x * d$y)

  
lattice::xyplot(w ~ x | Dose, d, cex = .1)

w = ID %>% 
  mutate(L0_i2 = L0_i^2, Linf_i2 = Linf_i^2, alpha_i2 = alpha_i^2, 
         L0Linf_i = L0_i * Linf_i, L0alpha_i = L0_i * alpha_i, 
         Linfalpha_i = Linf_i * alpha_i) %>% 
  mutate(L0_i_sd = as.numeric(scale(L0_i)),
         Linf_i_sd = as.numeric(scale(Linf_i)),
         alpha_i_sd = as.numeric(scale(alpha_i)),
         L0_i2_sd = as.numeric(scale(L0_i2)),
         Linf_i2_sd = as.numeric(scale(Linf_i2)),
         alpha_i2_sd = as.numeric(scale(alpha_i2)),
         L0Linf_i_sd = as.numeric(scale(L0Linf_i)),
         L0alpha_i_sd = as.numeric(scale(L0alpha_i)),
         Linfalpha_i_sd = as.numeric(scale(Linfalpha_i))) %>% 
  mutate(w_i = mu_w + a_w * Dose +
           L0_i_sd   * (a0_s + a1_s * Dose)  + 
           Linf_i_sd * (a0_s + a1_s * Dose)  +
           alpha_i_sd * (a0_s + a1_s * Dose)  +
           L0_i2_sd * (b0_s + a1_s * Dose)  + 
           Linf_i2_sd * (b0_s + a1_s * Dose)  +
           alpha_i2_sd * (b0_s + a1_s * Dose) +
           L0Linf_i_sd * (c0_s + a1_s * Dose)  + 
           L0alpha_i_sd * (c0_s + a1_s * Dose)  +
           Linfalpha_i_sd * (c0_s + a1_s * Dose)) %>% 
  mutate(w_rel = w_i / mean(w_i))
hist(w$w_rel)
w %>% 
  ggplot(aes(x = L0_i_sd, y = w_rel)) +
  geom_point() +
  geom_smooth(method = "lm", 
              formula=y ~ poly(x, 2, raw=TRUE)) +
  facet_wrap(~Dose)

w %>% 
  ggplot(aes(x = Linf_i_sd, y = w_rel)) +
  geom_point() +
  geom_smooth(method = "lm", 
              formula=y ~ poly(x, 2, raw=TRUE)) +
  facet_wrap(~Dose)

w %>% 
  ggplot(aes(x = alpha_i_sd, y = w_rel)) +
  geom_point() +
  geom_smooth(method = "lm", 
              formula=y ~ poly(x, 2, raw=TRUE)) +
  facet_wrap(~Dose)



## Growth-fitness-selection joint model ----
# merge fitness with growth dataframe
df = left_join(df, w)
df$l = with(df, L/Linf_mu)
dlist = list(l = df$l,
             w = df$w_rel,
             Dose = df$Dose,
             t = df$t,
             ID = df$ID)
set.seed(42)
inits = list(mu_a0 = rexp(1, 1),
             mu_b0 = rexp(1, 1),
             mu_c0 = rexp(1, 1))
constraints = list(mu = "lower=0",
                   mu_w = "lower=0",
                   beta0 = "lower=0")

form = alist(
  # Growth model
  l ~ normal(mu, sigma_obs),
  mu <- exp(linf) * exp( log( exp(l0) / exp(linf) ) *
                           exp( - exp(alpha) * t ) ),
  l0 <- a0[ID] + a1 * Dose,
  linf <- b0[ID] + b1 * Dose,
  alpha <- c0[ID] + c1 * Dose,
  a0[ID] ~ normal(mu_a0, sigma_a0),
  b0[ID] ~ normal(mu_b0, sigma_b0),
  c0[ID] ~ normal(mu_c0, sigma_c0),
  
  # Fitness error-in-variable
  w[ID] ~ normal(mu_w, sigma_w),
  mu_w <- beta0 + 
    beta1 * a0[ID] + beta2 * b0[ID]  + beta3 * c0[ID], #+
    # gamma1 * l0[ID]^2 + gamma2 * linf[ID]^2  + gamma3 * alpha[ID]^2 +
    # gamma4 * (l0[ID] * linf[ID]) + 
    # gamma5 * (l0[ID] * alpha[ID]) +
    # gamma6 * (linf[ID] * alpha[ID]),
  
  # Dose-response for fitness
  beta0 <- beta00 + d1 * Dose,
  
  # Dose-response for selection
  beta1 <- beta10 + beta11 * Dose,
  beta2 <- beta20 + beta21 * Dose,
  beta3 <- beta30 + beta31 * Dose,
  # gamma1 <- gamma10 + gamma11 * Dose,
  # gamma2 <- gamma20 + gamma21 * Dose,
  # gamma3 <- gamma30 + gamma31 * Dose,
  # gamma4 <- gamma40 + gamma41 * Dose,
  # gamma5 <- gamma50 + gamma51 * Dose,
  # gamma6 <- gamma60 + gamma61 * Dose,
  
  # Priors for growth
  mu_a0 ~ normal(log(.15), .1),
  mu_b0 ~ normal(log(1), .001),
  mu_c0 ~ normal(log(.03), .1),
  a1 ~ normal(0, 1),
  b1 ~ normal(0, 1),
  c1 ~ normal(0, 1),
  sigma_a0 ~ exponential(1),
  sigma_b0 ~ exponential(1),
  sigma_c0 ~ exponential(1),
  sigma_obs ~ exponential(1),
  
  # Priors for fitness
  beta00 ~ normal(1000, 200),
  d1 ~ normal(0, 1),
  beta10 ~ normal(0, 1),
  beta20 ~ normal(0, 1),
  beta30 ~ normal(0, 1),
  # gamma10 ~ normal(0, 1),
  # gamma20 ~ normal(0, 1),
  # gamma30 ~ normal(0, 1),
  # gamma40 ~ normal(0, 1),
  # gamma50 ~ normal(0, 1),
  # gamma60 ~ normal(0, 1),
  beta11 ~ normal(0, 1),
  beta21 ~ normal(0, 1),
  beta31 ~ normal(0, 1),
  # gamma11 ~ normal(0, 1),
  # gamma21 ~ normal(0, 1),
  # gamma31 ~ normal(0, 1),
  # gamma41 ~ normal(0, 1),
  # gamma51 ~ normal(0, 1),
  # gamma61 ~ normal(0, 1),
  sigma_w ~ exponential(1)
  )

nlmm.5 = ulam(
  form, 
  data = dlist,
  start = inits,
  constraints = constraints,
  control = list(adapt_delta = .999,
                 max_treedepth = 15),
  chains = 1,
  cores = 1, 
  # threads = 4,
  seed = 42)
precis(nlmm.5)
stancode(nlmm.5)