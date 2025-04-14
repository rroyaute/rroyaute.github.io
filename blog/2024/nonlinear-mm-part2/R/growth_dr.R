library(tidyverse); library(patchwork); library(tidybayes)
library(brms); library(marginaleffects); library(viridis)
library(easystats); library(kableExtra)

# Setting theme
theme_set(theme_bw(16))

# Gompertz functions ----
Gomp.fun = function(t, Linf, L0, alpha){
  Lhat = Linf * exp(log(L0/Linf)*exp(-alpha*t)) # predicted growth
  return(Lhat)
}

Gomp.fun.log = function(t, Linf, L0, alpha){
  logLhat = exp(-alpha * t) * (log(L0) - log(Linf)) + log(Linf) # predicted growth
  return(logLhat)
}

Gomp.fun.dose = function(t, a0_L0, a0_Linf, a0_alpha, 
                         a1_L0, a1_Linf, a1_alpha, D){
  L0 = a0_L0 + a1_L0 * D
  Linf = a0_Linf + a1_Linf * D
  alpha = a0_alpha + a1_alpha * D
  
  Lhat = Linf * exp(log(L0/Linf)*exp(-alpha*t)) # predicted growth
  return(Lhat)
}


Gomp.fun.log.dose = function(t, a0_L0, a0_Linf, a0_alpha, 
                             a1_L0, a1_Linf, a1_alpha, D){
  L0 = a0_L0 + a1_L0 * D
  Linf = a0_Linf + a1_Linf * D
  alpha = a0_alpha + a1_alpha * D
  
  logLhat = exp(-alpha * t) * (log(L0) - log(Linf)) + log(Linf) # predicted growth
  return(logLhat)}

# # lognormal functions ----
# # lognormal to data scale
# ln_mean = function(mu, sigma) {
#   exp(mu + (0.5 * sigma^2))
# }
# ln_sd = function(mu, sigma) {
#   v <- (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)
#   return(sqrt(v))
# }
# 
# #  data scale to lognormal
# mean_ln = function(m, s) {
#   log(m / sqrt(s^2 / m^2 + 1))
# }
# sd_ln = function(m, s) {
#   sqrt(log(s^2 / m^2 + 1))
# }
# 
# Simulate growth in controls ----
n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
CV = 0.2 # random noise around mean (coefficient of variation)
a0_L0 = 182 # initial length (micrometers)
a0_Linf = 1370 # maximal length (micrometers)
a0_alpha = 0.028 # Growth rate (hour^-1)

set.seed(42)
df.0 = crossing(times, Id = 1:n_id) %>% 
  mutate(logLhat = Gomp.fun.log(times, Linf = a0_Linf, L0 = a0_L0, alpha = a0_alpha)) %>% 
  # mutate(L = rnorm(n(), exp(logLhat), sigma)) %>% 
  mutate(L = rlnorm(n(), logLhat, CV))
df.0 %>% 
  ggplot(aes(y = L, x = times)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = Gomp.fun, 
                args = list(Linf = a0_Linf, L0 = a0_L0, alpha = a0_alpha),
                color = "red", size = 1) +
  labs(x = "Hours since hatching", 
       y = expression(paste("Body-length ", "(", mu, "m", ")")))

# Simulate growth by dose ----
n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
CV = 0.1 # random noise around mean (coefficient of variation)
a0_L0 = 182 # initial length (micrometers)
a0_Linf = 1370 # maximal length (micrometers)
a0_alpha = 0.028 # Growth rate (hour^-1)
a1_L0 = -.25 * a0_L0 # Change in mean trait per unit of contaminant concentration
a1_Linf = -.25 * a0_Linf
a1_alpha = -.25 * a0_alpha

D = c(0, .1, .3, .5, .9, 1.1, 1.2)

set.seed(42)
df.dr = crossing(times, Id = 1:n_id, Dose = D,
                 a0_L0, a0_Linf, 
                 a0_alpha, a1_L0,
                 a1_Linf, a1_alpha) %>% 
  mutate(logLhat = 
           Gomp.fun.log.dose(times, a0_L0 = a0_L0, a0_Linf = a0_Linf, 
                             a0_alpha = a0_alpha, a1_L0 = a1_L0,
                             a1_Linf = a1_Linf, a1_alpha = a1_alpha, D = Dose)) %>% 
  mutate(L = rlnorm(n(), logLhat, CV))
# mutate(L = rnorm(n(), exp(logLhat), sigma))
# 
df.dr %>% 
  filter(times == 120) %>%
  group_by(Dose) %>% 
  summarise(mu_L_log = mean(logLhat), 
            sd_L_log = sd(logLhat), 
            mu_L = mean(L),
            sd_L = sd(L)) %>% 
  mutate(CV = sd_L / mu_L)

p = df.dr %>% 
  ggplot(aes(y = L, x = times)) +
  geom_point(alpha = .2, size = 2) +
  facet_wrap(~Dose) +
  # geom_function(fun = Gomp.fun.dose,
  #               args = list(t = times, a0_L0 = a0_L0, a0_Linf = a0_Linf,
  #                           a0_alpha = a0_alpha, a1_L0 = a1_L0,
  #                           a1_Linf = a1_Linf, a1_alpha = a1_alpha, 
  #                           D = Dose),
  #               color = "red", size = 1) +
  # ylim(0, 1600) +
  labs(x = "Hours since hatching", 
       y = expression(paste("Body-length ", "(", mu, "m", ")")))
p

# Fit Gompertz growth without vi ----
bf.l = bf(L ~ 
            exp(-alpha * times) * (log(L0) - log(Linf)) + log(Linf), 
          L0 + Linf + alpha ~ 1, # This trick allows to put a lower bound on the intercept parameter
          
          # nlf(L0 ~ aL0), # L0 depends on a unique value (intercept)
          # nlf(Linf ~ a0Linf), # L0 depends on a unique value (intercept)
          # nlf(alpha ~ a0alpha), # L0 depends on a unique value (intercept)
          # a0L0 + a0Linf + a0alpha ~ 1, # This trick allows to put a lower bound on the intercept parameter
          family = lognormal,
          nl = T, 
          center = F)

default_prior(bf.l, df.0)

priors = 
  # Intercept priors (control values)
  prior(normal(200, .2 * 200), nlpar = L0, class = b, lb = 0) + # constrain to positive values with lb = 0
  prior(normal(1500, .2 * 1500), nlpar = Linf, class = b, lb = 0) +
  prior(normal(.03, .2 * .03), nlpar = alpha, class = b, lb = 0) +
  # Residual prior
  prior(exponential(10), class = sigma)

gomp.prior = brm(bf.l, 
                 df.0,
                 backend = "cmdstan",
                 prior = priors, 
                 iter = 1000,
                 sample_prior = "only",
                 file = "blog/2024/nonlinear-mm-part2/mods/gomp.prior",
                 seed = 42, 
                 cores = 4,
                 threads = threading(3))
model_parameters(gomp.prior, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.prior))
plot(conditional_effects(gomp.prior, method = "posterior_predict"))
plot(conditional_effects(gomp.prior, 
                         ndraws = 100, spaghetti = T), points = T)
pp_check(gomp.prior, ndraws = 100) + xlim(0, 3000)

gomp.post = brm(bf.l, 
                df.0,
                backend = "cmdstan",
                prior = priors,
                sample_prior = "yes",
                warmup = 3000,
                iter = 4000,
                file = "blog/2024/nonlinear-mm-part2/mods/gomp.post",
                chains = 4,
                seed = 42, 
                cores = 4,
                threads = 3)
model_parameters(gomp.post, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.post, method = "posterior_predict"), 
     points = T)
plot(conditional_effects(gomp.post, 
                         ndraws = 100, 
                         spaghetti = T), 
     points = T)
pp_check(gomp.post, ndraws = 100)

# Fit dose-response without vi ----
bf.l = bf(L ~ 
            exp(-alpha * times) * (log(L0) - log(Linf)) + log(Linf), 
          nlf(L0 ~ a0L0 + a1L0 * Dose), # dose-response formula on L0
          nlf(Linf~ a0Linf + a1Linf * Dose),
          nlf(alpha ~ a0alpha + a1alpha * Dose),
          a0L0 + a0Linf + a0alpha + 
            a1L0 + a1Linf + a1alpha ~ 1, # dose-response parameters have a single value
          family = lognormal,
          nl = T, 
          center = F)

default_prior(bf.l, df.dr)

priors = 
  # Intercept priors (control values)
  prior(normal(200, .2 * 200), nlpar = a0L0, class = b, lb = 0) + # constrain to positive values with lb = 0
  prior(normal(1500, .2 * 1500), nlpar = a0Linf, class = b, lb = 0) +
  prior(normal(.03, .2 * .03), nlpar = a0alpha, class = b, lb = 0) +
  # Regression priors
  prior(normal(0, .2 * 200), nlpar = a1L0, class = b) + # slope SD is +/- 20 % around control value
  prior(normal(0, .2 * 1500), nlpar = a1Linf, class = b) +
  prior(normal(0, .2 * .03), nlpar = a1alpha, class = b) +
  # Residual prior
  prior(exponential(10), class = sigma)

gomp.dr.prior = brm(bf.l, 
                    df.dr,
                    backend = "cmdstan",
                    prior = priors, 
                    iter = 1000,
                    sample_prior = "only",
                    file = "blog/2024/nonlinear-mm-part2/mods/gomp.dr.prior",
                    seed = 42, 
                    cores = 4,
                    threads = threading(3))
model_parameters(gomp.dr.prior, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.dr.prior))
plot(conditional_effects(gomp.dr.prior, method = "posterior_predict"))
plot(conditional_effects(gomp.dr.prior, 
                         ndraws = 100, spaghetti = T), points = T)
pp_check(gomp.dr.prior, ndraws = 100) + xlim(0, 3000)

gomp.dr.post = brm(bf.l, 
                   df.dr,
                   backend = "cmdstan",
                   prior = priors,
                   sample_prior = "yes",
                   warmup = 3000,
                   iter = 4000,
                   file = "blog/2024/nonlinear-mm-part2/mods/gomp.dr.post",
                   chains = 4,
                   seed = 42, 
                   cores = 4,
                   threads = 3)
model_parameters(gomp.dr.post, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.dr.post, method = "posterior_predict"), 
     points = T)
plot(conditional_effects(gomp.dr.post, 
                         ndraws = 100, 
                         spaghetti = T), 
     points = T)
pp_check(gomp.dr.post, ndraws = 100)

# Simulate individual differences in growth ----
n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
CVi = 0.25 # among-individual differences (coefficient of variation)
CVw = 0.05 # random noise around mean (coefficient of variation)
a0_L0 = 182 # initial length (micrometers)
a0_Linf = 1370 # maximal length (micrometers)
a0_alpha = 0.028 # Growth rate (hour^-1)

# Setup matrix for generating individual parameters from MVN distribution
rho = -.7 # Assume strong negative correlation between Linf and alpha
Mu     = c(0, 0, 0)
sigmas = c(a0_L0 * CVi, a0_Linf * CVi, a0_alpha * CVi) # 20 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = 1:n_id,
         L0_i = L0_i + a0_L0,
         Linf_i = Linf_i + a0_Linf,
         alpha_i = alpha_i + a0_alpha)
ID %>% 
  select(-ID) %>% 
  GGally::ggpairs() +
  theme_bw() 


df.0.vi = 
  crossing(times, ID) %>% 
  mutate(logLhat = Gomp.fun.log(times, Linf = Linf_i, L0 = L0_i, alpha = alpha_i)) %>% 
  mutate(L = rlnorm(n(), logLhat, CVw))
df.0.vi %>% 
  ggplot(aes(y = L, x = times, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .2, size = .5) +
  geom_function(fun = Gomp.fun, 
                args = list(Linf = a0_Linf, L0 = a0_L0, alpha = a0_alpha),
                color = "red", size = 1) +
  labs(x = "Hours since hatching", 
       y = expression(paste("Body-length ", "(", mu, "m", ")")))


# Fit Gompertz growth with vi ----
bf.l = bf(L ~ 
            exp(-alpha * times) * (log(L0) - log(Linf)) + log(Linf), 
          L0 + Linf + alpha ~ 1 + (1|c|ID), # This trick allows to put a lower bound on the intercept parameter
          family = lognormal,
          nl = T, 
          center = F)

default_prior(bf.l, df.0.vi)



priors = 
  # Intercept priors (control values)
  prior(normal(200, .2 * 200), nlpar = L0, class = b, lb = 0) + # constrain to positive values with lb = 0
  prior(normal(1500, .2 * 1500), nlpar = Linf, class = b, lb = 0) +
  prior(normal(.03, .2 * .03), nlpar = alpha, class = b, lb = 0) +
  # among-individual variance prior (20 % of variation around mean)
  prior(exponential(0.025), class = sd, nlpar = L0) +
  prior(exponential(0.0033), class = sd, nlpar = Linf) +
  prior(exponential(166.66), class = sd, nlpar = alpha) +
  # Residual prior
  prior(exponential(10), class = sigma) +
  # Correlation prior
  prior(lkj(2), class = cor)

# priors = 
#   # Intercept priors (control values)
#   prior(normal(200, .2 * 200), nlpar = L0, class = b, lb = 0) + # constrain to positive values with lb = 0
#   prior(normal(1500, .2 * 1500), nlpar = Linf, class = b, lb = 0) +
#   prior(normal(.03, .2 * .03), nlpar = alpha, class = b, lb = 0) +
#   # R2D2 sd prior
#   prior(R2D2(mean_R2 = .8, prec_R2 = 10, 
#                     cons_D2 = .5, main = T), class = "sd") +
#   # Correlation prior
#   prior(lkj(2), class = cor)
#   


gomp.prior.vi = brm(bf.l, 
                    df.0.vi,
                    backend = "cmdstan",
                    prior = priors, 
                    iter = 1000,
                    sample_prior = "only",
                    # file = "blog/2024/nonlinear-mm-part2/mods/gomp.prior.vi",
                    seed = 42, 
                    cores = 4,
                    threads = threading(3))
model_parameters(gomp.prior.vi, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.prior.vi, re_formula = NA))
plot(conditional_effects(gomp.prior.vi, re_formula = NA, method = "posterior_predict"))
plot(conditional_effects(gomp.prior.vi, 
                         ndraws = 100, spaghetti = T), points = T)
pp_check(gomp.prior.vi, ndraws = 100) + xlim(0, 3000)

gomp.post.vi = brm(bf.l, 
                   df.0.vi,
                   backend = "cmdstan",
                   prior = priors,
                   sample_prior = "yes",
                   warmup = 3000,
                   iter = 4000,
                   file = "blog/2024/nonlinear-mm-part2/mods/gomp.post.vi",
                   chains = 4,
                   seed = 42, 
                   cores = 4,
                   threads = 3)
model_parameters(gomp.post.vi, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.post.vi, re_formula = NA), 
     points = T)
plot(conditional_effects(gomp.post.vi, re_formula = NA, method = "posterior_predict"), 
     points = T)
plot(conditional_effects(gomp.post.vi, 
                         ndraws = 100, 
                         spaghetti = T), 
     points = T)
pp_check(gomp.post.vi, ndraws = 100)

# Simulate dose-response with vi ----
n_id = 30 # 30 individuals
times = seq(0, 126, by = 24) # One observation every day
CVi = 0.25 # among-individual differences (coefficient of variation)
CVw = 0.05 # random noise around mean (coefficient of variation)
a0_L0 = 182 # initial length (micrometers)
a0_Linf = 1370 # maximal length (micrometers)
a0_alpha = 0.028 # Growth rate (hour^-1)
a1_L0 = -.25 * a0_L0 # Change in mean trait per unit of contaminant concentration
a1_Linf = -.25 * a0_Linf
a1_alpha = -.25 * a0_alpha

D = c(0, .1, .3, .5, .9, 1.1, 1.2)

# Setup matrix for generating individual parameters from MVN distribution
rho = -.7 # Assume strong negative correlation between Linf and alpha
Mu     = c(0, 0, 0)
sigmas = c(a0_L0 * CVi, a0_Linf * CVi, a0_alpha * CVi) # 20 % CV around the mean
rho_mat = matrix(c(1, 0, 0,
                   0, 1, rho,
                   0, rho, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID = MASS::mvrnorm(n_id * length(D), Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("L0_i", "Linf_i", "alpha_i") %>% 
  mutate(ID = 1:(n_id * length(D)),
         Dose = rep(D, each = n_id)) %>% 
  mutate(ID = paste(Dose, ID, sep = "_"),
         L0_i = a0_L0 + a1_L0 * Dose + L0_i,
         Linf_i = a0_Linf + a1_Linf * Dose + Linf_i,
         alpha_i = a0_alpha + a1_alpha * Dose + alpha_i)

ID %>% 
  select(-ID) %>% 
  GGally::ggpairs() +
  theme_bw() 

df.dr.vi = 
  crossing(times, ID) %>% 
  mutate(logLhat = 
           Gomp.fun.log.dose(times, a0_L0 = a0_L0, a0_Linf = a0_Linf, 
                             a0_alpha = a0_alpha, a1_L0 = a1_L0,
                             a1_Linf = a1_Linf, a1_alpha = a1_alpha, D = Dose)) %>% 
  mutate(L = rlnorm(n(), logLhat, CVw))

p = df.dr.vi %>% 
  ggplot(aes(y = L, x = times)) +
  geom_point(alpha = .2, size = 2) +
  facet_wrap(~Dose) +
  # geom_function(fun = Gomp.fun.dose,
  #               args = list(t = times, a0_L0 = a0_L0, a0_Linf = a0_Linf,
  #                           a0_alpha = a0_alpha, a1_L0 = a1_L0,
  #                           a1_Linf = a1_Linf, a1_alpha = a1_alpha, 
  #                           D = Dose),
  #               color = "red", size = 1) +
  # ylim(0, 1600) +
  labs(x = "Hours since hatching", 
       y = expression(paste("Body-length ", "(", mu, "m", ")")))
p


# Fit dose-response with vi ----
## Define model ----
bf.l = bf(L ~ 
            exp(-alpha * times) * (log(L0) - log(Linf)) + log(Linf), 
          nlf(L0 ~ a0L0 + a1L0 * Dose), # dose-response formula on L0
          nlf(Linf~ a0Linf + a1Linf * Dose),
          nlf(alpha ~ a0alpha + a1alpha * Dose),
          a0L0 + a0Linf + a0alpha ~  1 + (1 |c| ID),
          a1L0 + a1Linf + a1alpha ~ 1, # dose-response parameters have a single value
          family = lognormal,
          nl = T, 
          center = F)

default_prior(bf.l, df.dr.vi)

## Define and plot priors ----
priors = 
  # Intercept priors (control values)
  prior(normal(200, .05 * 200), nlpar = a0L0, class = b, lb = 0) + # constrain to positive values with lb = 0
  prior(normal(1500, .05 * 1500), nlpar = a0Linf, class = b, lb = 0) +
  prior(normal(.03, .05 * .03), nlpar = a0alpha, class = b, lb = 0) +
  # among-individual variance prior (20 % of variation around mean)
  prior(normal(.1 * 200, .05 * 200), class = sd, nlpar = a0L0) +
  prior(normal(.1 * 1500, .05 * 1500), class = sd, nlpar = a0Linf) +
  prior(normal(.1 * .03, .05 * .03), class = sd, nlpar = a0alpha) +
  # prior(exponential(0.025), class = sd, nlpar = a0L0) +
  # prior(exponential(0.0033), class = sd, nlpar = a0Linf) +
  # prior(exponential(166.66), class = sd, nlpar = a0alpha) +
  # Regression priors
  prior(normal(0, .1 * 200), nlpar = a1L0, class = b) + # slope SD is +/- 10 % around control value
  prior(normal(0, .1 * 1500), nlpar = a1Linf, class = b) +
  prior(normal(0, .1 * .03), nlpar = a1alpha, class = b) +
  # Residual prior
  prior(normal(0, .1), class = sigma) +
  # prior(exponential(10), class = sigma) +
  # Correlation prior
  prior(lkj(1), class = cor)

p1 = priors %>% 
  parse_dist() %>% 
  filter(class == "b") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Intercepts") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

p2 = priors %>% 
  parse_dist() %>% 
  filter(class == "sd") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Among-individual variance") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

p3 = priors %>% 
  parse_dist() %>% 
  filter(class == "sigma") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Residual variance") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

(p1 + p2 + p3) + plot_layout(ncol = 1)


## Fit models ----
gomp.prior.vi.dr = brm(bf.l, 
                       df.dr.vi,
                       backend = "cmdstan",
                       prior = priors, 
                       iter = 1000,
                       sample_prior = "only",
                       # file = "blog/2024/nonlinear-mm-part2/mods/gomp.prior.vi.dr",
                       seed = 42, 
                       cores = 4,
                       threads = threading(3))
model_parameters(gomp.prior.vi.dr, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.prior.vi.dr, re_formula = NA))
plot(conditional_effects(gomp.prior.vi.dr, re_formula = NA, method = "posterior_predict"))
plot(conditional_effects(gomp.prior.vi.dr, 
                         ndraws = 100, spaghetti = T), points = T)
pp_check(gomp.prior.vi.dr, ndraws = 100) 

gomp.post.vi.dr = brm(bf.l, 
                      df.dr.vi,
                      backend = "cmdstan",
                      prior = priors,
                      sample_prior = "yes",
                      # warmup = 3000,
                      # iter = 4000,
                      control = list(adapt_delta = .995,
                                     max_treedepth = 20),
                      # file = "blog/2024/nonlinear-mm-part2/mods/gomp.post.vi.dr",
                      chains = 4,
                      seed = 42, 
                      cores = 4,
                      threads = 3)
model_parameters(gomp.post.vi.dr, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.post.vi.dr, re_formula = NA), 
     points = T)
plot(conditional_effects(gomp.post.vi.dr, re_formula = NA, method = "posterior_predict"), 
     points = T)
plot(conditional_effects(gomp.post.vi.dr, 
                         ndraws = 100, 
                         spaghetti = T), 
     points = T)
pp_check(gomp.post.vi.dr, ndraws = 100)
shinystan::launch_shinystan(gomp.post.vi.dr)

# Model reparametrization: divide by Linf ----
df.dr.vi$l = df.dr.vi$L / a0_Linf

p = df.dr.vi %>% 
  ggplot(aes(y = l, x = times)) +
  geom_point(alpha = .2, size = 2) +
  facet_wrap(~Dose) +
  # geom_function(fun = Gomp.fun.dose,
  #               args = list(t = times, a0_L0 = a0_L0, a0_Linf = a0_Linf,
  #                           a0_alpha = a0_alpha, a1_L0 = a1_L0,
  #                           a1_Linf = a1_Linf, a1_alpha = a1_alpha, 
  #                           D = Dose),
  #               color = "red", size = 1) +
  # ylim(0, 1600) +
  labs(x = "Hours since hatching", 
       y = expression(paste("Body-length ", "(", mu, "m", ")")))
p

bf.l = bf(l ~ 
            exp(-alpha * times) * (log(l0) - log(linf)) + log(linf), 
          nlf(l0 ~ a0l0 + a1l0 * Dose), # dose-response formula on L0
          nlf(linf ~ a0linf + a1linf * Dose),
          nlf(alpha ~ a0alpha + a1alpha * Dose),
          a0l0 + a0linf + a0alpha ~  1 + (1 |c| ID),
          a1l0 + a1linf + a1alpha ~ 1, # dose-response parameters have a single value
          family = lognormal,
          nl = T, 
          center = F)

default_prior(bf.l, df.dr.vi)

priors = 
  # Intercept priors (control values)
  prior(normal(0.13, .05 * 0.13), nlpar = a0l0, class = b, lb = 0) + # constrain to positive values with lb = 0
  prior(constant(1), nlpar = a0linf, class = b) +
  prior(normal(.03, .05 * .03), nlpar = a0alpha, class = b, lb = 0) +
  # among-individual variance prior (20 % of variation around mean)
  prior(normal(.1 * 0.13, .05 * 0.13), class = sd, nlpar = a0l0) +
  prior(normal(.1 , .05), class = sd, nlpar = a0linf) +
  prior(normal(.1 * .03, .05 * .03), class = sd, nlpar = a0alpha) +
  
  # prior(exponential(38), class = sd, nlpar = a0l0) +
  # prior(exponential(5), class = sd, nlpar = a0linf) +
  # prior(exponential(166), class = sd, nlpar = a0alpha) +
  # Regression priors
  prior(normal(0, .15 * 0.13), nlpar = a1l0, class = b) + 
  prior(normal(0, .15), nlpar = a1linf, class = b) +
  prior(normal(0, .15 * .03), nlpar = a1alpha, class = b) +
  # Residual prior
  prior(exponential(10), class = sigma) +
  # Correlation prior
  prior(lkj(1), class = cor)

gomp.prior.vi.dr.sc = brm(bf.l, 
                       df.dr.vi,
                       backend = "cmdstan",
                       prior = priors, 
                       iter = 1000,
                       sample_prior = "only",
                       # file = "blog/2024/nonlinear-mm-part2/mods/gomp.prior.vi.dr.sc",
                       seed = 42, 
                       cores = 4,
                       threads = threading(3))
model_parameters(gomp.prior.vi.dr.sc, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.prior.vi.dr.sc, re_formula = NA))
plot(conditional_effects(gomp.prior.vi.dr.sc, re_formula = NA, method = "posterior_predict"))
plot(conditional_effects(gomp.prior.vi.dr.sc, 
                         ndraws = 100, spaghetti = T), points = T)
pp_check(gomp.prior.vi.dr.sc, ndraws = 100) 

gomp.post.vi.dr.sc = brm(bf.l, 
                      df.dr.vi,
                      backend = "cmdstan",
                      prior = priors,
                      sample_prior = "yes",
                      # warmup = 3000,
                      # iter = 4000,
                      # control = list(adapt_delta = .995,
                      #                max_treedepth = 20),
                      # file = "blog/2024/nonlinear-mm-part2/mods/gomp.post.vi.dr.sc",
                      chains = 4,
                      seed = 42, 
                      cores = 4,
                      threads = 3)
model_parameters(gomp.post.vi.dr.sc, effects = "all") %>%
  kable(digits = 2)
plot(conditional_effects(gomp.post.vi.dr.sc, re_formula = NA), 
     points = T)
plot(conditional_effects(gomp.post.vi.dr.sc, re_formula = NA, method = "posterior_predict"), 
     points = T)
plot(conditional_effects(gomp.post.vi.dr.sc, 
                         ndraws = 100, 
                         spaghetti = T), 
     points = T)
pp_check(gomp.post.vi.dr.sc, ndraws = 100)
shinystan::launch_shinystan(gomp.post.vi.dr.sc)
