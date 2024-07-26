library(tidyverse); library(mvgam); library(gratia)
library(rethinking)
theme_set(theme_bw(14))

mu_w = 1 # Mean relative fitness

b1 = -.05 # Negative directional selection for x
b2 = .05 # Positive directional selection for y (small value to avoid 0 fitness values)
g11 = .05 # Disruptive selection for x
g22 = -.05 # Stabilizing selection for y (small value to avoid 0 fitness values)
g12 = .01 # Correlational selection for x,y

# Plot x,y,w, z = 0
x = seq(-3, 3, by = .5)
y = seq(-3, 3, by = .5)

fit_f = function(x, y){
  w = mu_w +
    b1 * x + b2 * y  +
    g11 * x^2 + g22 * y^2 +
    g12 * (x * y)
  return(w)
}
# w,x with y fixed at the mean value (0)
plot(x, fit_f(x, y=0), "l")
# w,y with x fixed at the mean value (0)
plot(y, fit_f(x=0, y), "l")

# Selection surface
z = outer(x, y, function(x,y) fit_f(x,y))
persp(x, y, z,
      col = "lightblue",
      xlim = c(-3, 3), 
      ylim = c(-3, 3),
      theta = 180+60)

df = crossing(x = x,
             y = y) %>% 
  mutate(w = mu_w +
           b1 * x + b2 * y  +
           g11 * x^2 + g22 * y^2 +
           g12 * (x * y))

gam.1 = gam(w ~ s(x) + s(y) + ti(x, y), data = df)
summary(gam.1)
plot.gam(gam.1)
draw(gam.1)
vis.gam(gam.1,
        color = "topo", 
        theta = 300, 
        n.grid = 50, 
        lwd = 0.4)
# Get standardized gradients
# TODO


# Simulate linear increase in selection with dose
mu_w = 1 # Mean relative fitness

b1_0 = -.05 # Negative directional selection for x at dose = 0
b2_0 = .05 # Positive directional selection for y  at dose = 0
g11_0 = .05 # Disruptive selection for x at dose = 0
g22_0 = -.05 # Stabilizing selection for y at dose = 0
g12_0 = -.01 # Correlational selection for x,y at dose = 0

d_w = -.5 # Average fitness decreases by 0.5 per dose unit
b1_1 = -.5 # Change in linear selection on x per dose unit
b2_1 = -.5 # Change in linear selection on y per dose unit
g11_1 = -.5 # Change in quadratic selection on x per dose unit
g22_1 = -.5 # Change in quadratic selection on y per dose unit
g12_1 = -.5 # # Change in correlational selection on x,y per dose unit


# Plot x,y,w, z = 0
x = seq(-3, 3, by = .5)
y = seq(-3, 3, by = .5)

fit_f_dose = function(x, y, dose){
  w = mu_w + d_w * dose +
    (b1_0 + b1_1 * dose) * x + 
    (b2_0 + b2_1 * dose) * y  +
    (g11_0 + g11_1 * dose) * x^2 + 
    (g22_0 + g22_1 * dose) * y^2 +
    (g12_0 + g12_1 * dose) * (x * y)
  return(w)
}
# w,x with y fixed at the mean value (0)
plot(x, fit_f_dose(x, y=0, dose = 0), "l")
plot(x, fit_f_dose(x, y=0, dose = .5), "l")
plot(x, fit_f_dose(x, y=0, dose = 1), "l")

# w,y with x fixed at the mean value (0)
plot(y, fit_f_dose(x=0, y, dose = 0), "l")
plot(y, fit_f_dose(x=0, y, dose = .5), "l")
plot(y, fit_f_dose(x=0, y, dose = 1), "l")


df = crossing(x = x,
              y = y,
              dose = seq(0, 1, .1)) %>% 
  mutate(w = mu_w + d_w * dose +
           (b1_0 + b1_1 * dose) * x + 
           (b2_0 + b2_1 * dose) * y  +
           (g11_0 + g11_1 * dose) * x^2 + 
           (g22_0 + g22_1 * dose) * y^2 +
           (g12_0 + g12_1 * dose) * (x * y))


# Estimate selection Lande & Arnold style ----
# Assume x,y uncorrelated for now

form = alist(
  # Growth model
  x ~ normal(mu_x, sigma_x),
  y ~ normal(mu_y, sigma_y),
  
  # Fitness error-in-variable
  w ~ normal(mu_w, sigma_w),
  mu_w <- b0 + 
    b1 * x + b2 * y +
    g11 * x^2 + g22 * y^2 +
    g12 * (x*y),
  # Dose-response for fitness
  b0 <- d0_w + d1_w * dose,
  
  # Dose-response for selection
  b1 <- d0_x + d1_x * dose,
  b2 <- d0_y + d1_y * dose,
  g11 <- d0_x2 + d1_x2 * dose, 
  g22 <- d0_y2 + d1_y2 * dose, 
  g12 <- d0_xy + d1_xy * dose, 
  
  # Priors for x,y,w
  mu_x ~ normal(0, 1),
  mu_y ~ normal(0, 1),
  mu_w ~ normal(1, .2),
  
  # Priors for dose-response
  # gradients at dose = 0
  d0_w ~ normal(0, 1),
  d0_x ~ normal(0, 1),
  d0_y ~ normal(0, 1),
  d0_x2 ~ normal(0, 1),
  d0_y2 ~ normal(0, 1),
  d0_xy ~ normal(0, 1),
  # gradients change with dose
  d1_w ~ normal(0, 1),
  d1_x ~ normal(0, 1),
  d1_y ~ normal(0, 1),
  d1_x2 ~ normal(0, 1),
  d1_y2 ~ normal(0, 1),
  d1_xy ~ normal(0, 1),
  
  sigma_x ~ exponential(1),
  sigma_y ~ exponential(1),
  sigma_w ~ exponential(1))

set.seed(42)
inits = list(mu_x = runif(1, -1, 1),
             mu_y = rnorm(1, -1, 1),
             mu_w = rnorm(1, -.4, .4))

stancode(form, data = df)

fit.1 = ulam(
  form, 
  data = df,
  start = inits,
  # constraints = constraints,
  control = list(adapt_delta = .999,
                 max_treedepth = 15),
  chains = 1,
  cores = 1, 
  # threads = 4,
  seed = 42)
precis(fit.1)
stancode(fit.1)


# GAM estimation ----
df %>% 
  ggplot(aes(x = x, y = w, level = dose)) +
  geom_point() +
  geom_smooth() + 
  facet_grid(y~dose)

gam.2 = gam(w ~ s(x, by = dose) + 
              s(y, by = dose) + 
              ti(x, y, by = dose), 
            data = df)
summary(gam.2)
plot.gam(gam.2)
draw(gam.2)
vis.gam(gam.2, 
        view = c("x", "y"),
        cond = list(dose = 0),
        color = "topo", 
        theta = 300, 
        n.grid = 50, 
        lwd = 0.4)
vis.gam(gam.2, 
        view = c("x", "y"),
        cond = list(dose = .5),
        theta = 300, 
        n.grid = 50, 
        lwd = 0.4)
vis.gam(gam.2, 
        view = c("x", "y"),
        cond = list(dose = 1),
        theta = 300, 
        n.grid = 50, 
        lwd = 0.4)


gam.3 = gam(w ~ s(x) + 
              s(y) + 
              ti(x, y) + 
              s(dose, bs = "re"), 
            data = df)
summary(gam.3)
plot.gam(gam.3)
draw(gam.3)
vis.gam(gam.3, 
        view = c("x", "y"),
        cond = list(dose = 0),
        color = "topo", 
        theta = 300, 
        n.grid = 50, 
        lwd = 0.4)
vis.gam(gam.3, 
        view = c("x", "y"),
        cond = list(dose = .5),
        theta = 300, 
        n.grid = 50, 
        lwd = 0.4)
vis.gam(gam.3, 
        view = c("x", "y"),
        cond = list(dose = 1),
        theta = 300, 
        n.grid = 50, 
        lwd = 0.4)

predictions(gam.3, 
                 newdata = datagrid(dose = df$dose,
                                    x = -3:3)) %>% 
  data.frame() %>% 
  ggplot(aes(x, estimate, level = dose)) +
  geom_line() + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  facet_grid(~dose) +
  labs(y = "Predicted selection", x = "x")

predictions(gam.3, 
            newdata = datagrid(dose = df$dose,
                               y = -3:3)) %>% 
  ggplot(aes(x = y, estimate, level = dose)) +
  geom_line() + 
  facet_grid(~dose) +
  labs(y = "Predicted fitness", x = "y")

