## General settings
set.seed(123) #for reproducibility
# days at which we assume outcome is measured
timevec <- c(0.1,1,3,5,7,10,14,21,28,35,42)

#different number of individuals per dose to make it clearer which is which
#also, that's the structure of the data which motivated the tutorial
Nlow = 7; Nmed = 8; Nhigh = 9; 
Nind = Nlow + Nmed + Nhigh; #total number of individuals

# Set values for dose
# since we only consider dose on a log scale
# we'll log transform right here and then always use it in those log units
high_dose = log(1000)
med_dose = log(100)
low_dose = log(10)
dosevec = c(rep(low_dose,Nlow),rep(med_dose,Nmed),rep(high_dose,Nhigh))
# we are also creating a version of the dose variable
# that consists of ordered categories instead of numeric values
# we'll use that mostly for plotting
dosevec_cat = ordered(c(rep("low", Nlow),
                        rep("medium",Nmed),
                        rep("high",Nhigh)),
                      levels=c("low","medium","high"))


## Setting parameter values
sigma = 1
a1 = 0.1
b1 = -0.1
m3_mua = 3
m3_mub = 1
m3_sigmaa = 0.1
m3_sigmab = 0.1
m3_a0 = rnorm(n=Nind, m3_mua, m3_sigmaa)
m3_b0 = rnorm(n=Nind, m3_mub, m3_sigmaa)

m3pars = c(sigma = sigma, a1 = a1, b1 = b1,
           a0_mu = m3_mua, b0_mu = m3_mub)
m3_alpha = m3_a0 + a1*(dosevec - med_dose)
m3_beta = m3_b0 + b1*(dosevec - med_dose)
m3_mu =  exp(m3_alpha) %*% t(log(timevec)) - exp(m3_beta) %*% t(timevec)
m3_y = rnorm(length(m3_mu),m3_mu, sigma)
m3_dat <- data.frame(id = rep(1:Nind,length(timevec)),
                     dose = rep(dosevec,length(timevec)),
                     dose_adj = rep(dosevec,length(timevec))-med_dose,
                     dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Nind)),
                     mu = as.vector(m3_mu),
                     outcome = as.vector(m3_y),
                     model = "m3")

m3_dat <- data.frame(id = rep(1:Nind,length(timevec)),
                     dose = rep(dosevec,length(timevec)),
                     dose_adj = rep(dosevec,length(timevec))-med_dose,
                     dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Nind)),
                     mu = as.vector(m3_mu),
                     outcome = as.vector(m3_y),
                     model = "m3")

ggplot(m3_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)") +
  theme_minimal()
#adaptive priors, partial-pooling model
#2x(N+2)+1 parameters
m4 <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0[id] + a1*dose_adj,
  beta <-  b0[id] + b1*dose_adj,
  a0[id] ~ dnorm(mu_a,  sigma_a),
  b0[id] ~ dnorm(mu_b, sigma_b),
  mu_a ~ dnorm(2, 1),
  mu_b ~ dnorm(0.5, 1),
  sigma_a ~ cauchy(0, 1),
  sigma_b ~ cauchy(0, 1),
  a1 ~ dnorm(0.3, 1),
  b1 ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0, 1))

fit <- ulam(m4,
            data = m3_dat,
            # start=startlist[[n]],
            # constraints=constraints,
            chains = 4, 
            cores = 4)
            # warmup = warmup, 
            # iter = iter,
            # seed = seed)
precis(fit, depth = 2)
stancode(fit)
