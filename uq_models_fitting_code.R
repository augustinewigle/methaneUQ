library(runjags)
library(bayesplot)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(HDInterval)
library(tidyr)
library(stringr)

second_trial_data <- read.csv("trial2_anon.csv")
first_trial_data <- read.csv("trial1_anon.csv")

source("helper_functions.R")

# QOGI provider C - uncertainty -------------------------------------------------

qogiC_dat <- data_prep(second_trial_data, technology_name = "qogiC")

qogiC_dat_external <- data_prep(first_trial_data, technology_name = "qogiC")

# Create model string
model_qogiC <- make_model(normalmean_str = "log(alpha0+alpha1*Q[i] + alpha2*Q[i]*Q[i]*(1-ind[i]) + (ind[i])*(beta0 + beta1*Q[i]))
                         ind[i] <- (Q[i] > t)",
                          normalprec_str = "tau_base + Q[i]/c",
                          priors_str = "

                       beta0 <- alpha2*t^2-beta1*t
                       tau_base <- 1/sigma_base^2
                       sigma_base ~ dnorm(0,0.01) T(0,) # Gelman 2006 - this is the SD of Qmeas when Q=0
                       c ~ dnorm(0, 0.01) T(0,)
                       alpha0 ~ dgamma(0.5,2)
                       alpha1 ~ dlnorm(0,1)
                       alpha2 ~ dnorm(0,1) T(0,)
                       beta1 ~ dlnorm(0, 1)
                       t ~ dunif(0,80)
                       ")
# Look at JAGS model code
cat(model_qogiC)

# Run the model in JAGS
samples_qogiC <- run.jags(model = model_qogiC,
                        monitor = c("alpha0", "alpha1","alpha2","beta1", "beta0","t", "tau_base","c"),
                        data = qogiC_dat,
                        n.chains = 2,
                        burnin = 10000,
                        inits = list(list(alpha0 = 0.001,
                                          alpha1 = rlnorm(1),
                                          alpha2 = 0.0001,
                                          sigma_base = runif(1, 0, 50),
                                          t = 1),
                                     list(alpha0 = 0.01,
                                          alpha1 = rlnorm(1),
                                          alpha2 = 0.0001,
                                          sigma_base = runif(1, 0, 50),
                                          t = 2)))

# Assess convergence
mcmc_trace(samples_qogiC$mcmc)
summary(samples_qogiC)

# Extracting the DIC
runjags::extract(samples_qogiC, what = "dic")

# Plotting prediction bands
# ignore warning
scatter_gg(data = qogiC_dat,
           newdata = qogiC_dat_external,
           samples_obj = samples_qogiC,
           Qtrue_vals = seq(0, 80, length.out = 100),
           plot_Qmeas_val_func = quad_linear,
           title = "QOGI C", hdi = T)

# NIR HSI Uncertainty -----------------------------------------------------------------
nirhsi_dat <- data_prep(all_data = second_trial_data,
                        technology_name = "nirhsi")

nirhsi_dat_external <- data_prep(all_data = first_trial_data,
                                 technology_name = "nirhsi")

model_nirhsi <- make_model(normalmean_str = "log(alpha0 + alpha1*Q[i])",
                        normalprec_str = "tau",
                        priors_str = "

                         tau <- 1/sigma^2
                         sigma ~ dnorm(0,0.01) T(0,) # Gelman 2006
                         alpha0 ~ dgamma(15,1.5)
                         alpha1 ~ dlnorm(0,1)
                         ")
cat(model_nirhsi)

samples_nirhsi <- run.jags(model = model_nirhsi,
                        monitor = c("alpha0","alpha1","sigma"),
                        data = nirhsi_dat,
                        n.chains = 2,
                        burnin = 10000,
                        inits = list(list(alpha0 = 0.001,
                                          alpha1 = rlnorm(1),
                                          sigma = runif(1, 0, 50)),
                                     list(alpha0 = 0.01,
                                          alpha1 = rlnorm(1),
                                          sigma = runif(1, 0, 50))))
# Assess convergence
mcmc_trace(samples_nirhsi$mcmc)
summary(samples_nirhsi)

# Extracting the DIC
runjags::extract(samples_nirhsi, what = "dic")

scatter_gg(data = nirhsi_dat,
           newdata = nirhsi_dat_external,
           samples_obj = samples_nirhsi,
           Qtrue_vals = seq(0, 80, length.out = 50),
           plot_Qmeas_val_func = simple_linear_func,
           title = "Aerial NIRHSI", hdi= T)

# Code for Application Section 4.2 -----------------------------------------------------

# Get posterior samples for QOGI C, thin = 2 means use half of the posterior samples
posterior <- combine.mcmc(samples_qogiC$mcmc, thin = 2)

M <- 5000 # Number of samples to take from the prior distribution
n_imp <- 4000 # number of samples to draw from desired dist in importance sampling

set.seed(23)
true_Q <- 25
Qmeas_new <- quad_linear(true_Q, posterior = as.matrix(posterior[sample(1:nrow(posterior), 5),]))

# Two different priors
Qprior1 <- runif(M, 0, 200)
Qprior2 <- rlnorm(M, meanlog = 2.6, sdlog = 1) # 2.6 comes from Zavala-Araiza et al 2018
p1 <- hist(Qprior1, freq = F, breaks = 50)
p2 <- hist(Qprior2, freq = F, breaks = 50)

# Get results - may take several mins to run
post_prior1_Qmeas1 <- predict_Q(Qmeas_new = Qmeas_new[1], Qnew_prior_sample = Qprior1,
                                posterior = posterior, n_imp = n_imp,
                                density_func = qogiC_density_func)
post_prior1_Qmeas2 <- predict_Q(Qmeas_new = Qmeas_new, Qnew_prior_sample = Qprior1,
                                posterior = posterior, n_imp = n_imp,
                                density_func = qogiC_density_func)


post_prior2_Qmeas1 <- predict_Q(Qmeas_new = Qmeas_new[1], Qnew_prior_sample = Qprior2,
                                posterior = posterior, n_imp = n_imp,
                                density_func = qogiC_density_func)
post_prior2_Qmeas2 <- predict_Q(Qmeas_new = Qmeas_new, Qnew_prior_sample = Qprior2,
                                posterior = posterior, n_imp = n_imp,
                                density_func = qogiC_density_func)
