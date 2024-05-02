

# JAGS Functions ----------------------------------------------------------------

#' Creates jags code string
#' @param logphi_str String of JAGS code that describes the specification of log(phi_i). indexing of data/parameters using "[i]"
#' @param sigmasq_inv_str String of JAGS code that describes the specification of the INVERSE of the variance parameter for observation i. indexing of data/parameters using "[i]"
#' @param priors_str String of JAGS code that specifies priors for all parameters in the model. See the SI and JAGS user manual for help.
#' @param logscale Should the likelihood be specified on the log of the measurements or on the raw measurements? Default is T
#' @param t Should the t distribution with 3 degrees of freedom be used instead of the normal distribution for the likelihood? Default is F. Can be used for sensitivity analysis.
make_model <- function(logphi_str,
                       sigmasq_inv_str,
                       priors_str,
                       logscale = T,
                       t = F) {

  if(logscale & t) {

    jags_str <- paste0("data {

  for(i in 1:nobs) {

    logQmeas[i] <- log(Qmeas[i])

  }

}

model {

  # Likelihood

  for(i in 1:nobs) {

    mean[i] <-", logphi_str, "

    logQmeas[i] ~ dt(mean[i],", sigmasq_inv_str,", 3)

  }

  # Priors

  ", priors_str, "

}")

  } else if(logscale & !t) {

    jags_str <- paste0("data {

  for(i in 1:nobs) {

    logQmeas[i] <- log(Qmeas[i])

  }

}

model {

  # Likelihood

  for(i in 1:nobs) {

    mean[i] <-", logphi_str, "

    logQmeas[i] ~ dnorm(mean[i],", sigmasq_inv_str,")

  }

  # Priors

  ", priors_str, "

}")

  } else {

    jags_str <- paste0("
model {

  # Likelihood

  for(i in 1:nobs) {

    mean[i] <-", normalmean_str, "

    Qmeas[i] ~ dlnorm(mean[i], ", sigmasq_inv_str, ")

  }

  # Priors

  tau <- 1/sigma^2

  ", priors_str, "

}")

  }

  return(jags_str)

}

#' selects data from only the desired technology, removes rows where the measured emission rate is zero
#' and returns a jags_data object
#' @param all_data data.frame of all data which must include columns `technology`, `actual`, and `estimate`
#' @param technology_name string of the technology name which must match the name in all_data
#' @return A list of vectors Qmeas, Q, and nobs
data_prep <- function(all_data, technology_name = "all", Qmeas_name = "estimate_kgh", Q_name = "actual_kgh") {

  if(technology_name == "all") {

    technology <- all_data

  } else {

    technology <- subset(all_data, technology == technology_name)

  }

  df <- setNames(data.frame(Qmeas = c(technology[,Qmeas_name]),
                            Q = c(technology[,Q_name])), c("Qmeas", "Q")) %>% subset(!is.na(Qmeas) & Qmeas > 0)

  return(list(Qmeas = df$Qmeas,
              Q = df$Q,
              nobs = nrow(df)))

}

# Plotting functions ---------------------------------------------------------------------------------

#' Create scatterplots with prediction bands, as shown in the publication
#' @param data The data used to fit the model. Should be a list with named elements Q and Qmeas, like from the helper function data_prep
#' @param newdata Optional; external data to compare to prediction bands. Same format as argument data
#' @param samples_obj a run.jags object which contains the MCMC samples from the posterior of Theta
#' @param Qtrue_vals Vector of x values to use when generating prediction bands. The longer the vector the longer it takes to run, but the smoother
#'                   the prediction bands will be
#' @param plot_Qmeas_val_func The name of an R function which draws from the likelihood of a measurement for fixed Theta and true emission rate. Arguments are Q_val and posterior
#' @param title_str A string to specify a custom title for the plot
#' @param hdi logical; Should highest density CIs be used? Default is F, and quantile-based CIs are calculated
#' @param fixed_yrange Either NA, where the range of the Y axis is determined automatically. Otherwise, vector of length 2 giving lower and upper limits for the Y axis, to zoom in or out of the plot
scatter_gg <- function(data,
                       newdata = NA,
                       samples_obj = NA,
                       Qtrue_vals,
                       plot_Qmeas_val_func,
                       title_str = "",
                       hdi = F,
                       fixed_yrange = NA) {


  any_new <- (length(newdata) == 3)



  if(any_new) {

    data_plot <- data.frame(Q = data$Q,
                            Qmeas = data$Qmeas,
                            data_type = rep("Used to fit model (second trial)"))

    data_plot <- rbind(data_plot,
                       data.frame(Q = newdata$Q,
                                  Qmeas = newdata$Qmeas,
                                  data_type = rep("External (first trial)")))

    # xlims <- range(c(range(data$Q), range(newdata$Q)))
    # ylims <- range(c(range(data$Qmeas), range(newdata$Qmeas)))

  } else {

    data_plot <- data.frame(Q = data$Q,
                            Qmeas = data$Qmeas,
                            data_type = rep("Used to fit model"))

  }

  if(length(samples_obj)> 1) {

    Qmeas_posterior <- generate_predictions(Qtrue_vals = Qtrue_vals,
                                            samples_obj = samples_obj,
                                            Qmeas_func = plot_Qmeas_val_func)

  } else {

    stop("samples_obj must be specified")

  }

  colnames(Qmeas_posterior) <- paste0("Q = ", Qtrue_vals) # each column corresponds to one Qtrue value

  Qmeas_posterior_long <- data.frame(Qmeas_posterior) %>% pivot_longer(everything()) %>% group_by(name) %>%
    reframe(point_est = median(value),
            quantile_df(value, c(0.50, 0.90, 0.95), hdi = hdi)) %>% pivot_wider(names_from = ci,
                                                                                values_from = val) %>%
    pivot_longer(cols = starts_with(c("upper", "lower")),
                 names_to = c("bound", "conf_level"),
                 names_sep = "0",
                 values_to = "value",
                 names_transform = list(conf_level = as.numeric)) %>%
    mutate(conf_level = paste0(conf_level*100, "%"),
           Q = as.numeric(str_remove(name,"Q..."))) %>% select(!name) %>%
    pivot_wider(names_from = bound,
                values_from = value) %>% mutate(data_type = NA, Qmeas = 1, line_dummy = "Posterior Median Prediction")

  line_df <- data.frame(slope = c(1, NA),
                        intercept = c(0,NA),
                        linetype = c(2, 1),
                        name = c("Perfect prediction line", "Posterior median prediction"))

  if(length(fixed_yrange) > 1) {

    g <- ggplot(data_plot, aes(x = Q, y = Qmeas, shape = data_type, fill = data_type)) +
      geom_ribbon(data = subset(Qmeas_posterior_long, conf_level == "50%"),
                  aes(x = Q, ymin = lower, ymax = upper, alpha = conf_level, col = conf_level), alpha = 0.4, linetype = 0, fill = "cornflowerblue") +
      geom_ribbon(data = subset(Qmeas_posterior_long, conf_level == "90%"),
                  aes(x = Q, ymin = lower, ymax = upper, alpha = conf_level, col = conf_level), alpha = 0.3, linetype = 0, fill = "cornflowerblue") +
      geom_ribbon(data = subset(Qmeas_posterior_long, conf_level == "95%", col = conf_level),
                  aes(x = Q, ymin = lower, ymax = upper, alpha = conf_level, col = conf_level), alpha = 0.3, linetype = 0, fill = "cornflowerblue") +
      scale_color_manual(values = rep("cornflowerblue", 3), breaks = c("95%", "90%", "50%")) +
      scale_alpha_identity() +
      guides(colour = guide_legend(reverse = T, override.aes = list(alpha = c(1/3, 0.6/3, 0.3/3)),
                                   order = 3,
                                   title = ifelse(hdi,
                                                  "Prediction band\n(HPD intervals)",
                                                  "Prediction band\n(Quantile-based intervals)"))) +
      geom_point(alpha = 1, col = "black") +
      geom_line(data = Qmeas_posterior_long, mapping = aes(x = Q, y = point_est), alpha = 1, col = "black", lwd = 1) +
      scale_shape_manual(values = c(24, 21, 24, 21),
                         breaks = c("External (first trial)", "Used to fit model (second trial)", "External", "Used to fit model")) +
      scale_fill_manual(values = c("firebrick", "goldenrod", "firebrick", "goldenrod"),
                        breaks =c("External (first trial)", "Used to fit model (second trial)", "External", "Used to fit model"))+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      geom_abline(data = line_df, aes(slope = slope, intercept = intercept, linetype = name, lwd = name)) +
      scale_linewidth_manual(values = c(1, 0.7), breaks = c("Posterior median prediction", "Perfect prediction line")) +
      scale_linetype_manual(values = c(1, 2), breaks = c("Posterior median prediction", "Perfect prediction line")) +
      guides(fill = guide_legend(name = "Data type", title = "Data type", order = 1),
             shape = guide_legend(name = "Data type", title = "Data type", order = 1),
             linetype = guide_legend(name = "name", title = "Prediction", order = 2),
             lwd = guide_legend(name = "name", title = "Prediction", order = 2)) +
      labs(title = title_str, x = "True Emission Rate (kg/h)", y = "Measured Emission Rate (kg/h)") +
      ylim(fixed_yrange) + #uncomment to zoom in on plot
      theme_bw()

  } else {

    g <- ggplot(data_plot, aes(x = Q, y = Qmeas, shape = data_type, fill = data_type)) +
      geom_ribbon(data = subset(Qmeas_posterior_long, conf_level == "50%"),
                  aes(x = Q, ymin = lower, ymax = upper, alpha = conf_level, col = conf_level), alpha = 0.4, linetype = 0, fill = "cornflowerblue") +
      geom_ribbon(data = subset(Qmeas_posterior_long, conf_level == "90%"),
                  aes(x = Q, ymin = lower, ymax = upper, alpha = conf_level, col = conf_level), alpha = 0.3, linetype = 0, fill = "cornflowerblue") +
      geom_ribbon(data = subset(Qmeas_posterior_long, conf_level == "95%", col = conf_level),
                  aes(x = Q, ymin = lower, ymax = upper, alpha = conf_level, col = conf_level), alpha = 0.3, linetype = 0, fill = "cornflowerblue") +
      scale_color_manual(values = rep("cornflowerblue", 3), breaks = c("95%", "90%", "50%")) +
      scale_alpha_identity() +
      guides(colour = guide_legend(reverse = T, override.aes = list(alpha = c(1/3, 0.6/3, 0.3/3)),
                                   order = 3,
                                   title = ifelse(hdi,
                                                  "Prediction band\n(HPD intervals)",
                                                  "Prediction band\n(Quantile-based intervals)"))) +
      geom_point(alpha = 1, col = "black") +
      geom_line(data = Qmeas_posterior_long, mapping = aes(x = Q, y = point_est), alpha = 1, col = "black", lwd = 1) +
      scale_shape_manual(values = c(24, 21, 24, 21),
                         breaks = c("External (first trial)", "Used to fit model (second trial)", "External", "Used to fit model")) +
      scale_fill_manual(values = c("firebrick", "goldenrod", "firebrick", "goldenrod"),
                        breaks =c("External (first trial)", "Used to fit model (second trial)", "External", "Used to fit model"))+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0)+
      geom_abline(data = line_df, aes(slope = slope, intercept = intercept, linetype = name, lwd = name)) +
      scale_linewidth_manual(values = c(1, 0.7), breaks = c("Posterior median prediction", "Perfect prediction line")) +
      scale_linetype_manual(values = c(1, 2), breaks = c("Posterior median prediction", "Perfect prediction line")) +
      guides(fill = guide_legend(name = "Data type", title = "Data type", order = 1),
             shape = guide_legend(name = "Data type", title = "Data type", order = 1),
             linetype = guide_legend(name = "name", title = "Prediction", order = 2),
             lwd = guide_legend(name = "name", title = "Prediction", order = 2)) +
      labs(title = title_str, x = "True Emission Rate (kg/h)", y = "Measured Emission Rate (kg/h)") +
      #xlim(c(0,15)) #uncomment to zoom in on plot
      theme_bw()

  }

  #return(Qmeas_posterior_long)
  return(g)

}

#' Helper function for scatter_gg
#' @param Qtrue_vals a vector of Qtrue_vals that you want to generature posterior predictions for
#' @param samples_obj a run.jags object contaitning posterior samples from model
#' @param Qmeas_func name of function which will draw predictions for Qmeas given a value of Q and a posterior of named parameters. Arguments are Q_val and posterior.
#' @return A matrix of dimensions (number of iterations in posterior) by (number of Qtrue values)

generate_predictions <- function(Qtrue_vals, samples_obj, Qmeas_func) {

  posterior <- combine.mcmc(samples_obj$mcmc)

  Qmeas_plot_post <- matrix(nrow = nrow(posterior), ncol = length(Qtrue_vals))

  for(i in 1:length(Qtrue_vals)) {

    Qmeas_plot_post[,i] <- Qmeas_func(Q_val = Qtrue_vals[i], posterior = posterior)

  }

  return(Qmeas_plot_post)

}

#' Helper function for scatter_gg
quantile_df <- function(x, prediction_levels, hdi = F) {

  if(hdi) {

    df <- data.frame(lower = rep(NA, length(prediction_levels)),
                     upper = rep(NA, length(prediction_levels)),
                     pred_level = prediction_levels)

    for(i in 1:length(prediction_levels)) {

      df[i,] <- c(hdi(x, prediction_levels[i]), prediction_levels[i])

    }

    pivot_longer(df, cols = c("lower", "upper"), names_to = "ci", values_to = "val") %>%
      transmute(val = val, ci = paste0(ci, pred_level))

  } else {

    alphas <- 0.5*(1-prediction_levels)

    probs <- c(alphas, 1-alphas)

    tibble(
      val = quantile(x, probs = probs),
      ci = paste0(c(rep("lower", length(prediction_levels)), rep("upper", length(prediction_levels))), rep(prediction_levels, 2))
    )

  }


}



# function which will draw predictions for Qmeas given a value of Q (Q_val) and a posterior of named parameters, for the QOGI C preferred model
quad_linear <- function(Q_val, posterior) {

  meani <- log(posterior[,"alpha0"] +
                 posterior[,"alpha1"]*Q_val +
                 posterior[,"alpha2"]*Q_val^2 +
                 (Q_val>posterior[,"gamma"])*
                 (posterior[,"alpha2"]*(posterior[,"gamma"]^2 - Q_val^2) +
                    posterior[,"beta1"]*(Q_val-posterior[,"gamma"])))
  epsi <- rnorm(nrow(posterior),
                mean = 0,
                sd = 1/sqrt(posterior[,"tau_base"] + Q_val/posterior[,"eta"]))
  exp(meani+epsi)
}

# function which will draw predictions for Qmeas given a value of Q (Q_val) and a posterior of named parameters, for the NIR HSI preferred model
simple_linear_func <- function(Q_val, posterior) {

  exp(log(posterior[,"alpha0"] + posterior[,"alpha1"]*Q_val) + rnorm(nrow(posterior), 0, posterior[,"sigma"]))

}


# Weighted bootstrap functions -----------------------------------------------
#' Get a sample from the posterior P(Q_new | Y_new, field data)
#' @param Qmeas_new vector of values of new observations
#' @param Qnew_prior_sample Vector of Qnew values drawn from the prior distribution,
#'        to be used as population to resample from in importance sampling
#' @param posterior matrix/dataframe of posterior samples to use in calculating conditional expectation
#' @param n_imp number of samples to take when doing importance sampling
#' @return the sample from the posterior
predict_Q <- function(Qmeas_new, Qnew_prior_sample, posterior, n_imp, density_func) {

  Y_new <- log(Qmeas_new)

  # Calculate conditional expectation for every sample from prior Qnew_prior_sample - this is the bottleneck
  means <- lapply(as.list(Qnew_prior_sample), e_theta, Y_new = Y_new, posterior = posterior, density_func = density_func)
  weights <- unlist(means)/sum(unlist(means)) # scale so weights sum to 1

  # Do important sampling
  imp_samp <- sample(Qnew_prior_sample, size = n_imp, replace = T, prob = weights)

  return(imp_samp)

}

#' Find the conditional expectation for one value of Qnew
#' @param Q_new is a numeric value of Q_new
#' @param Y_new is numeric value of Y_new (log(Qmeas_new))
#' @param posterior samples from the posterior distribution
#' @param density_func
e_theta <- function(Q_new, Y_new, posterior, density_func) {

  npts <- length(Y_new)
  ps <- matrix(nrow = nrow(posterior), ncol = npts)
  for(i in 1:npts){

    ps[,i] <- density_func(y = Y_new[i], q = Q_new, posterior = posterior)
  }

  p_Ynews <- apply(ps, 1, prod)

  return(mean(p_Ynews))

}

qogiC_density_func <- function(y, q, posterior) {

  small <- if_else(q > posterior[,"gamma"], 0, 1)

  # return the normal density
  dnorm(x=y,
        mean = log(small*(posterior[,"alpha0"] + posterior[,"alpha1"]*q + posterior[,"alpha2"]*q^2) +
                     (1-small)*(posterior[,"alpha0"] + posterior[,"beta0"] + (posterior[,"alpha1"] + posterior[,"beta1"])*q)),
        sd = sqrt(1/(posterior[,"tau_base"] + q/posterior[,"eta"])))

}
