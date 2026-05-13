#GAS model
#Filter sequance 
#The scores
#Plotting the scores (scaled)  against the innovation of kalman filter
load("fMRI-ROI-time-series.RData")
#install.packages("pracma", type ="binary")
#install.packages("gasmodel", type ="binary")
library(gasmodel)

ts.gaus <- Y[59,,16,1]
ts.nong <- Y[59,,13,1]
distr(filter_type = "duration", filter_dim = "uni")
distr()
ts.plot(ts.gaus, ylim = c(-15,15))


library(ggplot2)

# Fit the model
est_norm <- gas(y = ts.gaus, distr = "t", scaling = "fisher_inv_sqrt")
est_norm

# Extract quantities
T_len     <- length(ts.gaus)
mu_t      <- est_norm$fit$par_tv[, 1]          # filtered mean
sigma_est <- exp(est_norm$fit$coef_est["log(var)_omega"]) # static variance (if estimated)

# Scaled score = innovation for normal with fisher_inv
scaled_score <- est_norm$fit$score_tv[, 1]     # already scaled
innovation   <- ts.gaus - mu_t                  # y_t - mu_t (Kalman innovation)

# Combine into data frame
df <- data.frame(
  t            = 1:T_len,
  scaled_score = scaled_score,
  innovation   = innovation
)

### Plot the scaled score vs the innovation from score driven model ####
ggplot(df, aes(x = t)) +
  geom_line(aes(y = scaled_score, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innovation, colour = "v"),
            linewidth = 0.7,
            linetype = "dashed") +
  labs(
    title    = "Scaled GAS Score vs Innovation",
    subtitle = "Normal distribution with Fisher-inverse scaling",
    x        = "Time",
    y        = "Value",
    colour   = NULL
  ) +
  theme_minimal() +
  scale_colour_manual(
    values = c("u" = "steelblue",
               "v" = "firebrick"),
    labels = c(
      "u" = expression(u[t]),
      "v" = expression(v[t])
    )
  )


# Fit the model
# For ts.nong (duration data)
est_dur <- gas(y = ts.nong, distr = "t", scaling = "fisher_inv_sqrt")
est_dur
# Extract quantities
T_len     <- length(ts.nong)
mu_tt      <- est_dur$fit$par_tv[, 1]          # filtered mean
sigma_est <- exp(est_dur$fit$coef_est["log(var)_omega"]) # static variance (if estimated)

# Scaled score = innovation for normal with fisher_inv
scaled_score_t <- est_dur$fit$score_tv[, 1]     # already scaled
innovation_t   <- ts.nong - mu_t                  # y_t - mu_t (Kalman innovation)

# Combine into data frame
df <- data.frame(t= 1:T_len, scaled_score = scaled_score_t, innovation   = innovation_t)


### Plot the scaled score vs the innovation from score driven model for non-gaussian model ####

ggplot(df, aes(x = t)) +
  geom_line(aes(y = scaled_score_t, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innovation_t, colour = "v"),
            linewidth = 0.7,
            linetype = "dashed") +
  labs(
    title    = "Scaled GAS Score vs Innovation",
    subtitle = "Student-t distribution with Fisher-inverse scaling",
    x        = "Time",
    y        = "Value",
    colour   = NULL
  ) +
  theme_minimal() +
  scale_colour_manual(
    values = c("u" = "steelblue",
               "v" = "firebrick"),
    labels = c(
      "u" = expression(u[t]),
      "v" = expression(v[t])
    )
  )



### Plot of the scaled scores from the both GAS models ####
#install.packages("latex2exp", type = "binary")
library("tidyverse")
df_st_v <- data.frame(time = c(seq_along(mu_t),
                               seq_along(mu_tt)),
                      v = c(scaled_score,
                            scaled_score_t),
                      name = c(rep("mu_t", 404),
                               rep("mu_tt", 404))) |> 
  mutate(name = factor(name),
         outside = factor(abs(v) > 1.96, labels = c("no", "yes")))
library(latex2exp)
ggplot(df_st_v, aes(x = time, y = v ,
                    color = name)) +
  geom_line() +
  geom_hline(yintercept = c(-1.96,1.96),
             linetype = "dashed",
             col = "blue") +
  theme(legend.position = "none") +
  ylab(TeX("$U_t$")) +
  xlab("Time") +
  theme(axis.title.y = element_text(angle = 270),
        axis.text.y = element_text(angle = 270)) +
  coord_flip()

### Plot of the scaled score vs the innovation from Kalman filter ####
# Loading Dataset #
setwd("C:/Users/comet/Desktop/ATS/ATS-coursework")
load("C:/Users/comet/Desktop/ATS/ATS-coursework/fMRI-ROI-time-series.RData")
l <- 404

#gaussian
y1 <- Y[59,,16,1]
#non-gaussian
y2 <- Y[59,,13,1]

is.na(y1) |> any()
is.na(y2) |> any()

## Plotting ##

ts.plot(y2, col = "red")
lines(y1)

## Loading Functions ###
source("ats_lab_2_functions_with_mean.R")

# Estimating Kalman Filter #
#Setting initial values
theta0 <- c(0,.5,1,1)

theta1_mle <- estimator(y1,theta0)
mod_1_coef <- data.frame(Estimate = unlist(theta1_mle$theta_list),
                         Std.Err = theta1_mle$vcov |> diag() |> sqrt(),
                         row.names = c("omega", "phi", "sigma_epsilon", "sigma_eta")) |> 
  mutate(z_stat = Estimate/Std.Err) |> 
  mutate(p_value = 2*(1-pnorm(abs(z_stat))))
mod_1_coef[3:4, 3:4] <- NA

theta2_mle <- estimator(y2, theta0)

mod_2_coef <- data.frame(Estimate = unlist(theta2_mle$theta_list),
                         Std.Err = theta2_mle$vcov |> diag() |> sqrt(),
                         row.names = c("omega", "phi", "sigma_epsilon", "sigma_eta")) |> 
  mutate(z_stat = Estimate/Std.Err) |> 
  mutate(p_value = 2*(1-pnorm(abs(z_stat))))
mod_2_coef[3:4, 3:4] <- NA

KF1 <- KF(y1, 
          theta1_mle$theta_list$hat_omega,
          theta1_mle$theta_list$hat_phi,
          theta1_mle$theta_list$hat_sigma_e,
          theta1_mle$theta_list$hat_sigma_eta)
KF2 <- KF(y2, 
          theta2_mle$theta_list$hat_omega,
          theta2_mle$theta_list$hat_phi,
          theta2_mle$theta_list$hat_sigma_e,
          theta2_mle$theta_list$hat_sigma_eta)
# Combine into data frame
df <- data.frame(t= 1:T_len, scaled_score = scaled_score, innovation   = KF1$innov)

# Plot
ggplot(df, aes(x = t)) +
  geom_line(aes(y = scaled_score, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innovation, colour = "v"),
            linewidth = 0.7,
            linetype = "dashed") +
  labs(
    title    = "Scaled GAS Score vs Kalman Filter Innovation",
    subtitle = "Gaussian distribution with Fisher-inverse scaling",
    x        = "Time",
    y        = "Value",
    colour   = NULL
  ) +
  theme_minimal() +
  scale_colour_manual(
    values = c("u" = "steelblue",
               "v" = "firebrick"),
    labels = c(
      "u" = expression(u[t]),
      "v" = expression(v[t])
    )
  )


### Plot of the scaled score vs the innovation from kalman filter for non-gaussian model ####
df <- data.frame(t= 1:T_len, scaled_score = scaled_score_t, innovation   = KF2$innov)

# Plot
ggplot(df, aes(x = t)) +
  geom_line(aes(y = scaled_score_t, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innovation, colour = "v"),
            linewidth = 0.7,
            linetype = "dashed") +
  labs(
    title    = "Scaled GAS Score vs Kalman Filter Innovation",
    subtitle = "Student-t distribution with Fisher-inverse scaling",
    x        = "Time",
    y        = "Value",
    colour   = NULL
  ) +
  theme_minimal() +
  scale_colour_manual(
    values = c("u" = "steelblue",
               "v" = "firebrick"),
    labels = c(
      "u" = expression(u[t]),
      "v" = expression(v[t])
    )
  )
