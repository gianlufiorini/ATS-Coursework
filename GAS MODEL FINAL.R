#GAS model
library(patchwork)

load("fMRI-ROI-time-series.RData")
#install.packages("pracma", type ="binary")
#install.packages("gasmodel", type ="binary")
library(gasmodel)

y1 <- Y[59,,16,1]
y2 <- Y[59,,13,1]
# distr(filter_type = "duration", filter_dim = "uni")
# distr()
ts.plot(y1, ylim = c(-15,15))


library(ggplot2)

# Fit the model
gas_y1 <- gas(y = y1, distr = "t", scaling = "fisher_inv_sqrt")
gas_y1

# Extract quantities
T_len     <- length(y1)
mu_t1      <- gas_y1$fit$par_tv[, 1]          # filtered mean
sigma_est <- (gas_y1$fit$coef_est["var"]) # static variance (if estimated)

# Scaled score = innovation for normal with fisher_inv
u_1 <- gas_y1$fit$score_tv[, 1]     # already scaled
innov_1   <- y1 - mu_t1                  # y_t - mu_t1 

# Combine into data frame
df <- data.frame(
  t            = 1:T_len,
  u_1 = u_1,
  innovation   = innov_1
)

### Plot the scaled score vs the innovation from score driven model ####
ggplot(df, aes(x = t)) +
  geom_line(aes(y = u_1, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innov_1, colour = "v"),
            linewidth = 0.7,
            linetype = "dashed") +
  labs(
    title    = "Scaled Score vs Innovation",
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
# For y2 (duration data)
gas_y2 <- gas(y = y2, distr = "t", scaling = "fisher_inv_sqrt")
gas_y2
# Extract quantities
T_len     <- length(y2)
mu_t2      <- gas_y2$fit$par_tv[, 1]          # filtered mean
sigma_est <- (gas_y2$fit$coef_est["var"]) # static variance (if estimated)

# Scaled score = innovation for normal with fisher_inv
u_2 <- gas_y2$fit$score_tv[, 1]     # already scaled
innov_2   <- y2 - mu_t2                 # y_t - mu_t2

# Combine into data frame
df <- data.frame(t= 1:T_len, u_1 = u_2, innovation   = innov_2)


### Plot the scaled score vs the innovation from score driven model for non-gaussian model ####

ggplot(df, aes(x = t)) +
  geom_line(aes(y = u_2, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innov_2, colour = "v"),
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
df_st_v <- data.frame(time = c(seq_along(mu_t1),
                               seq_along(mu_t2)),
                      v = c(u_1,
                            u_2),
                      name = c(rep("u_1", 404),
                               rep("u_2", 404))) |> 
  mutate(name = factor(name),
         outside = factor(abs(v) > 1.96, labels = c("no", "yes")))
library(latex2exp)
g1 <- ggplot(df_st_v, aes(x = time, y = v ,
                          color = name)) +
  geom_line() +
  theme(legend.position = "none") +
  ylab(TeX("$U_t$")) +
  xlab("Time") +
  theme(axis.title.y = element_text(angle = 270),
        axis.text.y = element_text(angle = 270)) +
  coord_flip()

g2 <- ggplot(df_st_v, aes(x = v, fill = name)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", 
                 alpha = 0.4, 
                 bins = 30) +
  theme(legend.position = "inside") +
  scale_fill_discrete(name = "", labels = (c("v1" = "std_innovations_y1",
                                             "v2" = "std_innovations_y2",
                                             "royalblue" = "N(0,1)"))) +
  stat_function(inherit.aes = F, aes(x = v), fun = dt, color = "royalblue", linewidth = .3,
                args = list(df = gas_y2$fit$coef_est["df"]))+
  xlab("") +
  theme_minimal()


g2/g1 + plot_layout(heights = c(1, 4.5))

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
df <- data.frame(t= 1:T_len, u_1 = u_1, innovation   = KF1$innov)

# Plot
ggplot(df, aes(x = t)) +
  geom_line(aes(y = u_1, colour = "u"), linewidth = 0.7) +
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
df <- data.frame(t= 1:T_len, u_1 = u_2, innovation   = KF2$innov)

# Plot
ggplot(df, aes(x = t)) +
  geom_line(aes(y = u_2, colour = "u"), linewidth = 0.7) +
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
