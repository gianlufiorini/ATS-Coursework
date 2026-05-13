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
  u = u_1,
  innovation   = innov_1
)

#IMPORTANT!!!! NOTE THAT THE INNOVATIONS HERE ARE *N O T* FROM THE KALMAN FILTER, THEY ARE THE RAW INNOVATIONS
#COMPUTED FROM THE GAS MODEL!!

### Plot the scaled score vs the innovation from score driven model FOR Y1 ####
ggplot(df, aes(x = t)) +
  geom_line(aes(y = u, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innovation, colour = "v"),
            linewidth = 0.7,
            linetype = "dashed") +
  labs(
    title    = "Scaled Score vs Innovation",
    subtitle = "t distribution with Fisher-inverse scaling",
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
df <- data.frame(t= 1:T_len, u = u_2, innovation  = innov_2)


### Plot the scaled score vs the innovation from score driven model for non-gaussian model ####

ggplot(df, aes(x = t)) +
  geom_line(aes(y = u, colour = "u"), linewidth = 0.7) +
  geom_line(aes(y = innovation, colour = "v"),
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

#DIAGNOSTICS ####
#Testing whiteness of score sequence and innovations

#Should be uncorrelated
acf(u_1)
acf(u_2)

##should be uncorrelated
acf(innov_1)
acf(innov_2)


#Should have kurtosis < 3
moments::kurtosis(u_1)
moments::kurtosis(u_2)

#Should be uncorrelated (reject H0)
Box.test(u_1, fitdf = 5, lag = 20)
Box.test(u_2, fitdf = 5, lag = 20)

Box.test(innov_1, fitdf = 5, lag = 20)
Box.test(innov_2, fitdf = 5, lag = 20)

### Plot of the scaled scores from the both GAS models ####
#install.packages("latex2exp", type = "binary")
library("tidyverse")
df_st_v <- data.frame(time = c(seq_along(mu_t1),
                               seq_along(mu_t2)),
                      u_t = c(u_1,
                              u_2),
                      v_t = c(innov_1,
                              innov_2),
                      name_ut = c(rep("u_1", 404),
                                  rep("u_2", 404)),
                      name_vt = c(rep("v_1", 404),
                                  rep("v_2", 404)))

#u_t is a winsorized verion of the innovation, so we expect the behavior we see in the histogram
#winsorized means a robust version, were the outlying values are downweighted, "censored" in a sense
library(latex2exp)
g1 <- ggplot(df_st_v, aes(x = time, y = u_t ,
                          color = name_ut)) +
  geom_line() +
  theme(legend.position = "none") +
  ylab(TeX("$u_t$")) +
  xlab("Time") +
  theme(axis.title.y = element_text(angle = 270),
        axis.text.y = element_text(angle = 270)) +
  coord_flip()

g2 <- ggplot(df_st_v, aes(x = u_t, fill = name_ut)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", 
                 alpha = 0.4, 
                 bins = 30) +
  theme(legend.position = "inside") +
  stat_function(inherit.aes = F, aes(x = u_t), fun = dt, color = "royalblue", linewidth = .3,
                args = list(df = gas_y2$fit$coef_est["df"]))+
  xlab("") +
  theme_minimal()


g2/g1 + plot_layout(heights = c(1, 4.5))


#Innovations should be t-distributed

g1 <- ggplot(df_st_v, aes(x = time, y = v_t ,
                          color = name_vt)) +
  geom_line() +
  theme(legend.position = "none") +
  ylab(TeX("$u_t$")) +
  xlab("Time") +
  theme(axis.title.y = element_text(angle = 270),
        axis.text.y = element_text(angle = 270)) +
  coord_flip()

g2 <- ggplot(df_st_v, aes(x = v_t, fill = name_vt)) +
  geom_histogram(aes(y = after_stat(density)), 
                 position = "identity", 
                 alpha = 0.4, 
                 bins = 30) +
  theme(legend.position = "inside") +
  stat_function(inherit.aes = F, aes(x = u_t), fun = dt, color = "turquoise", linewidth = .3,
                args = list(df = gas_y2$fit$coef_est["df"]))+
  stat_function(inherit.aes = F, aes(x = u_t), fun = dt, color = "coral", linewidth = .3,
                args = list(df = gas_y1$fit$coef_est["df"]))+
  xlab("") +
  theme_minimal()


g2/g1 + plot_layout(heights = c(1, 4.5))

#qqplots not that useful
# qqnorm(innov_1)
# qqnorm(innov_2)


# NEEDED FOR SECTION 5  (Ignore) ####
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
df <- data.frame(t= 1:T_len, u = u_1, innovation   = KF1$innov)

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
df <- data.frame(t= 1:T_len, u = u_2, innovation = KF2$innov)


# Plot
ggplot(df, aes(x = t)) +
  geom_line(aes(y = u, colour = "u"), linewidth = 0.7) +
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

#Comparing estimates states #GOOD
df_mupred <- data.frame(t = rep(1:T_len, 4),
                        y = c(y1,y2, y1, y2),
                        mu_gas = c(mu_t1, mu_t2, KF1$mu_pred, KF2$mu_pred),
                        name = factor(rep(c("y1","y2", "y1", "y2"), each = T_len)),
                        method = factor(c(rep("GAS Estimate", T_len * 2),
                                          rep("Kalman Filter Estimate", T_len * 2))))
ggplot(data = df_mupred,
       aes(x = t, colour = method, y = mu_gas)) +
  geom_line() +
  geom_line(aes(y = y), col = "black", alpha = .3, linewidth = 1) + 
  geom_point(size = .1, pch = 5) +
  facet_grid(rows = vars(name), scale = "free_y") +
  ylab("Estimated States") +
  xlab("Time") +
  scale_color_manual(values = c("Kalman Filter Estimate" = "red",
                                "GAS Estimate" = "green"),
                     position = "top", name = "Method") +
  guides(alpha = "none")


#Comparing state estimates TO CHECK IF CORRECT (DONT DELETE)
df_mupred <- data.frame(t = rep(1:T_len, 2),
                        mu_gas = c(mu_t1, mu_t2),
                        mu_kalman = c(KF1$mu_pred, KF2$mu_pred),
                        name = factor(c(rep("y1", T_len),
                                        rep("y2", T_len))))

