# ##### ATS COURSEWORK #### #
# Riccardo Visentini
# Nikola Stevanovic
# Gianluca Fiorini
# Tommaso Di Bene


#Loading Libraries  & data#####
rm(list=ls())

library(ggplot2)
library(gasmodel)
library(moments)
library(tseries)
library(knitr)
library(latex2exp)
library(kableExtra)
library(patchwork)
library(carData)
library(tidyverse)
library(car)
library(MASS)
library(tidyr)
library(ggpmisc)

#Loading the dataset
load("fMRI-ROI-time-series.RData")
l <- 404
#Selection of the 2 time series
y1 <- Y[59,,16,1] #gaussian
y2 <- Y[59,,13,1] #non-gaussian

#Plotting time series

bold_df <- data.frame(time = seq_along(y1), y1 = as.numeric(y1), y2 = as.numeric(y2))
  #Create a data frame with both series
df1 <- data.frame(time = seq_along(y1), value = as.numeric(y1))
df2 <- data.frame(time = seq_along(y2), value = as.numeric(y2))
bold_df_long <- pivot_longer(bold_df, cols = c("y1", "y2"), names_to = "Variable", 
                             values_to = "Value")

#plot of the first time series with local-mean and variance
#followed by summary statistics, graphs for ACF and PACF, and tests for normality
fit1 <- loess( y1 ~ time, data = data.frame(y1, time =  1:length(y1)), family = "symmetric",
               degree = 1)$fitted
df_sd_y1<-data.frame(time = seq_along(y1), fit = fit1, value=sqrt(loess(sqrt(devs) ~ time, 
          data = data.frame(devs = (y1 - fit1)^2, time = 1:length(y1)), degree = 1, 
          family = "symmetric")$fitted)) %>% mutate(ub1 = fit1 + 2*value, lb1 = fit1 - 2*value,
          sdev = fit1 + value)

g1<-ggplot(df1, aes(x = time, y = value)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_line(data=df_sd_y1, aes(x = time, y = fit1), color = "green", linewidth = 0.75)+
  geom_line(data=df_sd_y1, aes(x = time, y = sdev), color = "purple", linewidth = 0.75)+
  theme_minimal() + labs(x = "Time", y = "y1")
g1

xbar <- mean(y1) #close to 0
var <- var(y1) #1.69
kur <- kurtosis(y1) #close to 3
skew <- skewness(y1) #close to 0

shapiro.test(y1) 

par(mfrow = c(1,2), mar = c(5,4,3,0.1), cex.main = 1, cex.axis = .7)
acf(y1)
pacf(y1)

#plot of the second time series with local-mean and variance
#followed by summary statistics and the graphs for ACF and PACF
fit2 <- loess( y2 ~ time, data = data.frame(y2, time =  1:length(y2)),
               family = "symmetric", degree = 1)$fitted

df_sd_y2<-data.frame(time = seq_along(y2), fit = fit2,value=sqrt(loess(sqrt(devs) ~ time,
          data = data.frame(devs = (y2 - fit2)^2, time = 1:length(y2)),  family = "symmetric",
          degree = 1)$fitted)) %>% mutate(ub2 = fit2 + 2*value, lb2 = fit2 - 2*value,
          sdev = fit2 + value)

g2<-ggplot(df2, aes(x = time, y = value)) +
  geom_line(color = "black", linewidth = 0.75) +
  geom_line(data=df_sd_y2, aes(x = time, y = fit2), color = "green", linewidth = 0.75)+
  geom_line(data=df_sd_y2, aes(x = time, y = sdev), color = "purple", linewidth = 0.75) +
  theme_minimal() + labs(x = "Time Index", y = "y2")
g2

xbar <- mean(y2) #close to 0
var <- var(y2) #11.82
kur <- kurtosis(y2) #close to 3
skew <- skewness(y2) #close to 0

shapiro.test(y2)

par(mfrow = c(1,2), mar = c(5,4,3,0.1), cex.main = 1, cex.axis = .7)
acf(y2)
pacf(y2)

#qqplot for both series
par(cex.lab = .7, cex.axis = .7, cex.main = .7, mfrow = c(1,2), mar = c(3,2,2,1))
qqnorm(y1, xlab = "Normal Quantiles", ylab = "Empirical Quantiles", cex = .7,
       main = "qqplot for y1")
qqline(y1)
qqnorm(y2, xlab = "Normal Quantiles", ylab = "", cex = .7, main = "qqplot for y2")
qqline(y2)

# Time-varying Location model with Gaussian Noise ####
## Model Fitting ####
source("ats_lab_2_functions_with_mean.R") #loading functions for Kalman filter


### AR(1) signal + gaussian noise model for y1 & y2 ####

#Estimating Kalman Filter
#Setting initial values
theta0 <- c(0,.5,2,2) #initial values

theta1_mle <- estimator(y1,theta0) #MLEs for y1
mod_1_coef <- data.frame(Estimate = unlist(theta1_mle$theta_list),
                         Std.Err = theta1_mle$vcov |> diag() |> sqrt(),
                         row.names = c("omega", "phi", "sigma_epsilon", "sigma_eta")) |> 
  mutate(z_stat = Estimate/Std.Err) |> mutate(p_value = 2*(1-pnorm(abs(z_stat)))) #preparing output
mod_1_coef[3:4, 3:4] <- NA #Wald test not appropriate for variance parameters...

theta2_mle <- estimator(y2, theta0) #MlEs for y2
mod_2_coef <- data.frame(Estimate = unlist(theta2_mle$theta_list),
                         Std.Err = theta2_mle$vcov |> diag() |> sqrt(),
                         row.names = c("omega", "phi", "sigma_epsilon", "sigma_eta")) |> 
  mutate(z_stat = Estimate/Std.Err) |> mutate(p_value = 2*(1-pnorm(abs(z_stat)))) #preparing output
mod_2_coef[3:4, 3:4] <- NA #Wald test not appropriate for variance parameters...

#Kalman filter for y1
KF1 <- KF(y1, theta1_mle$theta_list$hat_omega, theta1_mle$theta_list$hat_phi,
          theta1_mle$theta_list$hat_sigma_e, theta1_mle$theta_list$hat_sigma_eta)

#computing log-lik and AIC for model1
ll_AIC_1 <- c(log.likelihood = KF1$llk-202*log(2*pi), AIC = -2*(KF1$llk-202*log(2*pi)) 
              + 2*nrow(mod_1_coef))

#Kalman filter for y2
KF2 <- KF(y2, theta2_mle$theta_list$hat_omega, theta2_mle$theta_list$hat_phi,
          theta2_mle$theta_list$hat_sigma_e, theta2_mle$theta_list$hat_sigma_eta)

#computing log-lik and AIC for model2
ll_AIC_2 <- c(log.likelihood = KF2$llk -202*log(2*pi) , AIC = -2*(KF2$llk-202*log(2*pi)) 
              + 2*nrow(mod_2_coef))

#Model output for AR(1) signal + Gaussian noise model for y1
print(mod_1_coef,digits = 2)
print(ll_AIC_1, digits = 3)

#Kalman Filter state estimates for y1
df1 <- data.frame(time = seq_along(y1), value = as.numeric(y1), name = "y1")
df2 <- data.frame(time = seq_along(y2), value = as.numeric(y2), name = "y2")

df1 <- rbind(df1, data.frame(time = seq_along(KF1$mu_pred), value = KF1$mu_pred,
        name = "mupred_1")) |> mutate(name = factor(name, levels = c("y1", "mupred_1")))
df2 <- rbind(df2, data.frame(time = seq_along(KF2$mu_pred), value = KF2$mu_pred,
        name = "mupred_2"))|> mutate(name = factor(name, levels = c("y2", "mupred_2")))

ggplot(df1, aes(x = time, y = value, color = name, linewidth = name)) +
  geom_line() + theme_minimal() + theme(legend.text = element_text(size = 10)) +
  labs(title = "Observations vs. Kalman Filter State Estimates", x = "Time") +
  scale_color_manual(values = c("mupred_1" = "green", "y1" = "grey3"),
  labels = c("mupred_1" = TeX("$\\mu_{t|t-1}$"), "y1" = TeX("$y_1$")), name = NULL,) +
  scale_linewidth_discrete(range = c(.5, .75)) + guides(linewidth = "none") 

#AR(1) signal + gaussian noise model for y2
#Model output for AR(1) signal + Gaussian noise model for y2
print(mod_2_coef, digits = 3)
print(ll_AIC_2, digits = 3)

#Kalman Filter state estimates for y2
ggplot(df2, aes(x = time, y = value, color = name, linewidth = name)) +
  geom_line() + theme_minimal() + theme(legend.text = element_text(size = 10)) +
  labs(title = "Observations vs. Kalman Filter State Estimates", x = "Time") +
  scale_color_manual(values = c("mupred_2" = "green", "y2" = "grey3"),
  labels = c("mupred_2" = TeX("$\\mu_{t|t-1}$"), "y2" = TeX("$y1$")), name = NULL,) +
  scale_linewidth_discrete(range = c(.5, .75)) + guides(linewidth = "none")

## DIAGNOSTICS ####
### Diagnostics for the AR(1) signal + gaussian noise model for y1 ####

#Computing innovations
v1 <- y1 - KF1$mu_pred
v2 <- y2 - KF2$mu_pred

#Standardized innovations. If assumptions hold they should be IID N(0,1)
st_v1 <- v1 / sqrt(KF1$Ft)
st_v2 <- v2 / sqrt(KF2$Ft)

df_st_v <- data.frame(time = c(seq_along(v1), seq_along(v2)), v = c(st_v1,st_v2),
            name = c(rep("v1", 404), rep("v2", 404))) |> mutate(name = factor(name),
            outside = factor(abs(v) > 1.96, labels = c("no", "yes")))

#Diagnostics for the AR(1) signal + gaussian noise model for y2
#Time series of standardized innovations and Histograms. The blue line in the 
#histogram is a standard normal density.
g1 <- ggplot(df_st_v, aes(x = time, y = v, color = name)) + geom_line() +
  theme(legend.position = "none") + ylab(TeX("$F_t^{-1/2}v_t$")) + xlab("Time") +
  theme(axis.title.y = element_text(angle = 270), axis.text.y = element_text(angle = 270)) +
  coord_flip() + theme_minimal() + scale_x_reverse(position = "top")

g2 <- ggplot(df_st_v, aes(x = v, fill = name)) + geom_histogram(aes(y = after_stat(density)), 
        position = "identity", alpha = 0.4, bins = 30) + theme(legend.position = "inside") +
        scale_fill_discrete(name = "", labels = (c("v1" = "std_innovations_y1",
        "v2" = "std_innovations_y2", "royalblue" = "N(0,1)"))) +
        stat_function(inherit.aes = F, aes(x = v), fun = dnorm, color = "royalblue", 
        linewidth = .3) + xlab("") + theme_minimal()

g2/g1+ plot_layout(heights = c(1, 4.5))

#Hypotheses tests and summary statistics
Box.test(st_v1, lag = 20, fitdf = 4, type = "Ljung-Box")
shapiro.test(st_v1)
jarque.bera.test(st_v1)

Box.test(st_v2, lag = 20, fitdf = 4, type = "Ljung-Box")
shapiro.test(st_v2)
jarque.bera.test(st_v2)

kurtosis(st_v1)
kurtosis(st_v2)

#ACF for standardized innovations for y1
acf(st_v1, lag= 40, main = NA)
#ACF for standardized innovations for y2
acf(st_v2, lag= 40, main = NA)

#QQplots for standardized innovations
par(cex.lab = .7, cex.axis = .7, cex.main = .7, mfrow = c(1,2), mar = c(3,2,2,1))
qqPlot(st_v1, xlab = "Normal Quantiles", ylab = "Empirical Quantiles", cex = .7,
       main = "qqplot for std. innovations of y1", id = F)
qqPlot(st_v2, xlab = "Normal Quantiles", ylab = "", cex = .7,
       main = "qqplot for std. innovations of y2", id = F)

## Time-Varying Location Model with Student-t noise ####
library(gasmodel)

### Signal + Student t noise for y1 ####
gas_y1 <- gas(y = y1, distr = "t", scaling = "fisher_inv_sqrt")

# Extract quantities
T_len <- length(y1)
mu_t1 <- gas_y1$fit$par_tv[, 1] # filtered mean
sigma_est <- (gas_y1$fit$coef_est["var"]) # static variance (if estimated)

with(gas_y1$fit, data.frame(Estimate = coef_est, Std.Err = coef_sd, 
    z_stat = coef_zstat, p_value = coef_pval)) |> print(digits = 3) #preparing output
data.frame(log.likelihood = sum(gas_y1$fit$loglik_tv), AIC = gas_y1$fit$aic,
           row.names = NULL) |> print(digits = 3) #computing loh-likelihood and AIC

#Scaled score = innovation for normal with fisher_inv
u_1 <- gas_y1$fit$score_tv[, 1] # already scaled
innov_1 <- y1 - mu_t1 # y_t - mu_t1 

#Combine into data frame
gf <- data.frame(t = rep(1:T_len, 2), y = c(y1, mu_t1),
      name = factor(rep(c("y1", "mu_t1"), each = T_len), levels = c("y1", "mu_t1")))

df <- data.frame(t = 1:T_len, u = u_1, innovation = innov_1)
#these innovations are not from the Kalman Filter, they are the raw innovations
#computed from the GAS model

g1 <- ggplot(gf,  aes(x = t, y = y, color = name, linewidth = name)) +
  geom_line() + theme_minimal() + theme(legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("mu_t1" = "green", "y1" = "grey3"),
  labels = c("mu_t1" = TeX("$\\mu_{t}$"), "y1" = TeX("$y_1$")), name = NULL,) +
  scale_linewidth_discrete(range = c(.3, .5)) + xlab("Time") + ylab("Value") + 
  scale_x_continuous(position = "top") + guides(linewidth = "none") 

#Plot the scaled score vs the innovation from score driven model FOR Y1
g2 <- ggplot(df, aes(x = t)) + geom_line(aes(y = u, colour = "u"), linewidth = 0.5) +
      geom_line(aes(y = innovation, colour = "v"), linewidth = 0.5, linetype = "dashed") +
      labs(x= "Time", y= "Value", colour   = NULL) + theme_minimal() +
      scale_colour_manual(values = c("u" = "steelblue", "v" = "firebrick"),
      labels = c("u" = expression(u[t]), "v" = expression(v[t])))

g1/g2

### Signal + Student t noise for y2 ####
gas_y2 <- gas(y = y2, distr = "t", scaling = "fisher_inv_sqrt")

#Extract quantities
T_len <- length(y2)
mu_t2 <- gas_y2$fit$par_tv[, 1] # filtered mean
sigma_est <- (gas_y2$fit$coef_est["var"]) # static variance (if estimated)

with(gas_y2$fit, data.frame(Estimate = coef_est, Std.Err = coef_sd, 
    z_stat = coef_zstat, p_value = coef_pval)) |> print(digits = 3)
data.frame(log.likelihood = sum(gas_y2$fit$loglik_tv), AIC = gas_y2$fit$aic,
           row.names = NULL) |> 
  print(digits = 3)

#Scaled score = innovation for normal with fisher_inv
u_2 <- gas_y2$fit$score_tv[, 1] # already scaled
innov_2 <- y2 - mu_t2 # y_t - mu_t2 

#Combine into data frame
gf <- data.frame(t = rep(1:T_len, 2), y = c(y2, mu_t2),
      name = factor(rep(c("y2", "mu_t2"), each = T_len), levels = c("y2", "mu_t2")))

df <- data.frame(t = 1:T_len, u = u_2, innovation = innov_2)
#these innovations are not from the Kalman Filter, they are the raw innovations
#computed from the GAS model

g1 <- ggplot(gf,  aes(x = t, y = y, color = name, linewidth = name)) +
  geom_line() + theme_minimal() + theme(legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("mu_t2" = "green", "y2" = "grey3"),
  labels = c("mu_t2" = TeX("$\\mu_{t}$"), "y2" = TeX("$y_2$")), name = NULL,) +
  scale_linewidth_discrete(range = c(.3, .5)) + xlab("Time") + ylab("Value") + 
  scale_x_continuous(position = "top") + guides(linewidth = "none") 

#Plot the scaled score vs the innovation from score driven model for y1
g2 <- ggplot(df, aes(x = t)) + geom_line(aes(y = u, colour = "u"), linewidth = 0.5) +
  geom_line(aes(y = innovation, colour = "v"), linewidth = 0.5, linetype = "dashed") +
  labs(x= "Time", y = "Value", colour = NULL) + theme_minimal() +
  scale_colour_manual(values = c("u" = "steelblue", "v" = "firebrick"),
  labels = c("u" = expression(u[t]), "v" = expression(v[t])))

g1/g2

## Diagnostics for the signal + Student t noise model for both y1 and y2 ####

df_st_v <- data.frame(time = c(seq_along(mu_t1), seq_along(mu_t2)),
            u_t = c(u_1, u_2), v_t = c(innov_1, innov_2),
            name_ut = c(rep("u_1", 404), rep("u_2", 404)),
            name_vt = c(rep("v_1", 404), rep("v_2", 404)))

g1 <- ggplot(df_st_v, aes(x = time, y = u_t , color = name_ut)) + geom_line() +
  theme(legend.position = "none") + ylab(TeX("$u_t$")) + xlab("Time") +
  theme(axis.title.y = element_text(angle = 270), axis.text.y = element_text(angle = 270)) +
  coord_flip() + theme_minimal()+ guides(color = "none") +
  scale_x_reverse(position = "bottom")

g2 <- ggplot(df_st_v, aes(x = u_t, fill = name_ut)) + 
    geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.4, 
    bins = 30) + theme(legend.position = "none",
    legend.title = element_text("Scaled Score")) + guides(fill = "none") +
    xlab("") + theme_minimal() + theme_minimal() + guides(color = "none")

#Innovations should be t-distributed
dt_locscale <- function(x, df, mu = 0, sigma = 1,...) {
  dt((x - mu)/sigma, df = df,...)/sigma
}

g3 <- ggplot(df_st_v, aes(x = time, y = v_t , color = name_vt)) + geom_line() +
  theme(legend.position = "none") + ylab(TeX("$v_t$")) + xlab("Time") +
  theme(axis.title.y = element_text(angle = 270), axis.text.y = element_text(angle = 270)) +
  coord_flip() + theme_minimal() + guides(color = "none") +
  scale_x_reverse(position = "top")

g4 <- ggplot(df_st_v, aes(x = v_t, fill = name_vt)) +
      geom_histogram(aes(y = after_stat(density)), position = "identity", alpha = 0.4, 
      bins = 30) + theme(legend.position = "right") +
      stat_function(inherit.aes = F, aes(x = v_t), fun = dt_locscale, color = "turquoise",
      linewidth = .3, args = list(df = gas_y2$fit$coef_est["df"],
      sigma = sqrt(gas_y2$fit$coef_est["var"])))+
      stat_function(inherit.aes = F, aes(x = v_t), fun = dt_locscale, color = "coral", 
      linewidth = .3, args = list(df = gas_y1$fit$coef_est["df"],
      sigma = sqrt(gas_y1$fit$coef_est["var"])))+ labs(fill = "Time Series") +
      xlab("") + theme_minimal() + scale_fill_discrete(labels = c("v_1" = "y_1", "v_2" = "y_2"))

(g2/g1)+ plot_layout(heights = c(1,4.5))|(g4/g3) + plot_layout(heights = c(1,4.5))

#Hypotheses tests and summary statistics
Box.test(u_1, fitdf = 5, lag = 20)
Box.test(u_2, fitdf = 5, lag = 20)

Box.test(innov_1, fitdf = 5, lag = 20)
Box.test(innov_2, fitdf = 5, lag = 20)

moments::kurtosis(u_1)
moments::kurtosis(u_2)

#ACF for Scaled Score and Raw Innovations for y1
par(mfrow = c(1,2), mar = c(5,4,3,0.1), cex.main = 1, cex.axis = .7)
acf(u_1, main = "Scaled Score for y1")
acf(innov_1, main = "Raw Innovations for y1")

#ACF for Scaled Score and Raw Innovations for y2
par(mfrow = c(1,2), mar = c(5,4,3,0.1), cex.main = 1, cex.axis = .7)
acf(u_2, main = "Scaled Score for y2")
acf(innov_2, main = "Raw Innovations for y2")

#Comparing estimates states
df_mupred <- data.frame(t = rep(1:T_len, 4), y = c(y1,y2, y1, y2),
              mu_gas = c(mu_t1, mu_t2, KF1$mu_pred, KF2$mu_pred),
              name = factor(rep(c("y1","y2", "y1", "y2"), each = T_len)),
              method = factor(c(rep("GAS Estimate", T_len * 2),
                                          rep("Kalman Filter Estimate", T_len * 2))))
ggplot(data = df_mupred, aes(x = t, colour = method, y = mu_gas)) + geom_line() +
  geom_line(aes(y = y), col = "black", alpha = .3, linewidth = 1) + 
  geom_point(size = .1, pch = 5) + coord_flip() + facet_grid(cols = vars(name), scale = "free_x") +
  ylab("Estimated States") + xlab("Time") +
  scale_color_manual(values = c("Kalman Filter Estimate" = "red", "GAS Estimate" = "green"),
  position = "top", name = "Method", labels = c("Kalman Filter Estimate" = "KF",
  "GAS Estimate" = "GAS")) + scale_x_reverse(position = "top") +
  guides(alpha = "none") + theme_minimal() + theme(axis.title.y = element_text(angle = 270)) 



