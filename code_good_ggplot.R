#loading libraries ####
library(ggplot2)
library(tseries)
library(carData)
library(car)
library(MASS)
library(tidyr)
library(ggpmisc)

# Loading Dataset ####
setwd("C:/Users/comet/Desktop/ATS/ATS-coursework")
load("C:/Users/comet/Desktop/ATS/ATS-coursework/fMRI-ROI-time-series.RData")
l <- 404

#gaussian
y1 <- Y[59,,16,1]
#non-gaussian
y2 <- Y[59,,13,1]

is.na(y1) |> any()
is.na(y2) |> any()

## Plotting ####
ts.plot(y2, col = "red")
lines(y1)

#This is if we want to represent them one at a time
df1 <- data.frame(
  time = seq_along(y1), 
  value = as.numeric(y1)
)

df2 <- data.frame(
  time = seq_along(y2), 
  value = as.numeric(y2)
)

ggplot(df1, aes(x = time, y = value)) +
  geom_line(color = "steelblue", size = 1) +
  #  geom_point(alpha = 0.5) + optional: adds points to see individual data entries
  theme_minimal() +
  labs(
    title = "Gaussian Time Series Plot",
    x = "Time Index",
    y = "Value"
  )

ggplot(df2, aes(x = time, y = value)) +
  geom_line(color = "steelblue", size = 1) +
#  geom_point(alpha = 0.5) + optional: adds points to see individual data entries
  theme_minimal() +
  labs(
    title = "Non-Gaussian Time Series Plot",
    x = "Time Index",
    y = "Value"
  ) 


#Plotting them together
#Create a data frame with both series
df <- data.frame(
  time = seq_along(y1),
  Gaussian = as.numeric(y1),
  Not_Gaussian = as.numeric(y2)
)

#Reshape from "wide" to "long"
df_long <- pivot_longer(df, cols = c("Gaussian", "Not_Gaussian"), 
                        names_to = "Variable", 
                        values_to = "Value")

#Plot
ggplot(df_long, aes(x = time, y = Value, color = Variable)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Overlapping Time Series", x = "Time", y = "Value")


## Loading Functions ####
source("ats_lab_2_functions.R")

# Estimating Kalman Filter ####
#Setting initial values
theta0 <- c(.5,2,2)

theta1_mle <- estimator(y1,theta0)
theta2_mle <- estimator(y2, theta0)

KF1 <- KF(y1, 
          theta1_mle$theta_list$hat_phi,
          theta1_mle$theta_list$hat_sigma_e,
          theta1_mle$theta_list$hat_sigma_eta)
KF1

KF2 <- KF(y2, 
          theta2_mle$theta_list$hat_phi,
          theta2_mle$theta_list$hat_sigma_e,
          theta2_mle$theta_list$hat_sigma_eta)
KF2

## Computing innovations ####
v1 <- y1 - KF1$mu_pred
v2 <- y2 - KF2$mu_pred

#Standardized innovations. If assumptions hold they should be IID N(0,1)
st_v1 <- v1 / KF1$Ft
st_v2 <- v2 / KF1$Ft

## Plotting ####
# Add the predictions to the data frame
df1$prediction <- as.numeric(KF1$mu_pred)
df2$prediction <- as.numeric(KF2$mu_pred)

par(mfrow = c(1,2))
ts.plot(y1, ylim = c(min(min(y1,y2)),
                     max(max(y1,y2))))
lines(KF1$mu_pred, col = "green", lwd = 2)


ts.plot(y2, ylim = c(min(min(y1,y2)),
                     max(max(y1,y2))))
lines(KF2$mu_pred, col = "green", lwd = 2)

#Now with ggplot
#Plot with the new layer
ggplot(df1, aes(x = time)) +
  geom_line(aes(y = value), color = "black", alpha = 0.9) + #original data
  geom_line(aes(y = prediction), color = "green", size = 0.8) + #KF1 prediction
  theme_minimal() +
  labs(
    title = "Actual Data vs. Kalman Filter Prediction",
    y = "Value",
    x = "Time"
  )

#Plot with the new layer
ggplot(df2, aes(x = time)) +
  geom_line(aes(y = value), color = "black", alpha = 0.9) + #original data
  geom_line(aes(y = prediction), color = "green",size = 0.8) + #KF1 prediction
  theme_minimal() +
  labs(
    title = "Actual Data vs. Kalman Filter Prediction",
    y = "Value",
    x = "Time"
  )

#One can play with alpha and size options in both geom_line commands to adjust 
#for opacity and size.


#Plot of standardized innovations
par(mfrow = c(1,1))
ts.plot(st_v1, ylim = c(min(min(st_v1,st_v2)),
                     max(max(st_v1,st_v2))))
lines(st_v2, ylim = c(min(min(st_v1,st_v2)),
                   max(max(st_v1,st_v2))), col = "red")
abline(h = 2, lty = 2, lwd = 1.5)
abline(h = -2, lty = 2, lwd = 1.5)

df_st <- data.frame(
  time = seq_along(st_v1),
  v1 = as.numeric(st_v1),
  v2 = as.numeric(st_v2)
)

#Reshape to long format for easy plotting
df_st <- pivot_longer(df_st, cols = c("v1", "v2"), 
                           names_to = "innovation_type", 
                           values_to = "value")

ggplot(df_st, aes(x = time, y = value)) +
  geom_line(color = "gray30") +
  facet_wrap(~innovation_type, ncol = 1) +  # Put one above the other
  geom_hline(yintercept = 0, color = "red") +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(title = "Standardized Innovations (Time Series)",
       subtitle = "Should look like white noise within blue dashed lines",
       y = "Standardized Value")

ggplot(df, aes(x = date, y = unemploy)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2007-09-15"),
             linetype = 2, color = 2, linewidth = 1)

KF1$Pt1_t
KF2$Pt1_t

alpha <- .05
alpha_bonf <- .05

#https://kalman-filter.com/normalized-innovation-squared/
sum(NIS1 < F_alpha_2 | NIS1 > F_1_alpha_2)
sum(NIS2 < F_alpha_2 | NIS2 > F_1_alpha_2)
403 *.05

KF1$llk
KF2$llk

acf(st_v1, lag = 400)
acf(st_v2, lag = 400)

Box.test(st_v1, lag = 20, fitdf = 3, type = "Ljung-Box")
Box.test(st_v2, lag = 20, fitdf = 3, type = "Ljung-Box")

shapiro.test(st_v1)
shapiro.test(st_v2)

jarque.bera.test(st_v1)
jarque.bera.test(st_v2)

ks.test(st_v1, y = "pnorm")
ks.test(st_v2, y = "pnorm")

qqPlot(st_v1)
qqPlot(st_v2)

