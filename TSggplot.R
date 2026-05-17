load("C:/Users/tommy/OneDrive/Desktop/Advanced Time Series/fMRI-ROI-time-series.RData")

library(ggplot2)

# Sample dataset
set.seed(123)
df <- data.frame(
  Date = seq(as.Date("2023-01-01"), by = "month", length.out = 12),
  Value = round(runif(12, min = 50, max = 100), 2)
)


# Time series plot
ggplot(df, aes(x = Date, y = Value)) +
  geom_line(color = "steelblue") +
  geom_point(color = "darkred") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
#  geom_area(fill = "gray", alpha = 0.9) + 
#  geom_smooth(se = F) +
  theme_minimal()



library(lubridate)

# Data
df <- economics[economics$date > as.Date("2000-01-01"), ]

# New column with the corresponding year for each date
df$year <- year(df$date)

ggplot(df, aes(x = date, y = unemploy)) +
  geom_line() +
  facet_wrap(~year, scales = "free") +
  scale_x_date(date_labels = "%Y-%m-%d")

library(dplyr)

# Data
df2 <- economics_long[economics_long$date > as.Date("2000-01-01"), ] %>%
  filter(variable == "pce" | variable == "unemploy")

ggplot(df2, aes(x = date, y = value, color = variable)) +
  geom_line() +
  theme(legend.position = "bottom")

install.packages("ggpmisc")
library(ggpmisc)

# Data
df <- economics[economics$date > as.Date("2000-01-01"), ]

ggplot(df, aes(x = date, y = unemploy)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2007-09-15"),
             linetype = 2, color = 2, linewidth = 1)

ggplot(df, aes(x = date, y = unemploy)) +
  geom_line() +
  stat_peaks(geom = "point", span = 15, color = "steelblue3", size = 2) +
  stat_peaks(geom = "label", span = 15, color = "steelblue3", angle = 0,
             hjust = -0.1, x.label.fmt = "%Y-%m-%d") +
  stat_peaks(geom = "rug", span = 15, color = "blue", sides = "b")

ggplot(df, aes(x = date, y = unemploy)) +
  geom_line() +
  stat_valleys(geom = "point", span = 11, color = "red", size = 2) +
  stat_valleys(geom = "label", span = 11, color = "red", angle = 0,
               hjust = -0.1, x.label.fmt = "%Y-%m-%d") +
  stat_valleys(geom = "rug", span = 11, color = "red", sides = "b")


