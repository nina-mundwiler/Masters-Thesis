# Station: NEU (Neuchâtel), precipitation dataset rre150d0
# Daily total precipitation, from 6 UTC to 6 UTC (next day)
# Temporal coverage: 31.12.1863 - 31.12.2023

# Packages for data processing

install.packages("dplyr")
library(dplyr)

install.packages("ggplot2")
library(ggplot2)

install.packages("zoo")
library(zoo)

# Preprocessing
read.table("NEU_rre150d0.txt",header=T,as.is=F)
NEU_P <- NEU_rre150d0

colnames(NEU_P)[2] <- "date"
colnames(NEU_P)[3] <- "precip"

NEU_P$date <- as.Date(as.character(NEU_P$date), format = "%Y%m%d")
NEU_P$year <- format(NEU_P$date, "%Y")
NEU_P$month <- format(NEU_P$date, "%m")

# Check for missing values
sum(is.na(NEU_P))

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
  
# R10mm
# Annual count of days with precipitation ≥ 10 mm
R10mm_days_NEU <- NEU_P[NEU_P$precip >= 10, ]
R10mm_summary_NEU <- aggregate(R10mm_days_NEU$precip, 
  by = list(Year = R10mm_days_NEU$year), FUN = length)
colnames(R10mm_summary_NEU)[2] <- "Count_R10mm_NEU"
R10mm_summary_NEU$Year <- as.numeric(R10mm_summary_NEU$Year)

# Linear trend analysis (based on all availible data)
model_R10mm_NEU <- lm(Count_R10mm_NEU ~ Year, data = R10mm_summary_NEU)
summary(model_R10mm_NEU)

# Data homogenization period (NEU) : 1864 to 2010
R10mm_summary_NEU$Dataset <- ifelse(R10mm_summary_NEU$Year <= 2010, 
                                    "homogenized","non-homogenized")
R10mm_homog_NEU <- R10mm_summary_NEU[R10mm_summary_NEU$Dataset == "homogenized", ]

# Graph with data homogenization period
ggplot(R10mm_summary_NEU, aes(x = Year, y = Count_R10mm_NEU)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "R10mm: Annual count of days with precipitation ≥ 10 mm (NEU)",
  x = "Year",y = "Number of days",
  caption = "p-value: non-significant") +
  scale_y_continuous(limits = c(10, 75), breaks = seq(10, 75, by = 10)) +
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) +
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal() + theme(legend.position = "none")

# Summary statistics (homogenized data only)
  #Mean
mean_R10mm_homog_NEU <- mean(R10mm_homog_NEU$Count_R10mm_NEU, na.rm = TRUE)

  #Usual range 
range_R10mm_homog_NEU <- quantile(R10mm_homog_NEU$Count_R10mm, 
  probs = c(0.05, 0.95), na.rm = TRUE)

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------

# R95pTOT
# Total annual precipitation from days exceeding the 95th percentile,
# Based on wet days (P ≥ 1 mm), during the 1961–1990 reference period  
  
# 95th percentile for wet days (1961–1990)
NEU_reference_period <- NEU_P[NEU_P$date >= as.Date("1961-01-01") & NEU_P$date <= as.Date("1990-12-31") &
  NEU_P$precip >= 1, ]
threshold_95_NEU <- quantile(NEU_reference_period$precip, 0.95, na.rm = TRUE)

# Days above the 95th percentile, and compute annual total
R95p_summary_NEU <- NEU_P[NEU_P$precip >= 1 & NEU_P$precip > threshold_95_NEU, ]
R95p_summary_NEU <- aggregate(R95p_summary_NEU$precip, 
  by = list(Year = R95p_summary_NEU$year),FUN = sum, na.rm = TRUE)
colnames(R95p_summary_NEU)[2] <- "Total_R95p_NEU"
R95p_summary_NEU$Year <- as.numeric(R95p_summary_NEU$Year)

# Linear trend analysis (based on all availible data)
model_R95p_NEU <- lm(Total_R95p_NEU ~ Year, data = R95p_summary_NEU)
summary(model_R95p_NEU)

# Data homogenization period (NEU) : 1864 to 2010
R95p_summary_NEU$Dataset <- ifelse(R95p_summary_NEU$Year <= 2010, 
  "homogenized","non-homogenized")
R95p_homog_NEU <- R95p_summary_NEU[R95p_summary_NEU$Dataset == "homogenized", ]

# Graph with data homogenization period
ggplot(R95p_summary_NEU, aes(x = Year, y = Total_R95p_NEU)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "R95p: Annual total precipitation from days > 95th percentile (NEU)",
  x = "Year", y = "Annual precipitation (mm)", caption = "p-value: non-significant") +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) +
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal() + theme(legend.position = "none") +
  annotate("text", x = 1860, y = 590, label = "95th percentile = 22.2 mm", hjust = 0, vjust = 1, size = 3, color = "black")

# Summary statistics (homogenized data only)
  # Mean
mean_R95p_homog_NEU <- mean(R95p_homog_NEU$Total_R95p_NEU, na.rm = TRUE)

 # Range
range_R95p_homog_NEU <- quantile(R95p_homog_NEU$Total_R95p_NEU, 
  probs = c(0.05, 0.95), na.rm = TRUE)
  
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
  
# Complementary approach : R95p days
# Annual count of days with precipitation > 95th percentile
# Based on wet days (P ≥ 1 mm), during the 1961–1990 reference period  
  
# 95th percentile for wet days (1961–1990)
NEU_reference_period <- NEU_P[NEU_P$date >= as.Date("1961-01-01") & NEU_P$date <= as.Date("1990-12-31") &
  NEU_P$precip >= 1, ]
threshold_95_days_NEU <- quantile(NEU_reference_period$precip, 0.95, na.rm = TRUE)

# Number of days per year exceeding this threshold
R95p_days_summary_NEU <- NEU_P[NEU_P$precip > threshold_95_days_NEU, ]
R95p_days_summary_NEU <- aggregate(R95p_days_summary_NEU$precip, 
  by = list(Year = R95p_days_summary_NEU$year),FUN = length)
colnames(R95p_days_summary_NEU)[2] <- "Days_R95p_NEU"
R95p_days_summary_NEU$Year <- as.numeric(R95p_days_summary_NEU$Year)  
  
# Linear trend analysis (based on all availible data)
model_R95p_days_NEU <- lm(Days_R95p_NEU ~ Year, data = R95p_days_summary_NEU)
summary(model_R95p_days_NEU)

# Data homogenization period (NEU) : 1864 to 2010
R95p_days_summary_NEU$Dataset <- ifelse(R95p_days_summary_NEU$Year <= 2010, 
  "homogenized","non-homogenized")
R95p_days_homog_NEU <- R95p_days_summary_NEU[R95p_days_summary_NEU$Dataset == "homogenized", ]

# Graph with data homogenization period
ggplot(R95p_days_summary_NEU, aes(x = Year, y = Days_R95p_NEU)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "Annual number of days with daily precipitation > 95th percentile (NEU)",
  x = "Year", y = "Number of days",
  caption = "p-value: non-significant") +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 2)) +
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) +
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal() + theme(legend.position = "none") +
  annotate("text", x = 1860, y = 19, 
  label = paste0("95th percentile = ", round(threshold_95_days_NEU, 1), " mm"), hjust = 0, vjust = 1, size = 3, color = "black")

# Summary statistics (homogenized data only)
  # Mean
mean_R95p_days_homog_NEU <- mean(R95p_days_homog_NEU$Days_R95p_NEU, na.rm = TRUE)

  # Range
range_R95p_days_homog_NEU <- quantile(R95p_days_homog_NEU$Days_R95p_NEU, 
 probs = c(0.05, 0.95),na.rm = TRUE)
  
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
  
# Rx1day
# Annual maximum of daily precipitation 
Rx1day_summary_NEU <- aggregate(NEU_P$precip, 
  by = list(Year = NEU_P$year), FUN = max, na.rm = TRUE)
colnames(Rx1day_summary_NEU)[2] <- "Yearly_Max_Rx1day_NEU"
Rx1day_summary_NEU$Year <- as.numeric(Rx1day_summary_NEU$Year)
  
# Linear trend analysis (based on all availible data)
model_Rx1day_NEU <- lm(Yearly_Max_Rx1day_NEU ~ Year, data = Rx1day_summary_NEU)
summary(model_Rx1day_NEU)
  
# Data homogenization period (NEU) : 1864 to 2010
Rx1day_summary_NEU$Dataset <- ifelse(Rx1day_summary_NEU$Year <= 2010, 
  "homogenized","non-homogenized")
Rx1day_homog_NEU <- Rx1day_summary_NEU[Rx1day_summary_NEU$Dataset == "homogenized", ]

# Graph with data homogenization period
ggplot(Rx1day_summary_NEU, aes(x = Year, y = Yearly_Max_Rx1day_NEU)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "Rx1day: Annual maximum of daily precipitation (NEU)",
  x = "Year",y = "Daily precipitation (mm)",
  caption = "p-value: non-significant") +
  scale_y_continuous(limits = c(15, 120), breaks = seq(20, 120, by = 10)) +
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) +
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal() + theme(legend.position = "none")
  
# Summary statistics (homogenized data only)
  # Mean 
mean_Rx1day_homog_NEU <- mean(Rx1day_homog_NEU$Yearly_Max_Rx1day_NEU, na.rm = TRUE)

  # Range 
range_Rx1day_homog_NEU <- quantile(Rx1day_homog_NEU$Yearly_Max_Rx1day_NEU, 
  probs = c(0.05, 0.95),na.rm = TRUE)

-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
  
# CWD
# Maximum number of consecutive wet days
  
# Binary wet-day indicator
NEU_P$Wet_Day <- ifelse(NEU_P$precip >= 1, 1, 0)

# Create the loop to compute the maximum wet spell length
CWD_summary_NEU <- data.frame(Year = unique(NEU_P$year), Max_CWD_NEU = NA)
for (i in seq_along(CWD_summary_NEU$Year)) 
  {year_data <- NEU_P[NEU_P$year == CWD_summary_NEU$Year[i], ]
  rle_wet <- rle(year_data$Wet_Day)
  CWD_summary_NEU$Max_CWD_NEU[i] <- if (any(rle_wet$values == 1)) 
  {max(rle_wet$lengths[rle_wet$values == 1], na.rm = TRUE)} 
  else
  {0}}  
CWD_summary_NEU$Year <- as.numeric(CWD_summary_NEU$Year)

# Linear trend analysis (based on all availible data)
model_CWD_NEU <- lm(Max_CWD_NEU ~ Year, data = CWD_summary_NEU)
summary(model_CWD_NEU)

# Data homogenization period (NEU) : 1864 to 2010
CWD_summary_NEU$Dataset <- ifelse(CWD_summary_NEU$Year <= 2010, 
  "homogenized","non-homogenized")
CWD_homog_NEU <- CWD_summary_NEU[CWD_summary_NEU$Dataset == "homogenized", ]

# Graph with data homogenization period
ggplot(CWD_summary_NEU, aes(x = Year, y = Max_CWD_NEU)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) +
  labs(title = "CWD: Maximum number of consecutive wet days (NEU)",
  x = "Year",y = "Number of days",caption = "p-value: non-significant") +
  scale_y_continuous(limits = c(3, 23), breaks = seq(4, 23, by = 2)) +
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) +
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal() + theme(legend.position = "none")

# Summary statistics (homogenized data only)
  # Mean 
mean_CWD_homog_NEU <- mean(CWD_homog_NEU$Max_CWD_NEU, na.rm = TRUE)

  # Range 
range_CWD_homog_NEU <- quantile(CWD_homog_NEU$Max_CWD_NEU, 
  probs = c(0.05, 0.95),  na.rm = TRUE)

---
  #fin 
