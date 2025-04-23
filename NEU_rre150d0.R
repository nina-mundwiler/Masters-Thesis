#P; somme journalière conventionnelle 6 UTC - 6 UTC du jour suivant;  18631231-20231231

read.table("NEU_rre150d0.txt",header=T,as.is=F)

NEU_sommeP <- NEU_rre150d0
colnames(NEU_sommeP)[2] <- "date"
colnames(NEU_sommeP)[3] <- "precip"

sum(is.na(NEU_sommeP))


NEU_sommeP$date <- as.Date(as.character(NEU_sommeP$date), format = "%Y%m%d")

#packages
install.packages("lubridate")
library(lubridate)
install.packages("dplyr")
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("zoo")
library(zoo)

NEU_sommeP$year <- format(NEU_sommeP$date, "%Y")
NEU_sommeP$month <- format(NEU_sommeP$date, "%m")


-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
#Rx1day : Précipitations maximales sur une journée pendant un mois donné 
yearly_max <- NEU_sommeP %>% group_by(year) %>%
  summarise(Yearly_Max_Rx1day = max(precip, na.rm = TRUE),   
    Date_Max = date[which.max(precip)]) %>%
  ungroup() %>%mutate(year = as.numeric(year)) 

#régression linéaire
modele_rx1day <- lm(Yearly_Max_Rx1day ~ year, data = yearly_max)
summary_modele_rx1day <- summary(modele_rx1day)

r_squared_rx1day <- summary_modele_rx1day$r.squared
p_value_rx1day <- summary_modele_rx1day$coefficients[2, 4] 
slope_estimate_rx1day <- summary_modele_rx1day$coefficients[2, 1]  
slope_error_rx1day <- summary_modele_rx1day$coefficients[2, 2]  

ggplot(yearly_max, aes(x = year, y = Yearly_Max_Rx1day)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) + 
  labs(title = "Rx1day : Annual maximum of daily precipitation (NEU)",
       x = "Year",
       y = "Daily precipitation (mm)") +
  scale_y_continuous(limits = c(20, 120), breaks = seq(20, 120, by = 10)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()+theme(plot.title = element_text(size = 12)) 


#avec homognéisation 

yearly_max$Dataset <- ifelse(yearly_max$year <= 2010, "homogenized", "non-homogenized")

ggplot(yearly_max, aes(x = year, y = Yearly_Max_Rx1day)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) + 
  labs(
    title = "Rx1day: Annual maximum of daily precipitation (NEU)",
    x = "Year",
    y = "Daily precipitation (mm)",
    color = "Dataset",
    caption = "p-value: non-significant"
  ) +
  scale_y_continuous(limits = c(15, 120), breaks = seq(20, 120, by = 10)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))



--
  
  # 1. Filtrer les données homogénéisées
  Rx1day_NEU_homog <- yearly_max %>%
  filter(Dataset == "homogenized")

mean_Rx1day_NEU <- mean(Rx1day_NEU_homog$Yearly_Max_Rx1day, na.rm = TRUE)
quantiles_Rx1day_NEU <- quantile(Rx1day_NEU_homog$Yearly_Max_Rx1day, probs = c(0.05, 0.95), na.rm = TRUE)

mean_Rx1day_NEU
quantiles_Rx1day_NEU


  


-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
#Rx5days : Précipitations maximales sur cinq jours pendant un mois donné 
NEU_sommeP <- NEU_sommeP %>% arrange(date) %>%  group_by(year) %>%
  mutate(Rx5day = rollapply(precip, width = 5, FUN = sum, align = "right", fill = NA, na.rm = TRUE)) %>%
  ungroup()

yearly_max_rx5day <- NEU_sommeP %>%group_by(year) %>%
  summarise(Yearly_Max_Rx5day = max(Rx5day, na.rm = TRUE)) %>%ungroup() %>%
  mutate(year = as.numeric(year))  

modele_rx5day <- lm(Yearly_Max_Rx5day ~ year, data = yearly_max_rx5day)
summary_modele_rx5day <- summary(modele_rx5day)

r_squared_rx5day <- summary_modele_rx5day$r.squared
p_value_rx5day <- summary_modele_rx5day$coefficients[2, 4]  
slope_estimate_rx5day <- summary_modele_rx5day$coefficients[2, 1] 
slope_error_rx5day <- summary_modele_rx5day$coefficients[2, 2] 

ggplot(yearly_max_rx5day, aes(x = year, y = Yearly_Max_Rx5day)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) + 
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
  labs(title = "Rx5day : Évolution des précipitations maximales annuelles sur 5 jours consécutifs (NEU)",
  x = "Année",
  y = "Rx5day (mm)") +
  scale_y_continuous(limits = c(40, 140), breaks = seq(40, 140, by = 10)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()+
  theme(plot.title = element_text(size = 12)) 



---------------------------------------------------------------------------
---------------------------------------------------------------------------
#sdii: simple precipitation intensity index 
#supprimer valeurs NA
NEU_sommeP <- NEU_sommeP[!is.na(NEU_sommeP$precip), ]

#stocker les valeurs SDII par année
SDII_yearly <- numeric(length(unique(NEU_sommeP$year)))
years <- unique(NEU_sommeP$year)

for (i in seq_along(years)) {year_data <- NEU_sommeP[NEU_sommeP$year == years[i], ]
  wet_days <- year_data[year_data$precip >= 1, ]
  total_precipitation <- sum(wet_days$precip, na.rm = TRUE)
  number_of_wet_days <- nrow(wet_days)
  if (number_of_wet_days > 0) {SDII_yearly[i] <- total_precipitation / number_of_wet_days} else {
    SDII_yearly[i] <- NA  }}

results_yearly_sdii <- data.frame( Year = years,SDII = SDII_yearly)

results_yearly_sdii$Year <- as.numeric(as.character(results_yearly_sdii$Year))

#modèle linéaire
modele_sdii <- lm(SDII ~ Year, data = results_yearly_sdii[results_yearly_sdii$Year != 1874, ])
summary_modele_sdii <- summary(modele_sdii)

r_squared_sdii <- summary_modele_sdii$r.squared
p_value_sdii <- summary_modele_sdii$coefficients[2, 4] 
slope_estimate_sdii <- summary_modele_sdii$coefficients[2, 1] 
slope_error_sdii <- summary_modele_sdii$coefficients[2, 2] 


#graph
ggplot(results_yearly_sdii[results_yearly_sdii$Year != 1874, ], aes(x = Year, y = SDII)) +
  geom_line(color = "blue", linewidth = 0.5, group = 1) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) + 
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.5) + 
  labs(title = "SDII : Indice simple de l’intensité des précipitations sur les jours de pluie ≥ 1mm (NEU)",
       x = "Année",
       y = "SDII (mm)") +
  scale_y_continuous(limits = c(5, 11), breaks = seq(5, 11, by = 1)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()+
  theme(plot.title = element_text(size = 12)) 



-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
#R10mm : Nombre annuel de jours où précipitations ≥ 10mm
R10mm <- NEU_sommeP %>%
  filter(precip >= 10) %>%  
  group_by(year) %>%
  summarise(Count_R10mm = n()) %>% 
  ungroup() %>%
  mutate(year = as.numeric(year)) 


mean(R10mm$Count_R10mm)

R10mm$Dataset <- ifelse(R10mm$year < 2009, "homogenized", "non-homogenized")
R10mm_homog <- R10mm %>%
  filter(Dataset == "homogenized")
quantile(R10mm_homog$Count_R10mm, probs = c(0.05, 0.95), na.rm = TRUE)



modele_r10mm <- lm(Count_R10mm ~ year, data = R10mm)
summary_modele_r10mm <- summary(modele_r10mm)

r_squared_r10mm <- summary_modele_r10mm$r.squared
p_value_r10mm <- summary_modele_r10mm$coefficients[2, 4]  
slope_estimate_r10mm <- summary_modele_r10mm$coefficients[2, 1]  
slope_error_r10mm <- summary_modele_r10mm$coefficients[2, 2] 
  
#graph
ggplot(R10mm, aes(x = year, y = Count_R10mm)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
  labs(title = "R10mm : Évolution du nombre de jours avec précipitations ≥ 10 mm (NEU)",
  x = "Année",
  y = "Nombre de jours") +
  scale_y_continuous(limits = c(10,60), breaks = seq(10,60, by = 10)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()

#avec couleurs différenciées selon homogénéisation des données 
R10mm$Dataset <- ifelse(R10mm$year < 2010, "homogenized", "non-homogenized")

p_value_r10mm <- 0.59

ggplot(R10mm, aes(x = year, y = Count_R10mm)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +
  labs(
    title = "R10mm: Annual count of days when daily precipitation is ≥ 10 mm (NEU)",
    x = "Year",
    y = "Number of days",
    caption = "p-value: non-significant"
  ) +
  scale_y_continuous(limits = c(10, 75), breaks = seq(10, 75, by = 10)) +
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) +
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal()

---
  
  
  # 1. Filtrer les données homogénéisées
  R10mm_NEU_homog <- R10mm %>%
  filter(Dataset == "homogenized")

# 2. Moyenne
mean_R10mm_NEU <- mean(R10mm_NEU_homog$Count_R10mm, na.rm = TRUE)

# 3. Quantiles 5 % et 95 %
quantiles_R10mm_NEU <- quantile(R10mm_NEU_homog$Count_R10mm, probs = c(0.05, 0.95), na.rm = TRUE)

# Résultats
mean_R10mm_NEU
quantiles_R10mm_NEU








-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
#R20mm : Nombre annuel de jours où précipitations ≥ 20mm
R20mm <- NEU_sommeP %>%
  filter(precip >= 20) %>% 
  group_by(year) %>%
  summarise(Count_R20mm = n()) %>%  
  ungroup() %>%
  mutate(year = as.numeric(year)) 

#modèle linéaire 
modele_r20mm <- lm(Count_R20mm ~ year, data = R20mm)
summary_modele_r20mm <- summary(modele_r20mm)

r_squared_r20mm <- summary_modele_r20mm$r.squared
p_value_r20mm <- summary_modele_r20mm$coefficients[2, 4]  
slope_estimate_r20mm <- summary_modele_r20mm$coefficients[2, 1] 
slope_error_r20mm <- summary_modele_r20mm$coefficients[2, 2] 

  
#graph
ggplot(R20mm, aes(x = year, y = Count_R20mm)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) + 
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) + 
  labs(title = "R20mm : Évolution du nombre de jours avec précipitations ≥ 20 mm (NEU)",
      x = "Année",
      y = "Nombre de jours") +
  scale_y_continuous(limits = c(0,25), breaks = seq(0,25, by = 5)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()
----------------------

#test 30 et 50mm/24h = degré de danger 2 et 3 (non sig)
R50mm <- NEU_sommeP %>%
  filter(precip >= 50) %>% 
  group_by(year) %>%
  summarise(Count_R50mm = n()) %>%  
  ungroup() %>%
  mutate(year = as.numeric(year)) 

#modèle linéaire 
modele_r50mm <- lm(Count_R50mm ~ year, data = R50mm)
summary_modele_r50mm <- summary(modele_r50mm)

r_squared_r50mm <- summary_modele_r50mm$r.squared
p_value_r50mm <- summary_modele_r50mm$coefficients[2, 4]  
slope_estimate_r50mm <- summary_modele_r50mm$coefficients[2, 1] 
slope_error_r50mm <- summary_modele_r50mm$coefficients[2, 2] 


#graph
ggplot(R50mm, aes(x = year, y = Count_R50mm)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) + 
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) + 
  labs(title = "R50mm : Évolution du nombre de jours avec précipitations ≥ 50 mm (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  scale_y_continuous(limits = c(0,3), breaks = seq(0,3, by = 1)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()


-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
#CDD max length of dry spell (max consec days with RR<1mm) : Durée maximale d’une période sèche, nombre maximal de jours consécutifs où RR < 1mm
cdd_data <- NEU_sommeP %>%
  mutate(Dry_Spell = ifelse(precip < 1, 1, 0)) %>%  
  group_by(year) %>%
  summarise(Max_CDD = max(rle(Dry_Spell)$lengths[rle(Dry_Spell)$values == 1], na.rm = TRUE)) %>% 
  ungroup() %>%
  mutate(year = as.numeric(year)) 

#modèle linéaire 
modele_cdd <- lm(Max_CDD ~ year, data = cdd_data)
summary_modele_cdd <- summary(modele_cdd)

r_squared_cdd <- summary_modele_cdd$r.squared
p_value_cdd <- summary_modele_cdd$coefficients[2, 4]  
slope_estimate_cdd <- summary_modele_cdd$coefficients[2, 1] 
slope_error_cdd <- summary_modele_cdd$coefficients[2, 2]  

ggplot(cdd_data, aes(x = year, y = Max_CDD)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
  labs(title = "CDD : Évolution du nombre de jours secs consécutifs (<1mm) (NEU)",
  x = "Année",
 y = "Nombre de jours") +
  scale_y_continuous(limits = c(10,60), breaks = seq(10,60, by = 5)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()

-------------------------------------------------------------------------------
  -------------------------------------------------------------------------------
#CWD max length of wet spell (max consec days with RR>=1mm) : Durée maximale d’une période pluvieuse, nombre maximal de jours consécutifs où RR ≥ 1mm
cwd_data <- NEU_sommeP %>%
  mutate(Wet_Spell = ifelse(precip >= 1, 1, 0)) %>%
  group_by(year) %>%
  summarise(Max_CWD = max(rle(Wet_Spell)$lengths[rle(Wet_Spell)$values == 1], na.rm = TRUE)) %>% 
  ungroup() %>%
  mutate(year = as.numeric(year)) 

#modèle linéaire 
modele_cwd <- lm(Max_CWD ~ year, data = cwd_data)
summary_modele_cwd <- summary(modele_cwd)

r_squared_cwd <- summary_modele_cwd$r.squared
p_value_cwd <- summary_modele_cwd$coefficients[2, 4] 
slope_estimate_cwd <- summary_modele_cwd$coefficients[2, 1] 
slope_error_cwd <- summary_modele_cwd$coefficients[2, 2] 

  
#graph
ggplot(cwd_data, aes(x = year, y = Max_CWD)) +
   geom_line(color = "blue", size = 0.5) +
    geom_point(color = "blue") +
   geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) + 
 geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
    labs(title = "CWD : Évolution du nombre de jours humides consécutifs (≥1mm) (NEU)",
     x = "Année",
     y = "Nombre de jours") +
  scale_y_continuous(limits = c(4,15), breaks = seq(4,15, by = 1)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()
  


#inclure homogénéisation 
cwd_data$Dataset <- ifelse(cwd_data$year <= 2010, "homogenized", "non-homogenized")

ggplot(cwd_data, aes(x = year, y = Max_CWD)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) + 
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
  labs(
    title = "CWD: Maximum number of consecutive wet days (daily precipitation ≥ 1 mm) (NEU)",
    x = "Year",
    y = "Number of days", 
    color = "Dataset",
    caption = "p-value: non-significant"
  ) +
  scale_y_continuous(limits = c(3, 23), breaks = seq(4, 23, by = 2)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  scale_color_manual(values = c("homogenized" = "darkblue", "non-homogenized" = "skyblue")) +
  theme_minimal()

--
  # 1. Filtrer les données homogénéisées
  cwd_NEU_homog <- cwd_data %>%
  filter(Dataset == "homogenized")

# 2. Moyenne
mean_cwd_NEU <- mean(cwd_NEU_homog$Max_CWD, na.rm = TRUE)

# 3. Quantiles 5 % et 95 %
quantiles_cwd_NEU <- quantile(cwd_NEU_homog$Max_CWD, probs = c(0.05, 0.95), na.rm = TRUE)

# Résultats
mean_cwd_NEU
quantiles_cwd_NEU

  
  
  
  
  
  --
  
  

-------------------------------------------------------------------------------
  -------------------------------------------------------------------------------
#R95pTOT : Total annuel de PRCP lorsque RR > 95e percentile

#95e centile pour les jours humides (RR >= 1 mm) ; 1961 et 1990
centile95_6190 <- NEU_sommeP %>%
  filter(date >= as.Date("1961-01-01") & date <= as.Date("1990-12-31") & precip >= 1) %>%
  summarise(RRwn95 = quantile(precip, 0.95, na.rm = TRUE)) %>%
  pull(RRwn95)


#jours avec P > 95e et somme annuelle
R95pTOT <- NEU_sommeP %>%
  filter(precip >= 1 & precip > centile95_6190) %>%  
  group_by(year) %>%
  summarise(Total_R95p = sum(precip, na.rm = TRUE)) %>% 
  ungroup() %>%
  mutate(year = as.numeric(year))  


#modèle linéaire 
modele_r95pTOT <- lm(Total_R95p ~ year, data = R95pTOT)
summary_modele_r95pTOT <- summary(modele_r95pTOT)

r_squared_r95pTOT <- summary_modele_r95pTOT$r.squared
p_value_r95pTOT <- summary_modele_r95pTOT$coefficients[2, 4]  # P-value de la pente
slope_estimate_r95pTOT <- summary_modele_r95pTOT$coefficients[2, 1]  # Estimation de la pente
slope_error_r95pTOT <- summary_modele_r95pTOT$coefficients[2, 2]  # Erreur de la pente

  
#graph
ggplot(R95pTOT, aes(x = year, y = Total_R95p)) +
  geom_line(color = "blue", size = 0.5) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
  labs(title = "R95pTOT : Évolution du total annuel des précipitations lorsque P > 95e centile (NEU)",
    x = "Année",y = "R95pTOT (mm)") +
  scale_y_continuous(limits = c(0,560), breaks = seq(0,560, by = 50)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal() +
  annotate("text", x = 2000, y = 550, label = "95e centile = 22,2 mm", 
           hjust = 0, vjust = 1, size = 3, color = "black")



#avec homogénéisation 
R95pTOT$Dataset <- ifelse(R95pTOT$year <= 2010, "homogenized", "non-homogenized")

ggplot(R95pTOT, aes(x = year, y = Total_R95p)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", 
              linetype = "dashed", size = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
  labs(
    title = "R95p: Annual total precipitation from days exceeding the 95th percentile (NEU)",
    x = "Year",
    y = "Annual precipitation (mm)",
    color = "Dataset",
    caption = "p-value: non-significant"
  ) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, by = 100)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  scale_color_manual(values = c("homogenized" = "darkblue", 
                                "non-homogenized" = "skyblue")) +
  theme_minimal() +
  annotate("text", x = 1860, y = 590, label = "95th percentile = 22,2 mm", 
           hjust = 0, vjust = 1, size = 3, color = "black")

--
  
  # 1. Filtrer les données homogénéisées
  R95pTOT_NEU_homog <- R95pTOT %>%
  filter(Dataset == "homogenized")

# 2. Moyenne
mean_R95pTOT_NEU <- mean(R95pTOT_NEU_homog$Total_R95p, na.rm = TRUE)

# 3. Quantiles 5 % et 95 %
quantiles_R95pTOT_NEU <- quantile(R95pTOT_NEU_homog$Total_R95p, probs = c(0.05, 0.95), na.rm = TRUE)

# Résultats
mean_R95pTOT_NEU
quantiles_R95pTOT_NEU


  
  
  ---
  








#et suite du test, voir the annual number of days exceeding the 95th percentile
centile95_6190 <- NEU_sommeP %>%
  filter(date >= as.Date("1961-01-01") & date <= as.Date("1990-12-31") & precip >= 1) %>%
  summarise(RRwn95 = quantile(precip, 0.95, na.rm = TRUE)) %>%
  pull(RRwn95)

R95p_days <- NEU_sommeP %>%
  filter(precip >= 1 & precip > centile95_6190) %>%  
  group_by(year) %>%
  summarise(Days_R95p = n()) %>%  
  ungroup() %>%
  mutate(year = as.numeric(year))

#régression linéaire
model_r95p_days <- lm(Days_R95p ~ year, data = R95p_days)
summary(model_r95p_days)
p_value_r95p_days <- summary(model_r95p_days)$coefficients["year", "Pr(>|t|)"]


R95p_days$Dataset <- ifelse(R95p_days$year <= 2010, "homogenized", "non-homogenized")

ggplot(R95p_days, aes(x = year, y = Days_R95p)) +
  geom_line(aes(color = Dataset), size = 0.5) +
  geom_point(aes(color = Dataset)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", 
              linetype = "dashed", size = 0.5) +  
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
  labs(
    title = "Annual number of days with daily precipitation > 95th percentile (NEU)",
    x = "Year",
    y = "Number of days",
    color = "Dataset",
    caption = "p-value: non-significant"
  ) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 2)) +
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) +
  scale_color_manual(values = c("homogenized" = "darkblue", 
                                "non-homogenized" = "skyblue")) +
  theme_minimal() +
  annotate("text", x = 1860, y = 19, 
           label = paste0("95th percentile = ", round(centile95_6190, 1), " mm"),
           hjust = 0, vjust = 1, size = 3, color = "black")




---
  
  # 1. Filtrer les données homogénéisées
  R95p_days_NEU_homog <- R95p_days %>%
  filter(Dataset == "homogenized")

# 2. Moyenne
mean_R95p_days_NEU <- mean(R95p_days_NEU_homog$Days_R95p, na.rm = TRUE)

# 3. Quantiles 5 % et 95 %
quantiles_R95p_days_NEU <- quantile(R95p_days_NEU_homog$Days_R95p, probs = c(0.05, 0.95), na.rm = TRUE)

# Résultats
mean_R95p_days_NEU
quantiles_R95p_days_NEU

  
  
  
  
  ---
  
  




-------------------------------------------------------------------------------
  -------------------------------------------------------------------------------
#R99pTOT : Total annuel de PRCP lorsque RR > 99e percentile
#99e centile pour les jours humides (RR >= 1 mm), selon norme 61-90
centile99_6190 <- NEU_sommeP %>%
  filter(date >= as.Date("1961-01-01") & date <= as.Date("1990-12-31") & precip >= 1) %>%
  summarise(RRwn99 = quantile(precip, 0.99, na.rm = TRUE)) %>%
  pull(RRwn99)

#jours avec des P > 99e
R99pTOT <- NEU_sommeP %>%
  filter(precip >= 1 & precip > centile99_6190) %>% 
  group_by(year) %>%
  summarise(Total_R99p = sum(precip, na.rm = TRUE)) %>%  
  ungroup()

#modèle linéaire 
modele_r99pTOT <- lm(Total_R99p ~ year, data = R99pTOT)
summary_modele_r99pTOT <- summary(modele_r99pTOT)

r_squared_r99pTOT <- summary_modele_r99pTOT$r.squared
p_value_r99pTOT <- summary_modele_r99pTOT$coefficients[2, 4]  
slope_estimate_r99pTOT <- summary_modele_r99pTOT$coefficients[2, 1] 
slope_error_r99pTOT <- summary_modele_r99pTOT$coefficients[2, 2] 

R99pTOT$year <- as.numeric(as.character(R99pTOT$year)) 

#graph
ggplot(R99pTOT, aes(x = year, y = Total_R99p)) +
   geom_line(color = "blue", size = 0.5) +
   geom_point(color = "blue") +
   geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) +  # Droite de régression
  geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  # Courbe LOESS
  labs(title = "R99pTOT : Évolution du total annuel des précipitations lorsque P > 99e centile (NEU)",
   x = "Année",      y = "R99pTOT (mm)") +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,250, by = 50)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal() +
  annotate("text", x = 2000, y = 250, label = "95e centile = 39,7 mm", 
           hjust = 0, vjust = 1, size = 3, color = "black")


-------------------------------------------------------------------------------
  -------------------------------------------------------------------------------
#PRCPTOT : Total annuel des précipitations les jours pluvieux
PRCPTOT <- NEU_sommeP %>%
  filter(precip >= 1) %>%  
  group_by(year) %>%
  summarise(Total_Precipitation = sum(precip, na.rm = TRUE)) %>% 
  ungroup() %>%
  mutate(year = as.numeric(year))  

#modèle linéaire
modele_prcptot <- lm(Total_Precipitation ~ year, data = PRCPTOT)
summary_modele_prcptot <- summary(modele_prcptot)
  
r_squared_prcptot <- summary_modele_prcptot$r.squared
p_value_prcptot <- summary_modele_prcptot$coefficients[2, 4]  
slope_estimate_prcptot <- summary_modele_prcptot$coefficients[2, 1] 
slope_error_prcptot <- summary_modele_prcptot$coefficients[2, 2] 
  
#graph
ggplot(PRCPTOT, aes(x = year, y = Total_Precipitation)) +
 geom_line(color = "blue", size = 0.5) +
 geom_point(color = "blue") +
   geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5) + 
 geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.5) +  
 labs(title = "PRCPTOT : Évolution du total annuel des précipitations sur les jours humides (NEU)",
     x = "Année",
      y = "PRCPTOT (mm)") +
  scale_y_continuous(limits = c(500,1500), breaks = seq(500,1500, by = 250)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()

----------------------------------------
----------------------------------------
  
#autres mesures perso (rapport climatologique)

#évolution des sommes annuelles de P 
somme_annuelle_P <- NEU_sommeP %>%
  group_by(year = format(date, "%Y")) %>% 
  summarise(total_precip = sum(precip, na.rm = TRUE)) %>% 
  ungroup()

somme_annuelle_P$year <- as.numeric(somme_annuelle_P$year)

ggplot(somme_annuelle_P, aes(x = year, y = total_precip)) +
  geom_line(color = "blue", size = 0.5) + 
  geom_smooth(method = "loess", color = "black") +  
  labs(title = "Évolution des sommes annuelles de précipitations (NEU)",
       x = "Année",
       y = "Somme annuelle de précipitations (mm)") +
  scale_y_continuous(limits = c(500,1500), breaks = seq(500,1500, by = 250)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()  




#évolution des sommes par saison 
NEU_sommeP <- NEU_sommeP %>%
  mutate(season = case_when(
    month(date) %in% c(12, 1, 2) ~ "Hiver",
    month(date) %in% c(3, 4, 5) ~ "Printemps",
    month(date) %in% c(6, 7, 8) ~ "Été",
    month(date) %in% c(9, 10, 11) ~ "Automne",
    TRUE ~ "Inconnu"  
  ))

somme_precip_saison <- NEU_sommeP %>%
  group_by(year = format(date, "%Y"), season) %>%
  summarise(total_precip = sum(precip, na.rm = TRUE)) %>%
  ungroup()

somme_precip_saison$year <- as.numeric(somme_precip_saison$year)

somme_precip_saison$season <- factor(somme_precip_saison$season, 
                                     levels = c("Printemps", "Été", "Automne", "Hiver"))

ggplot(somme_precip_saison, aes(x = year, y = total_precip, color = season)) +
  geom_line(size = 0.25) +
  geom_smooth(method = "loess", se = FALSE, aes(group = season), size = 0.5) + 
  labs(title = "Évolution des sommes annuelles de précipitations par saison (NEU)",
       x = "Année",
       y = "Précipitations totales (mm)") +
  scale_color_manual(values = c("Hiver" = "cornflowerblue", 
                                "Printemps" = "darkolivegreen2", 
                                "Été" = "gold", 
                                "Automne" = "darkorange")) +
  scale_y_continuous(limits = c(50,650), breaks = seq(100,650, by = 100)) +  
  scale_x_continuous(breaks = seq(1875, 2025, by = 25)) + 
  theme_minimal()


#nombre de jours de P>=20mm par saison, et évolution
#Printemps 
printemps_data <- NEU_sommeP %>%
  filter(month(date) %in% c(3, 4, 5)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

printemps_data$year <- as.numeric(printemps_data$year)

ggplot(printemps_data, aes(x = year, y = count)) +
  geom_line(color = "darkolivegreen2") + 
  geom_point(color = "darkolivegreen2") +  
  labs(title = "Nombre de jours avec précipitations ≥ 20 mm - Printemps (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  theme_minimal()

    #calcul moyenne saisonnière 2000-2023
printemps_data_2000_2023 <- NEU_sommeP %>%
  filter(month(date) %in% c(3, 4, 5) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

moyenne_printemps <- mean(printemps_data_2000_2023$count, na.rm = TRUE)


#été 
ete_data <- NEU_sommeP %>%
  filter(month(date) %in% c(6, 7, 8)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

ete_data$year <- as.numeric(ete_data$year)

ggplot(ete_data, aes(x = year, y = count)) +
  geom_line(color = "gold") + 
  geom_point(color = "gold") +  
  labs(title = "Nombre de jours avec précipitations ≥ 20 mm - Été (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +  
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  theme_minimal()

ete_data_2000_2023 <- NEU_sommeP %>%
  filter(month(date) %in% c(6, 7, 8) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

moyenne_ete <- mean(ete_data_2000_2023$count, na.rm = TRUE)

#automne 
automne_data <- NEU_sommeP %>%
  filter(month(date) %in% c(9, 10, 11)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

automne_data$year <- as.numeric(automne_data$year)

ggplot(automne_data, aes(x = year, y = count)) +
  geom_line(color = "darkorange") + 
  geom_point(color = "darkorange") +  
  labs(title = "Nombre de jours avec précipitations ≥ 20 mm - Automne (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  scale_y_continuous(limits = c(0, 13), breaks = seq(0, 13, by = 2)) +  
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  theme_minimal()

automne_data_2000_2023 <- NEU_sommeP %>%
  filter(month(date) %in% c(9, 10, 11) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

moyenne_automne <- mean(automne_data_2000_2023$count, na.rm = TRUE)


#hiver
hiver_data <- NEU_sommeP %>%
  filter(month(date) %in% c(12, 1, 2)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

hiver_data$year <- as.numeric(hiver_data$year)

ggplot(hiver_data, aes(x = year, y = count)) +
  geom_line(color = "cornflowerblue") + 
  geom_point(color = "cornflowerblue") +  
  labs(title = "Nombre de jours avec précipitations ≥ 20 mm - Hiver (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2)) +  
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  theme_minimal()

hiver_data_2000_2023 <- NEU_sommeP %>%
  filter(month(date) %in% c(12, 1, 2) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 20, na.rm = TRUE)) %>%  
  ungroup()

moyenne_hiver <- mean(hiver_data_2000_2023$count, na.rm = TRUE)

#nombre de jours de P>=99th centile (39.71mm) par saison, et évolution
#printemps 
printemps_data_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(3, 4, 5)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

printemps_data_99th$year <- as.numeric(printemps_data_99th$year)

ggplot(printemps_data_99th, aes(x = year, y = count)) +
  geom_line(color = "darkolivegreen2") + 
  geom_point(color = "darkolivegreen2") +  
  labs(title = "Nombre de jours avec précipitations ≥ 39.718 mm - Printemps (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  theme_minimal()

printemps_data_2000_2023_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(3, 4, 5) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

moyenne_printemps_99th <- mean(printemps_data_2000_2023_99th$count, na.rm = TRUE)


#été 
ete_data_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(6, 7, 8)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

ete_data_99th$year <- as.numeric(ete_data_99th$year)

ggplot(ete_data_99th, aes(x = year, y = count)) +
  geom_line(color = "gold") + 
  geom_point(color = "gold") +  
  labs(title = "Nombre de jours avec précipitations ≥ 39.718 mm - Été (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  theme_minimal()

ete_data_2000_2023_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(6, 7, 8) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

moyenne_ete_99th <- mean(ete_data_2000_2023_99th$count, na.rm = TRUE)


#automne 
automne_data_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(9, 10, 11)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

automne_data_99th$year <- as.numeric(automne_data_99th$year)

ggplot(automne_data_99th, aes(x = year, y = count)) +
  geom_line(color = "darkorange") + 
  geom_point(color = "darkorange") +  
  labs(title = "Nombre de jours avec précipitations ≥ 39.718 mm - Automne (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  theme_minimal()

automne_data_2000_2023_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(9, 10, 11) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

moyenne_automne_99th <- mean(automne_data_2000_2023_99th$count, na.rm = TRUE)


hiver_data_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(12, 1, 2)) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

hiver_data_99th$year <- as.numeric(hiver_data_99th$year)

ggplot(hiver_data_99th, aes(x = year, y = count)) +
  geom_line(color = "cornflowerblue") + 
  geom_point(color = "cornflowerblue") +  
  labs(title = "Nombre de jours avec précipitations ≥ 39.718 mm - Hiver (NEU)",
       x = "Année",
       y = "Nombre de jours") +
  theme_minimal()

hiver_data_2000_2023_99th <- NEU_sommeP %>%
  filter(month(date) %in% c(12, 1, 2) & year(date) >= 2000 & year(date) <= 2023) %>%  
  group_by(year = format(date, "%Y")) %>% 
  summarise(count = sum(precip >= 39.718, na.rm = TRUE)) %>%  
  ungroup()

moyenne_hiver_99th <- mean(hiver_data_2000_2023_99th$count, na.rm = TRUE)




  
  
# écart à la norme des saisons
#Printemps 
norme6190_printemps <- 70 + 64 + 84
printemps_data <- NEU_sommeP %>% filter(month(date) %in% c(3, 4, 5))  

somme_printemps_annuelle <- printemps_data %>%
  group_by(year = format(date, "%Y")) %>% 
  summarise(total_precip_printemps = sum(precip, na.rm = TRUE)) %>%  
  ungroup()

somme_printemps_annuelle$year <- as.numeric(somme_printemps_annuelle$year)

somme_printemps_annuelle$ecart_norme <- somme_printemps_annuelle$total_precip_printemps - norme6190_printemps

ggplot(somme_printemps_annuelle, aes(x = year)) +
  geom_bar(aes(y = ecart_norme, fill = ecart_norme >= 0), stat = "identity") +  
  geom_smooth(aes(y = ecart_norme), method = "loess", color = "black", size = 0.5) +  
  scale_fill_manual(values = c("red", "blue"), labels = c("Écart positif", "Écart négatif")) + 
  labs(title = "Écart des précipitations printanières par rapport à la norme (1961-1990) (NEU)",
       x = "Année",
       y = "Écart à la norme (mm)") +
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  scale_y_continuous(limits = c(-160, 300), breaks = seq(-150,300, by = 50)) +  
  theme_minimal() +  
  theme(legend.position = "none") +
  annotate("text", x = 1864, y = 296, label = "norme printanière 61-90 = 218 mm", 
           hjust = 0, vjust = 1, size = 3, color = "black")

#Eté
norme6190_ete <- 94 + 81 + 99
ete_data <- NEU_sommeP %>% filter(month(date) %in% c(6, 7, 8))

somme_ete_annuelle <- ete_data %>%
  group_by(year = format(date, "%Y")) %>%
  summarise(total_precip_ete = sum(precip, na.rm = TRUE)) %>%
  ungroup()

somme_ete_annuelle$year <- as.numeric(somme_ete_annuelle$year)
somme_ete_annuelle$ecart_norme <- somme_ete_annuelle$total_precip_ete - norme6190_ete

ggplot(somme_ete_annuelle, aes(x = year)) +
  geom_bar(aes(y = ecart_norme, fill = ecart_norme >= 0), stat = "identity") +
  geom_smooth(aes(y = ecart_norme), method = "loess", color = "black", size = 0.5) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Écart positif", "Écart négatif")) +
  labs(title = "Écart des précipitations estivales par rapport à la norme (1961-1990) (NEU)",
       x = "Année",
       y = "Écart à la norme (mm)") +
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  scale_y_continuous(limits = c(-200, 350), breaks = seq(-200,350, by = 50)) +  
  theme_minimal() +
  theme(legend.position = "none")+
  annotate("text", x = 1864, y = 325, label = "norme estivale 61-90 = 274 mm", 
           hjust = 0, vjust = 1, size = 3, color = "black")

# Automne
norme6190_automne <- 82+68+81
automne_data <- NEU_sommeP %>% filter(month(date) %in% c(9, 10, 11))

somme_automne_annuelle <- automne_data %>%
  group_by(year = format(date, "%Y")) %>%
  summarise(total_precip_automne = sum(precip, na.rm = TRUE)) %>%
  ungroup()

somme_automne_annuelle$year <- as.numeric(somme_automne_annuelle$year)
somme_automne_annuelle$ecart_norme <- somme_automne_annuelle$total_precip_automne - norme6190_automne

ggplot(somme_automne_annuelle, aes(x = year)) +
  geom_bar(aes(y = ecart_norme, fill = ecart_norme >= 0), stat = "identity") +
  geom_smooth(aes(y = ecart_norme), method = "loess", color = "black", size = 0.5) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Écart positif", "Écart négatif")) +
  labs(title = "Écart des précipitations automnales par rapport à la norme (1961-1990) (NEU)",
       x = "Année",
       y = "Écart à la norme (mm)") +
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  scale_y_continuous(limits = c(-200, 370), breaks = seq(-200,350, by = 50)) +  
  theme_minimal() +
  theme(legend.position = "none")+
  annotate("text", x = 1864, y = 370, label = "norme automnale 61-90 = 231 mm", 
           hjust = 0, vjust = 1, size = 3, color = "black")

# Hiver
norme6190_hiver <- 79+71+69
hiver_data <- NEU_sommeP %>% filter(month(date) %in% c(12, 1, 2))

somme_hiver_annuelle <- hiver_data %>%
  group_by(year = format(date, "%Y")) %>%
  summarise(total_precip_hiver = sum(precip, na.rm = TRUE)) %>%
  ungroup()

somme_hiver_annuelle$year <- as.numeric(somme_hiver_annuelle$year)
somme_hiver_annuelle$ecart_norme <- somme_hiver_annuelle$total_precip_hiver - norme6190_hiver

ggplot(somme_hiver_annuelle, aes(x = year)) +
  geom_bar(aes(y = ecart_norme, fill = ecart_norme >= 0), stat = "identity") +
  geom_smooth(aes(y = ecart_norme), method = "loess", color = "black", size = 0.5) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Écart positif", "Écart négatif")) +
  labs(title = "Écart des précipitations hivernales par rapport à la norme (1961-1990) (NEU)",
       x = "Année",
       y = "Écart à la norme (mm)") +
  scale_x_continuous(breaks = c(1850,1875, 1900, 1925, 1950, 1975, 2000, 2025)) +  
  scale_y_continuous(limits = c(-200, 250), breaks = seq(-200,250, by = 50)) +  
  theme_minimal() +
  theme(legend.position = "none")+
  annotate("text", x = 1864, y = 245
           , label = "norme hivernale 61-90 = 219 mm", 
           hjust = 0, vjust = 1, size = 3, color = "black")

  

  
  
----------------------------------------------------------------------------------
  ----------------------------------------------------------------------------------
  #pas fait sur les autres stations 
  ----------------------------------------------------------------------------------
  ----------------------------------------------------------------------------------
  ##tests supplémentaires rapport climatologique 
  
  
#juin 23 le + sec depuis début des mesures ?
P_juin <- subset(NEU_sommeP, format(date, "%m") == "06")
sommeP_juin <- aggregate(precip ~ format(date, "%Y"), data = P_juin, sum)
colnames(sommeP_juin) <- c("Année", "Somme_Precipitations")
barplot(sommeP_juin$Somme_Precipitations,names.arg = sommeP_juin$Année,
        xlab = "Année", ylab = "Précipitations (mm)",main = "Somme des précipitations en juin (NEU)",
        col = "blue",las = 1) 
  #classement des juin les plus secs 
  classement_juin <- sommeP_juin[order(sommeP_juin$Somme_Precipitations), ]
  juin_secs <- head(classement_juin, 15)
  print(juin_secs)
  

#novembre 23 le + arrosé ?
P_novembre <- subset(NEU_sommeP, format(date, "%m") == "11")
sommeP_novembre <- aggregate(precip ~ format(date, "%Y"), data = P_novembre, sum)
colnames(sommeP_novembre) <- c("Année", "Somme_Precipitations")
classement_novembre <- sommeP_novembre[order(-sommeP_novembre$Somme_Precipitations), ]
barplot(sommeP_novembre$Somme_Precipitations,names.arg = sommeP_novembre$Année,
        xlab = "Année",ylab = "Précipitations (mm)",main = "Somme des précipitations en novembre (NEU)",
        col = "blue",las = 1)
  #classement des novembre les plus arrosés 
  classement_novembre <- head(classement_novembre, 10)
  print(classement_novembre)
 
  
#automne plus arrosé ?
automne_data <- subset(NEU_sommeP, format(date, "%m") %in% c("09", "10", "11"))
somme_precip_automne <- aggregate(precip ~ format(date, "%Y"), data = automne_data, sum)
colnames(somme_precip_automne) <- c("Année", "Somme_Precipitations")
print(somme_precip_automne)
classement_automne <- somme_precip_automne[order(-somme_precip_automne$Somme_Precipitations), ]
automne_plus_pluvieux <- head(classement_automne, 15)
print(automne_plus_pluvieux)
barplot(automne_plus_pluvieux$Somme_Precipitations,
        names.arg = automne_plus_pluvieux$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Automnes les plus pluvieux (NEU)",
        col = ifelse(automne_plus_pluvieux$Année == "2023", "red", "blue"),
        ylim = c(0, 600),
        las = 3)
abline(h = 231, col = "green", lty = 2)
abline(h = 241, col = "purple", lty = 2)
legend("topright", legend = c("Norme de 1961-1990", "Norme de 1991-2020"),
       col = c("green", "purple"), lty = 2, bty = "n")


#été sec ? 
ete_data <- subset(NEU_sommeP, format(date, "%m") %in% c("06", "07", "08"))
somme_precip_ete <- aggregate(precip ~ format(date, "%Y"), data = ete_data, sum)
colnames(somme_precip_ete) <- c("Année", "Somme_Precipitations")
print(somme_precip_ete)
classement_ete <- somme_precip_ete[order(somme_precip_ete$Somme_Precipitations), ]
ete_plus_sec <- head(classement_ete, 20)
print(ete_plus_sec)
barplot_heights_ete <- barplot(ete_plus_sec$Somme_Precipitations,
                               names.arg = ete_plus_sec$Année,
                               xlab = "Année",
                               ylab = "Précipitations (mm)",
                               main = "Étés les plus secs (NEU)",
                               col = ifelse(ete_plus_sec$Année == "2023", "red", "blue"),
                               ylim = c(0, 400),  
                               las = 3)  
abline(h = 274, col = "green", lty = 2)  
abline(h = 278, col = "purple", lty = 2)  
legend("topright", legend = c("Norme de 1961-1990", "Norme de 1991-2020"),
       col = c("green", "purple"), lty = 2, bty = "n") 


#nb jours >=20mm
jours_plusou_egal20mm <- aggregate(precip ~ format(date, "%Y"), data = NEU_sommeP[NEU_sommeP$precip >= 20, ], length)
colnames(jours_plusou_egal20mm) <- c("Année", "Nombre_Jours_20mm")
print(jours_plusou_egal20mm)
barplot(jours_plusou_egal20mm$Nombre_Jours_20mm,
        names.arg = jours_plusou_egal20mm$Année,
        xlab = "Année",
        ylab = "Nombre de Jours (≥ 20 mm)",
        main = "Nombre de jours où P ≥ 20 mm (NEU)",
        col = "blue",
        las = 1)


#nb jours >99th
percentile_99 <- quantile(NEU_sommeP$precip, 0.99)
print(percentile_99)
jours_plus_99 <- NEU_sommeP[NEU_sommeP$precip > percentile_99, ]
nombre_jours_99 <- aggregate(precip ~ format(date, "%Y"), data = jours_plus_99, length)
colnames(nombre_jours_99) <- c("Année", "Nombre_Jours_99th")
print(nombre_jours_99)
barplot(nombre_jours_99$Nombre_Jours_99th,
        names.arg = nombre_jours_99$Année,
        xlab = "Année",
        ylab = "Nombre de Jours (> 99e percentile)",
        main = "Fréquence des Jours avec Précipitations > 99e Percentile par Année",
        col = "blue",
        las = 1) 

  #test que sur j >1mm 
jours_pluie <- NEU_sommeP[NEU_sommeP$precip > 1, ]
percentile_99 <- quantile(jours_pluie$precip, 0.99)
print(percentile_99)
jours_plus_99 <- NEU_sommeP[NEU_sommeP$precip > percentile_99, ]
nombre_jours_99 <- aggregate(precip ~ format(date, "%Y"), data = jours_plus_99, length)
colnames(nombre_jours_99) <- c("Année", "Nombre_Jours_99th")
print(nombre_jours_99)
barplot(nombre_jours_99$Nombre_Jours_99th,
        names.arg = nombre_jours_99$Année,
        xlab = "Année",
        ylab = "Nombre de Jours (> 99e percentile)",
        main = "Nombre de jours avec P > 99e percentile",
        col = "blue",
        las = 1) 
  
  #test 95th 
percentile_95 <- quantile(jours_pluie$precip, 0.95)
print(percentile_95)
jours_plus_95 <- NEU_sommeP[NEU_sommeP$precip > percentile_95, ]
nombre_jours_95 <- aggregate(precip ~ format(date, "%Y"), data = jours_plus_95, length)
colnames(nombre_jours_95) <- c("Année", "Nombre_Jours_95th")
print(nombre_jours_95)
barplot(nombre_jours_95$Nombre_Jours_95th,
        names.arg = nombre_jours_95$Année,
        xlab = "Année",
        ylab = "Nombre de Jours (> 95e percentile)",
        main = "Fréquence des Jours avec Précipitations > 95e Percentile par Année",
        col = "blue",
        las = 1) 
  

#somme P pour l'année 2023 
precip_2023 <- NEU_sommeP[format(NEU_sommeP$date, "%Y") == "2023", ]
somme_precip_2023 <- sum(precip_2023$precip, na.rm = TRUE)
print(somme_precip_2023)

  
  
  
  
  
#somme de précipitations pour chaque année 
#pour voir si dans les 20 dernière années il y a eu + d'années pluvieuses 
somme_precip_annuelle <- aggregate(precip ~ format(date, "%Y"), data = NEU_sommeP, sum)
colnames(somme_precip_annuelle) <- c("Année", "Somme_Precipitations")
print(somme_precip_annuelle)
classement_annuel <- somme_precip_annuelle[order(-somme_precip_annuelle$Somme_Precipitations), ]
vingt_plus_pluvieuses <- head(classement_annuel, 20)
vingt_plus_seches <- tail(classement_annuel, 20)
print(vingt_plus_pluvieuses)
print(vingt_plus_seches)

  #pluvieuses 
couleurs <- ifelse(as.numeric(vingt_plus_pluvieuses$Année) < 1990, "deepskyblue", 
                   "blue")
barplot(vingt_plus_pluvieuses$Somme_Precipitations,
        names.arg = vingt_plus_pluvieuses$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Années les plus pluvieuses (NEU)",
        col = couleurs,
        las = 2.5, 
        ylim = c(0, 1600))

legend(x = 20, y = 1600, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue")) 
#sèches 
couleurs_seches <- ifelse(as.numeric(vingt_plus_seches$Année) < 1990, "deepskyblue", 
                          "blue")
barplot(vingt_plus_seches$Somme_Precipitations,
        names.arg = vingt_plus_seches$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Années les plus sèches (NEU)",
        col = couleurs_seches,
        las = 2.5, 
        ylim = c(0, 900)) 
legend(x = 20, y = 900, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))

#somme de précipitations pour chaque printemps
printemps_data <- subset(NEU_sommeP, format(date, "%m") %in% c("03", "04", "05"))
somme_precip_printemps_annuelle <- aggregate(precip ~ format(date, "%Y"), data = printemps_data, sum)
colnames(somme_precip_printemps_annuelle) <- c("Année", "Somme_Precipitations")
print(somme_precip_printemps_annuelle)
classement_printemps_annuel <- somme_precip_printemps_annuelle[order(-somme_precip_printemps_annuelle$Somme_Precipitations), ]
vingt_plus_pluvieuses_printemps <- head(classement_printemps_annuel, 20)
vingt_plus_seches_printemps <- tail(classement_printemps_annuel, 20)
print(vingt_plus_seches_printemps)
print(vingt_plus_pluvieuses_printemps)
  #pluvieuses 
couleurs_printemps <- ifelse(as.numeric(vingt_plus_pluvieuses_printemps$Année) 
                             < 1990, "deepskyblue", "blue") 
barplot(vingt_plus_pluvieuses_printemps$Somme_Precipitations,
        names.arg = vingt_plus_pluvieuses_printemps$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Printemps les plus pluvieux (NEU)",
        col = couleurs_printemps,
        las = 2.5, 
        ylim = c(0, 600))
legend(x = 20, y = 600, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))
  #sèches 
couleurs_printemps_sec <- ifelse(as.numeric(vingt_plus_seches_printemps$Année) < 1990, "deepskyblue", 
                                 "blue") 

barplot(vingt_plus_seches_printemps$Somme_Precipitations,
        names.arg = vingt_plus_seches_printemps$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Printemps les plus secs (NEU)",
        col = couleurs_printemps_sec,
        las = 2.5, 
        ylim = c(0, 150))

legend(x = 20, y = 150, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))


#somme de précipitations pour chaque été
ete_data <- subset(NEU_sommeP, format(date, "%m") %in% c("06", "07", "08"))
somme_precip_ete_annuelle <- aggregate(precip ~ format(date, "%Y"), data = ete_data, sum)
colnames(somme_precip_ete_annuelle) <- c("Année", "Somme_Precipitations")
print(somme_precip_ete_annuelle)
classement_ete_annuel <- somme_precip_ete_annuelle[order(-somme_precip_ete_annuelle$Somme_Precipitations), ]
vingt_plus_pluvieuses_ete <- head(classement_ete_annuel, 20)
vingt_plus_seches_ete <- tail(classement_ete_annuel, 20)
print(vingt_plus_pluvieuses_ete)
print(vingt_plus_seches_ete)
  #pluvieuses 
couleurs_ete <- ifelse(as.numeric(vingt_plus_pluvieuses_ete$Année) < 1990, "deepskyblue",  
                       "blue")  
barplot(vingt_plus_pluvieuses_ete$Somme_Precipitations,
        names.arg = vingt_plus_pluvieuses_ete$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Étés les plus pluvieux (NEU)",
        col = couleurs_ete,
        las = 2.5, 
        ylim = c(0, 700))
legend(x = 20, y = 700, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))
  #sèches 
couleurs_ete_sec <- ifelse(as.numeric(vingt_plus_seches_ete$Année) < 1990, "deepskyblue", 
                           "blue") 

barplot(vingt_plus_seches_ete$Somme_Precipitations,
        names.arg = vingt_plus_seches_ete$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Étés les plus secs (NEU)",
        col = couleurs_ete_sec,
        las = 2.5, 
        ylim = c(0, 200))

legend(x = 20, y = 200, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))


#somme P pour chaque automne 
automne_data <- subset(NEU_sommeP, format(date, "%m") %in% c("09", "10", "11"))
somme_precip_automne_annuelle <- aggregate(precip ~ format(date, "%Y"), data = automne_data, sum)
colnames(somme_precip_automne_annuelle) <- c("Année", "Somme_Precipitations")
print(somme_precip_automne_annuelle)
classement_automne_annuel <- somme_precip_automne_annuelle[order(-somme_precip_automne_annuelle$Somme_Precipitations), ]
vingt_plus_pluvieuses_automne <- head(classement_automne_annuel, 20)
vingt_plus_seches_automne <- tail(classement_automne_annuel, 20)
print(vingt_plus_pluvieuses_automne)
print(vingt_plus_seches_automne)

  #pluvieuses
couleurs_automne <- ifelse(as.numeric(vingt_plus_pluvieuses_automne$Année) < 1990, "deepskyblue",  
                           "blue")  
barplot(vingt_plus_pluvieuses_automne$Somme_Precipitations,
        names.arg = vingt_plus_pluvieuses_automne$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Automnes les plus pluvieux (NEU)",
        col = couleurs_automne,
        las = 2.5, 
        ylim = c(0, 600))

legend(x = 20, y = 600, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))

  # secs
couleurs_automne_sec <- ifelse(as.numeric(vingt_plus_seches_automne$Année) < 1990, "deepskyblue",  
                               "blue") 
barplot(vingt_plus_seches_automne$Somme_Precipitations,
        names.arg = vingt_plus_seches_automne$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Automnes les plus secs (NEU)",
        col = couleurs_automne_sec,
        las = 2.5, 
        ylim = c(0, 180))

legend(x = 20, y = 180, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))



#somme des P hiver 
hiver_data <- subset(NEU_sommeP, format(date, "%m") %in% c("12", "01", "02"))
somme_precip_hiver_annuelle <- aggregate(precip ~ format(date, "%Y"), data = hiver_data, sum)
colnames(somme_precip_hiver_annuelle) <- c("Année", "Somme_Precipitations")
print(somme_precip_hiver_annuelle)
classement_hiver_annuel <- somme_precip_hiver_annuelle[order(-somme_precip_hiver_annuelle$Somme_Precipitations), ]
vingt_plus_pluvieuses_hiver <- head(classement_hiver_annuel, 20)
vingt_plus_seches_hiver <- tail(classement_hiver_annuel, 20)
print(vingt_plus_pluvieuses_hiver)
print(vingt_plus_seches_hiver)

  # pluvieux
couleurs_hiver <- ifelse(as.numeric(vingt_plus_pluvieuses_hiver$Année) < 1990, "deepskyblue",  
                         "blue")  

barplot(vingt_plus_pluvieuses_hiver$Somme_Precipitations,
        names.arg = vingt_plus_pluvieuses_hiver$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Hivers les plus pluvieux (NEU)",
        col = couleurs_hiver,
        las = 2.5, 
        ylim = c(0, 500))

legend(x = 20, y = 500, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))

  #secs
couleurs_hiver_sec <- ifelse(as.numeric(vingt_plus_seches_hiver$Année) < 1990, "deepskyblue",  
                             "blue") 

barplot(vingt_plus_seches_hiver$Somme_Precipitations,
        names.arg = vingt_plus_seches_hiver$Année,
        xlab = "Année",
        ylab = "Précipitations (mm)",
        main = "Hivers les plus secs (NEU)",
        col = couleurs_hiver_sec,
        las = 2.5, 
        ylim = c(0, 150))

legend(x = 20, y = 150, legend = c("1864-1989", "1990-2023"),
       fill = c("deepskyblue", "blue"))











