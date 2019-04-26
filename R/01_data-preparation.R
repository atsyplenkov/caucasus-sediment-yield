################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part 1. Data preparation
#
# Anatolii Tsyplenkov, Valentin Golosov, Matthias Vanmaercke
# atsyplenkov@gmail.com
#                                                                                                    
################################################################################
library(tidyverse)
library(extrafont)
library(ggpubr)
library(outliers)
library(elevatr)
library(sf)
library(raster)

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Set working directory
source("R/00_own-functions.R")

# 1) Load database ------------------------------------------------------------
sy <- read.csv("data/raw/sy_caucasus.csv", header = T, sep = ",")

sy %>% 
  mutate(Country = as.character(Country),
         River_name = as.character(River_name),
         Station_name = as.character(Station_name),
         Catchment_area = as.numeric(Catchment_area),
         Measuring_period = as.character(Measuring_period),
         Source_data = as.character(Source_data),
         Method = as.character(Method)) %>% 
  rename(sy = SSY_average) %>% 
  as_tibble() -> sy

################
# Keep stations with mp > 10 years
sy %>% 
  mutate(MP_length = ifelse(is.na(MP_length), 10, MP_length)) %>% 
  filter(MP_length >= 10) -> sy

################

# Catchment area analysis
sy %>% 
  mutate(area_g = cut(Catchment_area,
                      breaks = c(0, 1000, 5000,
                                 10000, 15000, Inf),
                      labels = c("<1000", "1000-5000", "5000-10000",
                                 "10000-15000", ">15000"))) %>% 
  group_by(area_g) %>% 
  summarise(n = n())

# Histogram of MP length
sy %>% 
  ggplot(aes(x = MP_length )) +
  geom_histogram(fill = "#EFC000",
                 color = "black",
                 binwidth = 5,
                 alpha = .7) +
  geom_vline(aes(xintercept = median(MP_length , na.rm = T)),
             size = 1,
             linetype = "dashed") +
  scale_x_continuous(expand = c (0, 0),
                     breaks = seq(10, 90, 10)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Продолжительность наблюдения, лет",
       # title = "Количество постов",
       subtitle = paste0("Медианная продолжительность ",
                         median(sy$MP_length, na.rm = T),
                         " лет, n = ",
                         tally(sy)),
       y = "Кол-во постов") +
  theme_pubr(base_family = "Ubuntu") -> sy_mp

ggsave("figures/3-6_caucasus-sy_mp.png", dpi = 500,
       width = 9, height = 4)


# Проверка на нормальность
# Гистограмма
sy %>% 
  ggplot(aes(x = (sy))) +
  geom_histogram(binwidth = 200, color = "dimgrey", fill = "white") +
  labs(y = "Количество",
       x =  expression("SSY,"*~"т"%.%"км"^-2)) +
  ggpubr::theme_pubclean(base_family = "Ubuntu") -> sy_hist

qplot(sample = log10(sy),
      data = sy) +
  stat_qq() +
  stat_qq_line() +
  labs(y = "Наблюдаемые квантили",
       x = "Квантили нормального распределения") +
  ggpubr::theme_pubclean(base_family = "Ubuntu") -> sy_qq

# Диаграмма размахов
sy %>%
  mutate(sy_log = log10(sy)) %>% 
  mutate(outlier.high = sy_log > quantile(sy_log, .75) + 1.50*IQR(sy_log),
         outlier.low = sy_log < quantile(sy_log, .25) - 1.50*IQR(sy_log),
         outlier.color = case_when(outlier.high ~ "red",
                                   outlier.low ~ "steelblue",
                                   outlier.low == F | outlier.high == F ~ "black")) %>% 
  ungroup() %>% 
  ggplot(aes(x = "",
             y = (sy))) +
  stat_boxplot(geom ='errorbar',
               width = .25) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = outlier.color),
              width = .1,
              alpha = .6,
              show.legend = F) +
  stat_summary(fun.data = n_fun,
               geom = "text",
               hjust = .5) +
  scale_y_log10(labels = fancyNumbers,
                breaks = prettyLogs) +
  ggsci::scale_color_lancet() +
  xlab(label = NULL) +
  ylab(expression(italic("log"[10])*"SSY,"*~"т"%.%"км"^-2)) +
  ggpubr::theme_pubclean(base_family = "Ubuntu") -> sy_box

ggsave("figures/3-7_caucasus-sy_norm.png", 
       plot = ggpubr::ggarrange(sy_hist,
                                sy_qq,
                                sy_box,
                                ncol = 3,
                                labels = "AUTO"),
       dpi = 500, width = 11, height = 4)

rm(sy_box, sy_hist, sy_qq)

# Remove outliers based on Grubbs test
# Note that the p-value threshold value is equal to 0.5
# This made for more strict removing
while (grubbs.test(log10(sy$sy))$p.value < 0.5) {
  if (grepl("lowest.*", grubbs.test(log10(sy$sy))$alternative) == T) {
    sy <- sy[-which.min(sy$sy),]  
  } else {
    sy <- sy[-which.max(sy$sy),]
  }
}

# Normality test
shapiro.test(log10(sy$sy))

# Save data
save(sy, file = "data/tidy/sy_10-caucasus.Rdata")
