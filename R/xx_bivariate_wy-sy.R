################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part X. Bivariate map of KKNN with elevation
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
library(kknn)
library(extrafont)
library(viridis)
library(sabre)
library(rnaturalearth)
library(elevatr)
library(cowplot)

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Load
source("R/00_own-functions.R")
load("data/spatial/sy_10-caucasus-spatial.Rdata")
load("data/spatial/kknn-h_raster.Rdata")

# 1) KKNN with water yield
# Download elevation -----------------------------------------------------
# For stations
sy %>% 
  mutate(Runoff_average = case_when(
    Source_data %in% c("Abduev, 2015", "Eyubova, 2015") ~ Runoff_average * 10^3 / Catchment_area,
    TRUE ~ Runoff_average
  )) %>% 
  # select only data with water yield
  filter(!is.na(Runoff_average)) %>% 
  get_aws_points() -> sy_wy

sy_wy <- sy_wy[[1]]

sy_wy %<>%
  mutate(Runoff_average = log10(Runoff_average))

wy_train <- data.frame(wy = sy_wy$Runoff_average,
                       h = sy_wy$elevation,
                       lon = st_coordinates(sy_wy)[,1], 
                       lat = st_coordinates(sy_wy)[,2])


# create empty result data frame
wy_result <- data.frame(wy = as.numeric(NA),
                        h = hc_raster$h,
                        lon = st_coordinates(hc_raster)[, 1], 
                        lat = st_coordinates(hc_raster)[, 2])
# 2) run KKNN
k <- 10
kernel <- "triangular"


wy_kknn <- kknn::kknn(wy ~ ., 
                      train = wy_train, 
                      test = wy_result, 
                      kernel = kernel,
                      # distance = 1.5,
                      k = k)

# bring back to result data frame
# only retain the probability of the dominant dialect at that grid cell
wy_result %<>%
  # extract the interpolated sy at each grid cell with the 
  # kknn::fitted function
  mutate(wy = fitted(wy_kknn))

wy_raster <- st_as_sf(wy_result, 
                      coords = c("lon", "lat"),
                      crs = projection(sy),
                      remove = F)

# 3) Calibrate
raster::rasterFromXYZ(dplyr::select(wy_result,
                                    x = lon,
                                    y = lat,
                                    wy),
                      crs = projection(sy)) %>% 
  raster::extract(., as(sy_wy, "Spatial")) -> wy_train$wy_pred

wy_train %>% 
  filter(!is.na(wy_pred)) %>% 
  ggplot(aes(y = wy_pred,
             x = wy)) +
  geom_abline(intercept = 0, alpha = .8, linetype = "dotted") +
  geom_point(color = "black",
             fill = "gray60",
             shape = 21,
             size = 2) +
  geom_smooth(method = "lm",
              formula = y ~ x, 
              se = F,
              color = "black",
              linetype = "dashed") +
  ggpmisc::stat_poly_eq(aes(label =  paste(stat(adj.rr.label))),
                        formula = y ~ x, 
                        parse = TRUE) +
  expand_limits(x = 1, y = 1) +
  theme_clean()

# 4) Plot map
wy_raster$wy <- 10^wy_raster$wy

# Define quantiles and labels
no_classes <- 6
quantiles <- quantile(wy_raster$wy, 
                      probs = seq(0, 1, length.out = no_classes + 1))

# here I define custom labels (the default ones would be ugly)
labels <- prettyNum(round(as.vector(quantiles[-1])), big.mark = " ")

wy_raster$wy_class <- cut(wy_raster$wy,
                          breaks = quantiles,
                          labels = labels,
                          include.lowest = T)



ggplot() +
  # add raster geom for the knn result
  geom_raster(data = wy_raster,
              aes(x = lon, y = lat, fill = wy_class)) +
  annotate("text", x = -9610, y = 4710000,
           label = expression(italic("Черное море")),
           color = "gray80", size = 6, angle = -30) +
  # add caucasus shape
  geom_sf(data = caucasus, color = "gray60", fill = NA, size = .6) +
  # add coastline
  geom_sf(data = coast, colour = "gray60", fill = NA, size = .6) +
  # add fancy legend
  # https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
  scale_fill_viridis(option = "E",
                     # name = expression("SSY"*","*~t%.%km^"-2"),
                     name = expression("q"*","*~"л"%.%"с"%.%"км"^"-2"),
                     discrete = T,
                     direction = -1,
                     guide = guide_legend(
                       direction = "horizontal",
                       keyheight = unit(2, units = "mm"),
                       keywidth = unit(80 / length(labels), units = "mm"),
                       title.position = 'top',
                       title.hjust = 0.5,
                       label.hjust = 1,
                       nrow = 1,
                       byrow = T,
                       reverse = F,
                       label.position = "bottom"
                     )) +
  theme_map() +
  theme(legend.position = "bottom")


# Bivariate MAP --------------------------------------------------------------
hc_raster$wy <- wy_raster$wy

quantiles_sy <- hc_raster %>%
  pull(sy) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create 3 buckets for Elevation (h)
quantiles_wy <- hc_raster%>%
  pull(wy) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create color scale that encodes two variables
# red for SY and blue for Elevation (h)
# the special notation with gather is due to readibility reasons
bivariate_color_scale <- tibble(
  "3 - 3" = "#3F2949", # high inequality, high income
  "2 - 3" = "#435786",
  "1 - 3" = "#4885C1", # low inequality, high income
  "3 - 2" = "#77324C",
  "2 - 2" = "#806A8A", # medium inequality, medium income
  "1 - 2" = "#89A1C8",
  "3 - 1" = "#AE3A4E", # high inequality, low income
  "2 - 1" = "#BC7C8F",
  "1 - 1" = "#CABED0" # low inequality, low income
) %>%
  gather("group", "fill")

# cut into groups defined above and join fill
hc_raster %>%
  mutate(sy_quantiles = cut(sy,
                            breaks = quantiles_sy,
                            include.lowest = TRUE),
         wy_quantiles = cut(wy,
                            breaks = quantiles_wy,
                            include.lowest = TRUE),
         # by pasting the factors together as numbers we match the groups defined
         # in the tibble bivariate_color_scale
         group = paste(as.numeric(sy_quantiles),
                       "-",
                       as.numeric(wy_quantiles))) %>% 
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group") -> tt

ggplot() +
  # add raster geom for the knn result
  geom_raster(data = tt,
              aes(x = lon, y = lat, fill = fill)) +
  # add caucasus shape
  geom_sf(data = caucasus, color = "gray60", fill = NA, size = .6) +
  # add coastline
  geom_sf(data = coast, colour = "gray60", fill = NA, size = .6) +
  # add annotations
  annotate("text", x = -64055, y = 4786418,
           label = expression(italic("Черное море")),
           # label = expression(italic("Black Sea")),
           color = "gray80", size = 6, angle = -30) +
  # annotate("text", x = 919866, y = 4609282,
  #          label = expression(italic("Caspian Sea")),
  #          color = "gray80", size = 6, angle = -55) +
  # labs
  # labs(title = "Среднегодовой модуль стока взвешенных наносов Кавказа",
  #      subtitle = expression(SSY*","~"т"%.%"км"^"-2"),
  #      caption = "На основании наблюдений на гидрометеопостах с ≈1930х по 2015\nАнатолий Цыпленков\natsyplenkov@gmail.com") +
  # labs(title = "Caucasus average annual Suspended Sediment Yield",
  #      subtitle = expression(SSY*","~t%.%km^"-2"),
  #      caption = "Based on gauging stations data from ≈1930s to 2015\nAnatolii Tsyplenkov\natsyplenkov@gmail.com") +
  scale_fill_identity() +
  theme_map() -> map
