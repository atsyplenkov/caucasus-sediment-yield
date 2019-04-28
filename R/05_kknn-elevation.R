################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part 5. KKNN with elevation
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

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Load
source("R/00_own-functions.R")
load("data/spatial/sy_10-caucasus-spatial.Rdata")

# 1) Download elevation -----------------------------------------------------
# For stations
sy %>% 
  get_aws_points() -> sy_h

sy_h <- sy_h[[1]]

# For grid
# may take couple of minutes
grid %>% 
  as("Spatial") %>% 
  get_aws_points() -> grid_h

grid_h <- grid_h[[1]]

grid_h <- st_as_sf(grid_h)

# 2) KKNN with elevation -----------------------------------------------------
set.seed(125)

sample_row <- sample(nrow(sy_h), 40)
hc_train <- data.frame(sy = sy_h$sy[-sample_row],
                       h = sy_h$elevation[-sample_row],
                       lon = st_coordinates(sy_h)[-sample_row, 1], 
                       lat = st_coordinates(sy_h)[-sample_row, 2])

hc_test <- data.frame(sy = sy_h$sy[sample_row],
                      h = sy_h$elevation[sample_row],
                      lon = st_coordinates(sy_h)[sample_row, 1], 
                      lat = st_coordinates(sy_h)[sample_row, 2])


# create empty result data frame
hc_result <- data.frame(sy = as.numeric(NA),
                        h = grid_h$elevation,
                        lon = st_coordinates(grid_h)[, 1], 
                        lat = st_coordinates(grid_h)[, 2])
# run KKNN
k <- 12
kernel <- "triangular"


hc_kknn <- kknn::kknn(sy ~ ., 
                      train = hc_train, 
                      test = hc_result, 
                      kernel = kernel,
                      # distance = 1.5,
                      k = k)

# bring back to result data frame
# only retain the probability of the dominant dialect at that grid cell
hc_result %<>%
  # extract the interpolated sy at each grid cell with the 
  # kknn::fitted function
  mutate(sy = fitted(hc_kknn))

hc_raster <- st_as_sf(hc_result, 
                      coords = c("lon", "lat"),
                      crs = projection(sy),
                      remove = F)
# 3) Validate ---------------------------------------------------------------
# Extract values
# Train
raster::rasterFromXYZ(dplyr::select(hc_result,
                                    x = lon,
                                    y = lat,
                                    sy),
                      crs = projection(sy_h)) %>% 
  raster::extract(., as(sy_h[-sample_row,], "Spatial")) -> hc_train$sy_pred

# Validate  
raster::rasterFromXYZ(dplyr::select(hc_result,
                                    x = lon,
                                    y = lat,
                                    sy),
                      crs = projection(sy_h)) %>% 
  raster::extract(., as(sy_h[sample_row,], "Spatial")) -> hc_test$sy_pred


bind_rows(hc_train %>% mutate(type = "train"),
          hc_test %>% mutate(type = "validate")) %>%
  ggplot(aes(y = sy_pred,
             x = sy)) +
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
  labs(x = expression(italic("log"[10])*"SSY"),
       y = expression(italic("log"[10])*"SSY"["мод"])) +
  expand_limits(x = 1, y = 1) +
  theme_clean() +
  facet_wrap(~type, labeller = labeller(
    type = c(train = "Модель",
             validate = "Валидация")
  )) -> caucasus_kknn_h_validate_graph

bind_rows(hc_train %>% mutate(type = "train"),
          hc_test %>% mutate(type = "validate")) %>% 
  mutate(diff = sy_pred - sy) %>% 
  group_by(type) %>% 
  summarise(RMSE = sqrt(mean(diff^2, na.rm = T)),
            ME = mean(diff, na.rm = T),
            NSE = hydroGOF::NSE(sy_pred, sy),
            KGE = hydroGOF::KGE(sy_pred, sy),
            R2 = hydroGOF::br2(sy_pred, sy),
            r = cor(sy_pred, sy,
                    method = "spearman",
                    use = "pairwise.complete.obs")) %>% 
  mutate_if(is.numeric, list(~signif(.,3))) -> caucasus_kknn_h_validate

# 4) Vizualise ---------------------------------------------------------------
# !
hc_raster$sy <- 10^hc_raster$sy

# Define quantiles and labels
no_classes <- 6
quantiles <- quantile(hc_raster$sy, 
                      probs = seq(0, 1, length.out = no_classes + 1))

# here I define custom labels (the default ones would be ugly)
labels <- prettyNum(round(as.vector(quantiles[-1])), big.mark = " ")

hc_raster$sy_class <- cut(hc_raster$sy,
                          breaks = quantiles,
                          labels = labels,
                          include.lowest = T)

ggplot() +
  # add raster geom for the knn result
  geom_raster(data = hc_raster,
              aes(x = lon, y = lat, fill = sy_class)) +
  annotate("text", x = -9610, y = 4710000,
           label = expression(italic("Черное море")),
           color = "gray80", size = 6, angle = -30) +
  # add caucasus shape
  geom_sf(data = caucasus, color = "gray60", fill = NA, size = .6) +
  # add coastline
  geom_sf(data = coast, colour = "gray60", fill = NA, size = .6) +
  # add training and validation points
  geom_point(data = bind_rows(hc_train %>% mutate(type = "Модель"),
                              hc_test %>% mutate(type = "Валидация")),
             aes(x = lon, y = lat, color = type),
             # color = "black",
             # fill = "white",
             # shape = 21,
             alpha = 0.8,
             size = 2) +
  # add fancy legend
  # https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
  scale_fill_viridis(option = "E",
                     # name = expression("SSY"*","*~t%.%km^"-2"),
                     name = expression("SSY"*","*~"т"%.%"км"^"-2"),
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
  ggsci::scale_color_d3(name = "Посты",
                        guide = guide_legend(title.position = 'top',
                                             title.hjust = 0.5,
                                             label.hjust = 1)) +
  theme_map() +
  theme(legend.position = "bottom") -> ssy_caucasus_kknn_h

# SAVE -----------------------------------------------------------------------
# Figures
ggsave("figures/3-15_ssy_caucasus_kknn_h.png",
       plot = ssy_caucasus_kknn_h,
       dpi = 500, w = 8, h = 6)

ggsave("figures/3-14_caucasus_kknn_h_validate_graph.png",
       plot = caucasus_kknn_h_validate_graph,
       dpi = 500, w = 6, h = 4)

# Tables
load("data/tidy/model_validate.Rdata")
save(caucasus_kknn_validate, caucasus_krige_validate,
     caucasus_kknn_h_validate,
     file = "data/tidy/model_validate.Rdata")

# Temporary data
save(hc_raster, sy_h,
     file = "data/spatial/kknn-h_raster.Rdata")
