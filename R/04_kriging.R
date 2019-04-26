################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part 4. Ordinary kriging
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
library(gstat)
library(automap)
library(extrafont)
library(viridis)
library(sabre)
library(rnaturalearth)

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Load
source("R/00_own-functions.R")
load("data/spatial/sy_10-caucasus-spatial.Rdata")

# 1) Sub-sample ----------------------------------------------------------------
# Create validation and modeling dataset 
set.seed(125)
sample_row <- sample(nrow(sy), 40)


sy_m <- data.frame(sy_log = sy$sy[-sample_row],
                   lon = st_coordinates(sy)[-sample_row, 1], 
                   lat = st_coordinates(sy)[-sample_row, 2])

sy_val <- data.frame(sy_log = sy$sy[sample_row],
                     lon = st_coordinates(sy)[sample_row, 1], 
                     lat = st_coordinates(sy)[sample_row, 2])

# Converting to spatial objects
coordinates(sy_m) <- ~ lon + lat
coordinates(sy_val) <- ~ lon + lat

projection(sy_m) <- projection(sy)
projection(sy_val) <- projection(sy)

grid <- as(grid, "Spatial")
# projection(grid) <- projection(sy)

# 2) Cross-validation ------------------------------------------------------------
# Cross-Validation and selection best fit model
exp_var <- autofitVariogram(sy_log ~ 1, input_data = sy_m, model = "Exp")
lin_var <- autofitVariogram(sy_log ~ 1, input_data = sy_m, model = "Lin")
ok_var <- autofitVariogram(sy_log ~ 1, input_data = sy_m, model = "Sph")
gau_var <- autofitVariogram(sy_log ~ 1, input_data = sy_m, model = "Gau")

# Function for model errors calculation
kriging_coval <- function(var_model){
  
  krige.cv(sy_log ~ 1, sy_m,
           model = var_model$var_model,
           nfold = nrow(sy_m),
           verbose = FALSE) %>% 
    as_tibble() %>% 
    summarise(ME = mean(residual),
              MAE = mean(abs(residual)),
              RMSE = sqrt(mean(residual^2)),
              MSDR = mean(residual^2 / var1.var))
}

# Diagnostics from cross-validation of variogrammodels
cbind(model = c("Exp", "Lin", "Sph", "Gau"),
      rbind(kriging_coval(exp_var),
            kriging_coval(lin_var),
            kriging_coval(ok_var),
            kriging_coval(gau_var)) %>% signif(3)) %>% 
  as_tibble() %>% 
  mutate(type = "Ordinary Kriging") -> krig_variograms

krig_variograms %<>%
  mutate_if(is.numeric, list(~signif(.,3)))

# Check the distribution of residuals
krige.cv(sy_log ~ 1, sy_m,
         model = ok_var$var_model,
         nfold = nrow(sy_m),
         verbose = FALSE) %>% 
  as_tibble() -> cv.o

cv.o %>% 
  ggplot(aes(x = residual)) +
  geom_histogram(fill = "white", color = "black",
                 binwidth = .15) +
  annotate("text", x= -.75, y = 30, hjust = 1,
           label = paste0("W = ",
                          signif(shapiro.test(as_tibble(cv.o)$residual)$statistic, 2),
                          "\n p-value = ",
                          signif(shapiro.test(as_tibble(cv.o)$residual)$p.value, 1))) +
  labs(x = "",
       y = "Количество") +
  ggpubr::theme_pubclean(base_family = "Ubuntu")

rm(exp_var, lin_var, gau_var)

# 3) Ordinary kriging -------------------------------------------------------
ok_auto <- autoKrige(sy_log ~ 1, input_data = sy_m,
                     new_data = grid, model = "Sph")

# Export Variogram
ok_var$exp_var %>% 
  ggplot(aes(x = dist, y = gamma)) +
  geom_point(size = 2) +
  geom_line(data = variogramLine(ok_var$var_model,
                                 maxdist = max(ok_var$exp_var$dist)),
            color = "dimgrey",
            linetype = "dashed",
            size = 1) +
  # Highlight range
  geom_vline(xintercept = ok_var$var_model[2,3],
             linetype = "dotted") +
  geom_curve(data = tibble(x = 3*10^5,
                           y = 0.11,
                           xend = ok_var$var_model[2,3],
                           yend = 0), 
             aes(x = x, y = y, xend = xend, yend = yend),
             curvature = -.3,
             color = "black", 
             arrow = arrow(type = "closed",
                           length = unit(.15, "cm"))) +
  annotate("text",
           x = 3*10^5,
           y = 0.11,
           vjust = -0.6,
           label = as.expression(bquote(italic("a") == .(round(ok_var$var_model[2,3])))),
           family = "Ubuntu") +
  # Highlight nugget
  geom_curve(data = tibble(x = 1*10^5,
                           y = 0.3,
                           xend = 0,
                           yend = ok_var$var_model[1,2]), 
             aes(x = x, y = y, xend = xend, yend = yend),
             curvature = .3,
             color = "black", 
             arrow = arrow(type = "closed",
                           length = unit(.15, "cm"))) +
  annotate("text",
           x = 1*10^5,
           y = .3,
           vjust = -0.6,
           label = as.expression(bquote(italic("c"[0]) == .(round(ok_var$var_model[1,2], 2)))),
           family = "Ubuntu") +
  labs(x = "Расстояние",
       y = "Вариация",
       caption = expression(italic("c"[0])*~"- самородок;"
                            *~italic("a")*~"- радиус влияния")) +
  expand_limits(y = 0,
                x = 0) +
  scale_y_continuous(expand = c(.005, .005)) +
  scale_x_continuous(expand = c(.005, .005)) +
  theme_clean() -> geom_variog

# 4) Validation ---------------------------------------------------------------
raster::rasterFromXYZ(as.data.frame(ok_auto$krige_output)[, 1:3]) %>% 
  raster::extract(., sy_m) -> sy_m$sy_pred

raster::rasterFromXYZ(as.data.frame(ok_auto$krige_output)[, 1:3]) %>% 
  raster::extract(., sy_val) -> sy_val$sy_pred

bind_rows(as_tibble(sy_m) %>% mutate(type = "train"),
          as_tibble(sy_val) %>% mutate(type = "validate")) %>%
  ggplot(aes(y = sy_pred,
             x = sy_log)) +
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
  )) -> caucasus_krige_validate_graph

bind_rows(as_tibble(sy_m) %>% mutate(type = "train"),
          as_tibble(sy_val) %>% mutate(type = "validate")) %>%
  mutate(diff = sy_pred - sy_log) %>% 
  group_by(type) %>% 
  summarise(RMSE = sqrt(mean(diff^2, na.rm = T)),
            ME = mean(diff, na.rm = T),
            NSE = hydroGOF::NSE(sy_pred, sy_log),
            KGE = hydroGOF::KGE(sy_pred, sy_log),
            R2 = hydroGOF::br2(sy_pred, sy_log),
            r = cor(sy_pred, sy_log,
                    method = "spearman",
                    use = "pairwise.complete.obs")) %>% 
  mutate_if(is.numeric, list(~signif(.,3))) -> caucasus_krige_validate

# 5) Plot --------------------------------------------------------------------
krig_result <- st_as_sf(ok_auto$krige_output, 
                        coords = c("coords.x1", "coords.x2"),
                        crs = projection(sy),
                        remove = F) 

krig_result %<>%
  mutate(sy = 10^var1.pred,
         lon = ok_auto$krige_output@coords[,1],
         lat = ok_auto$krige_output@coords[,2])


# Define quantiles and labels
no_classes <- 6
quantiles <- quantile(krig_result$sy, 
                      probs = seq(0, 1, length.out = no_classes + 1))

# here I define custom labels (the default ones would be ugly)
labels <- prettyNum(round(as.vector(quantiles[-1])), big.mark = " ")

krig_result$sy_class <- cut(krig_result$sy,
                            breaks = quantiles,
                            labels = labels,
                            include.lowest = T)

ggplot() +
  # add raster geom for the knn result
  geom_raster(data = krig_result,
              aes(x = lon, y = lat, fill = sy_class)) +
  annotate("text", x = -9610, y = 4710000,
           label = expression(italic("Черное море")),
           color = "gray80", size = 6, angle = -30) +
  # add coastline
  geom_sf(data = coast, colour = "gray60", fill = NA, size = .6) +
  # add training and validation points
  geom_point(data = bind_rows(as_tibble(sy_m) %>% mutate(type = "Модель"),
                              as_tibble(sy_val) %>% mutate(type = "Валидация")),
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
  theme(legend.position = "bottom") -> ssy_caucasus_krig

# SAVE -----------------------------------------------------------------------
# Figures
ggsave("figures/3-13_ssy_caucasus_krig.png",
       plot = ssy_caucasus_krig,
       dpi = 500, w = 8, h = 6)

ggsave("figures/3-12_caucasus_krige_validate_graph.png",
       plot = caucasus_krige_validate_graph,
       dpi = 500, w = 6, h = 4)

ggsave("figures/3-11_caucasus_krige_variogram.png",
       plot = geom_variog,
       dpi = 500, w = 7, h = 4)

# Tables
load("data/tidy/model_validate.Rdata")
save(caucasus_kknn_validate, caucasus_krige_validate,
     file = "data/tidy/model_validate.Rdata")

# Export to EXCEL
caucasus_book <- XLConnect::loadWorkbook("analysis/summary_caucasus.xlsx")
XLConnect::createSheet(caucasus_book, "Kriging variograms")
XLConnect::writeWorksheet(object = caucasus_book,
                          data = krig_variograms,
                          sheet = "Kriging variograms")
XLConnect::saveWorkbook(object = caucasus_book, file = "analysis/summary_caucasus.xlsx")

