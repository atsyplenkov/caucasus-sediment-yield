################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part X. Classification comparison
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

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Load
source("R/00_own-functions.R")
load("data/spatial/sy_10-caucasus-spatial.Rdata")

# Read Gabrielyan [1971]
# gab <- st_read("/Users/ATSYPLENKOV/Dropbox/02_Cтатьи/2018/SY EUROPE/Mountains/publ/00_raw-data/gabrielyan_1971.shp") %>% 
#   # reproject to UTM
#   st_transform(., 32638) %>% 
#   # repair topology
#   st_buffer(., dist = 0)

# Mozzherin V.V., Sharifullin A.G. [2014]
mozh <- st_read("/Users/ATSYPLENKOV/Dropbox/02_Cтатьи/2018/SY EUROPE/Mountains/publ/00_raw-data/mozherin_2014.shp") %>% 
  # reproject to UTM
  st_transform(., 32638) %>% 
  # repair topology
  st_buffer(., dist = 0)



# 4) Classification ----------------------------------------------------------
hc_train <- data.frame(class = sy$class_m, 
                       lon = st_coordinates(sy)[, 1], 
                       lat = st_coordinates(sy)[, 2])

# create empty result data frame
hc_result <- data.frame(class = as.factor(NA), 
                        lon = st_coordinates(grid)[, 1], 
                        lat = st_coordinates(grid)[, 2])
# run KKNN
k <- 12
kernel <- "triangular"

hc_kknn <- kknn::kknn(class ~ ., 
                      train = hc_train, 
                      test = hc_result, 
                      kernel = kernel,
                      distance = 1.5,
                      k = k)

# bring back to result data frame
# only retain the probability of the dominant dialect at that grid cell
hc_result %<>%
  # extract the interpolated dialect at each grid cell with the 
  # kknn::fitted function
  mutate(class = fitted(hc_kknn),
         # only retain the probability of the interpolated dialect,
         # discard the other 7
         prob = apply(hc_kknn$prob,
                      1,
                      function(x) max(x)))

hc_raster <- st_as_sf(hc_result, 
                      coords = c("lon", "lat"),
                      crs = projection(sy),
                      remove = F)
# 5) Vizualise ---------------------------------------------------------------
# download countries borders
countries10 <- rnaturalearth::ne_download(scale = 10,
                                          "countries",
                                          returnclass = "sf")
countries10 %>% 
  st_cast(., "MULTILINESTRING") %>% 
  # crop to plot boundaries
  st_crop(., sf::st_bbox(st_transform(caucasus, 4326)) - .42) %>% 
  # reproject to UTM
  st_transform(., 32638) -> coast

# Specify plot labels
labels <- levels(hc_raster$class)

ggplot() +
  # add raster geom for the knn result
  geom_raster(data = hc_raster,
              aes(x = lon, y = lat, fill = class, alpha = prob)) +
  # add coastline
  geom_sf(data = coast, colour = "gray60", fill = NA, size = .6) +
  annotate("text", x = -9610, y = 4710000,
           label = expression(italic("Черное море")),
           color = "gray80", size = 6, angle = -30) +
  geom_sf(data = mozh, aes(color = "black"), size = .4, fill = NA) +
  # add training points
  geom_point(data = hc_train,
             aes(x = lon, y = lat),
             color = "black",
             fill = "white",
             shape = 21,
             alpha = 0.8,
             size = 2) +
  # remove propability legend
  scale_alpha_continuous(guide = F) +
  # add fancy legend
  # https://timogrossenbacher.ch/2016/12/beautiful-thematic-maps-with-ggplot2-only/
  scale_fill_viridis(option = "E",
                     # name = expression("Denudation rate"*","*~mm%.%y^"-1"),
                     name = expression("Годичный слой денудации"*","*~"мм"%.%"год"^"-1"),
                     discrete = T,
                     direction = -1,
                     guide = guide_legend(
                       direction = "horizontal",
                       keyheight = unit(2, units = "mm"),
                       keywidth = unit(80 / length(labels), units = "mm"),
                       title.position = 'top',
                       # title.hjust = 0.5,
                       label.hjust = 1,
                       nrow = 1,
                       byrow = T,
                       reverse = F,
                       label.position = "bottom"
                     )) +
  scale_color_manual(name = "",
                     values = "black",
                     labels = "границы областей равной денудации по \n [Мозжерин, Шарифуллин, 2014]",
                     guide = guide_legend(title.position = 'top',
                                          # title.hjust = 0.5,
                                          label.hjust = 1)) +
  theme_map() +
  theme(legend.position = "bottom") -> hc_caucasus_classes

# 6) Compare ------------------------------------------------------------------
# Reclass previous maps
# Gabrielyan [1971]
# gab %<>%
#   mutate(class = case_when(
#     id == 1 ~ "> 2000",
#     id == 2 ~ labels[7],
#     id == 3 ~ labels[6],
#     id == 4 ~ labels[5],
#     id == 5 ~ labels[4],
#     id == 6 ~ labels[3],
#     id == 7 ~ labels[2],
#     id == 8 ~ labels[1]
#   ))

# Mozzherin V.V., Sharifullin A.G. [2014]
mozh %<>% 
  mutate(class = case_when(
    id == 6 ~ labels[1],
    id == 7 ~ labels[2],
    id == 8 ~ labels[3],
    id == 9 ~ labels[4]
  ))

# Convert raster to polygons
coordinates(hc_result) <- ~lon+lat
projection(hc_result) <- projection(sy)
hc_polygons <- raster::rasterFromXYZ(hc_result[, 1])
projection(hc_polygons) <- projection(sy)
hc_polygons <- raster::rasterToPolygons(hc_polygons, dissolve = T)

# Convert to sf object
hc_polygons <- st_as_sf(hc_polygons) %>% 
  mutate(class = case_when(
    layer == 1 ~ labels[1],
    layer == 2 ~ labels[2],
    layer == 3 ~ labels[3],
    layer == 4 ~ labels[4]
  ))

# SABRE!
vmeasure_calc(x = hc_polygons, x_name = class,
              y = mozh, y_name = class)

hc_polygons$area <-  as.numeric(st_area(hc_polygons) / 10^6)
mozh$area <-  as.numeric(st_area(mozh) / 10^6)


mozh %>% 
  group_by(class) %>% 
  summarise(area = sum(area)) %>%
  mutate(area_p = area*100/sum(area)) %>% 
  as_tibble() %>% 
  left_join(., as_tibble(hc_polygons), by = "class") %>%  
  dplyr::select(class, area_p, area.y) %>% 
  mutate(area_hc = area.y*100/sum(area.y)) %>% 
  dplyr::select(-area.y) %>% 
  gather(layer, area, -class) %>% 
  ggplot(aes(x = class, y = area)) +
  geom_col(aes(fill = layer),
           color = "white",
           alpha = .8,
           position = "dodge") +
  geom_text(aes(label = glue::glue("{round(area)}%"), group = layer),
            position = position_dodge(0.9), vjust = -.3) +
  expand_limits(y = c(0,35)) +
  labs(x = expression("Темпы денудации"*","*~"мм"%.%"год"^"-1"),
       y = "Площадь, %") +
  scale_fill_manual(name = "",
                    values = c("#0288b7", "#a90010"),
                    labels = c("KKNN", "[Мозжерин, Шарифуллин, 2014]")) +
  scale_x_discrete(limits = labels) +
  theme_clean() -> class_area

# SAVE -----------------------------------------------------------------------
ggsave("figures/20_sy-caucasus-classes.png",
       plot = hc_caucasus_classes,
       dpi = 500, w = 8, h = 6)

ggsave("figures/21_sy-caucasus-classes-area.png",
       plot = class_area,
       dpi = 500, w = 7, h = 4)

# Summarize ------------------------------------------------------------------
caucasus_total_sy <- sum(sy$Catchment_area * 10^sy$sy)

sy %>% 
  as_tibble() %>% 
  mutate(sy = 10^sy) %>% 
  summarise(n = n(),
            mean = mean(sy),
            sd = sd(sy),
            med = median(sy),
            CV = sd/mean,
            tot = sum(sy*Catchment_area),
            tot_p = 100 * tot / caucasus_total_sy) %>% 
  mutate_if(is.numeric, list(~ signif(., 3))) %>% 
  mutate_if(is.numeric, list(~ prettyNum(., big.mark = " "))) -> sy_summary_tot

sy %>% 
  as_tibble() %>% 
  mutate(group = cut(Catchment_area,
                     breaks = c(0, 1000, Inf),
                     labels = c("<1000", ">1000")),
         sy = 10^sy) %>% 
  group_by(group) %>% 
  summarise(n = n(),
            mean = mean(sy),
            sd = sd(sy),
            med = median(sy),
            CV = sd/mean,
            tot = sum(sy*Catchment_area),
            tot_p = 100 * tot / caucasus_total_sy) %>% 
  mutate_if(is.numeric, list(~ signif(., 3))) %>% 
  mutate_if(is.numeric, list(~ prettyNum(., big.mark = " "))) -> sy_summary_group

# Save to workbook
caucasus_book <- loadWorkbook("analysis/summary_caucasus.xlsx")
createSheet(caucasus_book, "SSY summary")
writeWorksheet(object = caucasus_book,
               data = bind_rows(sy_summary_group, sy_summary_tot),
               sheet = "SSY summary")
saveWorkbook(object = caucasus_book, file = "analysis/summary_caucasus.xlsx")


sy %>% 
  as_tibble() %>% 
  mutate(group = cut(Catchment_area,
                     breaks = c(0, 1000, Inf),
                     labels = c("<1000", ">1000")),
         sy = 10^sy) %>% 
  group_by(group) %>% 
  ggplot(aes(x = group, y = sy)) +
  geom_boxplot(aes(fill = group)) +
  stat_boxplot(geom ='errorbar', width = 0.1) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
  labs(x = expression("Площадь бассейна, "*"км"^2),
       y = expression(italic(log[10])*SSY)) +
  scale_y_log10(breaks = prettyLogs,
                labels = fancyNumbers) +
  scale_fill_manual(name = "",
                    values = c("#0288b7", "#a90010"),
                    guide = F) +
  theme_clean()

sy %>% 
  as_tibble() %>% 
  mutate(group = cut(Catchment_area,
                     breaks = c(0, 1000, Inf),
                     labels = c("small", "big")),
         sy = 10^sy) %>% 
  kruskal.test(sy ~ group, data = .)

sy %>% 
  as_tibble() %>% 
  mutate(aspect = case_when(
    Country %in% c("Russia") ~ "north",
    TRUE ~ "south"
  ),
  sy = 10^sy) %>% 
  # filter(Catchment_area <= 1000) %>% 
  group_by(aspect) %>% 
  summarise(n = n(),
            mean = mean(sy),
            hc = mean/2.65/1000,
            sd = sd(sy),
            med = median(sy),
            CV = sd/mean,
            tot = sum(sy*Catchment_area),
            tot_p = 100 * tot / caucasus_total_sy) %>% 
  mutate_if(is.numeric, list(~ signif(., 3))) %>% 
  mutate_if(is.numeric, list(~ prettyNum(., big.mark = " ")))
