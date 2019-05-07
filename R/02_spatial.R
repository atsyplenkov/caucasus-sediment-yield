################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part 2. Spatial analysis
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
library(XLConnect)

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Load
source("R/00_own-functions.R")
load("data/tidy/sy_10-caucasus.Rdata")

# 2) Convert sy to hc --------------------------------------------------------
sy %<>% 
  mutate(hc = (sy / 2.65) / 1000,
         class = cut(hc,
                     breaks = c(0, 25, 50, 125, 250, 500, 1000, 2000),
                     labels = c("0-25", "25-50", "50-125", "125-250",
                                "250-500", "500-1000", "1000-2000")),
         class_m = cut(hc,
                       breaks = c(0, 0.025, 0.1, 0.25, Inf),
                       labels = c("<0.025", "0.025 - 0.1",
                                  "0.1 - 0.25", "> 0.25")),
         sy = log10(sy))

# Convert Point Data into Geodata
# sy %<>%
#   st_as_sf(coords = c("Longitude", "Lattitude"),
#            crs = "+proj=longlat +ellps=WGS84") %>% 
#   # reproject to UTM
#   st_transform(., 32638)

# 3) Load spatial data -------------------------------------------------------
caucasus <- sf::read_sf("data/raw/caucasus.shp") %>% 
  # reproject to UTM
  st_transform(., 32638)

# Create a fishnet of points
grid <- st_make_grid(caucasus, cellsize = 2000, what = "centers") %>% 
  st_intersection(., caucasus)

# Intersect Caucasus and SY-points
sy <- st_intersection(sy, caucasus)

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

# 4) Subset and explore
# Create validation and modeling dataset 
set.seed(125)
sample_row <- sample(nrow(sy), 40)


sy_m <- data.frame(sy = sy$sy[-sample_row],
                   lon = st_coordinates(sy)[-sample_row, 1], 
                   lat = st_coordinates(sy)[-sample_row, 2])

sy_val <- data.frame(sy = sy$sy[sample_row],
                      lon = st_coordinates(sy)[sample_row, 1], 
                      lat = st_coordinates(sy)[sample_row, 2])

# Compare statistics of subsample
rbind(
  sy_m %>% 
    as_tibble() %>% 
    mutate(sy = 10^sy) %>% 
    summarise(subset = "model",
              n = n(),
              min = min(sy) %>% signif(3),
              max = max(sy),
              mean = mean(sy) %>% signif(3),
              median = median(sy),
              sd = signif(sd(sy), 2),
              se = (sd(sy) / sqrt(length(sy))) %>% signif(3)),
  sy_val %>% 
    as_tibble() %>% 
    mutate(sy = 10^sy) %>% 
    summarise(subset = "val",
              n = n(),
              min = min(sy) %>% signif(3),
              max = max(sy),
              mean = mean(sy) %>% signif(3),
              median = median(sy),
              sd = signif(sd(sy), 2),
              se = (sd(sy) / sqrt(length(sy))) %>% signif(3))) %>% 
  mutate_if(is.numeric, list(~signif(.,3))) %>% 
  mutate_if(is.numeric, list(~prettyNum(., big.mark = " "))) -> sy_caucasus_subsamples

# Describe subsets
rbind(sy_val %>%
        as_tibble() %>% 
        mutate(type = "Test"),
      sy_m %>%
        as_tibble() %>% 
        mutate(type = "Train")) -> sy_clean

eda_boxplot <- ggplot(sy_clean, aes(x = type,
                                    y = 10^sy,
                                    fill = type)) +
  geom_boxplot() +
  stat_boxplot(geom ='errorbar',
               width = .15) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = 0.5) +
  scale_fill_manual(values = c("#0288b7", "#a90010"), guide = FALSE) + 
  scale_y_log10(labels = fancyNumbers,
                breaks = prettyLogs) +
  labs(x = NULL,
       y = expression(italic("log"[10])*"SSY")) +
  theme_clean()

eda_histogram <- ggplot(sy_clean,
                        aes(x = 10^sy,
                            fill = type)) +
  geom_histogram(bins = 11,
                 color = "white") +
  scale_fill_manual(values = c("#0288b7", "#a90010"), guide = FALSE) + 
  scale_x_log10(labels = fancyNumbers,
                breaks = prettyLogs) +
  # labs(y = "Количество",
  labs(y = "Count",
       x = expression(italic("log"[10])*"SSY")) +
  facet_wrap(~ type,
             nrow = 2,
             scales = "free_y") +
  theme_clean() +
  theme(panel.grid.major.x = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())

eda_ridges <- ggplot(sy_clean,
                     aes(x = 10^sy,
                         y = fct_rev(type),
                         fill = type)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE,
                                quantiles = 2,
                                scale = 3,
                                color = "white") + 
  scale_fill_manual(values = c("#0288b7", "#a90010"),
                    guide = FALSE) + 
  scale_x_log10(labels = fancyNumbers,
                breaks = prettyLogs) +
  labs(x = expression(italic("log"[10])*"SSY"),
       y = NULL) +
  theme_clean()


ggarrange(ggarrange(eda_boxplot, eda_histogram,
                    labels = c("A", "B")),
          eda_ridges,
          labels = c("", "C"),
          nrow = 2) -> sy_subsets

# Levene test
car::leveneTest(sy_clean$sy, as_factor(sy_clean$type))

# Save ----------------------------------------------------------------------
save(sy, grid, caucasus, coast,
     file = "data/spatial/sy_10-caucasus-spatial.Rdata")

ggsave("figures/3-08_caucasus-subset.png",
       plot = sy_subsets,
       dpi = 500,
       height = 5,
       width = 8)

# Export to EXCEL
caucasus_book <- loadWorkbook("analysis/summary_caucasus.xlsx", create = T)

# Subset summary
createSheet(caucasus_book, "Subset summary")
writeWorksheet(object = caucasus_book,
               data = sy_caucasus_subsamples,
               sheet = "Subset summary")

saveWorkbook(object = caucasus_book, file = "analysis/summary_caucasus.xlsx")
