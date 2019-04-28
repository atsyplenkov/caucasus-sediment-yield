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

# 5) Bivariate map ------------------------------------------------------------
# https://timogrossenbacher.ch/2019/04/bivariate-maps-with-ggplot2-and-sf/

# create 3 buckets for SY
quantiles_sy <- hc_raster %>%
  pull(sy) %>%
  quantile(probs = seq(0, 1, length.out = 4))

# create 3 buckets for Elevation (h)
quantiles_h <- hc_raster%>%
  pull(h) %>%
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
         h_quantiles = cut(h,
                           breaks = quantiles_h,
                           include.lowest = TRUE),
         # by pasting the factors together as numbers we match the groups defined
         # in the tibble bivariate_color_scale
         group = paste(as.numeric(sy_quantiles),
                       "-",
                       as.numeric(h_quantiles))) %>% 
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
  labs(title = "Среднегодовой модуль стока взвешенных наносов Кавказа",
       subtitle = expression(SSY*","~"т"%.%"км"^"-2"),
       caption = "На основании наблюдений на гидрометеопостах с ≈1930х по 2015\nАнатолий Цыпленков\natsyplenkov@gmail.com") +
  # labs(title = "Caucasus average annual Suspended Sediment Yield",
  #      subtitle = expression(SSY*","~t%.%km^"-2"),
  #      caption = "Based on gauging stations data from ≈1930s to 2015\nAnatolii Tsyplenkov\natsyplenkov@gmail.com") +
  scale_fill_identity() +
  theme_map() -> map

# ANNOTATE
# Get coordinates manually
# http://epsg.io/map#srs=32638
annotations <- tibble(
  # label = c(
  #   "Grey areas mean\nlow elevation and\nlow SSY",
  #   "Blue areas mean\nhigh elevation and\nlow SSY",
  #   "Violet areas mean\nhigh elevation and\nhigh SSY",
  #   "Red areas mean\nlow elevation and\nhigh SSY"
  # ),
  label = c(
    "Серый означает\nмалые высоты и\nмалый сток наносов",
    "Синий означает\nбольшие высоты и\nмалый сток наносов",
    "Фиолетовый означает\nбольшие высоты и\nбольшой сток наносов",
    "Красный означает\nмалые высоты и\nбольшой сток\nнаносов"
  ),
  arrow_from = c(
    "732188,4200904", # grey
    "264049,4351782", # blue
    "267077,4963366", # violet
    "791291,4864285" # red
  ),
  arrow_to = c(
    "822521,4415178", # grey
    "534624,4443962", # blue
    "258729,4804781", # violet
    "597891,4803845" # red
  ),
  curvature = c(
    0.2, # grey
    -0.2, # blue
    -0.1, # violet
    0.2 # red
  ),
  nudge = c(
    "-3000, -2000", # grey
    "0, 5000", # blue
    "0, 15000", # violet
    "0, 3000" # red
  ),
  just = c(
    "1,0", # grey
    "1,1", # blue
    "0.5,0", # violet
    "0,1" # red
  )
) %>%
  separate(arrow_from, into = c("x", "y")) %>%
  separate(arrow_to, into = c("xend", "yend")) %>%
  separate(nudge, into = c("nudge_x", "nudge_y"), sep = "\\,") %>%
  separate(just, into = c("hjust", "vjust"), sep = "\\,")  

# add annotations one by one by walking over the annotations data frame
# this is necessary because we cannot define nudge_x, nudge_y and curvature
# in the aes in a data-driven way like as with x and y, for example
annotations %>%
  pwalk(function(...) {
    # collect all values in the row in a one-rowed data frame
    current <- tibble(...)
    
    # convert all columns from x to vjust to numeric
    # as pwalk apparently turns everything into a character (why???)
    current %<>%
      mutate_at(vars(x:vjust), as.numeric)
    
    # update the plot object with global assignment
    map <<- map +
      # for each annotation, add an arrow
      geom_curve(
        data = current,
        aes(
          x = x,
          xend = xend,
          y = y,
          yend = yend
        ),
        # that's the whole point of doing this loop:
        curvature = current %>% pull(curvature),
        size = 0.2,
        arrow = arrow(
          length = unit(0.005, "npc")
        )
      ) +
      # for each annotation, add a label
      geom_text(
        data = current,
        aes(
          x = x,
          y = y,
          label = label,
          hjust = hjust,
          vjust = vjust
        ),
        # that's the whole point of doing this loop:
        nudge_x = current %>% pull(nudge_x),
        nudge_y = current %>% pull(nudge_y),
        # other styles
        family = "Ubuntu",
        color = "dimgrey",
        size = 3
      )
  })


# Draw the legend
bivariate_color_scale %<>%
  separate(group, into = c("sy", "h"), sep = " - ") %>%
  mutate(sy = as.integer(sy),
         h = as.integer(h))

ggplot() +
  geom_tile(
    data = bivariate_color_scale,
    mapping = aes(
      x = sy,
      y = h,
      fill = fill)
  ) +
  scale_fill_identity() +
  labs(x = paste0("Сток наносов ⟶️\nот ",round(range(quantiles_sy))[1],
                  " до ", prettyNum(round(range(quantiles_sy))[2], big.mark = " "),
                  " т/км"),
       y = paste0("Высота ⟶️\nот 0 до 5642 мБС")) +
  # labs(x = "Higher SSY ⟶️",
  #      y = "Higher Elevation ⟶️") +
  theme_map() +
  # make font small enough
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8,
                                    angle = 90),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.background = element_blank()) +
  # quadratic tiles
  coord_fixed() -> legend

# Combine plots
ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.07, 0.075, 0.27, 0.27) -> ssy_caucasus_bivariate

# SAVE -----------------------------------------------------------------------
ggsave("figures/3-XX_ssy_caucasus_bivariate_ru.png",
       plot = ssy_caucasus_bivariate,
       dpi = 500, w = 8, h = 6)
