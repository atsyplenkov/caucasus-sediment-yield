################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part 7. Select best model
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
library(rnaturalearth)
library(elevatr)
library(gstat)
library(automap)
library(purrr)

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Load
source("R/00_own-functions.R")
load("data/spatial/sy_10-caucasus-spatial.Rdata")

# Load models results
load("data/tidy/model_validate.Rdata")

# 1) Choose the best model
rbind(
  caucasus_kknn_validate %>% mutate(model = "KNN"),
  caucasus_krige_validate %>% mutate(model = "Кригинг"),
  caucasus_kknn_h_validate %>% mutate(model = "KNN-h"),
  caucasus_cokrige_validate %>% mutate(model = "Ко-кригинг")
) %>% 
  mutate(type = case_when(
    type == "train" ~ "модель",
    TRUE ~ "валидация"
  )) -> caucasus_ssy_models

caucasus_ssy_models %<>%
  dplyr::select(8, 1:ncol(.))

# Export to EXCEL
caucasus_book <- XLConnect::loadWorkbook("analysis/summary_caucasus.xlsx")
XLConnect::createSheet(caucasus_book, "Best model")
XLConnect::writeWorksheet(object = caucasus_book,
                          data = caucasus_ssy_models,
                          sheet = "Best model")
XLConnect::saveWorkbook(object = caucasus_book, file = "analysis/summary_caucasus.xlsx")


# 2) Graphically
caucasus_ssy_models %>% 
  dplyr::select(-ME) %>% 
  gather(coef, value, -type, -model) %>%
  ggplot(aes(y = value, x = model, fill = type)) +
  geom_col(position = position_dodge2()) +
  theme_clean() +
  coord_flip() +
  facet_wrap(~coef, scales = "free_x")
