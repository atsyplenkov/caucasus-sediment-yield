################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part X. Caucasus SSY assessment
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
load("data/spatial/kknn-h_raster.Rdata")

# Calculate total SY of Caucasus
sy_down <- st_read("data/tidy/caucasus_sy_downstream.shp")


sy_down %>% 
  as_tibble() %>% 
  mutate(sy = 10^sy) %>% 
  summarise(sum(sy * Ctchmn_)) %>% pull() -> caucasus_total_sy

# Read hydrobasins
hydrobasins <- st_read("/WORK/00_GLOBAL/ASIA/HYDROBASINS/hybas_eu_lev04_v1c.shp")

hydrobasins %>% 
  st_buffer(.,0) %>% 
  st_intersection(., st_transform(caucasus, projection(hydrobasins))) %>% 
  st_transform(., projection(caucasus)) -> tt

st_write(tt, "data/raw/hydrobasins.shp", delete_dsn = T)
