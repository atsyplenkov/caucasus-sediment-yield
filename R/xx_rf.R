################################################################################
#                                                                                                        
# Caucasus Sediment Yield
# Part xx. RF
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
library(ranger)

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Load
source("R/00_own-functions.R")
load("data/spatial/sy_10-caucasus-spatial.Rdata")

# 1) Random forest --------------------------------------------------------------
# Create validation and modeling dataset 
set.seed(125)
sample_row <- sample(nrow(sy), 40)


rf_train <- sy[-sample_row,] %>% 
  as_Spatial()

rf_test <- sy[sample_row,] %>% 
  as_Spatial()

# create pixel dataframe
grid %>% 
  as_Spatial() %>%
  sp::SpatialPixelsDataFrame(., data = as_tibble(.)) -> rf_pixel


# Derive buffer distance
rf_train <- rf_train[-c(134, 136),]
buff_dist <- GSIF::buffer.dist(rf_train["sy"],
                               rf_pixel,
                               as.factor(1:nrow(rf_train)))

# The value of the target variable SSY can be now modeled as
# a function of buffer distances:

dn0 <- paste(names(buff_dist), collapse="+")
fm0 <- as.formula(paste("sy ~ ", dn0))

# First we overlay points and grids to create a regression matrix
ov.sy <- over(rf_train["sy"], buff_dist)
rm.sy <- cbind(rf_train@data["sy"], ov.sy)

# to estimate also the prediction error variance i.e.
# prediction intervals we set quantreg=TRUE which initiates
# the Quantile Regression RF approach (Meinshausen, 2006):
m.sy <- ranger(fm0, rm.sy,
               num.trees = 600,
               splitrule = "extratrees",
               min.node.size = 7,
               num.random.splits = 3,
               quantreg = T,
               seed = 1)

hydroGOF::gof(predictions(m.sy), rf_train$sy)

tt <- ranger(sy ~ .,
             hc_train,
             num.trees = 600,
             quantreg = T,
             splitrule = "extratrees",
             num.random.splits = 3, seed = 1)

tt
