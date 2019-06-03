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
library(XLConnect)

Sys.setlocale("LC_ALL", "Russian_Russia")

# Delete all data in memory
rm(list = ls())

# Set working directory
source("R/00_own-functions.R")

# 1) Load database ------------------------------------------------------------
sy <- read.csv("data/raw/sy_caucasus.csv", header = T, sep = ",", encoding = "UTF-8")

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

# Convert Point Data into Geodata
sy %<>%
  st_as_sf(coords = c("Longitude", "Lattitude"),
           crs = "+proj=longlat +ellps=WGS84") %>% 
  # reproject to UTM
  st_transform(., 32638)

# 3) Load spatial data -------------------------------------------------------
caucasus <- sf::read_sf("data/raw/caucasus.shp") %>% 
  # reproject to UTM
  st_transform(., 32638)

# Intersect Caucasus and SY-points
sy <- st_intersection(sy, caucasus)

# Catchment area analysis
sy %>% 
  as_tibble() %>% 
  mutate(area_g = cut(Catchment_area,
                      breaks = c(0, 100, 1000,
                                 10000, Inf),
                      labels = c("<100", "100-1000", "1000-10000",
                                  ">10000"))) %>% 
  group_by(area_g) %>% 
  summarise(n = n()) -> sy_area_table

skimr::skim(sy)

sy %>% 
  as_tibble() %>% 
  group_by(Country) %>% 
  summarise(n = n(),
            Source = paste(unique(Source_data), collapse = ";")) -> sy_source_table

sy %>% 
  filter(Country == "Russia") %>% 
  get_years(.,"Measuring_period") %>% nrow()

# 2) MP length analysis -----------------------------------------------------
# Histogram of MP length
sy %>%
  ggplot(aes(x = MP_length )) +
  geom_histogram(fill = "#0188B7",
                 color = "white",
                 binwidth = 5) +
  geom_vline(aes(xintercept = median(MP_length , na.rm = T)),
             size = 1,
             linetype = "dashed") +
  scale_x_continuous(expand = c (0, 0),
                     breaks = seq(10, 90, 10)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = expression(italic(MP)~", years"),
       # title = "Количество постов",
       subtitle = paste0("Median length ",
                         median(sy$MP_length, na.rm = T),
                         " years, n = ",
                         tally(sy)),
       y = "Count") +
  theme_clean() -> sy_mp

# Plot year availibility
get_years <- function(z, MP = NULL){
  
  z %>% 
    as_tibble() %>% 
    separate(MP, c("year1", "year2", "year3"), ",") %>% 
    mutate(
      year1 = str_replace_all(year1, pattern = "-", replacement = ":"),
      year2 = str_replace_all(year2, pattern = "-", replacement = ":"),
      year3 = str_replace_all(year3, pattern = "-", replacement = ":")
    ) -> year
  
  year %>% 
    filter(str_detect(year1, "<")) %>% 
    dplyr::select(year1, MP_length) %>% 
    mutate(year1 = str_replace_all(year1, pattern = "<", replacement = ""),
           year1 = as.numeric(year1),
           year4 = year1 - MP_length,
           year4 = glue::glue("{year4}:{year1}")) %>%
    dplyr::select(year4) -> year4
  
  year %>% 
    filter(!is.na(year2)) %>% 
    dplyr::select(year2) -> year2
  
  year %>% 
    filter(!is.na(year3)) %>% 
    dplyr::select(year3) -> year3
  
  year %>% 
    filter(!str_detect(year1, "<")) %>% 
    dplyr::select(year1) -> year1
  
  year <- c(as.vector(t(year1)),
            as.vector(t(year2)),
            as.vector(t(year3)),
            as.vector(t(year4))) %>% 
    strsplit(., ":") %>% 
    lapply(as.numeric) %>% 
    lapply(function(x) {
      if (length(x) == 2) {
        seq(x[1], x[2], 1)
      } else {
        x
      }
    }) %>% 
    unlist() %>%
    tibble() %>% 
    magrittr::set_colnames(c("year"))
  
  return(year)
}

get_years(sy, "Measuring_period") %>%
  ggplot(aes(x = year)) +
  geom_histogram(fill = "#A91511",
                 color = "white",
                 binwidth = 5) +
  geom_vline(aes(xintercept = median(year , na.rm = T)),
             size = 1,
             linetype = "dashed") +
  scale_x_continuous(expand = c (0, 0),
                     breaks = seq(1925, 2015, 10)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Year",
       y = "Count") +
  theme_clean() -> sy_year

ggsave("figures/3-6_caucasus-sy_mp.png",
       plot = ggarrange(sy_mp,
                        sy_year,
                        nrow = 2,
                        labels = "AUTO"),
       dpi = 500,
       width = 5, height = 5)


# Проверка на нормальность ---------------------------------------------------
# Гистограмма
sy %>% 
  as_tibble() %>% 
  ggplot(aes(x = (sy))) +
  geom_histogram(binwidth = 200, color = "dimgrey", fill = "white") +
  labs(y = "Количество",
       x =  expression("SSY,"*~"т"%.%"км"^-2)) +
  ggpubr::theme_pubclean(base_family = "Ubuntu") -> sy_hist

qplot(sample = log10(sy),
      data = as_tibble(sy)) +
  stat_qq() +
  stat_qq_line() +
  labs(y = "Наблюдаемые квантили",
       x = "Квантили нормального распределения") +
  ggpubr::theme_pubclean(base_family = "Ubuntu") -> sy_qq

# Диаграмма размахов
sy %>%
  as_tibble() %>% 
  mutate(sy_log = log10(sy)) %>% 
  mutate(outlier.high = sy_log > quantile(sy_log, .75) + 1.50*IQR(sy_log),
         outlier.low = sy_log < quantile(sy_log, .25) - 1.50*IQR(sy_log),
         outlier.color = case_when(outlier.high ~ "red",
                                   outlier.low ~ "red",
                                   outlier.low == F | outlier.high == F ~ "black")) %>% 
  ungroup() %>% 
  ggplot(aes(x = "",
             y = sy)) +
  stat_boxplot(geom ='errorbar',
               width = .25) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = outlier.color),
              width = .1,
              alpha = .6,
              show.legend = F) +
  # stat_summary(fun.data = n_fun,
  #              geom = "text",
  #              hjust = .5) +
  scale_y_log10(labels = fancyNumbers,
                breaks = prettyLogs) +
  ggsci::scale_color_lancet() +
  xlab(label = NULL) +
  # ylab(expression(italic("log"[10])*"SSY,"*~"т"%.%"км"^-2)) +
  ylab(expression(italic(SSY)*", t"%.%"km"^"-2"%.%"year"^"-1")) +
  theme_clean() -> sy_box

sy %>% 
  elevatr::get_aws_points() -> tt

tt %>% 
  pluck(1) %>% 
  as_tibble() %>% 
  dplyr::select(Station_name, A = Catchment_area, h = elevation) %>% 
  mutate(A_log = log10(A)) %>% 
  mutate(outlier.high = A_log > quantile(A_log, .75) + 1.50*IQR(A_log),
         outlier.low = A_log < quantile(A_log, .25) - 1.50*IQR(A_log),
         outlier.color = case_when(outlier.high ~ "red",
                                   outlier.low ~ "red",
                                   outlier.low == F | outlier.high == F ~ "black")) %>% 
  ungroup() %>% 
  ggplot(aes(x = "",
             y = A)) +
  stat_boxplot(geom ='errorbar',
               width = .25) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = outlier.color),
              width = .1,
              alpha = .6,
              show.legend = F) +
  scale_y_log10(labels = fancyNumbers,
                breaks = prettyLogs) +
  ggsci::scale_color_lancet(name = "",
                            labels = c()) +
  xlab(label = NULL) +
  ylab(expression(italic(A)*", km"^"2")) +
  theme_clean() -> A_box

tt %>% 
  pluck(1) %>% 
  as_tibble() %>% 
  dplyr::select(Station_name, A = Catchment_area, h = elevation) %>% 
  mutate(h_log = (h)) %>%
  mutate(outlier.high = h_log > quantile(h_log, .75, na.rm = T) + 1.50*IQR(h_log, na.rm = T),
         outlier.low = h_log < quantile(h_log, .25, na.rm = T) - 1.50*IQR(h_log, na.rm = T),
         outlier.color = case_when(outlier.high ~ "red",
                                   outlier.low ~ "steelblue",
                                   outlier.low == F | outlier.high == F ~ "black")) %>%
  ungroup() %>%
  ggplot(aes(x = "",
             y = h)) +
  stat_boxplot(geom ='errorbar',
               width = .25) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = outlier.color),
              width = .1,
              alpha = .6,
              show.legend = F) +
  ggsci::scale_color_lancet() +
  scale_y_continuous(labels = function(...){prettyNum(..., big.mark = " ")}) +
  xlab(label = NULL) +
  ylab(expression(italic(H)*", m")) +
  theme_clean()  -> H_box

ggsave("figures/3-7_caucasus-sy_norm.png", 
       plot = ggpubr::ggarrange(sy_hist,
                                sy_qq,
                                sy_box,
                                ncol = 3,
                                labels = "AUTO"),
       dpi = 500, width = 11, height = 4)

ggsave("figures/3-7_caucasus-sy_explore.png", 
       plot = ggpubr::ggarrange(sy_box,
                                H_box,
                                A_box,
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

# SAVE --------------------------------------------------------------------
save(sy, file = "data/tidy/sy_10-caucasus.Rdata")

# Export to EXCEL
caucasus_book <- loadWorkbook("analysis/summary_caucasus.xlsx", create = T)

# Subset summary
createSheet(caucasus_book, "Database summary")
writeWorksheet(object = caucasus_book,
               data = sy_source_table,
               sheet = "Database summary")

saveWorkbook(object = caucasus_book, file = "analysis/summary_caucasus.xlsx")
