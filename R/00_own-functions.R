# For boxplots:
# SOURCE: https://owi.usgs.gov/blog/boxplots/

n_fun <- function(x){
  return(data.frame(y = 0.99*log10(10^4),
                    label = paste0("n = ", length(x))))
}

prettyLogs <- function(x){
  pretty_range <- range(x[x > 0])
  pretty_logs <- 10^(-10:10)
  log_index <- which(pretty_logs < pretty_range[2] & 
                       pretty_logs > pretty_range[1])
  log_index <- c(log_index[1]-1,log_index, log_index[length(log_index)]+1)
  pretty_logs_new <-  pretty_logs[log_index] 
  return(pretty_logs_new)
}

fancyNumbers <- function(n){
  nNoNA <- n[!is.na(n)]
  x <-gsub(pattern = "1e",replacement = "10^",
           x = format(nNoNA, scientific = TRUE))
  exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))
  
  base <- ifelse(exponents == 0, "1", ifelse(exponents == 1, "10","10^"))
  exponents[base == "1" | base == "10"] <- ""
  textNums <- rep(NA, length(n))  
  textNums[!is.na(n)] <- paste0(base,exponents)
  
  textReturn <- parse(text=textNums)
  return(textReturn)
}

# Map theme
# SOURCE: https://timogrossenbacher.ch/2018/03/categorical-spatial-interpolation-with-r/
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family = "Ubuntu", color = "#22211d"),
      # remove all axes
      axis.line = element_blank(),
      # axis.text.x = element_blank(),
      # axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # add a subtle grid
      panel.grid.major = element_line(color = "#dbdbd6", size = 0.6),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#ffffff", color = NA), 
      plot.margin = unit(c(.5, .5, .2, .5), "cm"),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "#ffffff", color = NA), 
      panel.spacing = unit(c(-.1, 0.2, .2, 0.2), "cm"),
      legend.background = element_rect(fill = "#ffffff", color = NA),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11, hjust = 0, color = "#4e4d47"),
      plot.title = element_text(size = 16, hjust = 0.5, color = "#4e4d47"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#4e4d47", 
                                   margin = margin(b = -0.1, 
                                                   t = -0.1, 
                                                   l = 2, 
                                                   unit = "cm"), 
                                   debug = F),
      plot.caption = element_text(size = 9, 
                                  hjust = .5, 
                                  margin = margin(t = 0.2, 
                                                  b = 0, 
                                                  unit = "cm"), 
                                  color = "#939184"),
      ...
    )
}

# GGPLOT2 THEME -----------------------------------------------------------
theme_clean <- function(base_font_family = "Ubuntu",
                        base_font_size = 12,
                        legend = "bottom") {
  
  ggpubr::theme_pubclean(base_family = base_font_family,
                         base_size = base_font_size) +
    theme(
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = legend,
      strip.background = element_blank()
    )
}