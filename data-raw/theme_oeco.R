## code to prepare `theme_oeco` dataset goes here
library(ggplot2)
theme_oeco <- theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.line.x = element_line(size = 0.35, colour = 'grey50'),
        axis.line.y = element_line(size = 0.35, colour = 'grey50'),
        axis.ticks = element_line(size = 0.25, colour = 'grey50'),
        legend.justification = c(1, 0.75), legend.position = "right",
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        legend.text.align = 0, legend.background = element_blank(),
        plot.subtitle = element_text(size = 10, vjust = 0),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
        strip.placement = 'outside', panel.spacing.x = unit(-0.5, 'cm'))

usethis::use_data(theme_oeco, overwrite = TRUE)
