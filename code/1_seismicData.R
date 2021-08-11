
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

dat_list <- readRDS("../data/dat_seismic_introFigure.rds")



# plot the observed curves ------------------------------------------------
gg1 <- ggplot(dat_list$data_obs, aes(index, value, group = id, col = hypocentral_distance)) +
  geom_line(size = 1.2) +
  scale_color_manual("hypocentral\ndistance", values = c("#002F6D","#3A6FB3","#9ECAE1")) +
  # ggtitle("Observed curves\n") +
  xlab(expression(t[]^{"*"}~"[observed]")) + ylab("abs. ground velocity [m/s]") +
  # ylab("absolute ground velocity") +
  theme_minimal(base_size = 30) +
  theme(legend.position    = "none",
        # plot.title         = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.title.y       = element_text(hjust = 1))



# plot the curves without leading zeros -----------------------------------
gg2 <- ggplot(dat_list$data_cut, aes(index, value, group = id, col = hypocentral_distance)) +
  geom_line(size = 1.2) +
  scale_color_manual("hypocentral\ndistance", values = c("#002F6D","#3A6FB3","#9ECAE1")) +
  # ggtitle("Observed curves\nwithout leading zeros") +
  xlab(expression(t[0]^{"*"}~"[observed]")) + ylab("abs. ground velocity [m/s]") +
  xlim(c(0,30)) +
  # ylab("absolute ground velocity") +
  theme_minimal(base_size = 30) +
  theme(legend.position    = "right",
        # plot.title         = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.title.y       = element_text(hjust = 1))



# final plot --------------------------------------------------------------
# extract the legend from one plot to plot it individually
grobs    <- ggplotGrob(gg2)$grobs
gglegend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# hide legend in gg2
gg2 <- gg2 + theme(legend.position = "none")

# plot grid
cowplot::plot_grid(gg1, gg2, gglegend, nrow = 1, rel_widths = c(0.4,0.4,0.2))
# ggsave("../figures/1_seismicData.pdf", width = 20, height = 5)
