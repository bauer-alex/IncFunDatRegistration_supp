
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

library(registr) # our registration approach
library(fdasrvf) # SRVF approach of Srivastava et al. (2011)



# simulate some curves ----------------------------------------------------
index_grid <- seq(0, 1, length.out = 100)
dat_long   <- data.frame(index = rep(index_grid, times = 3),
                         value = c(dnorm(index_grid, mean = 0.4, sd = 0.09),
                                   dnorm(index_grid, mean = 0.57, sd = 0.07),
                                   dnorm(index_grid, mean = 0.93, sd = 0.15)),
                         id    = factor(rep(1:3, each = length(index_grid))))

gg1 <- ggplot(dat_long, aes(index, value, group = id, col = id)) +
  geom_line(size = 1.2) +
  # ggtitle("Observed curves") +
  xlab("t* [observed]") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5),
        axis.text       = element_blank(),
        axis.title.y    = element_blank())



# complete curve registration with fdasrvf --------------------------------
curve_matrix <- matrix(dat_long$value,
                       ncol = length(unique(dat_long$id)), byrow = F)

reg_srvf <- time_warping(f = curve_matrix, time = index_grid, method = "mean",
                         MaxItr = 100, showplot = F)

# plot registered curves
reg_dat <- data.frame(id    = factor(rep(1:3, each = nrow(reg_srvf$fn))),
                      index = rep(index_grid, times = 3),
                      value = as.vector(reg_srvf$fn))
gg_srvf1 <- ggplot(reg_dat, aes(index, value, col = id)) +
  geom_line(size = 1.2) +
  ggtitle("Complete curve registration") + xlab("t [registered]") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5),
        axis.text       = element_blank(),
        axis.title.y    = element_blank())

# plot warping functions
warp_dat <- dat_long %>%
  mutate(index_reg = as.vector(reg_srvf$gam))
gg_srvf2 <- ggplot(warp_dat, aes(index_reg, index, col = id)) +
  geom_abline(intercept = 0, slope = 1, col = "gray90") +
  geom_line(size = 1.2) +
  ggtitle("Complete curve registration") + xlab("t* [observed]") + ylab("t [registered]") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5),
        axis.text       = element_blank())



# incomplete curve registration with registr ------------------------------
reg_varEM <- registr(Y = dat_long, incompleteness = "trailing", lambda = 0)

# plot registered curves
gg_varEM1 <- ggplot(reg_varEM$Y, aes(index, value, col = id)) +
  geom_line(size = 1.2) +
  ggtitle("Incomplete curve registration") + xlab("t [registered]") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5),
        axis.text       = element_blank(),
        axis.title.y    = element_blank())

# plot warping functions
gg_varEM2 <- ggplot(reg_varEM$Y, aes(tstar, index, col = id)) +
  geom_abline(intercept = 0, slope = 1, col = "gray90") +
  geom_line(size = 1.2) +
  ggtitle("Incomplete curve registration") + xlab("t* [observed]") + ylab("t [registered]") +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5),
        axis.text       = element_blank())



# main plot ---------------------------------------------------------------
gg_empty <- ggplot() + theme_minimal(base_size = 30) + theme(plot.title = element_text(hjust = 0.5))

plot_grid(plot_grid(gg_empty + ggtitle("Observed curves"),
                    gg_empty + ggtitle("Registered curves"),
                    gg_empty + ggtitle("Warping functions"),
                    nrow = 1),
          plot_grid(plot_grid(gg_empty, gg1, gg_empty,
                              ncol = 1, rel_heights = c(.25, .5, .25)),
                    plot_grid(gg_srvf1,  gg_srvf2,
                              gg_varEM1, gg_varEM2,
                              nrow = 2, ncol = 2, byrow = TRUE),
                    ncol = 2, rel_widths = c(1/3, 2/3)),
          nrow = 2, rel_heights = c(.1, .9))
# ggsave("../figures/1_incompleteness.pdf", width = 24, height = 10)
