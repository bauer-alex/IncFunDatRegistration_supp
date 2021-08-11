
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

library(fda)     # contains raw Berkeley data
library(mgcv)    # data smoothing
library(registr) # contains Berkeley data with simulated incompleteness
data("growth_incomplete")



# prepare the smoothed data before incompleteness was simulated -----------
dat_raw <- fda::growth

# transform to long dataset (and ensure appropriate sorting)
dat <- data.frame(id    = factor(rep(c(colnames(dat_raw$hgtm), colnames(dat_raw$hgtf)),
                                     each = length(dat_raw$age))),
                  index = rep(dat_raw$age, times = ncol(dat_raw$hgtm) + ncol(dat_raw$hgtf)),
                  value = c(as.vector(dat_raw$hgtm), as.vector(dat_raw$hgtf))) %>%
  arrange(id, index)

# slightly smooth the raw curves
smooth_list <- lapply(unique(dat$id), function(curve_id) {
  d = dat %>% filter(id == curve_id)
  m = mgcv::gam(value ~ s(index, bs = "cr", k = 15), data = d)
  d$value = unname(mgcv::predict.gam(m))
  return(d)
})
dat_smooth <- dplyr::bind_rows(smooth_list)

# take the first derivative
deriv_list <- lapply(levels(dat$id), function(curve_id) {
  d = dat_smooth %>% filter(id == curve_id)
  data.frame(id               = curve_id,
             index            = d$index[2:nrow(d)],
             value            = diff(d$value) / diff(d$index),
             stringsAsFactors = FALSE)
})
dat_deriv <- dplyr::bind_rows(deriv_list) %>%
  mutate(id = factor(id))



# prepare the plots -------------------------------------------------------
# spaghetti plot of first derivative with incompleteness
gg1 <- ggplot(growth_incomplete, aes(x = index, y = value, group = id)) +
  geom_line(alpha = .1, col = "black", lwd = 1) +
  # ggtitle("First derivative of incomplete growth curves") +
  xlab("t* [observed]") + ylab("Derivative") +
  scale_x_continuous(limits = c(0, 18)) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none",
        plot.title      = element_text(hjust = 0.5))

# lasagna plot of first derivative with incompleteness
ids <- levels(growth_incomplete$id)
growth_incomplete$id <- factor(growth_incomplete$id, levels = ids[order(sapply(ids, function(curve_id) max(growth_incomplete$index[growth_incomplete$id == curve_id])))])
gg3 <- growth_incomplete %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = index, y = id, col = value)) +
  geom_line(lwd = 2.5) +
  scale_color_continuous("Derivative", high = "midnightblue", low = "lightskyblue1") +
  xlab("t* [observed]") + ylab("curve") +
  # ggtitle("Observed curves\nwith simulated incompleteness") +
  theme_minimal(base_size = 30) +
  theme(panel.grid      = element_blank(),
        axis.text.y     = element_blank(),
        # plot.title      = element_text(hjust = 0.5),
        legend.position = "none")

# lasagna plot of curves without incompleteness in same ordering
dat_deriv_sorted <- dat_deriv
dat_deriv_sorted$id <- factor(dat_deriv$id, levels = ids[order(sapply(ids, function(curve_id) max(growth_incomplete$index[growth_incomplete$id == curve_id])))])
gg2 <- ggplot(dat_deriv_sorted, aes(x = index, y = id, col = value)) +
  geom_line(lwd = 2.5) +
  scale_color_continuous("Derivative", high = "midnightblue", low = "lightskyblue1") +
  xlab("t* [observed]") + ylab("curve") +
  # ggtitle("Observed curves\n") +
  theme_minimal(base_size = 30) +
  theme(panel.grid  = element_blank(),
        axis.text.y = element_blank())
        # plot.title  = element_text(hjust = 0.5))



# final plot --------------------------------------------------------------
# extract the legend from one plot to plot it individually
grobs    <- ggplotGrob(gg2)$grobs
gglegend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# hide legend in gg2
gg2 <- gg2 + theme(legend.position = "none")

# plot grid
cowplot::plot_grid(gglegend, gg2, gg3, gg1, nrow = 1, rel_widths = c(.1,.25,.25,.4))
# ggsave("../figures/A4_1_BerkeleyData.pdf", width = 20, height = 6)

