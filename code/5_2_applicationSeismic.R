
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

library(registr) # our accompanying package
library(mgcv)    # for smoothing

library(IncRegHelpers) # helper functions



# data preparation --------------------------------------------------------
dat <- readRDS("../data/dat_seismic.rds")


# model estimation --------------------------------------------------------
# load the pre-estimated model
res <- readRDS("../results/seismicApplication_finalModel.rds")
# code for model estimation (takes ~3.5 hours on 10 cores)
# res <- register_fpca(Y              = dat,
#                      family         = "gamma",
#                      npc_criterion  = c(.95,.02),
#                      Kt             = 10,
#                      Kh             = 4,
#                      max_iterations = 100,
#                      fpca_type      = "two-step",
#                      gradient       = FALSE,
#                      incompleteness = "trailing",
#                      cores          = 10,
#                      lambda_inc     = 0.004)
# saveRDS(res, file = "../results/seismicApplication_finalModel.rds")

# flip the first FPC by multiplication with -1, just cause it's nicer in the paper
res$fpca_obj$efunctions[,1] <- -1 * res$fpca_obj$efunctions[,1]
res$fpca_obj$scores[,1]     <- -1 * res$fpca_obj$scores[,1]

# add the estimated warping functions to the data, based on the observed 'index'
# identical(dat$index, res$Y$tstar) # -> TRUE
dat$t_hat <- res$Y$t_hat

# plot all inverse warping functions
ggplot(dat, aes(x = index, y = t_hat, group = id)) + geom_line(alpha = 0.05)



# main plot 1: observed + registered curves + FPCs ------------------------
# Observed curves, spaghetti plot
ylim_obs <- range(dat$value)
gg1 <- dat %>%
  plot_spaghetti(log10_scale = TRUE, y_breaks = c(.1,1,10), y_minorBreaks = NULL,
                 title = "Observed curves", xlab = expression(t[0]^{"*"}~"[observed]"),
                 ylim = ylim_obs, ylab = "ground velocity [m/s]\non log10 scale", alpha = .01) +
  geom_smooth(data = dat, aes(x = index, y = value), method = "gam",
              method.args = list(family = Tweedie(p = 2)), col = "dodgerblue3", se = FALSE)

# Registered curves, spaghetti plot
dat_regMean <- data.frame(index = res$fpca_obj$t_vec,
                          value = exp(as.vector(res$fpca_obj$mu)))
gg2 <- res$Y %>%
  plot_spaghetti(log10_scale = TRUE, y_breaks = c(.1,1,10), y_minorBreaks = NULL,
                 title = "Registered curves", index_var = "t_hat", xlab = expression(t[0]^{""}~"[registered]"),
                 ylim = ylim_obs, alpha = .01, hide_yLabels = TRUE) +
  geom_line(data = dat_regMean, aes(x = index, y = value), col = "dodgerblue3", size = 1.1)

# Represented curves
gg3      <- plot_representedCurves_spaghetti(fpca_obj = res$fpca_obj, exp_transformation = TRUE, log10_scale = TRUE,
                                             y_breaks = c(.1, 1, 10), y_minorBreaks = NULL, ylim = ylim_obs,
                                             title = "Represented curves", xlab = expression(t[0]^{""}~"[registered]"),
                                             ylab = "ground velocity [m/s]\non log10 scale", alpha = .01, hide_yLabels = TRUE)

# Estimated FPCs
ylim_FPCs <- c(0.1,6)
gg4       <- plot_FPC(res, FPC_index = 1, log10_scale = TRUE, title = "1.FPC", ylim = ylim_FPCs, xlab = expression(t[0]^{""}~"[registered]"), ylab = "ground velocity [m/s]\non log10 scale")
gg5       <- plot_FPC(res, FPC_index = 2, log10_scale = TRUE, title = "2.FPC", ylim = ylim_FPCs, xlab = expression(t[0]^{""}~"[registered]"), hide_yLabels = TRUE)

cowplot::plot_grid(gg1, gg2, gg3, gg4, gg5, nrow = 1, rel_widths = c(.28,.22,.22,.28,.22))
# ggsave("../figures/5_2_FPCs.png", width = 14.5, height = 3)



# lasagna plots of observed and registered curves -------------------------
xlim <- range(c(dat$index, res$Y$t_hat))

ids_newLevels <- dat %>% group_by(id) %>% summarize(max_value = max(value)) %>% 
  arrange(desc(max_value)) %>% pull(id) %>% as.character()

legend_limits <- log10(range(dat$value))

gg3_obs  <- dat %>% mutate(id = factor(id, levels = ids_newLevels)) %>% 
  plot_lasagna(log10_scale = TRUE, title = "Observed curves", xlab = expression(t[0]^{"*"}~"[observed]"), xlim = xlim,
               legend.position = "left", legend_limits = legend_limits, ylab = "curve", legend_title = "log10\nvalue") +
  scale_color_viridis_c("log10\nvalue", limits = legend_limits)
gg3_reg  <- res$Y %>% mutate(id = factor(id, levels = ids_newLevels)) %>% 
  plot_lasagna(log10_scale = TRUE, title = "Registered curves", index_var = "t_hat", xlab = expression(t[0]^{"*"}~"[registered]"),
               xlim = xlim, ylab = NULL, legend.position = "none", legend_limits = legend_limits) +
  scale_color_viridis_c(limits = legend_limits)
gg3_repr <- plot_representedCurves_lasagna(fpca_obj = res$fpca_obj, ids_sorted = ids_newLevels, 
                                           exp_transformation = TRUE,
                                           log10_scale = TRUE, title = "Represented curves",
                                           xlab = expression(t[0]^{"*"}~"[registered]"), xlim = xlim, ylab = NULL,
                                           legend.position = "none", legend_limits = legend_limits) +
  scale_color_viridis_c(limits = legend_limits)

# Estimated FPCs
ylim_FPCs <- c(0.1,6)
gg4       <- plot_FPC(res, FPC_index = 1, log10_scale = TRUE, title = "1.FPC", ylim = ylim_FPCs, xlab = expression(t[0]^{""}~"[registered]"), ylab = "abs. ground velocity [m/s]\non log10 scale")
gg5       <- plot_FPC(res, FPC_index = 2, log10_scale = TRUE, title = "2.FPC", ylim = ylim_FPCs, xlab = expression(t[0]^{""}~"[registered]"), hide_yLabels = TRUE)


# joint plot
cowplot::plot_grid(gg3_obs, gg3_reg, gg3_repr, gg4, gg5, nrow = 1,
                   rel_widths = c(.30, .19, .19, .28, .22))
# ggsave("../figures/5_2_lasagna_FPCs.png", width = 14.5, height = 3.5)



# main plot 2: amplitude and phase variation by hypo.dis and md -----------
# phase variation: overall absolute domain dilation after first 10 / 20 seconds
gg_distortion5  <- plot_phaseVar(dat, index_ref = 5)
gg_distortion20 <- plot_phaseVar(dat, index_ref = 20, hide_yLabels = TRUE)

gg_distortion5_spatial  <- plot_phaseVar_overSpace(dat, index_ref = 5) +
  theme(legend.key.width = unit(25, units = "points"))
gg_distortion20_spatial <- plot_phaseVar_overSpace(dat, index_ref = 20, hide_yLabels = TRUE) +
  theme(legend.key.width = unit(25, units = "points"))

# amplitude variation: FPC scores
# identical(dat$index, res$Y$tstar) # -> TRUE
dat_scores <- res$Y %>% 
  mutate(hypo.dis   = dat$hypo.dis,
         receiver.x = dat$receiver.x,
         receiver.y = dat$receiver.y,
         md         = dat$md) %>% 
  group_by(id) %>% slice(1) %>% ungroup() %>% 
  mutate(score_FPC1 = res$fpca_obj$scores[,1],
         score_FPC2 = res$fpca_obj$scores[,2])

gg_scores1 <- plot_ampVar(dat_scores, FPC_index = 1)
gg_scores2 <- plot_ampVar(dat_scores, FPC_index = 2, hide_yLabels = TRUE)

gg_scores1_spatial <- plot_ampVar_overSpace(dat_scores, FPC_index = 1) +
  theme(legend.key.width = unit(25, units = "points"))
gg_scores2_spatial <- plot_ampVar_overSpace(dat_scores, FPC_index = 2, hide_yLabels = TRUE) +
  theme(legend.key.width = unit(25, units = "points"))

# joint plots
cowplot::plot_grid(gg_distortion5, gg_distortion20, gg_scores1, gg_scores2, nrow = 1,
                   rel_widths = c(1, .9, 1, .9))
# ggsave("../figures/5_2_scoresAndWarps.pdf", width = 12, height = 3.2, device = cairo_pdf)

cowplot::plot_grid(gg_distortion5_spatial, gg_distortion20_spatial,
                   gg_scores1_spatial, gg_scores2_spatial, nrow = 1,
                   rel_widths = c(1, .8, 1, .8))
# ggsave("../figures/A5_2_scoresAndWarps_spatial.pdf", width = 15, height = 4, device = cairo_pdf)


# check the height profile of the locations of the measurement stations
# ggplot(dat, aes(receiver.x, receiver.y, color = receiver.z)) +
#   geom_point(size = 5)



# plot the topographic map over the evaluated region ----------------------
topo <- readRDS("../data/dat_seismic_topography.rds")

topo_regional   <- topo %>% 
  filter(x >= min(dat$receiver.x), x <= max(dat$receiver.x),
         y >= min(dat$receiver.y), y <= max(dat$receiver.y))

# data for the auxiliary lines for depicting different epicentral distances
dat_circle <- get_finalCircleDat(dat)

# mark the epicenter, located at seismometer "62"
epi_seismometer <- "62"
dat_epi <- dat %>% 
  filter(grepl("_62", id)) %>% slice(1) %>% select(receiver.x, receiver.y) %>% 
  dplyr::rename(x = receiver.x, y = receiver.y) %>% 
  mutate(x = x / 1000, y = y / 1000)

topo_regional %>% 
  mutate(x = x / 1000, y = y / 1000) %>% 
  ggplot(aes(x=x,y=y)) + 
  geom_point(aes(color = z), shape=15, size=1) +
  geom_point(data = dat_epi, color = "red", size = 2) +
  geom_line(data = dat_circle, aes(x, y, group = radius), lty = 2, col = gray(.5)) +
  scale_colour_gradientn(colours = terrain.colors(20), name="Height [m]") +
  xlab("Easting [km]") + ylab("Northing [km]") +
  ggtitle("Topography of the subregion evaluated in the seismic application") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5))
# ggsave("../figures/5_2_subregionTopography.pdf", width = 9, height = 6)


# joint plot of phase and amplitude variation and of the topography
gg_topo <- topo_regional %>% 
  mutate(x = x / 1000, y = y / 1000) %>% 
  ggplot(aes(x=x,y=y)) + 
  geom_point(aes(color = z), shape=15, size=1) +
  geom_point(data = dat_epi, color = "red", size = 4) +
  geom_line(data = dat_circle, aes(x, y, group = radius), lty = 2, col = gray(.5)) +
  scale_colour_gradientn(colours = terrain.colors(20), name="Height [m]") +
  xlab("Easting [km]") + ylab("Northing [km]") +
  ggtitle("Topography") +
  theme_minimal(base_size = 16) +
  theme(plot.title      = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(25, units = "points"))

cowplot::plot_grid(gg_distortion5_spatial, gg_distortion20_spatial,
                   gg_scores1_spatial, gg_scores2_spatial, gg_topo, nrow = 1,
                   rel_widths = c(1, .8, 1, .8, 1))
# ggsave("../figures/A5_2_scoresAndWarps_spatial_topo.pdf", width = 18, height = 4, device = cairo_pdf)



# Appendix plot: warping functions by hypo.dis and md ---------------------
# all warping functions
dat %>% 
  mutate(hypo.dis_cat = cut(hypo.dis, breaks = c(15000, 20000, 25000, 30000, 40000),
                            labels = c("< 20.000km","20.000 - 25.000 km","25.000 - 30.000 km","> 30.000 km")),
         md_cat       = cut(md, breaks = c(0, 0.2, 0.3, 0.4, 0.5),
                            labels = c("[0.1,0.2]","(0.2,0.3]","(0.3,0.4]","(0.4,0.5]"))) %>% 
  ggplot(aes(x = index, y = t_hat)) +
  geom_line(aes(group = id), alpha = 0.05) +
  geom_smooth(se = FALSE, col = "dodgerblue3") +
  facet_grid(md_cat ~ hypo.dis_cat) +
  xlab("t* [observed]") + ylab("t [registered]") +
  ggtitle("Inverse warping functions") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5))
# ggsave("../figures/A5_1_fullWarpings.pdf", width = 10, height = 7.5, device = cairo_pdf)
