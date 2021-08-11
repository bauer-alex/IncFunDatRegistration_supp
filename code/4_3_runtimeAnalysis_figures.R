
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots



# data preparation --------------------------------------------------------
dat <- readRDS("../results/results_runtimeAnalysis.rds")

# check if an error occurred in some setting
table(dat$error_message)
# -> all good


cols_runtimes <- grep("runtimes.run", colnames(dat))

dat_long_list <- lapply(1:nrow(dat), function(i) {
  data.frame(method         = dat$method[i],
             distribution   = dat$distribution[i],
             incompleteness = dat$incompleteness[i],
             n_curves       = dat$n_curves[i],
             n_timeGrid     = dat$n_timeGrid[i],
             n_cores        = dat$n_cores[i],
             run            = 1:length(cols_runtimes),
             runtime        = unlist(dat[i, cols_runtimes], use.names = FALSE) / 1000000000) # microbenchmark$time is given in nanoseconds
})
# Note: In the end, all methods were run with a complete curve assumption.
#       A complete curve setting was needed sincce fdasrvf doesn't work for incomplete curves.
dat_long <- dplyr::bind_rows(dat_long_list) %>%
  mutate(setting = paste0("N = ",n_curves,"\nD_i = ",n_timeGrid)) %>%
  mutate(method_label = case_when(method == "fdasrvf"                              ~ "SRVF",
                                  method == "registr old variationalEM complete"   ~ "varEM 1.0",
                                  method == "registr new two-step incomplete"      ~ "FGAMM",
                                  method == "registr new variationalEM complete"   ~ "varEM 2.1"),
                                  # method == "registr new variationalEM incomplete" ~ "varEM 2.1 [inc.]"),
         method_label = factor(method_label, levels = c("SRVF","varEM 1.0","varEM 2.1","FGAMM")))



# version 1: line plots ---------------------------------------------------
# only plot specific settings
settings_focus <- c("N = 1000\nD_i = 50","N = 1000\nD_i = 100","N = 1000\nD_i = 150","N = 2000\nD_i = 50","N = 3000\nD_i = 50")
dat_plot <- dat_long %>%
  filter(setting %in% settings_focus) %>%
  mutate(setting = factor(setting, levels = settings_focus))

# focus on the median runtimes
dat_plot <- dat_plot %>%
  mutate(method_label = as.character(method_label),
         method_label = case_when(distribution == "gamma" & method_label == "FGAMM" ~ "FGAMM (Gamma)",
                                  TRUE ~ method_label),
         method_label = factor(method_label, levels = c("SRVF","varEM 1.0","varEM 2.1","FGAMM","FGAMM (Gamma)"))) %>%
  group_by(method_label, setting) %>%
  summarize(median_runtime = median(runtime))

dat_plot_varyingN <- dat_plot %>%
  filter(grepl("D_i = 50", setting)) %>%
  mutate(issue   = "D_i = 50, varying N",
         setting = gsub("\nD_i = 50", "", setting))
dat_plot_varyingDi <- dat_plot %>%
  filter(grepl("N = 1000", setting)) %>%
  mutate(issue   = "N = 1000, varying D_i",
         setting = gsub("N = 1000\n", "", setting))
dat_plot <- dplyr::bind_rows(dat_plot_varyingN, dat_plot_varyingDi) %>%
  mutate(setting = factor(setting, levels = c("N = 1000","N = 2000","N = 3000","D_i = 50","D_i = 100","D_i = 150")))

colors_algorithm <- scales::hue_pal()(5) # get default ggplot colors, in accordance to the sim study plots
# add one more color for the 'Wrobel 1.0' algorithm
colors_algorithm <- c(colors_algorithm, gray(0.2))

ggplot(dat_plot, aes(x = setting, y = median_runtime / 60, group = method_label, col = method_label, linetype = method_label)) +
  geom_line() + geom_point() +
  facet_grid(~ issue, scales = "free_x") +
  scale_color_manual(values = colors_algorithm[c(1,6,5,2,2)]) +
  scale_linetype_manual(values = c(1,1,1,1,2)) +
  ylab("Runtime [min]") +
  theme_minimal(base_size = 16) +
  theme(legend.title    = element_blank(),
        axis.title.x    = element_blank())
# ggsave("../figures/4_3_runtimes.pdf", width = 10, height = 5)



# version 2: boxplots -----------------------------------------------------
# only plot specific settings
settings_focus <- c("N = 1000\nD_i = 50","N = 1000\nD_i = 100","N = 1000\nD_i = 150","N = 2000\nD_i = 50","N = 3000\nD_i = 50")
dat_plot <- dat_long %>%
  filter(setting %in% settings_focus) %>%
  mutate(setting = factor(setting, levels = settings_focus))

gg_gaussian <- dat_plot %>%
  filter(distribution == "gaussian") %>%
  ggplot(aes(x = method_label, y = runtime, col = method_label)) +
  geom_hline(yintercept = c(0,800), lty = 2, alpha = 0.2) +
  geom_boxplot() +
  facet_grid(. ~ setting, drop = FALSE) +
  ggtitle("Gaussian") + ylab("Runtime [s]") +
  theme_minimal(base_size = 14) +
  theme(legend.position  = "none",
        plot.title       = element_text(hjust = 0.5),
        # panel.grid.minor = element_blank(),
        axis.title.x     = element_blank(),
        axis.text.x      = element_blank())

gg_dummy <- ggplot(dat_plot, aes(y = runtime, col = method_label)) + geom_boxplot() +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        legend.title    = element_blank())
gg_legend <- cowplot::get_legend(gg_dummy)

col_twoStep <- tail(scales::hue_pal()(length(unique(dat_plot$method))), 1) # retrieve ggplot color
gg_gamma <- dat_plot %>%
  filter(distribution == "gamma") %>%
  ggplot(aes(x = method_label, y = runtime, col = method_label)) +
  geom_hline(yintercept = c(0,800), lty = 2, alpha = 0.2) +
  geom_boxplot() +
  facet_grid(. ~ setting, drop = FALSE) +
  scale_color_manual(values = col_twoStep) +
  ggtitle("Gamma") + ylab("Runtime [s]") + ylim(c(0,NA)) +
  theme_minimal(base_size = 14) +
  theme(legend.position  = "none",
        plot.title       = element_text(hjust = 0.5),
        # panel.grid.minor = element_blank(),
        axis.title.x     = element_blank(),
        axis.text.x      = element_blank())

cowplot::plot_grid(cowplot::plot_grid(gg_gaussian, gg_gamma, nrow = 1, rel_widths = c(.55,.45)),
                   gg_legend, nrow = 2, rel_heights = c(.9,.1))
# ggsave("../figures/4_3_runtimes_boxplots.pdf", width = 10, height = 3)
