
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

# set a global theme for ggplots
theme_set(theme_minimal(base_size = 16) +
            theme(plot.title = element_text(hjust = 0.5)))

library(IncRegHelpers) # helper functions




# 1) Distribution of the data ---------------------------------------------
### Gaussian distribution
dat_gaussian <- simulate_curves(N              = 30,
                                random_warping = TRUE,
                                distribution   = "gaussian",
                                seed           = 2020)

gg_gaussian1 <- ggplot(dat_gaussian, aes(x = index_raw, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nGaussian structure")
gg_gaussian2 <- ggplot(dat_gaussian, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nGaussian structure")
gg_gaussian3 <- ggplot(dat_gaussian, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nGaussian structure")
plot_grid(gg_gaussian1, gg_gaussian2, gg_gaussian3, nrow = 1)
# ggsave("../figures/A3_simSettings_1_GaussianDist.pdf", width = 15, height = 4)


### Gamma distribution
dat_gamma <- simulate_curves(N              = 30,
                             random_warping = TRUE,
                             distribution   = "gamma",
                             seed           = 2020)

gg_gamma1 <- ggplot(dat_gamma, aes(x = index_raw, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nGamma structure")
gg_gamma2 <- ggplot(dat_gamma, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nGamma structure")
gg_gamma3 <- ggplot(dat_gamma, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nGammma structure")
plot_grid(gg_gamma1, gg_gamma2, gg_gamma3, nrow = 1)
# ggsave("../figures/A3_simSettings_1_GammaDist.pdf", width = 15, height = 4)



# 2) Rank of amplitude variation ------------------------------------------
### Fixed FPCs as basis
time_grid <- seq(0, 1, length.out = 50)
dat_mean  <- data.frame(index = time_grid,
                        value = dnorm(x = time_grid, mean = 0.45, sd = 0.2))
# create orthogonal eigenfunctions as orthogonal polynomials
poly_matrix <- poly(time_grid, degree = 5)
dat_FPC <- data.frame(id    = rep(1:4, each = length(time_grid)) %>% factor(),
                      index = rep(time_grid, times = 4),
                      value = as.vector(poly_matrix[,2:5]))
groups <- paste0(rep(1:4, times = 2), "_", rep(c("-2ef","+2ef"), each = 4))
dat_FPCvar <- data.frame(id    = rep(dat_FPC$id, times = 2),
                         group = rep(c("mean - 2*FPC","mean + 2*FPC"), each = 4*length(time_grid)) %>% factor(),
                         index = rep(time_grid, times = 2*4),
                         value = c(rep(dat_mean$value, times = 4) - 2*dat_FPC$value,
                                   rep(dat_mean$value, times = 4) + 2*dat_FPC$value)) %>%
  mutate(curve_id = paste(as.character(id), as.character(group)) %>% factor())
ggplot(mapping = aes(x = index, y = value)) +
  geom_line(data = dat_mean) +
  geom_line(data = dat_FPCvar, aes(group = curve_id, col = group)) +
  xlab("t [registered]") + ggtitle("Fixed FPCs used for all simulations") +
  scale_color_manual("", values = c("firebrick2","dodgerblue2")) +
  facet_wrap(~ id) +
  theme(legend.position = "bottom")
# ggsave("../figures/A3_simSettings_2_FPCs.pdf", width = 10, height = 5)


### Rank 1
dat_rank1 <- simulate_curves(N              = 30,
                             random_warping = TRUE,
                             distribution   = "gaussian",
                             FPCA_structure = list(n_FPCs = 1, eigenvalues = 1),
                             seed           = 2020)
dat_rank1 <- dat_rank1$data

gg_rank1_1 <- ggplot(dat_rank1, aes(x = index_raw, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nAmplitude rank 1")
gg_rank1_2 <- ggplot(dat_rank1, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nAmplitude rank 1")
gg_rank1_3 <- ggplot(dat_rank1, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nAmplitude rank 1")
plot_grid(gg_rank1_1, gg_rank1_2, gg_rank1_3, nrow = 1)
# ggsave("../figures/A3_simSettings_2_rank1.pdf", width = 15, height = 4)


### Rank 2-3 (eigenvalues fastly descending)
dat_rank3 <- simulate_curves(N              = 30,
                             random_warping = TRUE,
                             distribution   = "gaussian",
                             FPCA_structure = list(n_FPCs = 3, eigenvalues = c(0.7,0.25,0.05)),
                             seed           = 2020)
dat_rank3 <- dat_rank3$data

gg_rank3_1 <- ggplot(dat_rank3, aes(x = index_raw, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nAmplitude rank 3\n(Eigenvalues: 0.7, 0.25, 0.05)")
gg_rank3_2 <- ggplot(dat_rank3, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nAmplitude rank 3\n")
gg_rank3_3 <- ggplot(dat_rank3, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nAmplitude rank 3\n")
plot_grid(gg_rank3_1, gg_rank3_2, gg_rank3_3, nrow = 1)
# ggsave("../figures/A3_simSettings_2_rank3.pdf", width = 15, height = 4)


### Rank 3-4 (eigenvalues slowly decreasing)
dat_rank4 <- simulate_curves(N              = 30,
                             random_warping = TRUE,
                             distribution   = "gaussian",
                             FPCA_structure = list(n_FPCs = 4, eigenvalues = c(0.4,0.3,0.2,0.1)),
                             seed           = 2020)
dat_rank4 <- dat_rank4$data

gg_rank4_1 <- ggplot(dat_rank4, aes(x = index_raw, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nAmplitude rank 4\n(Eigenvalues: 0.4, 0.3, 0.2, 0.1)")
gg_rank4_2 <- ggplot(dat_rank4, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nAmplitude rank 4\n")
gg_rank4_3 <- ggplot(dat_rank4, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nAmplitude rank 4\n")
plot_grid(gg_rank4_1, gg_rank4_2, gg_rank4_3, nrow = 1)
# ggsave("../figures/A3_simSettings_2_rank4.pdf", width = 15, height = 4)



# 3) Strength of incompleteness -------------------------------------------
### Complete curves (CC)
dat_cc <- simulate_curves(N              = 30,
                          random_warping = TRUE,
                          distribution   = "gaussian",
                          seed           = 2020)

gg_cc1 <- dat_cc %>%
  group_by(id) %>%
  mutate(max_index = max(index_raw)) %>%
  ungroup() %>%
  arrange(max_index, id) %>%
  mutate(id = factor(id, levels = unique(as.character(id)))) %>%
  ggplot(aes(x = index_raw, y = id, col = value)) +
  geom_line(size = 3) +
  scale_color_continuous("Value", high = "midnightblue", low = "lightskyblue1") +
  theme(axis.text.y = element_blank()) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nComplete curves")
gg_cc2 <- ggplot(dat_cc, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nComplete curves")
gg_cc3 <- ggplot(dat_cc, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nComplete curves")
plot_grid(gg_cc1, gg_cc2, gg_cc3, nrow = 1)
# ggsave("../figures/A3_simSettings_3_CC.pdf", width = 15, height = 4)


### Weak incompleteness (WIC)
dat_wic <- simulate_curves(N                   = 30,
                           random_warping      = TRUE,
                           incompleteness      = TRUE,
                           incompleteness_rate = 0.3,
                           distribution        = "gaussian",
                           seed                = 2020)

gg_wic1 <- dat_wic %>%
  group_by(id) %>%
  mutate(max_index = max(index_raw)) %>%
  ungroup() %>%
  arrange(max_index, id) %>%
  mutate(id = factor(id, levels = unique(as.character(id)))) %>%
  ggplot(aes(x = index_raw, y = id, col = value)) +
  geom_line(size = 3) +
  scale_color_continuous("Value", high = "midnightblue", low = "lightskyblue1") +
  theme(axis.text.y = element_blank()) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nWeak incompleteness")
gg_wic2 <- ggplot(dat_wic, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nWeak incompleteness")
gg_wic3 <- ggplot(dat_wic, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nWeak incompleteness")
plot_grid(gg_wic1, gg_wic2, gg_wic3, nrow = 1)
# ggsave("../figures/A3_simSettings_3_WIC.pdf", width = 15, height = 4)


### Strong incompleteness (SIC)
dat_sic <- simulate_curves(N                   = 30,
                           random_warping      = TRUE,
                           incompleteness      = TRUE,
                           incompleteness_rate = 0.6,
                           distribution        = "gaussian",
                           seed                = 2020)

gg_sic1 <- dat_sic %>%
  group_by(id) %>%
  mutate(max_index = max(index_raw)) %>%
  ungroup() %>%
  arrange(max_index, id) %>%
  mutate(id = factor(id, levels = unique(as.character(id)))) %>%
  ggplot(aes(x = index_raw, y = id, col = value)) +
  geom_line(size = 3) +
  scale_color_continuous("Value", high = "midnightblue", low = "lightskyblue1") +
  theme(axis.text.y = element_blank()) +
  xlab("t [registered]") + ggtitle("Raw curves before random warping\nStrong incompleteness")
gg_sic2 <- ggplot(dat_sic, aes(x = index, y = value, group = id)) +
  geom_line(col = "dodgerblue4", alpha = 0.2) + xlim(c(0,1)) +
  xlab("t* [observed]") + ggtitle("Randomly warped curves\nStrong incompleteness")
gg_sic3 <- ggplot(dat_sic, aes(x = index_raw, y = index)) +
  geom_line(aes(group = id), col = "dodgerblue4", alpha = 0.2) + geom_smooth(se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + xlab("t [registered]") + ylab("t* [observed]") +
  ggtitle("Warping functions (and their mean)\nStrong incompleteness")
plot_grid(gg_sic1, gg_sic2, gg_sic3, nrow = 1)
# ggsave("../figures/A3_simSettings_3_SIC.pdf", width = 15, height = 4)



# 4) Correlation structure ------------------------------------------------
### Complete independence of amplitude, phase and incompleteness (ID)
dat_id <- simulate_curves(N                             = 200,
                          random_warping                = TRUE,
                          incompleteness                = TRUE,
                          incompleteness_rate           = 0.6,
                          distribution                  = "gaussian",
                          FPCA_structure                = list(n_FPCs = 1, eigenvalues = c(1)),
                          seed                          = 2021)
dat_id <- dat_id$data

peaks <- sapply(unique(dat_id$id), function(y) { dat_id %>% filter(id == y) %>% pull(value) %>% max() })
dat_id <- dat_id %>%
  mutate(peak     = peaks[as.numeric(id)],
         peak_cat = cut(peak, breaks = c(0, quantile(peak, probs = c(1/3, 2/3)), 2.5)))

gg_id1 <- ggplot(dat_id, aes(x = index_raw, y = value, group = id, col = peak_cat)) +
  geom_line(alpha = 0.2) +
  xlim(c(0,1)) + xlab("t [registered]") +
  ggtitle("No correlation between amplitude, phase and incompleteness\nRaw curves before random warping\n(colored by amplitude to evaluate the correlation)") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_id2 <- ggplot(dat_id, aes(x = index, y = value, group = id, col = peak_cat)) +
  geom_line(alpha = 0.2) +
  xlim(c(0,1)) + xlab("t* [observed]") + ggtitle("Randomly warped curves") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_id3 <- ggplot(dat_id, aes(x = index_raw, y = index, group = id, col = peak_cat)) +
  geom_line(alpha = 0.1) + geom_smooth(aes(group = peak_cat), se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + ggtitle("Warping functions (and their means)") +
  xlab("t [registered]") + ylab("t* [observed]") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_id4 <- dat_id %>%
  group_by(id) %>%
  mutate(index_length = length(index)) %>%
  ggplot(aes(x = peak_cat, y = index_length, col = peak_cat)) + geom_violin() +
  ylim(c(0,50)) + xlab("category based on peak size") + ylab("# observed measurements") +
  ggtitle("Level of (in)completeness by amplitude size") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        plot.title      = element_text(hjust = 0.5))
plot_grid(gg_id1, gg_id2, gg_id3, gg_id4, ncol = 1)
# ggsave("../figures/A3_simSettings_4_ID.pdf", width = 8, height = 10)


### Correlation of amplitude and phase (AP)
dat_ap <- simulate_curves(N                    = 200,
                          random_warping       = TRUE,
                          incompleteness       = TRUE,
                          incompleteness_rate  = 0.6,
                          distribution         = "gaussian",
                          FPCA_structure       = list(n_FPCs = 1, eigenvalues = c(1)),
                          corr_amplitude_phase = -1,
                          seed                 = 2021)
dat_ap <- dat_ap$data

peaks <- sapply(unique(dat_ap$id), function(y) { dat_ap %>% filter(id == y) %>% pull(value) %>% max() })
dat_ap <- dat_ap %>%
  mutate(peak     = peaks[as.numeric(id)],
         peak_cat = cut(peak, breaks = c(0, quantile(peak, probs = c(1/3, 2/3)), 2.5)))

gg_ap1 <- ggplot(dat_ap, aes(x = index_raw, y = value, group = id, col = peak_cat)) +
  geom_line(alpha = 0.2) +
  xlim(c(0,1)) + xlab("t [registered]") +
  ggtitle("Correlation between amplitude and phase\nRaw curves before random warping\n(colored by amplitude to evaluate the correlation)") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_ap2 <- ggplot(dat_ap, aes(x = index, y = value, group = id, col = peak_cat)) +
  geom_line(alpha = 0.2) +
  xlim(c(0,1)) + xlab("t* [observed]") + ggtitle("Randomly warped curves") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_ap3 <- ggplot(dat_ap, aes(x = index_raw, y = index, group = id, col = peak_cat)) +
  geom_line(alpha = 0.1) + geom_smooth(aes(group = peak_cat), se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + ggtitle("Warping functions (and their means)") +
  xlab("t [registered]") + ylab("t* [observed]") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
# gg_ap4 <- dat_ap %>%
#   group_by(id) %>%
# mutate(index_length = length(index)) %>%
#   ggplot(aes(x = peak_cat, y = index_length, col = peak_cat)) + geom_violin() +
#   ylim(c(0,50)) + ggtitle("Level of (in)completeness by amplitude size") +
#   theme_minimal(base_size = 12) +
#   theme(legend.position = "none",
#         axis.text.x     = element_blank(),
#         plot.title      = element_text(hjust = 0.5))
plot_grid(gg_ap1, gg_ap2, gg_ap3, ncol = 1)
# ggsave("../figures/A3_simSettings_4_AP.pdf", width = 8, height = 7.5)


### Correlation of amplitude and incompleteness (AI)
dat_ai <- simulate_curves(N                             = 200,
                          random_warping                = TRUE,
                          incompleteness                = TRUE,
                          incompleteness_rate           = 0.6,
                          distribution                  = "gaussian",
                          FPCA_structure                = list(n_FPCs = 1, eigenvalues = c(1)),
                          corr_amplitude_incompleteness = 0.8,
                          seed                          = 2021)
dat_ai <- dat_ai$data

peaks <- sapply(unique(dat_ai$id), function(y) { dat_ai %>% filter(id == y) %>% pull(value) %>% max() })
dat_ai <- dat_ai %>%
  mutate(peak     = peaks[as.numeric(id)],
         peak_cat = cut(peak, breaks = c(0, quantile(peak, probs = c(1/3, 2/3)), 2.5)))

gg_ai1 <- ggplot(dat_ai, aes(x = index_raw, y = value, group = id, col = peak_cat)) +
  geom_line(alpha = 0.2) +
  xlim(c(0,1)) + xlab("t [registered]") +
  ggtitle("Correlation between amplitude and incompleteness\nRaw curves before random warping\n(colored by amplitude to evaluate the correlation)") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_ai2 <- ggplot(dat_ai, aes(x = index, y = value, group = id, col = peak_cat)) +
  geom_line(alpha = 0.2) +
  xlim(c(0,1)) + xlab("t* [observed]") + ggtitle("Randomly warped curves") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_ai3 <- ggplot(dat_ai, aes(x = index_raw, y = index, group = id, col = peak_cat)) +
  geom_line(alpha = 0.1) + geom_smooth(aes(group = peak_cat), se = FALSE) +
  xlim(c(0,1)) + ylim(c(0,1)) + ggtitle("Warping functions (and their means)") +
  xlab("t [registered]") + ylab("t* [observed]") +
  facet_wrap(~ peak_cat) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        strip.text      = element_blank(),
        plot.title      = element_text(hjust = 0.5))
gg_ai4 <- dat_ai %>%
  group_by(id) %>%
  mutate(index_length = length(index)) %>%
  ggplot(aes(x = peak_cat, y = index_length, col = peak_cat)) + geom_violin() +
  ylim(c(0,50)) + xlab("category based on peak size") + ylab("# observed measurements") +
  ggtitle("Level of (in)completeness by amplitude size") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x     = element_blank(),
        plot.title      = element_text(hjust = 0.5))
plot_grid(gg_ai1, gg_ai2, gg_ai3, gg_ai4, ncol = 1)
# ggsave("../figures/A3_simSettings_4_AI.pdf", width = 8, height = 10)
