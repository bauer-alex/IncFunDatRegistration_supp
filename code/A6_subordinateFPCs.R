
library(dplyr)   # general data handling
library(ggplot2) # data visualization

library(registr)       # GFPCA - FGAMM approach
library(IncRegHelpers) # helper functions




# data simulation ---------------------------------------------------------
# Take a simple simulation setting with
# - data with Gaussian noise
# - 100 complete curves with 100 measurements each
# - no random warping
# - three simulated eigenfunctions which explain 70%, 25%, 5% of the variation

sim_list <- simulate_curves(N              = 100,
                            n_timeGrid     = 100,
                            random_warping = FALSE,
                            FPCA_structure = list(n_FPCs = 3,
                                                  eigenvalues = c(.7,.25,.05)),
                            seed = 2021)
dat <- sim_list$data



# data visualization ------------------------------------------------------
# spaghetti plot
ggplot(dat, aes(x = index, y = value, group = id)) +
  geom_line() + ggtitle("Simulated curves")

# lasagna plot
ggplot(dat, aes(x = index, y = id, col = value)) +
  geom_line(size = 2) +
  theme_minimal(base_size = 16) +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())
# ggsave("../figures/A6_subordinateFPCs_simCurves.pdf", width = 6, height = 5)

# simulated FPCs, similar to the simulation study
ggplot(sim_list$FPC_data, aes(x = index, y = value, col = id)) +
  geom_line() + ggtitle("Simulated FPCs")



# FPCA estimation ---------------------------------------------------------
# extract the first 20 functional principal components
x <- gfpca_twoStep(Y = dat, npc = 20, Kc = 10)



# visualize the shares of explained variance ------------------------------
dat_evalues <- data.frame(index = 1:length(x$evalues),
                          pve   = x$evalues / x$evalues_sum)

# only visualize the explained variance shares starting from the fourth
# eigenfunction to focus on the slow decline of the eigenvalues
dat_evalues %>% 
  filter(index >= 4) %>% 
  ggplot(aes(x = index, y = pve)) +
  geom_line() +
  geom_point() +
  xlim(c(1,20)) +
  xlab("FPC index") + ylab("share of explained variance")

# extract the (cumulative) explained variance shares
round(x$evalues / x$evalues_sum, 3)
round(cumsum(x$evalues / x$evalues_sum), 3)



# visualize the functional principal components ---------------------------
dat_FPC <- data.frame(FPC   = rep(1:x$npc, each = length(x$t_vec)),
                      pve   = rep(x$evalues / x$evalues_sum, each = length(x$t_vec)),
                      index = rep(x$t_vec, times = x$npc),
                      value = as.vector(x$efunctions)) %>% 
  mutate(FPC_title = paste0("FPC ", FPC, "\nPVE: ", round(100*pve, 1), "%")) %>% 
  mutate(FPC_title = factor(FPC_title, levels = unique(FPC_title)))

dat_FPC %>% 
  ggplot(aes(x = index, y = value)) +
  geom_line(col = gray(0.3)) +
  facet_wrap(~ FPC_title) +
  scale_x_continuous("t", breaks = seq(0, 1, by = .25), labels = c("0","0.25","0.5","0.75","1")) +
  theme_minimal(base_size = 16)
# ggsave("../figures/A6_subordinateFPCs.pdf", width = 10, height = 7)
