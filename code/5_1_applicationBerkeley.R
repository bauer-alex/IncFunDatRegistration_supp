
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

library(registr) # our accompanying package

library(IncRegHelpers) # helper functions



# data prep ---------------------------------------------------------------
dat <- registr::growth_incomplete

ggplot(dat, aes(index, value, group = id)) + geom_line()



# model estimation --------------------------------------------------------
# 1) varEM with assumed completeness
reg_varEM_CC <- register_fpca(Y             = dat,
                              family        = "gaussian",
                              fpca_type     = "variationalEM",
                              npc_criterion = .9)

# 2) FGAMM with assumed completeness
reg_FGAMM_CC <- register_fpca(Y             = dat,
                              family        = "gaussian",
                              fpca_type     = "two-step",
                              npc_criterion = c(.9,.02))

# 3) varEM with assumed incompleteness
reg_varEM_IC <- register_fpca(Y              = dat,
                              family         = "gaussian",
                              fpca_type      = "variationalEM",
                              incompleteness = "full",
                              lambda_inc     = .025,
                              npc_criterion  = .9)

# 4) FGAMM with assumed incompleteness
reg_FGAMM_IC <- register_fpca(Y              = dat,
                              family         = "gaussian",
                              fpca_type      = "two-step",
                              incompleteness = "full",
                              lambda_inc     = .025,
                              npc_criterion  = c(.9,.02))

# gather all results in one list
model_list  <- list(reg_varEM_CC, reg_FGAMM_CC, reg_varEM_IC, reg_FGAMM_IC)
model_names <- c("varEM","FGAMM","varEM [inc.]","FGAMM [inc.]")



# main plot: registered curves + first two FPCs ---------------------------
# Row 1: Observed + registered curves
row1_1 <- plot_spaghetti(dat, title = "Observed curves",         xlab = "t* [observed]", ylab = "Derivative")
row1_2 <- plot_spaghetti(reg_FGAMM_CC$Y, title = "FGAMM",        index_var = "t_hat", xlab = "t [registered]", ylab = "Derivative")
row1_3 <- plot_spaghetti(reg_FGAMM_IC$Y, title = "FGAMM [inc.]", index_var = "t_hat", xlab = "t [registered]", hide_yLabels = TRUE)
row1_4 <- plot_spaghetti(reg_varEM_IC$Y, title = "varEM [inc.]", index_var = "t_hat", xlab = "t [registered]", hide_yLabels = TRUE)
row1_5 <- plot_spaghetti(reg_varEM_CC$Y, title = "varEM",        index_var = "t_hat", xlab = "t [registered]", hide_yLabels = TRUE)

row1_2to5 <- cowplot::plot_grid(row1_2, row1_3, row1_4, row1_5,
                                nrow = 1, rel_widths = c(1,.85,.85,.85))

# Row 2+3: FPCs 1 & 2
ylim      <- c(-2,29)
row2_2 <- plot_FPC(reg_FGAMM_CC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, ylab = "Derivative", hide_xLabels = TRUE)
row2_3 <- plot_FPC(reg_FGAMM_IC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row2_4 <- plot_FPC(reg_varEM_IC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row2_5 <- plot_FPC(reg_varEM_CC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row3_2 <- plot_FPC(reg_FGAMM_CC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, ylab = "Derivative")
row3_3 <- plot_FPC(reg_FGAMM_IC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, hide_yLabels = TRUE)
row3_4 <- plot_FPC(reg_varEM_IC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, hide_yLabels = TRUE)
row3_5 <- plot_FPC(reg_varEM_CC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, hide_yLabels = TRUE)

row2_1    <- ggplot() + theme_minimal()
row2_2to5 <- cowplot::plot_grid(row2_2, row2_3, row2_4, row2_5,
                                row3_2, row3_3, row3_4, row3_5,
                                nrow = 2, rel_heights = c(.45,.55),
                                rel_widths = c(1,.85,.85,.85))

# joint plot
cowplot::plot_grid(row1_1, row1_2to5,
                   row2_1, row2_2to5,
                   nrow = 2, rel_heights = c(.3,.7), rel_widths = c(.2,.8))
# ggsave("../figures/5_1_results_FPCs.pdf", width = 12, height = 7)



# Appendix plot: same as above, including warpings and four FPCs ----------
# Row 1: Observed + registered curves
# use the objects from above

# Row 2: Estimated warping functions
row2_1 <- ggplot() + theme_minimal()
row2_2 <- plot_warpings(reg_FGAMM_CC$Y)
row2_3 <- plot_warpings(reg_FGAMM_IC$Y, hide_yLabels = TRUE)
row2_4 <- plot_warpings(reg_varEM_IC$Y, hide_yLabels = TRUE)
row2_5 <- plot_warpings(reg_varEM_CC$Y, hide_yLabels = TRUE)

row2_2to5 <- cowplot::plot_grid(row2_2, row2_3, row2_4, row2_5,
                                nrow = 1, rel_widths = c(1,.85,.85,.85))

# Row 3-6: FPCs 1-4
ylim      <- c(-3,29)
row3_2 <- plot_FPC(reg_FGAMM_CC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, ylab = "Derivative", hide_xLabels = TRUE)
row3_3 <- plot_FPC(reg_FGAMM_IC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row3_4 <- plot_FPC(reg_varEM_IC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row3_5 <- plot_FPC(reg_varEM_CC, FPC_index = 1, sd_factor = 2, title = "FPC 1", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row4_2 <- plot_FPC(reg_FGAMM_CC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, ylab = "Derivative", hide_xLabels = TRUE)
row4_3 <- plot_FPC(reg_FGAMM_IC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row4_4 <- plot_FPC(reg_varEM_IC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row4_5 <- plot_FPC(reg_varEM_CC, FPC_index = 2, sd_factor = 2, title = "FPC 2", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row5_2 <- plot_FPC(reg_FGAMM_CC, FPC_index = 3, sd_factor = 2, title = "FPC 3", ylim = ylim, ylab = "Derivative", hide_xLabels = TRUE)
row5_3 <- plot_FPC(reg_FGAMM_IC, FPC_index = 3, sd_factor = 2, title = "FPC 3", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row5_4 <- plot_FPC(reg_varEM_IC, FPC_index = 3, sd_factor = 2, title = "FPC 3", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row5_5 <- plot_FPC(reg_varEM_CC, FPC_index = 3, sd_factor = 2, title = "FPC 3", ylim = ylim, hide_xLabels = TRUE, hide_yLabels = TRUE)
row6_2 <- plot_FPC(reg_FGAMM_CC, FPC_index = 4, sd_factor = 2, title = "FPC 4", ylim = ylim, ylab = "Derivative")
row6_3 <- plot_FPC(reg_FGAMM_IC, FPC_index = 4, sd_factor = 2, title = "FPC 4", ylim = ylim, hide_yLabels = TRUE)
row6_4 <- plot_FPC(reg_varEM_IC, FPC_index = 4, sd_factor = 2, title = "FPC 4", ylim = ylim, hide_yLabels = TRUE)
row6_5 <- plot_FPC(reg_varEM_CC, FPC_index = 4, sd_factor = 2, title = "FPC 4", ylim = ylim, hide_yLabels = TRUE)

row3_1    <- ggplot() + theme_minimal()
row3_2to5 <- cowplot::plot_grid(row3_2, row3_3, row3_4, row3_5,
                                row4_2, row4_3, row4_4, row4_5,
                                row5_2, row5_3, row5_4, row5_5,
                                row6_2, row6_3, row6_4, row6_5,
                                nrow = 4, rel_heights = c(.45,.45,.45,.55),
                                rel_widths = c(1,.85,.85,.85))

# joint plot
cowplot::plot_grid(row1_1, row1_2to5,
                   row2_1, row2_2to5,
                   row3_1, row3_2to5,
                   nrow = 3, rel_heights = c(1.1,1,5), rel_widths = c(.2,.8))
# ggsave("../figures/A4_2_results_FPCs_full.pdf", width = 12, height = 14)



# model estimation based on varying lambda values -------------------------
# 1) FGAMM with assumed incompleteness, lambda = 0
reg1 <- register_fpca(Y              = dat,
                      family         = "gaussian",
                      fpca_type      = "two-step",
                      incompleteness = "full",
                      lambda_inc     = 0,
                      npc_criterion  = c(.9,.02))

# 2) FGAMM with assumed incompleteness, lambda = 0.025
reg2 <- register_fpca(Y              = dat,
                      family         = "gaussian",
                      fpca_type      = "two-step",
                      incompleteness = "full",
                      lambda_inc     = .025,
                      npc_criterion  = c(.9,.02))

# 3) FGAMM with assumed incompleteness, lambda = 1
reg3 <- register_fpca(Y              = dat,
                      family         = "gaussian",
                      fpca_type      = "two-step",
                      incompleteness = "full",
                      lambda_inc     = 1,
                      npc_criterion  = c(.9,.02))



# Appendix plot: Results based on varying lambda values -------------------
# Row 1: registered curves
row1_1 <- plot_spaghetti(reg1$Y, title = expression(paste(lambda , " = 0")),     index_var = "t_hat", xlab = "t [registered]", ylab = "Derivative")
row1_2 <- plot_spaghetti(reg2$Y, title = expression(paste(lambda , " = 0.025")), index_var = "t_hat", xlab = "t [registered]", hide_yLabels = TRUE)
row1_3 <- plot_spaghetti(reg3$Y, title = expression(paste(lambda , " = 1")),     index_var = "t_hat", xlab = "t [registered]", hide_yLabels = TRUE)

# Row 2: estimated inverse warping functions
row2_1 <- plot_warpings(reg1$Y)
row2_2 <- plot_warpings(reg2$Y, hide_yLabels = TRUE)
row2_3 <- plot_warpings(reg3$Y, hide_yLabels = TRUE)

# Row 3: estimated overall domain dilations
row3_1 <- plot_domainDilation(reg1$Y)
row3_2 <- plot_domainDilation(reg2$Y, hide_yLabels = TRUE)
row3_3 <- plot_domainDilation(reg3$Y, hide_yLabels = TRUE)

# joint plot
cowplot::plot_grid(row1_1, row1_2, row1_3,
                   row2_1, row2_2, row2_3,
                   row3_1, row3_2, row3_3,
                   nrow = 3, byrow = TRUE, rel_widths = c(1,.85,.85))
# ggsave("../figures/A4_3_varyLambda.pdf", height = 8, width = 11)
