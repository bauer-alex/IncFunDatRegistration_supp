
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

library(registr) # contains the Berkeley data with simulated incompleteness
dat <- growth_incomplete

library(IncRegHelpers) # helper functions

# general settings
npc_criterion  <- c(.9, .02) # npc criterion
lambda_inc     <- .025       # penalization parameter
max_iterations <- 100        # max. number of joint iterations
n_cores        <- 1          # number of cores to be used for the computations



# define different template functions -------------------------------------
# template 1: mean over all curves
Y_template1 <- dat

# template 2: very early and small main peak
ids2        <- c("girl12","girl13","girl14")
Y_template2 <- dat %>% filter(id %in% ids2)

# template 3: main peak a bit later on and a bit more salient
ids3        <- c("boy02","boy03","boy05","boy08","boy19")
Y_template3 <- dat %>% filter(id %in% ids3)

# template 4: main peak even later
ids4        <- c("boy05","boy30")
Y_template4 <- dat %>% filter(id %in% ids4)



# main estimation step ----------------------------------------------------
# template 1
reg1 <- register_fpca(dat,
                      fpca_type      = "two-step", 
                      npc_criterion  = npc_criterion,
                      incompleteness = "full",
                      lambda_inc     = lambda_inc,
                      max_iterations = max_iterations,
                      cores          = n_cores)

# template 2
reg2 <- register_fpca(dat,
                      fpca_type      = "two-step",
                      npc_criterion  = npc_criterion,
                      incompleteness = "full",
                      lambda_inc     = lambda_inc,
                      Y_template     = Y_template2,
                      max_iterations = max_iterations,
                      cores          = n_cores)

# template 3
reg3 <- register_fpca(dat,
                      fpca_type      = "two-step", 
                      npc_criterion  = npc_criterion,
                      incompleteness = "full",
                      lambda_inc     = lambda_inc,
                      Y_template     = Y_template3,
                      max_iterations = max_iterations,
                      cores          = n_cores)

# template 4
reg4 <- register_fpca(dat,
                      fpca_type      = "two-step",
                      npc_criterion  = npc_criterion,
                      incompleteness = "full",
                      lambda_inc     = lambda_inc,
                      Y_template     = Y_template4,
                      max_iterations = max_iterations,
                      cores          = n_cores)



# joint plot of all relevant information ----------------------------------
base_size <- 12

# 1.row: template function
gg_row1_1 <- plot_template(dat, Y_template1, title = "Observed curves\nwith template 1", base_size = base_size)
gg_row1_2 <- plot_template(dat, Y_template2, title = "Observed curves\nwith template 2", ylab = NULL, base_size = base_size)
gg_row1_3 <- plot_template(dat, Y_template3, title = "Observed curves\nwith template 3", ylab = NULL, base_size = base_size)
gg_row1_4 <- plot_template(dat, Y_template4, title = "Observed curves\nwith template 4", ylab = NULL, base_size = base_size)
gg_row1 <- cowplot::plot_grid(gg_row1_1, gg_row1_2, gg_row1_3, gg_row1_4,
                              nrow = 1, rel_widths = c(1.08, 1, 1, 1))

# 2.row: spaghetti plot of registered curves
gg_row2_1 <- plot_spaghetti(reg1$Y, index_var = "t_hat", title = "Registered curves", ylab = "Derivative", xlab = "t [registered]", alpha = .15, base_size = base_size)
gg_row2_2 <- plot_spaghetti(reg2$Y, index_var = "t_hat", title = "Registered curves", xlab = "t [registered]", alpha = .15, base_size = base_size)
gg_row2_3 <- plot_spaghetti(reg3$Y, index_var = "t_hat", title = "Registered curves", xlab = "t [registered]", alpha = .15, base_size = base_size)
gg_row2_4 <- plot_spaghetti(reg4$Y, index_var = "t_hat", title = "Registered curves", xlab = "t [registered]", alpha = .15, base_size = base_size)
gg_row2 <- cowplot::plot_grid(gg_row2_1, gg_row2_2, gg_row2_3, gg_row2_4,
                              nrow = 1, rel_widths = c(1.08, 1, 1, 1))

# 3.row: inverse warping functions
gg_row3_1 <- plot_warpings(reg1$Y, title = "Inverse warping functions", alpha = .15, base_size = base_size)
gg_row3_2 <- plot_warpings(reg2$Y, title = "Inverse warping functions", alpha = .15, ylab = NULL, base_size = base_size)
gg_row3_3 <- plot_warpings(reg3$Y, title = "Inverse warping functions", alpha = .15, ylab = NULL, base_size = base_size)
gg_row3_4 <- plot_warpings(reg4$Y, title = "Inverse warping functions", alpha = .15, ylab = NULL, base_size = base_size)
gg_row3 <- cowplot::plot_grid(gg_row3_1, gg_row3_2, gg_row3_3, gg_row3_4,
                              nrow = 1, rel_widths = c(1.08, 1, 1, 1))

# 4.row: number of joint iterations
# dat_iter <- data.frame(template   = c("Template 1", "Template 2", "Template 3", "Template 4"),
#                        iterations = c(reg1$convergence$iterations, reg2$convergence$iterations,
#                                       reg3$convergence$iterations, reg4$convergence$iterations))
# gg_row4 <- ggplot(dat_iter, aes(x = template, weight = iterations)) + geom_bar() + ylab("# Iterations") +
#   ggtitle("Joint iterations") +
#   theme_minimal(base_size = base_size) +
#   theme(plot.title   = element_text(hjust = 0.5),
#         axis.title.x = element_blank())

# 5.row: 1.FPC
ylim_FPC1 <- c(-5, 37)
gg_row5_1 <- plot_FPC(reg1, FPC_index = 1, title = "FPC 1", ylab = "Derivative", ylim = ylim_FPC1, base_size = base_size)
gg_row5_2 <- plot_FPC(reg2, FPC_index = 1, title = "FPC 1", ylim = ylim_FPC1, base_size = base_size)
gg_row5_3 <- plot_FPC(reg3, FPC_index = 1, title = "FPC 1", ylim = ylim_FPC1, base_size = base_size)
gg_row5_4 <- plot_FPC(reg4, FPC_index = 1, title = "FPC 1", ylim = ylim_FPC1, base_size = base_size)
gg_row5 <- cowplot::plot_grid(gg_row5_1, gg_row5_2, gg_row5_3, gg_row5_4,
                              nrow = 1, rel_widths = c(1.08, 1, 1, 1))

# 6.row: 2.FPC
ylim_FPC2 <- c(-5.5, 29)
gg_row6_1 <- plot_FPC(reg1, FPC_index = 2, title = "FPC 2", ylab = "Derivative", ylim = ylim_FPC2, base_size = base_size)
gg_row6_2 <- plot_FPC(reg2, FPC_index = 2, title = "FPC 2", ylim = ylim_FPC2, base_size = base_size)
gg_row6_3 <- plot_FPC(reg3, FPC_index = 2, title = "FPC 2", ylim = ylim_FPC2, base_size = base_size)
gg_row6_4 <- plot_FPC(reg4, FPC_index = 2, title = "FPC 2", ylim = ylim_FPC2, base_size = base_size)
gg_row6 <- cowplot::plot_grid(gg_row6_1, gg_row6_2, gg_row6_3, gg_row6_4,
                              nrow = 1, rel_widths = c(1.08, 1, 1, 1))

# joint plot
cowplot::plot_grid(gg_row1, gg_row2, gg_row3, #gg_row4,
                   gg_row5, gg_row6,
                   ncol = 1, rel_heights = c(1.2, 1, 1, 1, 1))
# ggsave("../figures/A7_choiceOfTemplateFunction.pdf", width = 11, height = 12)
