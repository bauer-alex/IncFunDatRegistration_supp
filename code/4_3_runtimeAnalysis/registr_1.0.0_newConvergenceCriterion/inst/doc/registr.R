## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----load_libraries, echo = FALSE---------------------------------------------
library(registr)
library(ggplot2)
library(dplyr)

## ----sim_data2----------------------------------------------------------------
registration_data = simulate_unregistered_curves(I = 50, D = 200, seed = 2018)

head(registration_data)


## ----plot_sim2, echo = FALSE, fig.show='hold'---------------------------------

registration_data %>%
	ggplot(aes(index, plogis(latent_mean), group = id)) + theme_bw() + 
	geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")

registration_data %>%
	ggplot(aes(t, plogis(latent_mean), group = id)) + theme_bw() + 
	geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")


## ----sim_data1----------------------------------------------------------------
fpca_data = simulate_functional_data(I = 100, D = 200)

ls(fpca_data)

head(fpca_data$Y)

## ----plot1_sim1, fig.show='hold', fig.width = 2, echo = FALSE-----------------

Y = fpca_data$Y
pc_df = data.frame(pop_mean = fpca_data$alpha, 
									 psi1 = fpca_data$psi1,
									 psi2 = fpca_data$psi2,
									 index = seq(0, 1, length.out = 200),
									 id = 1)

ggplot(Y, aes(index, latent_mean, group = id)) + theme_bw() +
	geom_line(alpha = 0.25) + geom_line(data = pc_df, aes(y = pop_mean), color = "red") 

ggplot(pc_df, aes(index, psi1)) + theme_bw() + geom_line(color = "blue") 
ggplot(pc_df, aes(index, psi2)) + theme_bw() + geom_line(color = "blue") 


## ----plot2_sim1, echo = FALSE-------------------------------------------------
Y %>%
	filter(id == 7) %>%
	ggplot(aes(index, value)) + theme_bw() +
	geom_point(alpha = 0.75, size = 0.25) + geom_line(aes(y = plogis(latent_mean))) +
	labs(y = "Pr(Y = 1)")


## ----register_binary, message = FALSE-----------------------------------------

registr_bin = register_fpca(Y = registration_data, family = "binomial", Kt = 8, Kh = 3, npc = 1)

## ----plot_reg_bin, echo = FALSE, fig.show='hold', fig.width=2-----------------
Y = registr_bin$Y

ggplot(Y, aes(tstar, plogis(latent_mean), group = id)) + theme_bw() + 
	geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")

ggplot(Y, aes(t, plogis(latent_mean), group = id)) + theme_bw() + 
	geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")

ggplot(Y, aes(t_hat, plogis(latent_mean), group = id)) + theme_bw() + 
	geom_line(alpha = 0.25) + labs(y = "Pr(Y = 1)")


## ----plot_reg_bin_warp, echo = FALSE, fig.show='hold'-------------------------

ggplot(Y, aes(tstar, t, group = id)) + theme_bw() + 
	geom_line(alpha = 0.25)

ggplot(Y, aes(tstar, t_hat, group = id)) + theme_bw() + 
	geom_line(alpha = 0.25)



## ----register_gaussian, message = FALSE---------------------------------------
Y$value = Y$latent_mean
registr_gauss = register_fpca(Y = registration_data, family = "gaussian", Kt = 10)

## ----bfpca--------------------------------------------------------------------
bfpca_object = bfpca(fpca_data$Y, npc = 2, Kt = 8, print.iter = TRUE)


## ----plot_bfpca, echo = FALSE, fig.show='hold', fig.width=2-------------------
pc_df = pc_df %>%
	mutate(psi1_est = bfpca_object$efunctions[,1],
	psi2_est = bfpca_object$efunctions[,2],
	alpha_est = bfpca_object$alpha %>% as.vector())

ggplot(pc_df, aes(index, pop_mean)) + theme_bw() + geom_line(color = "blue") +
	geom_line(aes(y = alpha_est), linetype = 2, color = "red")

ggplot(pc_df, aes(index, psi1)) + theme_bw() + geom_line(color = "blue") +
	geom_line(aes(y = psi2_est), linetype = 2, color = "red")

ggplot(pc_df, aes(index, psi2)) + theme_bw() + geom_line(color = "blue") +
	geom_line(aes(y = psi1_est), linetype = 2, color = "red")


## ----registr_function---------------------------------------------------------
data_test_gradient = simulate_unregistered_curves(I = 50, D = 100)

start_time = Sys.time()
reg_analytic = registr(Y = data_test_gradient, family = "binomial", gradient = TRUE)
end_time = Sys.time()

analytic_gradient = as.numeric(round((end_time - start_time), 2))

start_time = Sys.time()
reg_numeric = registr(Y = data_test_gradient, family = "binomial", gradient = FALSE)
end_time = Sys.time()

numeric_gradient = as.numeric(round((end_time - start_time), 2))


