
library(dplyr)   # general data handling
library(ggplot2) # data visualization
library(cowplot) # grid of ggplots

library(IncRegHelpers) # helper functions

# Note: Use 'device = cairo_pdf' in ggsave s.t. greek characters written in
#       unicode are displayed correctly.



# gaussian results for fixed npc ------------------------------------------
pars_gaussian <- readRDS("../results/simStudy_npcFixed_gaussian_pars.rds")
res_gaussian  <- readRDS("../results/simStudy_npcFixed_gaussian.rds")

plot_results_simStudy(pars_gaussian, res_gaussian, corr_structure = "ID", log_scale = TRUE, base_size = 16)
# ggsave("../figures/4_2_gaussian_corrID.pdf", width = 11, height = 8, device = cairo_pdf)
plot_results_simStudy(pars_gaussian, res_gaussian, corr_structure = "AP", log_scale = TRUE, base_size = 16)
# ggsave("../figures/A3_gaussian_corrAP.pdf", width = 11, height = 8, device = cairo_pdf)
plot_results_simStudy(pars_gaussian, res_gaussian, corr_structure = "AI", log_scale = TRUE, base_size = 16)
# ggsave("../figures/A3_gaussian_corrAI.pdf", width = 10, height = 8, device = cairo_pdf)



# gaussian results for adaptively estimated npc ---------------------------
pars_gaussian <- readRDS("../results/simStudy_npcSelection_gaussian_pars.rds")
res_gaussian  <- readRDS("../results/simStudy_npcSelection_gaussian.rds")

plot_results_simStudy(pars_gaussian, res_gaussian, corr_structure = "ID",
                      log_scale = TRUE, show_allMeasures = TRUE, base_size = 16)
# ggsave("../figures/A3_npcSelection_gaussian_corrID.pdf", width = 11, height = 9, device = cairo_pdf)
plot_results_simStudy(pars_gaussian, res_gaussian, corr_structure = "AP",
                      log_scale = TRUE, show_allMeasures = TRUE, base_size = 16)
# ggsave("../figures/A3_npcSelection_gaussian_corrAP.pdf", width = 11, height = 9, device = cairo_pdf)
plot_results_simStudy(pars_gaussian, res_gaussian, corr_structure = "AI",
                      log_scale = TRUE, show_allMeasures = TRUE, base_size = 16)
# ggsave("../figures/A3_npcSelection_gaussian_corrAI.pdf", width = 10, height = 9, device = cairo_pdf)



# gamma results for fixed npc ---------------------------------------------
pars_gamma <- readRDS("../results/simStudy_npcFixed_gamma_pars.rds")
res_gamma  <- readRDS("../results/simStudy_npcFixed_gamma.rds")

plot_results_simStudy(pars_gamma, res_gamma, corr_structure = "ID", log_scale = TRUE, base_size = 16)
# ggsave("../figures/4_2_gamma_corrID.pdf", width = 11, height = 8, device = cairo_pdf)
plot_results_simStudy(pars_gamma, res_gamma, corr_structure = "AP", log_scale = TRUE, base_size = 16)
# ggsave("../figures/A3_gamma_corrAP.pdf", width = 11, height = 8, device = cairo_pdf)
plot_results_simStudy(pars_gamma, res_gamma, corr_structure = "AI", log_scale = TRUE, base_size = 16)
# ggsave("../figures/A3_gamma_corrAI.pdf", width = 10, height = 8, device = cairo_pdf)



# gamma results for adaptively estimated npc ------------------------------
pars_gamma <- readRDS("../results/simStudy_npcSelection_gamma_pars.rds")
res_gamma  <- readRDS("../results/simStudy_npcSelection_gamma.rds")

plot_results_simStudy(pars_gamma, res_gamma, corr_structure = "ID",
                      log_scale = TRUE, show_allMeasures = TRUE, base_size = 16)
# ggsave("../figures/A3_npcSelection_gamma_corrID.pdf", width = 11, height = 9, device = cairo_pdf)
plot_results_simStudy(pars_gamma, res_gamma, corr_structure = "AP",
                      log_scale = TRUE, show_allMeasures = TRUE, base_size = 16)
# ggsave("../figures/A3_npcSelection_gamma_corrAP.pdf", width = 11, height = 9, device = cairo_pdf)
plot_results_simStudy(pars_gamma, res_gamma, corr_structure = "AI",
                      log_scale = TRUE, show_allMeasures = TRUE, base_size = 16)
# ggsave("../figures/A3_npcSelection_gamma_corrAI.pdf", width = 10, height = 9, device = cairo_pdf)
