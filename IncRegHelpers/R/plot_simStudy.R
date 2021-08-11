
#' Plot the results of the simulation study
#' 
#' Plot function for the results table of the simulation study.
#' 
#' @param dat_pars,dat_res Main results tables from the simulation study.
#' @param corr_structure One of \code{c("ID","AP","AI")}.
#' @param log_scale Indicator if all y-axes should be log transformed. Defaults
#' to \code{FALSE}.
#' @param show_allMeasures Indicator if all measures should be plotted (for the
#' appendix plots) or not. Defaults to \code{FALSE}.
#' @param base_size Base text size as passed to \code{\link[ggplot2]{theme_minimal}}.
#' Defaults to 12.
#' 
#' @import cowplot dplyr ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_results_simStudy <- function(dat_pars, dat_res, corr_structure, log_scale = FALSE,
                                  show_allMeasures = FALSE, base_size = 12) {
  
  # remove potential ids that ran in an error
  ids <- dat_pars$job.id
  if (!all(ids %in% dat_res$job.id)) {
    ids_to_delete <- ids[!(ids %in% dat_res$job.id)]
    dat_pars      <- dat_pars[!(dat_pars$job.id %in% ids_to_delete),]
    message("Removed ",length(ids_to_delete)," of ",length(ids)," jobs since they ran into an error.")
  }
  
  # filter the data
  ids <- dat_pars %>%
    filter(correlation_structure == corr_structure) %>%
    pull(job.id)
  
  res_list <- dat_res %>%
    filter(job.id %in% ids) %>%
    pull(result)
  dat <- dat_pars %>%
    filter(job.id %in% ids) %>%
    select(job.id, algorithm, incompleteness_strength, amplitude_rank)
  
  
  # add the performance measures to dat
  if (show_allMeasures) {
    relevant_vars        <- c("MISE_curves","FPCs_nExtracted","LV_FPCA",#"wMISE_FPCA",
                              "MISE_warpings","MSE_domainLengths")
    relevant_vars_labels <- c("MISE y","n FPCs","1 - LV \u03C8",#"wMISE \u03C8",
                              "MISE h","MSE d")
  } else { 
    relevant_vars        <- c("MISE_curves","LV_FPCA","MISE_warpings","MSE_domainLengths")
    relevant_vars_labels <- c("MISE y","1 - LV \u03C8","MISE h","MSE d")
  }
  
  dat_list <- lapply(relevant_vars, function(var) {
    dat_new       <- dat
    dat_new$value <- sapply(res_list, function(x) x[[var]])
    dat_new$param <- var
    return(dat_new)
  })
  dat <- dplyr::bind_rows(dat_list)
  
  # remove FPC estimation rows of SRVF, since the approach doesn't return orthonormal FPCs in the original space
  dat <- dat %>% filter(!(algorithm == "srvf" & param == "wMISE_FPCA"))
  
  # final preparations
  if ("registrGaussianInc" %in% unique(dat$algorithm)) {
    dat <- dat %>% filter(!(algorithm %in% c("registrGammaInc","registrGamma")))
    dat$algorithm <- factor(dat$algorithm, levels = c("srvf","twoStep","twoStepInc","registrGaussianInc","registrGaussian"),
                            labels = c("SRVF","FGAMM","FGAMM [inc.]","varEM [inc.]","varEM"))
  } else {
    dat$algorithm <- factor(dat$algorithm, levels = c("srvf","twoStep","twoStepInc","registrInc","registr"),
                            labels = c("SRVF","FGAMM","FGAMM [inc.]","varEM [inc.]","varEM"))
  }
  dat$incompleteness_strength <- droplevels(factor(dat$incompleteness_strength, levels = c("CC","WIC","SIC"),
                                                   labels = c("complete curves","weak incompleteness","strong incompleteness")))
  # if (show_allMeasures) # turn the Larsson-Villani measure, s.t. all measures are interpreted as 'lower = better'
  dat$value[dat$param == "LV_FPCA"] <- 1 - dat$value[dat$param == "LV_FPCA"]
  dat$param <- factor(dat$param, levels = relevant_vars, labels = relevant_vars_labels)
  dat$amplitude_rank <- factor(dat$amplitude_rank, levels = c("Rank 1","Rank 2-3","Rank 3-4"),
                               labels = c("1","2-3","3-4"))
  
  colors_algorithm     <- scales::hue_pal()(5) # get default ggplot colors
  colors_mainPlot      <- colors_algorithm[-1]
  colors_completePlot  <- colors_algorithm[c(1,2,5)]
  # label function for 'scale_y_continuous' to control the number of printed decimals
  scaleFUN <- function(x) sprintf("%.4f", x)
  
  dat_mainPlot    <- dat
  
  if ("complete curves" %in% dat_mainPlot$incompleteness_strength) {
    gg_complete <- dat_mainPlot %>%
      filter(incompleteness_strength == "complete curves",
             param                   != "MSE d") %>%
      mutate(incompleteness_strength = droplevels(incompleteness_strength)) %>% 
      ggplot(aes(x = amplitude_rank, y = value, col = algorithm)) +
      geom_boxplot(fill = "gray95") +
      facet_grid(param ~ incompleteness_strength, scales = "free", switch = "y", drop = FALSE) +
      scale_color_manual(values = colors_completePlot) +
      xlab("\namplitude rank") +
      theme_minimal(base_size = base_size) +
      theme(legend.position  = "none",
            legend.title     = element_blank(),
            legend.text      = element_text(size = base_size),
            strip.background = element_rect(color = "gray90", fill = "gray90"),
            axis.title.y     = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "gray96"))
    if (log_scale) {
      gg_complete <- gg_complete + scale_y_continuous(trans = "log10", labels = scaleFUN)
    } else { # no log scale
      gg_complete <- gg_complete + ylim(c(0,NA))
    }
  }
  
  gg_incomplete <- dat_mainPlot %>%
    filter(incompleteness_strength != "complete curves") %>%
    mutate(incompleteness_strength = droplevels(incompleteness_strength)) %>%
    ggplot(aes(x = amplitude_rank, y = value, col = algorithm)) +
    geom_boxplot(fill = "gray95") +
    facet_grid(param ~ incompleteness_strength, scales = "free", switch = "y", drop = FALSE) +
    scale_color_manual(values = colors_mainPlot) +
    xlab("\namplitude rank") +
    theme_minimal(base_size = base_size) +
    theme(legend.position  = "none",
          legend.title     = element_blank(),
          legend.text      = element_text(size = base_size),
          strip.background = element_rect(color = "gray90", fill = "gray90"),
          axis.title.y     = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray96"))
  if (log_scale) {
    gg_incomplete <- gg_incomplete + scale_y_continuous(trans = "log10", labels = scaleFUN)
  } else { # no log scale
    gg_incomplete <- gg_incomplete + ylim(c(0,NA))
  }
  
  if (!("complete curves" %in% dat_mainPlot$incompleteness_strength)) {
    gg_legend <- cowplot::get_legend(ggplot(dat) + geom_boxplot(aes(y = value, col = algorithm)) +
                                       scale_color_manual(values = colors_mainPlot) +
                                       theme(legend.position = "bottom", legend.title = element_blank(),
                                             legend.text = element_text(size = base_size)))
    
    cowplot::plot_grid(gg_incomplete, gg_legend, ncol = 1, rel_heights = c(.9,.1))
  } else {
    gg_legend <- cowplot::get_legend(ggplot(dat) + geom_boxplot(aes(y = value, col = algorithm)) +
                                       theme(legend.position = "bottom", legend.title = element_blank(),
                                             legend.text = element_text(size = base_size)))
    
    cowplot::plot_grid(cowplot::plot_grid(gg_complete, gg_incomplete, nrow = 1, rel_widths = c(.35,.65)),
                       gg_legend, ncol = 1, rel_heights = c(.9,.1))
  }
  
}