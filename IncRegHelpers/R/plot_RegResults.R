
#' Plot inverse warping functions
#' 
#' Plot the inverse warping functions of a model estimated with
#' \code{registr::register_fpca}.
#' 
#' @inheritParams plot_spaghetti
#' @param Y Element \code{obj$Y} where \code{obj} is an object as returned
#' by \code{\link[registr]{register_fpca}}.
#' @param ylab Character title for the y-axis. Defaults to \code{"t [registered]"}.
#' @param title Optional character title for the plot.
#' 
#' @import ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_warpings <- function(Y, ylab = "t [registered]", title = NULL, alpha = .2,
                          base_size = 16, hide_yLabels = FALSE) {
  
  gg <- ggplot(Y, aes(x = tstar, y = t_hat, group = id)) +
    geom_line(alpha = alpha) +
    ggtitle(title) +
    scale_x_continuous("t* [observed]", limits = range(Y$tstar)) +
    scale_y_continuous(ylab, limits = range(Y$tstar) + c(-.01,.01)) +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  return(gg)
}



#' Plot one estimated functional principal component
#' 
#' Plot a specific functional principal component (FPC) of a model estimated with
#' \code{registr::register_fpca}.
#' 
#' @inheritParams plot_spaghetti
#' @param reg Object estimated with \code{\link[registr]{register_fpca}}.
#' @param FPC_index Numeric index specifying which FPC should be plotted.
#' @param sd_factor Numeric passed to \code{registr:::plot.fpca}.
#' @param hide_xLabels Indicator if the labels on the x-axis should be suppressed.
#' Defaults to \code{FALSE}.
#' 
#' @import ggplot2 registr
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_FPC <- function(reg, FPC_index, title = "<my_title>", sd_factor = 2, log10_scale = FALSE, ylim = NULL,
                     ylab = "", base_size = 16, xlab = "t [registered]", hide_xLabels = FALSE, hide_yLabels = FALSE) {
  
  gg <- plot(reg$fpca_obj, xlab = xlab, ylim = ylim, plot_FPCs = FPC_index, sd_factor = sd_factor,
             add_symbols = FALSE, subtitle = FALSE, ylab = ylab) +
    ggtitle(title)  +
    scale_color_manual(values = c(1, "dodgerblue3", "firebrick2")) +
    scale_linetype_manual(values = c(1,1,1)) +
    theme_minimal(base_size = base_size) +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "none")
  if (log10_scale) {
    gg <- gg + scale_y_log10(limits = ylim)
  }
  if (hide_xLabels) {
    gg <- gg +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())
  }
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  return(gg)
}


#' Visualize the estimated domain dilation
#' 
#' Plot the overall domain dilation of a model estimated with
#' \code{registr::register_fpca}.
#' 
#' @inheritParams plot_warpings
#' 
#' @import dplyr ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_domainDilation <- function(Y, base_size = 16, hide_yLabels = FALSE) {
  
  domainLenghts_dat <- Y %>%
    group_by(id) %>%
    summarize(length_observed   = diff(range(tstar)),
              length_registered = diff(range(t_hat)))
  diagonal_dat <- data.frame(x = range(Y$tstar),
                             y = range(Y$tstar))
  gg <- ggplot(domainLenghts_dat) +
    geom_line(data = diagonal_dat, aes(x = x, y = y), col = "gray70") +
    geom_point(aes(x = length_observed, y = length_registered, group = id),
               alpha = 0.5) +
    scale_x_continuous("domain length [observed]", limits = range(Y$tstar)) +
    scale_y_continuous("domain length [registered]", limits = range(Y$tstar)) +
    theme_minimal(base_size = base_size)
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  return(gg)
}



#' Retrieve datasets of the represented curves and the mean curve
#' 
#' Internal helper function based on the output of \code{\link[registr]{register_fpca}}.
#' The function returns a list with (i) a dataset with the represented version
#' of each observed curve (i.e. mean + individual FPC structure), and (ii) a 
#' dataset of the estimated mean curve.
#' 
#' @inheritParams plot_representedCurves_spaghetti
#' 
#' @import dplyr ggplot2
#' 
#' @return List of two \code{data.frame}s for (i) the represented curves and
#' (ii) the estimated mean curve.
#' 
get_dataList_representedCurves <- function(fpca_obj, exp_transformation = FALSE) {
  
  # dataset for the estimated overall mean curve
  mean_dat <- data.frame(index = fpca_obj$t_vec,
                         value = as.vector(fpca_obj$mu))
  
  # dataset for the represented curves
  dat <- fpca_obj$Yhat
  
  # Instead of the very dense time grid in Yhat (when the FGAMM approach was used),
  # only plot 200 points per curve
  if (fpca_obj$fpca_type == "two-step") {
    dat <- dat %>% 
      group_by(id) %>%
      slice(c(1, n(), sample(2:(n() - 1), size = 198))) %>% 
      ungroup()
  }
  
  if (exp_transformation) {
    dat$value      <- exp(dat$value)
    mean_dat$value <- exp(mean_dat$value)
  }
  
  return(list(dat      = dat,
              mean_dat = mean_dat))
}



#' Spaghetti plot of represented curves
#' 
#' Create a spaghetti plot of the represented curves (overall mean + individual
#' FPC structure) based on a model estimated with \code{\link[registr]{register_fpca}},
#' including the estimated overall mean curve.
#' 
#' @inheritParams plot_spaghetti
#' @param fpca_obj Element \code{fpca_obj} of an object estimated with
#' \code{\link[registr]{register_fpca}}.
#' @param exp_transformation Indicator if a log link was used and data should
#' be exp transformed. Defaults to \code{FALSE}.
#' 
#' @import ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_representedCurves_spaghetti <- function(fpca_obj, exp_transformation = FALSE, log10_scale = FALSE,
                                             title = "<my_title>", xlab = NULL, ylab = NULL, ylim = NULL,
                                             y_breaks = waiver(), y_minorBreaks = waiver(),
                                             base_size = 16, alpha = .2, hide_yLabels = FALSE) {
  
  dat_list <- get_dataList_representedCurves(fpca_obj, exp_transformation)
  
  gg <- plot_spaghetti(dat_list$dat, log10_scale = log10_scale, y_breaks = y_breaks, y_minorBreaks = y_minorBreaks,
                       title = title, xlab = xlab, ylab = ylab, ylim = ylim, alpha = alpha) +
    geom_line(data = dat_list$mean_dat, mapping = aes(index, value), col = "dodgerblue3", size = 1.2)
  
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  return(gg)
}



#' Lasagna plot of represented curves
#' 
#' Create a lasagna plot of the represented curves (overall mean + individual
#' FPC structure) based on a model estimated with \code{\link[registr]{register_fpca}},
#' including the estimated overall mean curve. By default, the plot is
#' sorted by \code{levels(fpca_obj$Yhat$id)} along the y-axis.
#' 
#' @inheritParams plot_lasagna
#' @param fpca_obj Element \code{fpca_obj} of an object estimated with
#' \code{\link[registr]{register_fpca}}.
#' @param ids_sorted Optional vector of the ids in the data as they should be
#' sorted along the y-axis of the plot.
#' @param exp_transformation Indicator if a log link was used and data should
#' be exp transformed. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link{plot_lasagna}}.
#' 
#' @import ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_representedCurves_lasagna <- function(fpca_obj, ids_sorted = NULL, exp_transformation = FALSE, log10_scale = FALSE,
                                           title = "<my_title>", base_size = 16, hide_yLabels = FALSE, ...) {
  
  dat_list <- get_dataList_representedCurves(fpca_obj, exp_transformation)
  
  if (!is.null(ids_sorted)) {
    dat_list$dat$id <- factor(dat_list$dat$id, levels = ids_sorted)
  }
  
  gg <- plot_lasagna(dat_list$dat, log10_scale = log10_scale, title = title, ...)
  
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  return(gg)
}
