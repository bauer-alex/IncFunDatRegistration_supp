
#' Spaghetti plots of observed / registered curves
#' 
#' Create a spaghetti plot of observed or registered curves.
#' 
#' @param Y \code{data.frame} with columns \code{c("id","index","value")}.
#' If the curves should be plotted over another domain as \code{"index"} this
#' can be specified by specifying argument \code{index_var}.
#' @param title Character title of the plot. Pass \code{NULL} to suppress the title.
#' @param index_var Character name of the index variable in \code{Y} which
#' specifies the domain over which the curves should be plotted. Defaults to
#' \code{"index"}.
#' @param log10_scale Indicator if the curves should be plotted on log10 scale.
#' Defaults to \code{FALSE}.
#' @param xlab,ylab Optional labels for the x-axis and y-axis.
#' @param ylim Optional limits for the y-axis.
#' @param y_breaks,y_minorBreaks Optional breaks for the y-axis.
#' @param base_size Base text size as passed to \code{\link[ggplot2]{theme_minimal}}.
#' Defaults to 16.
#' @param alpha Alpha value defining the transparency of the individual lines.
#' Defaults to 0.2.
#' @param hide_yLabels Indicator if the labels on the y-axis should be suppressed.
#' Defaults to \code{FALSE}.
#' 
#' @import ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_spaghetti <- function(Y, title = "<my_title>", index_var = "index", log10_scale = FALSE, 
                           xlab = NULL, ylab = NULL, ylim = NULL, y_breaks = waiver(), y_minorBreaks = waiver(),
                           base_size = 16, alpha = .2, hide_yLabels = FALSE) {
  
  colnames(Y)[match(index_var, colnames(Y))] <- "x"
  gg <- ggplot(Y) +
    geom_line(aes(x = x, y = value, group = id), alpha = alpha, col = "black", lwd = .5) +
    ggtitle(title) + ylab(ylab) + xlab(xlab) +
    scale_y_continuous(limits = ylim, breaks = y_breaks, minor_breaks = y_minorBreaks) + 
    theme_minimal(base_size = base_size) +
    theme(legend.position = "none",
          plot.title      = element_text(hjust = 0.5))
  if (log10_scale) {
    gg <- gg + coord_trans(y = "log10")
  }
  if (hide_yLabels) {
    gg <- gg + 
      theme(axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  return(gg)
}



#' Lasagna plots of observed / registered curves
#' 
#' Create a lasagna plot of observed or registered curves. The plot is
#' sorted by \code{levels(Y$id)} along the y-axis.
#' 
#' @inheritParams plot_spaghetti
#' @param xlim Option limits for the x-axis.
#' @param legend.position Character defining the position of the legend, as
#' passed to \code{\link[ggplot2]{theme}}. Defaults to \code{"bottom"}.
#' @param legend_limits Optional limits for the color scale.
#' @param legend_title Character title for the legend. Pass \code{NULL} to
#' suppress the title.
#' @param line_size Line size passed to \code{\link[ggplot2]{geom_line}}.
#' Defaults to 0.5.
#' 
#' @import ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_lasagna <- function(Y, title, index_var = "index", log10_scale = FALSE, 
                         xlim = NULL, xlab = index_var, ylab = "id",
                         base_size = 16, alpha = .2, hide_yLabels = FALSE,
                         legend.position = "bottom", legend_limits = NULL, legend_title = "value",
                         line_size = .5) {
  
  colnames(Y)[match(index_var, colnames(Y))] <- "x"
  
  if (log10_scale) {
    Y$value <- log10(Y$value)
  }
  
  gg <- ggplot(Y, aes(x = x, y = id, col = value)) +
    geom_line(size = line_size) +
    scale_x_continuous(limits = xlim) +
    scale_color_continuous(name = legend_title, limits = legend_limits) +
    ggtitle(title) + xlab(xlab) + ylab(ylab) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = legend.position,
          plot.title      = element_text(hjust = 0.5),
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank())
    
  return(gg)
}


#' Plot a template function
#' 
#' Plot the template function used for the initial registration step in
#' \code{\link[registr]{register_fpca}} together with a spaghetti plot of the
#' raw data.
#' 
#' @inheritParams plot_spaghetti
#' @param dat,dat_template \code{data.frame}s in the form of arguments \code{Y}
#' and \code{Y_template} as requested by function \code{\link[registr]{register_fpca}}.
#' 
#' @import ggplot2 mgcv
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_template <- function(dat, dat_template, title = "<my_title>", ylab = "Derivative", base_size = 16) {
  
  m        <- gam(value ~ s(index, sp = 0, k = 8), data = dat_template)
  mean_dat <- data.frame(index = seq(min(dat_template$index), max(dat_template$index), length.out = 100)) %>% 
    mutate(mean = predict(m, newdata = .))
  
  plot_spaghetti(dat, title = title, ylab = ylab, alpha = .15, 
                 xlab = "t* [observed]", base_size = base_size) +
    geom_line(data = mean_dat, aes(x = index, y = mean), col = "dodgerblue3", size = 1.2)
}
