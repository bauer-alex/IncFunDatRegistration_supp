
#' Plot the estimated amplitude variation vs hypocentral distance and the dynamic friction
#' 
#' Hexagonal binning plot function to visualize the estimated scores for a
#' specific functional principal component (FPC) versus the hypocentral distance
#' of the virtual seismogram and the dynamic coefficient of friction of the
#' seismic simulation.
#' 
#' @inheritParams plot_FPC
#' @param dat_scores \code{data.frame} with one row per observed curve and
#' columns containing the hypocentral distance of the virtual seismometer in
#' meters (column \code{"hypo.dis"}), the Easting and Northing coordinate of
#' the seismometer in meters (\code{"receiver.x"} and \code{"receiver.y"}),
#' the dynamic coefficient of friction of the seismic simulation (\code{"md"})
#' and one column \code{"score_FPCx"} for each estimated FPC \code{"x"}.
#' 
#' @import ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_ampVar <- function(dat_scores, FPC_index, hide_yLabels = FALSE, base_size = 16) {
  
  colnames(dat_scores)[colnames(dat_scores) == paste0("score_FPC",FPC_index)] <- "z"
  
  gg <- dat_scores %>% 
    mutate(hypo.dis = hypo.dis / 1000) %>% 
    ggplot(aes(x = hypo.dis, y = md, z = z)) +
    stat_summary_hex(fun = function(a) { mean(a) }, binwidth = c(1.1, 0.055)) +
    xlab("hypocentral distance [km]") + ylab(expression("\u03BC"[d])) +
    ggtitle(paste0("Scores for ",FPC_index,".FPC")) +
    scale_fill_gradient2("mean\nscore", high = "dodgerblue3", low = "firebrick2", mid = "gray95",
                         limits = function(x) { c(-1,1) * max(abs(x)) }) +
    theme_minimal(base_size = base_size) +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  return(gg)
}



#' Internal helper to create a dataset for drawing radial auxiliary lines
#' 
#' Helper function for visualizing amplitude or phase variation over space.
#' Based on the epicenter of the earthquake, this function creates a dataset
#' with n points that all lie on the circumference of a circle with given radius,
#' focusing on a section of the circle between a given interval of angles that
#' specify a window of cardinal directions.
#' 
#' @param center_x,center_y Easting and Northing coordinate of the epicenter.
#' @param radius Radius of the circle.
#' @param angle_min,angle_max Angles that specify the window of cardinal
#' directions. Must be specified as radians, i.e. between 0 and \code{2*pi}.
#' @param n_points Number of evaluation points.
#' 
#' @export
#' 
#' @return \code{data.frame} with \code{n_points} points on the specified part
#' of the circumference of the circle.
#' 
get_circle_points <- function(center_x, center_y, radius, angle_min, angle_max, n_points) {
  
  angles <- seq(angle_min, angle_max, length.out = n_points)
  data.frame(x      = center_x + radius * cos(angles),
             y      = center_y + radius * sin(angles),
             radius = radius)
}


#' Create a dataset of radial auxiliary lines with different radii
#' 
#' This is a wrapper function which calls \code{\link{get_circle_points}}
#' multiple times to create a dataset for plotting radial auxiliary lines
#' with different radii specifically for the seismic dataset.
#' Values 5km, 10km, 15km, 20km, 25km are used as radii.¸
#' 
#' @param dat Seismic dataset with columns \code{c("id","receiver.x","receiver.y")},
#' where the latter define the Easting and Northing coordinate of the virtual
#' seismometers in meters. The dataset must contain information on the seismometer
#' which is located directly at the epicenter (i.e. curve id's ending on \code{"_62"}).
#' 
#' @import dplyr
#' @export
#' 
#' @return \code{data.frame} of points on the different radial auxiliary lines.
#' 
get_finalCircleDat <- function(dat) {
  
  epi_seismometer <- "62"
  epi_coords <- dat %>% 
    filter(grepl("_62", id)) %>% slice(1) %>% select(receiver.x, receiver.y) %>% 
    mutate(receiver.x = receiver.x / 1000, receiver.y = receiver.y / 1000)
  epi_dists <- seq(5, 25, by = 5)
  dat_circleList <- lapply(epi_dists, function(epi_dist) {
    get_circle_points(center_x  = epi_coords$receiver.x,
                      center_y  = epi_coords$receiver.y,
                      radius    = epi_dist,
                      angle_min = .5*pi,
                      angle_max = 1*pi,
                      n_points  = 100)
  })
  dat_circle <- dplyr::bind_rows(dat_circleList)
  return(dat_circle)
}



#' Plot the estimated amplitude variation over space
#' 
#' Hexagonal binning plot function to visualize the estimated scores for a
#' specific functional principal component (FPC) over space.
#' 
#' @inheritParams plot_ampVar
#' @param plot_auxiliaryLines Indicator if radial auxiliary lines should be
#' drawn on the plot for highlighting 5km wide steps from the epicenter.
#' Defaults to \code{TRUE}.
#' 
#' @import dplyr ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_ampVar_overSpace <- function(dat_scores, FPC_index, plot_auxiliaryLines = TRUE,
                                  hide_yLabels = FALSE, base_size = 16) {
  
  colnames(dat_scores)[colnames(dat_scores) == paste0("score_FPC",FPC_index)] <- "z"
  
  gg <- dat_scores %>% 
    mutate(receiver.x = receiver.x / 1000,
           receiver.y = receiver.y / 1000) %>% 
    ggplot(aes(x = receiver.x, y = receiver.y)) +
    stat_summary_hex(aes(z = z), fun = function(a) { mean(a) }, binwidth = c(2, 3)) +
    ggtitle(paste0("Scores for ",FPC_index,".FPC")) +
    xlab("Easting [km]") + ylab("Northing [km]") +
    scale_fill_gradient2("mean\nscore", high = "dodgerblue3", low = "firebrick2", mid = "gray95",
                         limits = function(x) { c(-1,1) * max(abs(x)) }) +
    theme_minimal(base_size = base_size) +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "bottom")

  # create a dataset with the auxiliary lines
  if (plot_auxiliaryLines) {
    dat_circle <- get_finalCircleDat(dat_scores)
    
    gg <- gg +
      geom_line(data = dat_circle, aes(x, y, group = radius), lty = 2, col = gray(.5))
  }
  
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  return(gg)
}



#' Plot the estimated phase variation vs hypocentral distance and the dynamic friction
#' 
#' Hexagonal binning plot function to visualize the estimated phase variation
#' versus the hypocentral distance of the virtual seismogram and the dynamic
#' coefficient of friction of the seismic simulation. Phase variation is
#' measured as the overall absolute domain dilation estimated for the first
#' \code{index_ref} seconds.
#' 
#' @inheritParams plot_ampVar
#' @param dat \code{data.frame} of observed curves with columns
#' \code{c("id","index","value")} and with an additional column \code{"t_hat"}
#' which specifies the estimated registered time domain.
#' @param index_ref index value at which the amount of overall domain dilation
#' is looked at. If it's specified as \code{NA} (default), the overall domain
#' dilations are used of the individual curves.
#' 
#' @import dplyr ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_phaseVar <- function(dat, index_ref = NA, hide_yLabels = FALSE, base_size = 16) {
  
  if (!is.na(index_ref)) {
    dat <- dat %>% filter(index <= index_ref)
  }
  
  gg <- dat %>% 
    group_by(id) %>% filter(index == max(index)) %>% ungroup() %>% 
    mutate(hypo.dis = hypo.dis / 1000,
           t_hat    = t_hat - index) %>% 
    ggplot(aes(x = hypo.dis, y = md, z = t_hat)) +
    stat_summary_hex(fun = function(a) { mean(a) }, binwidth = c(1.1, 0.055)) +
    xlab("hypocentral distance [km]") + ylab(expression("\u03BC"[d])) +
    ggtitle(ifelse(!is.na(index_ref), 
                   paste0("Distortion after ",index_ref,"s"),
                   "Overall domain distortion")) +
    scale_fill_gradient2("mean\ndistortion [s]", high = "#542F85", low = "#88521C", mid = "gray95",
                         limits = function(x) { c(-1,1) * max(abs(x)) }) +
    theme_minimal(base_size = 16) +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  return(gg)
}



#' Plot the estimated phase variation over space
#' 
#' Hexagonal binning plot function to visualize the estimated phase variation
#' over space. Phase variation is measured as the overall absolute domain
#' dilation estimated for the first \code{index_ref} seconds.
#' 
#' @inheritParams plot_phaseVar
#' @param plot_auxiliaryLines Indicator if radial auxiliary lines should be
#' drawn on the plot for highlighting 5km wide steps from the epicenter.
#' Defaults to \code{TRUE}.
#' 
#' @import dplyr ggplot2
#' @export
#' 
#' @return Returns the plot in form of a \code{ggplot2} object.
#' 
plot_phaseVar_overSpace <- function(dat, index_ref = NA, plot_auxiliaryLines = TRUE,
                                    hide_yLabels = FALSE, base_size = 16) {
  if (!is.na(index_ref)) {
    dat <- dat %>% filter(index <= index_ref)
  }
  gg <- dat %>% 
    group_by(id) %>% filter(index == max(index)) %>% ungroup() %>% 
    mutate(receiver.x = receiver.x / 1000,
           receiver.y = receiver.y / 1000,
           t_hat      = t_hat - index) %>% 
    ggplot(aes(x = receiver.x, y = receiver.y)) +
    stat_summary_hex(aes(z = t_hat), fun = function(a) { mean(a) }, binwidth = c(2, 3)) +
    ggtitle(ifelse(!is.na(index_ref), 
                   paste0("Distortion after ",index_ref,"s"),
                   "Overall domain distortion")) +
    xlab("Easting [km]") + ylab("Northing [km]") +
    scale_fill_gradient2("mean\ndistortion [s]", high = "#542F85", low = "#88521C", mid = "gray95",
                         limits = function(x) { c(-1,1) * max(abs(x)) }) +
    theme_minimal(base_size = 16) +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  # create a dataset with the auxiliary lines
  if (plot_auxiliaryLines) {
    dat_circle <- get_finalCircleDat(dat)
    
    gg <- gg +
      geom_line(data = dat_circle, aes(x, y, group = radius), lty = 2, col = gray(.5))
  }
  
  if (hide_yLabels) {
    gg <- gg +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  }
  
  return(gg)
}
