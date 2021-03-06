% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_seismicApplication.R
\name{get_circle_points}
\alias{get_circle_points}
\title{Internal helper to create a dataset for drawing radial auxiliary lines}
\usage{
get_circle_points(center_x, center_y, radius, angle_min, angle_max, n_points)
}
\arguments{
\item{center_x, center_y}{Easting and Northing coordinate of the epicenter.}

\item{radius}{Radius of the circle.}

\item{angle_min, angle_max}{Angles that specify the window of cardinal
directions. Must be specified as radians, i.e. between 0 and \code{2*pi}.}

\item{n_points}{Number of evaluation points.}
}
\value{
\code{data.frame} with \code{n_points} points on the specified part
of the circumference of the circle.
}
\description{
Helper function for visualizing amplitude or phase variation over space.
Based on the epicenter of the earthquake, this function creates a dataset
with n points that all lie on the circumference of a circle with given radius,
focusing on a section of the circle between a given interval of angles that
specify a window of cardinal directions.
}
