
#' Simulate (in)complete curves with random warping and an FPC structure
#'
#' Simulate potentially incomplete curves with some random warping and a
#' Generalized Functional Principal Component structure.
#' 
#' @param N Number of curves to simulate. Defaults to 20.
#' @param n_timeGrid Number of (regular) evaluation points over the complete
#' time grid. Defaults to 50.
#' @param distribution One of \code{c("gaussian","t","gamma")}, specifying the
#' distribution of the data. Defaults to \code{"gaussian"}.
#' @param random_warping Indicator if the regular time grid should be randomly
#' warped. Defaults to \code{FALSE}.
#' @param incompleteness,incompleteness_rate Indicator if simulated curves
#' should be incomplete. Incompleteness is simulated by sampling individual
#' random cut-off indices in the latter part of the domain. This 'latter part'
#' are the last \code{<incompleteness_rate * 100>\%} of the domain.
#' Defaults to \code{FALSE} and the last 50\code{\%} of the domain.
#' @param FPCA_structure Optional named list with elements \code{"n_FPCs"}
#' (an integer between 1 and 4) and \code{"eigenvalues"} (a non-negative numeric
#' vector of length \code{n_FPCs} that sums to 1. The square roots are used as standard deviation
#' of a normal distribution to sample the individual FPC scores.). If given, structural
#' variation is added around the mean based on given functional principal
#' components.
#' @param corr_amplitude_phase Optional correlation between amplitude variation
#' (specifically the maximum peak) and phase variation (specifically how much
#' the initial part of the time domain is expanded). This parameter doesn't lead
#' to the two measures having exactly the specified correlation, but should be
#' seen as a general parameter, higher (absolute) values of which lead to
#' a clearer association. Can have values between -1 and 1. Defaults to 0.
#' @param corr_amplitude_incompleteness Optional correlation between amplitude
#' variation (specifically the maximum peak) and the curve's incompleteness
#' (specifically the length of the curve). This parameter doesn't lead to the
#' two measures having exactly the specified correlation, but should be seen as
#' a general parameter, higher (absolute) values of which lead to a clearer
#' association. Can have values between -1 and 1. Defaults to 0.
#' @param seed Optional numeric seed for \code{\link{set.seed}}
#' 
#' @import checkmate dplyr registr
#' @export
#' 
#' @return If \code{is.null(FPCA_structure)}, then only the simulated dataset
#' is returned (including the underlying mean structure).
#' If \code{!is.null(FPCA_structure)}, then a list is returned
#' with both the simulated data and a dataset with information on the raw,
#' underlying FPCs.
simulate_curves <- function(N                             = 20,
                            n_timeGrid                    = 50,
                            distribution                  = "gaussian",
                            random_warping                = FALSE,
                            incompleteness                = FALSE,
                            incompleteness_rate           = 0.5,
                            FPCA_structure                = NULL,
                            corr_amplitude_phase          = 0,
                            corr_amplitude_incompleteness = 0,
                            seed                          = NULL) {
  
  checkmate::assert_integerish(N)
  checkmate::assert_integerish(n_timeGrid)
  checkmate::assert_choice(distribution, choices = c("gaussian","t","gamma"))
  checkmate::assert_logical(random_warping)
  checkmate::assert_logical(incompleteness)
  checkmate::assert_numeric(incompleteness_rate, lower = 0, upper = 1)
  checkmate::assert_list(FPCA_structure, len = 2, null.ok = TRUE)
  if (!is.null(FPCA_structure)) {
    checkmate::assert_integerish(FPCA_structure$n_FPCs, lower = 1, upper = 4)
    checkmate::assert_numeric(FPCA_structure$eigenvalues, lower = 0, len = FPCA_structure$n_FPCs)
  }
  checkmate::assert_numeric(corr_amplitude_phase, lower = -1, upper = 1)
  checkmate::assert_numeric(corr_amplitude_incompleteness, lower = -1, upper = 1)
  checkmate::assert_integerish(seed, null.ok = TRUE)
  
  
  if (!is.null(seed))
    set.seed(seed)
  
  ### create a regular time grid as the basis for simulating the curves
  time_grid <- seq(from = 0, to = 1, length.out = n_timeGrid)
  
  ### simulate the curves
  # simulate the mean structure
  list_mean_raw <- lapply(1:N, function(i) {
    x <- dnorm(x = time_grid, mean = 0.45, sd = 0.2)
    x <- x / max(x)
    if (distribution == "gamma")
      x <- 3.5*x - 3 # ensure a mean and variation structure similar to the seismic data close to the hypocenter
    return(x)
  })
  # add FPCA-based variation around the mean
  if (is.null(FPCA_structure)) {
    list_mean <- list_mean_raw
  } else if (!is.null(FPCA_structure)) {
    # create orthogonal eigenfunctions as orthogonal polynomials
    poly_matrix <- poly(time_grid, degree = 5)
    dat_FPC <- data.frame(id    = rep(1:4, each = length(time_grid)) %>% factor(),
                          index = rep(time_grid, times = 4),
                          value = as.vector(poly_matrix[,2:5])) %>%
      filter(id %in% 1:FPCA_structure$n_FPCs)
    # add the eigenfunctions to the data
    list_mean <- lapply(list_mean_raw, function(x) {
      # simulate weights based on the eigenvalues
      weights <- rnorm(FPCA_structure$n_FPCs, mean = 0, sd = sqrt(FPCA_structure$eigenvalues))
      weighted_FPCsum <- dat_FPC %>%
        mutate(weight         = weights[id],
               value_weighted = value * weight) %>%
        group_by(index) %>%
        summarize(weighted_FPCsum = sum(value_weighted)) %>%
        pull(weighted_FPCsum)
      return(x + weighted_FPCsum)
    })
    
  }
  # simulate from the data distribution
  if (distribution == "gaussian") {
    sd <- 0.03 # fix the standard deviation to some value
    list_y <- lapply(1:N, function(i) { rnorm(n = n_timeGrid, mean = list_mean[[i]], sd = sd) })
    
  } else if (distribution == "t") {
    df <- 3 # fix the degrees of freedom to some value
    gaussian_sd <- 0.03 # make the simulated t values comparable to the utilized gaussian distribution
    list_y <- lapply(1:N, function(i) { list_mean[[i]] + gaussian_sd * rt(n = n_timeGrid, df = df) })
    
  } else if (distribution == "gamma") {
    # apply exp() as response function to ensure strict positivity
    list_mean <- lapply(list_mean, function(x) { exp(x) })
    shape  <- 5 # fix the shape to some value
    list_y <- lapply(1:N, function(i) { rgamma(n = n_timeGrid, shape = shape, scale = list_mean[[i]] / shape) })
  }
  
  max_peaks        <- sapply(list_y, max) %>% unlist()
  max_peaks_scaled <- (max_peaks - min(max_peaks)) / diff(range(max_peaks))
  
  ### randomly warp the regular time grid
  if (!random_warping) { # no random warping
    list_timeGrid <- lapply(1:N, function(i) time_grid)
  } else {
    list_timeGrid_unwarped <- lapply(1:N, function(i) time_grid)
    scaling_factors <- max_peaks_scaled
    if (corr_amplitude_phase < 0)
      scaling_factors <- 1 - scaling_factors
    list_timeGrid <- lapply(1:N, function(i) { 
      if (corr_amplitude_phase == 0) {
        basis_coefs <- runif(3, 0, 1) # c(.5,.5,.5) would lead to a diagonal warping function
        
      } else {
        # idea for if the correlation is positive (negative):
        #   higher (smaller) peaks should lead to a stronger (weaker) initial
        #   expansion of the time domain.
        #   For (absolute) correlation 1, fully specify the basis_coefs.
        #   The smaller the (absolute) correlation, the more noise is added to
        #   this initially given solution.
        basis_coefs <- c(0.5 + 0.5 * (scaling_factors[i] - 0.5) + rnorm(n = 1, sd = 0.5 * (1 - abs(corr_amplitude_phase))),
                         1 - (1 / (1 + exp(-6*(scaling_factors[i] - 0.5)))) + rnorm(n = 1, sd = 0.5 * (1 - abs(corr_amplitude_phase))),
                         1 - (1 / (1 + exp(-6*(scaling_factors[i] - 0.5)))) + rnorm(n = 1, sd = 0.5 * (1 - abs(corr_amplitude_phase))))
        basis_coefs <- sapply(basis_coefs, function(x) min(1,max(0,x)))
        if (all(basis_coefs == 0)) # can happen if abs(corr_amplitude_phase) is small
          basis_coefs <- c(0.5 + 0.5 * (scaling_factors[i] - 0.5),
                           rep(1 - (1 / (1 + exp(-6*(scaling_factors[i] - 0.5)))), 2))
      }
      registr:::grid_subj_create(basis_coefs, D = n_timeGrid) %>% as.vector()
    })
  }
  
  ### randomly introduce incompleteness
  if (incompleteness) {
    # simulate random cut-off indices in the latter part of the domain
    min_length_timeGrid <- ceiling((1 - incompleteness_rate) * n_timeGrid)
    cutOff_domain       <- min_length_timeGrid:n_timeGrid
    if (corr_amplitude_incompleteness == 0) {
      indices_cutOff <- sample(x = cutOff_domain, size = N, replace = TRUE)
      
    } else { # correlation between amplitude variation and incompleteness
      # idea for if the correlation is positive (negative):
      #   if the peak of curve i is higher (smaller), make the window of potential
      #   cut-off values smaller by some random number, s.t. higher (smaller)
      #   peaks lead to more complete curves. The random number is drawn with
      #   less variation the higher the absolute specified correlation value.
      # draw the scaling factor as random number
      scaling_factors <- sapply(1:N, function(i) {
        min_value_i <- max_peaks_scaled[i] * abs(corr_amplitude_incompleteness)
        max_value_i <- min_value_i + (1 - abs(corr_amplitude_incompleteness)) * (1 - min_value_i)
        runif(1, min = min_value_i, max = max_value_i)
      })
      scaling_factors <- (scaling_factors - min(scaling_factors)) / diff(range(scaling_factors))
      if (corr_amplitude_incompleteness < 0)
        scaling_factors <- 1 - scaling_factors
      indices_cutOff <- sapply(1:N, function(i) {
        cutOff_domain_i <- cutOff_domain[ceiling(scaling_factors[i] * length(cutOff_domain)):length(cutOff_domain)]
        cutOff_domain_i <- cutOff_domain_i[1:(floor((1 - abs(corr_amplitude_incompleteness)) * length(cutOff_domain_i)))]
        if (length(cutOff_domain_i) == 1)
          cutOff_domain_i <- rep(cutOff_domain_i, 2) # repeat the value to prevent sample() from sampling from 1:<cutOff_domain_i>
        sample(x = cutOff_domain_i, size = 1)
      })
    }
    list_timeGrid  <- lapply(1:N, function(i) { list_timeGrid[[i]][1:indices_cutOff[i]] })
    list_y         <- lapply(1:N, function(i) { list_y[[i]][1:indices_cutOff[i]] })
    list_mean_raw  <- lapply(1:N, function(i) { list_mean_raw[[i]][1:indices_cutOff[i]] })
    list_mean      <- lapply(1:N, function(i) { list_mean[[i]][1:indices_cutOff[i]] })
    if (random_warping)
      list_timeGrid_unwarped <- lapply(1:N, function(i) { list_timeGrid_unwarped[[i]][1:indices_cutOff[i]] })
  }
  
  ### final data preparation
  dat <- data.frame(id    = sapply(1:N, function(i) { rep(i, length(list_timeGrid[[i]])) }) %>% unlist() %>% factor(),
                    index = list_timeGrid %>% unlist(),
                    value = list_y %>% unlist(),
                    underlying_mean         = list_mean_raw %>% unlist(),
                    underlying_meanWithFPCs = list_mean %>% unlist())
  if (random_warping) # add the unwarped time to the data
    dat$index_raw <- list_timeGrid_unwarped %>% unlist()
  
  if (is.null(FPCA_structure)) { # no underlying FPC structure
    return(dat)
  } else { # return data and FPCs as a list
    return(list(data     = dat,
                FPC_data = dat_FPC))
  }
}
