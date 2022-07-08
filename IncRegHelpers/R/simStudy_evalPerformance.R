
#' Evaluation of the curve-specific mean structure
#' 
#' Function to evaluate the curve-specific mean structure by calculating the
#' mean integrated squared error (MISE) between (a) the simulated curves
#' (before adding random noise) and (b) the represented curves by the estimated
#' solution (i.e. estimated global mean structure + curve-specific functional
#' principal component structure).
#' 
#' @param dat_sim Simulated dataset.
#' @param reg_obj Registration object, as returned by \code{\link[registr]{register_fpca}}
#' or \code{\link[fdasrvf]{align_fPCA}}.
#' @param registr_fpcaType One of \code{c("variationalEM","two-step")}. Only
#' needed if \code{reg_obj} was estimated with \code{register_fpca}.
#' @param srvf_regularGrid Vector of the regular time grid. Only needed if
#' \code{reg_obj} was estimated with \code{align_fPCA}.
#' @param return_datasets Indicator if datasets for the true and estimated versions
#' of the (represented) curves without noise should be returned additionally to
#' the MISE value. Defaults to \code{FALSE}.
#' 
#' @import dplyr
#' @export
#' 
#' @return If \code{return_datasets == FALSE} only returns the calculated MISE
#' value. If \code{return_datasets == TRUE} returns a list with the MISE value
#' and two datasets containing the true and estimated versions of the
#' (represented) curves without noise.
#' 
eval_meanStructure <- function(dat_sim, reg_obj, registr_fpcaType, srvf_regularGrid, return_datasets = FALSE) {
  
  # truth: the true raw underlying global mean structure + the true individual FPCs evaluated on the internal time domain
  meanFPC_dat_true <- dat_sim %>%
    select(id, index_raw, index, underlying_meanWithFPCs) %>%
    dplyr::rename(index_internal = index_raw,
                  index_observed = index,
                  value          = underlying_meanWithFPCs)
  # ggplot(meanFPC_dat_true, aes(index_internal, value, group = id)) + geom_line(alpha = 0.2)
  
  # estimation: extract estimated global means + FPCs (curve-specific)
  meanFPC_dat_reg <- add_estMeanStructure(meanFPC_dat_true, reg_obj, srvf_regularGrid = srvf_regularGrid) %>%
    select(id, index_internal_est, index_observed, value_est) %>%
    dplyr::rename(index_internal = index_internal_est, value = value_est)
  if (class(reg_obj) == "registration") {
    if (reg_obj$family == "gamma" & registr_fpcaType == "two-step")
      meanFPC_dat_reg$value <- exp(meanFPC_dat_reg$value)
  }
  
  meanFPC_dat_true$squared_error <- (meanFPC_dat_true$value - meanFPC_dat_reg$value)^2
  MISE_curves <- mean(sapply(unique(meanFPC_dat_true$id), function(x) {
    rows <- which(meanFPC_dat_true$id == x)
    dat_model <- meanFPC_dat_true[rows,]
    f <- function(x) {
      y_approx <- approx(dat_model$index_observed, y = dat_model$squared_error, xout = x, method = "linear", rule = 2)
      return(y_approx$y)
    }
    
    integrate(f, lower = min(meanFPC_dat_true$index_observed[rows]), upper = max(meanFPC_dat_true$index_observed[rows]))$value
  }))
  
  if (return_datasets) {
    return(list(MISE_curves = MISE_curves,
                dat_true    = meanFPC_dat_true,
                dat_est     = meanFPC_dat_reg))
  } else {
    return(MISE_curves)
  }
}



#' Internal helper function for eval_meanStructure
#' 
#' Internal helper function to be called inside \code{\link{eval_meanStructure}},
#' for the comparison of a curve's mean structure
#' (global mean + curve-specific functional principal component structure) with
#' its true (simulated) mean structure.
#' This function takes a dataset of the true (simulated) mean structure with
#' columns \code{"index_internal"} and \code{"value"} (both containing the
#' underlying mean structure) as well as \code{"index_observed"} (containing the
#' warped version of \code{"index_internal"} based on the simulated warping).
#' Then, based on an object estimated either with \code{\link[registr]{register_fpca}}
#' or \code{\link[fdasrvf]{align_fPCA}}, the function adds two new columns to
#' this dataset: \cr
#' 1. \code{"index_internal_est"}: These are the index values on the registered
#' time domain that correspond to 'index_observed' based on the estimated
#' warping functions. \cr
#' 2. \code{"value_est"}: These are the values of the estimated mean structure
#' (global mean + curve-specific FPC structure) evaluated at
#' \code{"index_internal_est"}. \cr
#' With this information \code{"value"} and \code{"value_est"} can be directly
#' compared since both refer to the (true and estimated) mean structures
#' evaluated on the observed domain. \cr
#' Note: If \code{reg_obj} was computed with \code{align_fPCA},
#' the \code{"value_est"} values are computed by first extracting the estimated
#' mean structure from the \code{reg_obj} object and then linearly
#' interpolating it.
#' 
#' @inheritParams eval_meanStructure
#' @param dat_trueMean Dataset containing the simulated mean structure.
#' 
#' @import dplyr fdasrvf
#' @export
#' 
#' @return \code{dat_trueMean} with additional columns, see above description.
#' 
add_estMeanStructure <- function(dat_trueMean, reg_obj, srvf_regularGrid) {
  
  ids <- as.character(unique(dat_trueMean$id))
  
  # basic data prep and ensure appropriate sorting of the data
  if (class(reg_obj) == "registration") { # result of registr::register_fpca()
    dat_estWarping <- reg_obj$Y %>%
      dplyr::rename(index_observed = index, index_internal_est = t_hat) %>%
      arrange(id, index_observed)
    dat_trueMean <- dat_trueMean %>% arrange(id, index_observed) %>%
      mutate(index_internal_est = dat_estWarping$index_internal_est)
    
  } else if (class(reg_obj) == "list") { # result of fdasrvf::align_fPCA()
    dat_estWarping <- data.frame(id                 = rep(levels(dat_trueMean$id), each = nrow(reg_obj$gam)),
                                 index_observed     = as.vector(reg_obj$gam),
                                 index_internal_est = rep(srvf_regularGrid, times = ncol(reg_obj$gam))) %>%
      arrange(id, index_observed)
    # for each id: linearly interpolate the estimated warping function to add it to dat_trueMean
    dat_trueMean <- dat_trueMean %>% arrange(id, index_observed)
    dat_trueMean$index_internal_est <- unlist(lapply(as.character(unique(dat_trueMean$id)), function(id_focus) {
      new_index_observed <- dat_trueMean %>% filter(id == id_focus) %>% pull(index_observed)
      dat_id <- dat_estWarping %>% filter(id == id_focus)
      y_approx  <- approx(dat_id$index_observed, y = dat_id$index_internal_est, xout = new_index_observed, method = "linear", rule = 2)
      return(y_approx$y)
    }))
  }
  
  
  # extract the individual 'global mean + individual FPC structure' curves
  if (class(reg_obj) == "registration") { # result of registr::register_fpca()
    dat_estMean  <- reg_obj$fpca_obj$Yhat %>%
      arrange(id, index) %>%
      dplyr::rename(index_internal_est = index, value_est = value)
    
  } else if (class(reg_obj) == "list") { # result of fdasrvf::align_fPCA()
    dat_estMean_list <- lapply(1:length(ids), function(i) {
      npc <- ncol(reg_obj$vfpca$coef)
      y_srvfSpace <- reg_obj$mqn
      for (j in 1:npc) {
        y_srvfSpace <- y_srvfSpace + reg_obj$vfpca$coef[i,j] * reg_obj$vfpca$U[1:length(reg_obj$mqn),j]
      }
      # Problem:  In srsf_to_f() some f0 value has to be chosen which sets the initial value of the transformed curve.
      # Solution: First use the default 'f0 = 0'. Then use the difference between the transformed curve's mean
      #           and the true curve's mean as the final f0 value.
      #           In this way, we use the best possible representation of a curve for the performance evaluation of the SRVF approach.
      y_origSpace <- fdasrvf::srsf_to_f(q    = matrix(y_srvfSpace, ncol = 1),
                                        time = srvf_regularGrid)
      # linearly interpolate the time grid to evaluate y_origSpace at dat_trueMean$index_internal
      new_index <- dat_trueMean$index_internal[dat_trueMean$id == ids[i]]
      y_approx  <- approx(srvf_regularGrid, y = y_origSpace, xout = new_index, method = "linear", rule = 2)
      y_eval <- y_approx$y
      
      f0_optimal <- mean(dat_trueMean$value[dat_trueMean$id == ids[i]]) - mean(y_eval)
      y_origSpace <- fdasrvf::srsf_to_f(q    = matrix(y_srvfSpace, ncol = 1),
                                        time = srvf_regularGrid,
                                        f0   = f0_optimal)
      
      data.frame(id                 = ids[i],
                 index_internal_est = srvf_regularGrid,
                 value_est          = y_origSpace)
    })
    dat_estMean <- dplyr::bind_rows(dat_estMean_list)
  }
  
  # loop over all ids
  dat_prep_list <- lapply(ids, function(id) {
    dat_id         <- dat_trueMean[dat_trueMean$id == id,]
    dat_estMean_id <- dat_estMean[dat_estMean$id == id,]
    
    if (identical(dat_id$index_internal_est, dat_estMean_id$index_internal_est)) { # no interpolation needed (for register_fpca)
      dat_id$value_est <- dat_estMean_id$value_est
    } else { # interpolation needed (for fdasrvf)
      # for the values of index_internal_est, retrieve the values of the estimated
      # mean structure based on linearly interpolating the mean structure
      y_approx  <- approx(x = dat_estMean_id$index_internal_est, y = dat_estMean_id$value_est,
                          xout = dat_id$index_internal_est, method = "linear", rule = 2)
      dat_id$value_est <- y_approx$y
    }
    
    return(dat_id[,c("value_est"),drop = FALSE])
  })
  dat_prep <- dplyr::bind_rows(dat_prep_list)
  
  dat_trueMean <- dplyr::bind_cols(dat_trueMean, dat_prep)
  
  return(dat_trueMean)
}



#' Evaluation of the estimated global mean
#' 
#' Function to evaluate the performance in estimating the global mean curve by
#' calculating the mean squared error (MSE) between the simulated and the
#' estimated global mean.
#' 
#' @inheritParams eval_meanStructure
#' @param dat_simRegular Only needed if \code{reg_obj} was estimated with
#' \code{align_fPCA}. Interpolated version of \code{dat_sim} s.t. all curves are
#' observed on a regular grid.
#' @param distribution One of \code{c("gaussian","gamma")}, specifying if the
#' model was estimated with Gamma family.
#' @param return_datasets Indicator if datasets for the true and estimated versions
#' of the global mean curve should be returned additionally to
#' the MSE value. Defaults to \code{FALSE}.
#' 
#' @import dplyr fdasrvf
#' @export
#' 
#' @return If \code{return_datasets == FALSE} only returns the calculated MSE
#' value. If \code{return_datasets == TRUE} returns a list with the MSE value
#' and two datasets containing the true and estimated versions of the
#' global mean curve.
#' 
eval_globalMean <- function(dat_sim, dat_simRegular, distribution, reg_obj,
                            registr_fpcaType, return_datasets = FALSE) {
  # truth: the true raw underlying global mean structure evaluated on the internal time domain
  # use the global mean based on some complete curve id (just for covering the whole domain; the global mean is identical for all ids)
  id_complete <- unname(which.max(table(dat_sim$id)))
  mean_dat_true <- dat_sim %>%
    filter(id == id_complete) %>%
    select(index_raw, underlying_mean) %>%
    dplyr::rename(index_internal = index_raw, 
                  value          = underlying_mean)
  if (distribution == "gamma")
    mean_dat_true$value <- exp(mean_dat_true$value)
  # ggplot(mean_dat_true, aes(index_internal, value)) + geom_line()
  
  # estimation: extract estimated global mean
  if (class(reg_obj) == "registration") {
    mean_dat_reg <- data.frame(index_internal = reg_obj$fpca_obj$t_vec,
                               value          = reg_obj$fpca_obj$alpha[,1])
    if (reg_obj$family == "gamma" & registr_fpcaType == "two-step")
      mean_dat_reg$value <- exp(mean_dat_reg$value)
    
  } else { # fdasrvf::align_fPCA
    observed_initialMean <- dat_simRegular %>% filter(index == min(index)) %>% pull(value) %>% mean()
    regular_grid <- unique(dat_simRegular$index)
    mean_reg     <- srsf_to_f(q = matrix(reg_obj$mqn, ncol = 1), time = regular_grid, f0 = observed_initialMean)
    mean_dat_reg <- data.frame(index_internal = regular_grid,
                               value          = mean_reg)
  }
  
  # comparison: calculate the MSE between both curves
  MSE_globalMean <- calcMSE_globalMean(mean_dat_true, mean_dat_reg)
  
  if (return_datasets) {
    return(list(MSE_globalMean = MSE_globalMean,
                dat_trueMean   = mean_dat_true,
                dat_estMean    = mean_dat_reg))
  } else {
    return(MSE_globalMean)
  }
}



#' Internal helper function for eval_globalMean
#' 
#' Internal helper function to be called from within \code{\link{eval_globalMean}},
#' to calculate the MSE between the simulated and estimated global mean curve.
#' The comparison is based on the internal (registered) time domain.
#' 
#' Since the index values are never exactly identical between \code{dat_trueMean}
#' and \code{dat_estMean}, the latter is linearly interpolated.
#' 
#' @param dat_trueMean,dat_estMean Datasets containing the simulated and estimated
#' global mean curve, respectively, with relevant columns \code{"index_internal"}
#' and \code{"value"}.
#' 
#' @import dplyr
#' @export
#' 
#' @return MSE value.
#' 
calcMSE_globalMean <- function(dat_trueMean, dat_estMean) {
  
  # extract the true index and value vectors
  dat_trueMean    <- dat_trueMean %>% arrange(index_internal)
  index_trueMean  <- dat_trueMean$index_internal
  values_trueMean <- dat_trueMean$value
  
  # linearly interpolate dat_estMean
  y_approx  <- approx(dat_estMean$index_internal, y = dat_estMean$value,
                      xout = index_trueMean, method = "linear", rule = 2)
  values_estMean <- y_approx$y
  
  # plot(index_trueMean, values_trueMean, main = "true mean (black), estimated mean (red)")
  # points(index_trueMean, values_estMean, col = "red")
  
  # calculate the MSE
  mse_globalMean <- mean( (values_trueMean - values_estMean)^2 )
  
  return(mse_globalMean)
}



#' Evaluation of the estimated FPC structure based on wMISE
#' 
#' Function to evaluate the performance in estimating the functional principal
#' component (FPC) structure by calculating the weighted mean integrated squared
#' error (wMISE) between the simulated and estimated FPCs. The weighted arithmetic
#' mean is taken of the FPC-specific ISE values by weighting them with the
#' square roots of the simulated eigenvalues. \cr
#' Note that this function is also applicable to models estimated with the SRVF
#' approach. However, this in not used in the paper since the SRVF approach
#' only estimates orthonormal FPCs in the SRVF space, but not in the original
#' function space.
#' 
#' For each FPC, its ISE value is calculated by taking the minimum of
#' \code{c(<ISE to +1*FPC>, <ISE to -1*FPC>)}.
#' 
#' Some notes regarding the output of \code{\link[fdasrvf]{align_fPCA}}: \cr
#' 1. The FPCs are available in \code{x$vfpca$U}, the scores in \code{x$vfpca$coef} \cr
#' 2. Weirdly, nrow of the eigenfunction matrix doesn't equal \code{length(regular_grid)},
#' but \code{length(regular_grid) + 1}... We simply use the rows \code{1:(nrow - 1)}
#' of the matrix to have something matching with \code{length(regular_grid)}.
#' 
#' @inheritParams eval_meanStructure
#' @param FPCdat_sim Dataset containing the simulated FPCs, with columns
#' \code{c("id","index","value")}.
#' @param eigenvalues_sim Numeric vector of simulated eigenvalues. Their square
#' roots are used as weights for the wMISE calculation.
#' @param return_datasets Indicator if datasets for the true and estimated
#' FPCs (the latter scaled with -1, if necessary) should be returned
#' additionally to the wMISE. Defaults to \code{FALSE}.
#' 
#' @import dplyr fdasrvf
#' @export
#' 
#' @return List of wMISE values, based on (i) all FPCs, (ii) only the first FPC,
#' (iii) only the first two FPCs. If \code{return_datasets == TRUE} the returned
#' list also contains two datasets containing the true and estimated versions of
#' the FPCs.
#' 
eval_FPCstructure <- function(FPCdat_sim, reg_obj, eigenvalues_sim, srvf_regularGrid,
                              return_datasets = FALSE) {
  
  # for fdasrvf::align_fPCA: first linearly interpolate the simulated FPCs
  # to the same regular grid that was used for the estimation
  if (class(reg_obj) == "list") {
    dat_fpc_list <- lapply(unique(FPCdat_sim$id), function(i) {
      dat_i     <- FPCdat_sim %>% filter(id == i)
      y_approx  <- approx(dat_i$index, y = dat_i$value, xout = srvf_regularGrid, method = "linear", rule = 2)
      y_regGrid <- y_approx$y
      return(data.frame(index = srvf_regularGrid,
                        value = y_regGrid,
                        id    = i,
                        row.names = NULL))
    })
    FPCdat_sim <- dplyr::bind_rows(dat_fpc_list)
  }
  
  # retrieve the estimated FPCs
  # Note: the number of extracted FPCs can vary from the number of simulated FPCs if npc was determined using the share of explained variance
  if (class(reg_obj) == "registration") {
    npc     <- min(ncol(reg_obj$fpca_obj$efunctions), length(unique(FPCdat_sim$id)))
    FPC_dat <- data.frame(index = rep(reg_obj$fpca_obj$t_vec, times = npc),
                          value = as.vector(reg_obj$fpca_obj$efunctions[,1:npc, drop=FALSE]),
                          id    = as.character(rep(1:npc, each = length(reg_obj$fpca_obj$t_vec))),
                          stringsAsFactors = FALSE)
    
  } else { # fdasrvf
    npc          <- min(ncol(reg_obj$vfpca$coef), length(unique(FPCdat_sim$id)))
    FPC_dat_list <- lapply(1:npc, function(i) {
      fpc_transformed <- srsf_to_f(reg_obj$vfpca$U[1:length(reg_obj$mqn),i], time = srvf_regularGrid)
      FPC_dat <- data.frame(index = srvf_regularGrid,
                            value = as.vector(fpc_transformed),
                            id    = as.character(i),
                            stringsAsFactors = FALSE)
    })
    FPC_dat <- dplyr::bind_rows(FPC_dat_list)
  }
  
  # if reg_obj estimated fewer FPCs than in FPCdat_sim, crop FPCdat_sim
  if (npc < length(unique(FPCdat_sim$id))) {
    FPCdat_sim      <- FPCdat_sim %>% filter(id %in% 1:npc)
    eigenvalues_sim <- eigenvalues_sim[1:npc]
  }
  
  # for registr::register_fpca:
  # ensure that both are comparable on the same index grid.
  # Problem:  The estimated FPCs are not evaluated on the original time grid.
  # Solution: Spline-based interpolation to evaluate FPCs on the original grid.
  if (class(reg_obj) == "registration") {
    index_origGrid <- unique(FPCdat_sim$index)
    FPC_dat_origGrid_list <- lapply(1:npc, function(i) {
      dat_i      <- FPC_dat %>% filter(id == i)
      y_approx  <- approx(dat_i$index, y = dat_i$value, xout = index_origGrid, method = "linear", rule = 2)
      y_origGrid <- y_approx$y
      return(data.frame(index = index_origGrid,
                        value = y_origGrid,
                        id    = as.character(i),
                        row.names = NULL,
                        stringsAsFactors = FALSE))
    })
    FPC_dat_origGrid <- dplyr::bind_rows(FPC_dat_origGrid_list)
    
  } else { # fdasrvf
    FPC_dat_origGrid <- FPC_dat
  }
  # message(identical(FPCdat_sim$index, FPC_dat_origGrid$index)) # -> TRUE
  
  # compute the evaluation metric
  wMISE_FPCA <- calcWeightedMISE_GFPCA(dat_trueFPCs    = FPCdat_sim,
                                       dat_estFPCs     = FPC_dat_origGrid,
                                       eigenvalues_sim = eigenvalues_sim,
                                       return_datasets = return_datasets)
  wMISE_FPCA_onlyFirstFPC     <- calcWeightedMISE_GFPCA(dat_trueFPCs    = FPCdat_sim[FPCdat_sim$id == 1,],
                                                        dat_estFPCs     = FPC_dat_origGrid[FPC_dat_origGrid$id == 1,],
                                                        eigenvalues_sim = eigenvalues_sim[1])
  wMISE_FPCA_onlyFirstTwoFPCs <- calcWeightedMISE_GFPCA(dat_trueFPCs    = FPCdat_sim[FPCdat_sim$id %in% c(1,2),],
                                                        dat_estFPCs     = FPC_dat_origGrid[FPC_dat_origGrid$id %in% c(1,2),],
                                                        eigenvalues_sim = eigenvalues_sim[1:2])
  
  if (return_datasets) {
    return(list(wMISE_FPCA                  = wMISE_FPCA$wMISE_FPCA,
                wMISE_FPCA_onlyFirstFPC     = wMISE_FPCA_onlyFirstFPC,
                wMISE_FPCA_onlyFirstTwoFPCs = wMISE_FPCA_onlyFirstTwoFPCs,
                dat_trueFPCs                = FPCdat_sim %>% mutate(id = as.character(id)),
                dat_estFPCs                 = wMISE_FPCA$dat_estFPCs))
  } else {
    return(list(wMISE_FPCA                  = wMISE_FPCA,
                wMISE_FPCA_onlyFirstFPC     = wMISE_FPCA_onlyFirstFPC,
                wMISE_FPCA_onlyFirstTwoFPCs = wMISE_FPCA_onlyFirstTwoFPCs))
  }
}



#' Internal helper function for eval_FPCstructure
#' 
#' Helper function to be called internally from \code{\link{eval_FPCstructure}},
#' to calculate the goodness-of-fit for a GFPCA solution as the weighted MISE
#' over all functional principal components (FPCs) between the true and the
#' estimated FPC values. The square roots of the true, simulated eigenvalues are
#' used as weights. \cr
#' For each FPC, its ISE value is calculated by taking the minimum of
#' \code{c(<ISE to +1*FPC>, <ISE to -1*FPC>)}.
#' 
#' @inheritParams eval_FPCstructure
#' @param dat_trueFPCs,dat_estFPCs Datasets of true and estimated FPCs, with columns
#' \code{c("id","index","value")}.
#' 
#' @return Numeric wMISE value. If \code{return_datasets == TRUE} returns a list
#' which additionally contains datasets containing the true and estimated FPCs.
#' 
calcWeightedMISE_GFPCA <- function(dat_trueFPCs, dat_estFPCs, eigenvalues_sim, return_datasets = FALSE) {
  
  # ensure that both datasets are sorted identically
  if (!identical(as.character(dat_trueFPCs$id), as.character(dat_estFPCs$id)) |
      !identical(dat_trueFPCs$index, dat_estFPCs$index))
    stop("The specified FPC datasets must be sorted identically.")
  
  # calculate the FPC-specific MISEs by comparison with +1*FPC and -1*FPC
  ids <- unique(as.character(dat_trueFPCs$id))
  ISE_perID <- lapply(ids, function(i) {
    dat_true <- dat_trueFPCs[dat_trueFPCs$id == i,]
    dat_est  <- dat_estFPCs[dat_estFPCs$id == i,]
    
    # compute the integrals to +1*FPC and -1*FPC
    dat_true$squared_error_pos <- (dat_true$value - +1*dat_est$value)^2
    dat_true$squared_error_neg <- (dat_true$value - -1*dat_est$value)^2
    model_interpolation_pos    <- mgcv::gam(squared_error_pos ~ s(index, k = 40, sp = 0), dat = dat_true)
    model_interpolation_neg    <- mgcv::gam(squared_error_neg ~ s(index, k = 40, sp = 0), dat = dat_true)
    f_pos <- function(x) { mgcv::predict.gam(model_interpolation_pos, type = "response", newdata = data.frame(index = x)) }
    f_neg <- function(x) { mgcv::predict.gam(model_interpolation_neg, type = "response", newdata = data.frame(index = x)) }
    integral_pos <- integrate(f_pos, lower = min(dat_true$index), upper = max(dat_true$index))$value
    integral_neg <- integrate(f_neg, lower = min(dat_true$index), upper = max(dat_true$index))$value
    
    # return the minimum ISE of the comparison to +1*FPC and -1*FPC
    return(list(invert_FPC = (integral_pos > integral_neg),
                ISE        = min(integral_pos, integral_neg)))
  })
  
  ISE_vec <- sapply(ISE_perID, function(x) x$ISE)
  
  # potentially scale some estimated FPCs with factor -1
  for (i in 1:length(ids)) {
    if (ISE_perID[[i]]$invert_FPC) {
      rows <- which(dat_estFPCs$id == ids[i])
      dat_estFPCs$value[rows] <- -1*dat_estFPCs$value[rows]
    }
  }
  
  # take the average over all ISEs
  wMISE <- sum(sqrt(eigenvalues_sim) * ISE_vec) / sum(sqrt(eigenvalues_sim))
  if (return_datasets) {
    return(list(wMISE_FPCA   = wMISE,
                dat_trueFPCs = dat_trueFPCs,
                dat_estFPCs  = dat_estFPCs))
  } else {
    return(wMISE)
  }
}



#' Evaluation of the estimated FPC structure based on the LV measure
#' 
#' Function to evaluate the performance in estimating the functional principal
#' component (FPC) structure by calculating the adapted version of the
#' Larsson-Villani measure as introduced in the paper to compare not the
#' individual FPCs, but the span of the overall FPC basis. \cr
#' Note that this function is also applicable to models estimated with the SRVF
#' approach. However, this in not used in the paper since the SRVF approach
#' only estimates orthonormal FPCs in the SRVF space, but not in the original
#' function space.
#' 
#' Some notes regarding the output of \code{\link[fdasrvf]{align_fPCA}}: \cr
#' 1. The FPCs are available in \code{x$vfpca$U}, the scores in \code{x$vfpca$coef} \cr
#' 2. Weirdly, nrow of the eigenfunction matrix doesn't equal \code{length(regular_grid)},
#' but \code{length(regular_grid) + 1}... We simply use the rows \code{1:(nrow - 1)}
#' of the matrix to have something matching with \code{length(regular_grid)}.
#' 
#' @inheritParams eval_FPCstructure
#' 
#' @import fdasrvf
#' @export
#' 
#' @return LV measure.
#' 
eval_FPCstructure_LV <- function(FPCdat_sim, reg_obj, srvf_regularGrid) {
  
  # save the true, simulated FPCs in a matrix with one column per FPC
  matrix_FPCtrue <- matrix(FPCdat_sim$value,
                           byrow = FALSE,
                           ncol  = length(unique(FPCdat_sim$id)))
  index_true     <- unique(FPCdat_sim$index) # respective index vector
  
  # save the estimated FPCs in a matrix
  if (class(reg_obj) == "registration") {
    
    # extract estimated FPCs and save them in a matrix
    matrix_FPCest <- reg_obj$fpca_obj$efunctions
    index_est     <- reg_obj$fpca_obj$t_vec # respective index vector
    
  } else { # fdasrvf
    
    # extract a dataset of all FPCs
    npc          <- ncol(reg_obj$vfpca$coef)
    FPC_dat_list <- lapply(1:npc, function(i) {
      fpc_transformed <- srsf_to_f(reg_obj$vfpca$U[1:length(reg_obj$mqn),i], time = srvf_regularGrid)
      FPC_dat <- data.frame(index = srvf_regularGrid,
                            value = as.vector(fpc_transformed),
                            id    = i)
    })
    FPC_dat_est <- dplyr::bind_rows(FPC_dat_list)
    
    # create a matrix
    matrix_FPCest <- matrix(FPC_dat_est$value,
                            byrow = FALSE,
                            ncol  = length(unique(FPC_dat_est$id)))
    index_est     <- srvf_regularGrid
  }
  
  # interpolate the estimated FPCs s.t. they are evaluated on index_true
  FPCest_ip <- lapply(1:ncol(matrix_FPCest), function(i) {
    y_FPC   <- matrix_FPCest[,i]
    x_index <- index_est
    model   <- gam(y_FPC ~ s(x_index, k = 40, sp = 0))
    y_interpolated <- predict.gam(model, newdata = data.frame(x_index = index_true))
    return(y_interpolated)
  })
  matrix_FPCest_ip <- matrix(unlist(FPCest_ip),
                             byrow = FALSE,
                             ncol  = length(FPCest_ip))
  
  # if necessary, flip the estimated FPCs by scaling them with -1
  for (i in 1:min(ncol(matrix_FPCtrue), ncol(matrix_FPCest_ip))) {
    y_sim <- matrix_FPCtrue[,i]
    y_est <- matrix_FPCest_ip[,i]
    y_est_flipped <- (-1 * (y_est - 0.5)) + 0.5
    if (mean( (y_est_flipped - y_sim)^2 ) < mean( (y_est - y_sim)^2 ))
      matrix_FPCest_ip[,i] <- y_est_flipped
  }
  
  # calculate the LV measure
  lv_overlap(matrix_FPCtrue, matrix_FPCest_ip, dim_refBase = ncol(matrix_FPCtrue))
}



#' Calculate the adapted version of the measure of Larsson and Villani (2001)
#' 
#' Function to calculate the adaption of the measure of Larsson and Villani (2001)
#' as outlined in the paper, to calculate the overlap between a simulated and
#' an estimated functional principal component (FPC) basis. The resulting measure
#' is between 0 and 1.
#' 
#' @param a,b Matrices with the FPCs as columns, evaluated on the same index grid.
#' @param dim_refBase Dimension of the reference basis. E.g., if \code{a} contains
#' the simulated FPC basis and \code{b} the estimated basis, \code{dim_refBase}
#' should be set to \code{ncol(a)}. Defaults to \code{min(NCOL(a), NCOL(b))}.
#' 
#' @import checkmate
#' @export
#' 
#' @return Calculated LV measure.
#' 
lv_overlap <- function(a, b, dim_refBase = min(NCOL(a), NCOL(b))) {
  
  ## overlap measure scaled to [0, 1], not [0, min(NCOL(A),NCOL(B))]
  
  #Rolf Larsson, Mattias Villani (2001)
  #"A distance measure between cointegration spaces"
  checkmate::assert_matrix(a, any.missing = FALSE, mode = "numeric", min.rows = NCOL(a))
  checkmate::assert_matrix(b, any.missing = FALSE, mode = "numeric", min.rows = NCOL(b))
  stopifnot(NROW(a) == NROW(b),  NCOL(a) + NCOL(b) != 0)
  
  if (xor(NCOL(a) == 0, NCOL(b) == 0)) {
    return(0.0)
  }
  
  sv_a <- svd(a, nv = 0)$u
  sv_b <- svd(b, nv = 0)$u
  
  trace <- if (NCOL(b) <= NCOL(a)) {
    sum(diag(t(sv_b) %*% sv_a %*% t(sv_a) %*% sv_b))
  } else {
    sum(diag(t(sv_a) %*% sv_b %*% t(sv_b) %*% sv_a))
  }
  
  trace / dim_refBase
}



#' Evaluation of the estimated domain lengths
#' 
#' Function to evaluate the performance in estimating phase variation by
#' calculating the mean squared error (MSE) between the simulated and estimated
#' domain lengths of the individual curves.
#' 
#' @param dat_trueDomainLengths Dataset of simulated curves, with relevant columns
#' \code{c("id","index","value")}.
#' @param dat_estDomainLengths Dataset of registered curves, with relevant columns
#' \code{c("id","t_hat","value")}. E.g. the object \code{obj$Y} where \code{obj}
#' was returned by \code{\link[registr]{register_fpca}}.
#' 
#' @export
#' 
#' @return Calculated MSE.
#' 
calcMSE_domainLengths <- function(dat_trueDomainLengths, dat_estDomainLengths) {
  
  true_domainLengths <- dat_trueDomainLengths %>%
    group_by(id) %>%
    summarize(domain_length = diff(range(index))) %>%
    mutate(id = as.character(id)) %>%
    arrange(id)
  est_domainLengths <- dat_estDomainLengths %>%
    group_by(id) %>%
    summarize(domain_length = diff(range(t_hat))) %>%
    mutate(id = as.character(id)) %>%
    arrange(id)
  
  # ensure that both datasets are sorted accordingly
  if (!identical(true_domainLengths$id, est_domainLengths$id))
    warning("Problem with the sorting of the data objects. Check the source code of 'calcMSE_domainLengths()'.")
  
  mse <- mean( (true_domainLengths$domain_length - est_domainLengths$domain_length)^2 )
  return(mse)
}



#' Evaluation of the estimated warping functions
#' 
#' Function to evaluate the performance in estimating phase variation by
#' calculating the mean integrated squared error (MISE) between the simulated
#' and estimated warping functions.
#' 
#' If \code{reg_obj} was estimated with \code{\link[fdasrvf]{align_fPCA}} the
#' simulated warping functions are first linearly interpolated to the same
#' regular grid \code{srvf_regularGrid} that was used in the estimation step.
#' 
#' @inheritParams eval_meanStructure
#' 
#' @import dplyr
#' @export
#' 
#' @return Calculated MISE.
#' 
eval_warpingFunctions <- function(dat_sim, reg_obj, srvf_regularGrid) {
  
  # for fdasrvf::align_fPCA: first interpolate the simulated warping functions
  # to the same regular grid that was used for the estimation
  if (class(reg_obj) == "list") {
    dat_warp_list <- lapply(unique(dat_sim$id), function(i) {
      dat_i     <- dat_sim %>% filter(id == i)
      y_approx  <- approx(dat_i$index, y = dat_i$index_raw, xout = srvf_regularGrid, method = "linear", rule = 2)
      y_regGrid <- y_approx$y
      return(data.frame(index     = srvf_regularGrid,
                        index_raw = y_regGrid,
                        id        = i,
                        row.names = NULL))
    })
    dat_sim <- dplyr::bind_rows(dat_warp_list)
  }
  
  # retrieve the true warping functions
  true_time_unwarped <- dat_sim$index_raw
  true_time_warped   <- dat_sim$index
  
  # retrieve the estimated warping functions
  if (class(reg_obj) == "registration") {
    time_unregistered <- reg_obj$Y$tstar
    time_registered   <- reg_obj$Y$t_hat
  } else { # fdasrvf
    time_unregistered <- as.vector(reg_obj$gam)
    time_registered   <- rep(srvf_regularGrid, times = length(unique(dat_sim$id)))
  }
  
  # ensure that both are comparable on the same index grid
  # identical(true_time_warped, time_unregistered) # -> TRUE
  
  # compute the MISE
  squared_errors <- (true_time_unwarped - time_registered)^2
  MISE_warpings <- mean(sapply(unique(dat_sim$id), function(x) {
    rows <- which(dat_sim$id == x)
    dat_model <- data.frame(y_var = squared_errors[rows], x_var = true_time_unwarped[rows])
    
    f <- function(x) {
      y_approx <- approx(dat_model$x_var, y = dat_model$y_var, xout = x, method = "linear", rule = 2)
      return(y_approx$y)
    }
    
    integrate(f, lower = min(dat_model$x_var), upper = max(dat_model$x_var))$value
  }))
}
