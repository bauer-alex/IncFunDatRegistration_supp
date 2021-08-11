
### TODO If the simulation study should be run without fixing the number of FPCs,
#        but estimating it adaptively, the three rows in this code file where
#        'npc_criterion' is defined have to be uncommented.


#' Wrapper to call FGAMM with potentially assumed incompleteness
#' 
#' Wrapper function to call FGAMM with Gaussian or Gamma family in the simulation
#' study. Written specifically for the use with the \code{batchtools} package.
#' 
#' @param data,job General arguments as requested by \code{batchtools}.
#' @param instance A specific simulated data setting as returned by
#' \code{\link{simData_wrapper}}.
#' @param lambda_inc Penalization parameter passed to \code{\link[registr]{register_fpca}}.
#' Defaults to 0.
#' @param max_iterations Maximum number of iterations of the joint algorithm
#' of \code{register_fpca}.
#' @param n_cores Number of cores to use for computation.
#' @param incompleteness Indicator if incompleteness should be assumed or not.
#' Defaults to \code{TRUE}.
#' 
#' @import registr
#' @export
#' 
#' @return List with performance criteria.
#' 
FGAMMWithIncompleteness_wrapper <- function(data, job, instance, lambda_inc = 0,
                                            max_iterations, n_cores, incompleteness = TRUE) {
  
  # estimation step ---------------------------------------------------------
  if (incompleteness) {
    incompleteness <- "trailing"
  } else {
    incompleteness <- NULL
  }
  
  # new registr code
  runtime <- system.time({
    reg <- register_fpca(Y                            = instance$data,
                         Kt                           = 8,
                         Kh                           = 4,
                         family                       = instance$sim_settings$distribution,
                         incompleteness               = incompleteness,
                         lambda_inc                   = lambda_inc,
                         fpca_type                    = "two-step",
                         fpca_index_significantDigits = 4,
                         npc                          = instance$sim_settings$npc,
                         # npc_criterion                = c(.9, .02),
                         max_iterations               = max_iterations,
                         gradient                     = ifelse(instance$sim_settings$distribution == "gaussian", TRUE, FALSE),
                         cores                        = n_cores)
  })
  
  # explained share of variance
  FPCs_varExpl <- reg$fpca_obj$evalues / reg$fpca_obj$evalues_sum
  
  
  
  # 1) evaluation of the convergence ----------------------------------------
  convergence <- reg$convergence
  

  # 2) evaluation of the overall performance --------------------------------
  # --- 2.1) Performance of the estimation of the curve-specific mean structure
  #          (global mean + curve-specific FPC structure)
  MISE_curves <- eval_meanStructure(dat_sim          = instance$data,
                                    reg_obj          = reg,
                                    registr_fpcaType = "two-step")
  
  # --- 2.2) Performance of the estimation of the global mean
  #          (without the curve-specific FPC structure)
  MSE_globalMean <- eval_globalMean(dat_sim          = instance$data,
                                    distribution     = instance$sim_settings$distribution,
                                    reg_obj          = reg,
                                    registr_fpcaType = "two-step")
  
  
  # 3) evaluation of the performance wrt. amplitude variation ---------------
  # ----- 3.1) weighted MISE for the estimation of the FPCs
  resList_FPC <- eval_FPCstructure(FPCdat_sim      = instance$FPC_data,
                                   reg_obj         = reg,
                                   eigenvalues_sim = instance$sim_settings$FPC_eigenvalues)
  
  # ----- 3.2) Larsson-Villani measure for evaluating the estimated FPC basis
  LV_FPCA <- eval_FPCstructure_LV(FPCdat_sim = instance$FPC_data,
                                  reg_obj    = reg)
  
  
  # 4) evaluation of the performance wrt. phase variation -------------------
  # ----- 4.1) Performance of the estimation of the domain lengths
  dat_trueDomainLengths <- instance$data %>%
    select(-index) %>% dplyr::rename(index = index_raw)
  
  MSE_domainLengths <- calcMSE_domainLengths(dat_trueDomainLengths, reg$Y)
  
  
  # ----- 4.2) Performance of the estimation of the warping functions
  MISE_warpings <- eval_warpingFunctions(dat_sim = instance$data,
                                         reg_obj = reg)
  
  
  # return results ----------------------------------------------------------
  return(list(convergence       = convergence,
              FPCs_nExtracted   = reg$fpca_obj$npc,
              FPCs_varExpl      = FPCs_varExpl,
              MISE_curves       = MISE_curves,
              MSE_globalMean    = MSE_globalMean,
              wMISE_FPCA                  = resList_FPC$wMISE_FPCA,
              wMISE_FPCA_onlyFirstFPC     = resList_FPC$wMISE_FPCA_onlyFirstFPC,
              wMISE_FPCA_onlyFirstTwoFPCs = resList_FPC$wMISE_FPCA_onlyFirstTwoFPCs,
              LV_FPCA           = LV_FPCA,
              MSE_domainLengths = MSE_domainLengths,
              MISE_warpings     = MISE_warpings,
              runtime           = runtime[[3]]))
}



#' Wrapper to call FGAMM with assumed completeness
#' 
#' Wrapper function to call FGAMM with Gaussian or Gamma family and assumed
#' completeness in the simulation study. Written specifically for the use with
#' the \code{batchtools} package.
#' 
#' @inheritParams FGAMMWithIncompleteness_wrapper
#' 
#' @export
#' 
#' @return List with performance criteria.
#' 
FGAMM_wrapper <- function(data, job, instance, max_iterations, n_cores) {
  
  FGAMMWithIncompleteness_wrapper(data           = data,
                                  job            = job,
                                  instance       = instance,
                                  incompleteness = FALSE,
                                  lambda_inc     = NULL,
                                  max_iterations = max_iterations,
                                  n_cores        = n_cores)
}



#' Wrapper to call varEM with potentially assumed incompleteness
#' 
#' Wrapper function to call varEM with Gaussian family in the simulation study.
#' Written specifically for the use with the \code{batchtools} package.
#' 
#' @inheritParams FGAMMWithIncompleteness_wrapper
#' @param distribution_reg Optional distribution (one of \code{c("gaussian","gamma")})
#' to be used in the registration step of \code{\link[registr]{register_fpca}}
#' instead of the original distribution of the simulated data found in
#' \code{instance$sim_settings$distribution}.
#' 
#' @import registr
#' @export
#' 
#' @return List with performance criteria.
#' 
varEMWithIncompleteness_wrapper <- function(data, job, instance, lambda_inc = 0,
                                            max_iterations, n_cores, incompleteness = TRUE,
                                            distribution_reg = NULL) {
  
  # estimation step ---------------------------------------------------------
  if (incompleteness) {
    # if (instance$sim_settings$incompleteness_strength != "CC")
    incompleteness <- "trailing"
  } else {
    incompleteness <- NULL
  }
  
  if (is.null(distribution_reg))
    distribution_reg <- instance$sim_settings$distribution
  
  runtime <- system.time({
    reg <- register_fpca(Y                         = instance$data,
                         Kt                        = 8,
                         Kh                        = 4,
                         family                    = distribution_reg,
                         incompleteness            = incompleteness,
                         lambda_inc                = lambda_inc,
                         fpca_type                 = "variationalEM",
                         npc                       = instance$sim_settings$npc,
                         # npc_criterion             = .9,
                         max_iterations            = max_iterations,
                         cores                     = n_cores)
  })
  
  # explained share of variance
  if ("evalues_sum" %in% names(reg$fpca_obj)) { # estimation with npc_criterion
    FPCs_varExpl <- reg$fpca_obj$evalues / reg$fpca_obj$evalues_sum
  } else { # estimation with npc
    FPCs_varExpl <- rep(-1, reg$fpca_obj$npc)
  }
  
  
  # 1) evaluation of the convergence ----------------------------------------
  convergence <- reg$convergence
  
  
  # 2) evaluation of the overall performance --------------------------------
  # --- 2.1) Performance of the estimation of the curve-specific mean structure
  #          (global mean + curve-specific FPC structure)
  MISE_curves <- eval_meanStructure(dat_sim          = instance$data,
                                    reg_obj          = reg,
                                    registr_fpcaType = "variationalEM")
  
  
  # --- 2.2) Performance of the estimation of the global mean
  #          (without the curve-specific FPC structure)
  MSE_globalMean <- eval_globalMean(dat_sim          = instance$data,
                                    distribution     = distribution_reg,
                                    reg_obj          = reg,
                                    registr_fpcaType = "variationalEM")
  
  
  # 3) evaluation of the performance wrt. amplitude variation ---------------
  # ----- 3.1) weighted MISE for the estimation of the FPCs
  resList_FPC <- eval_FPCstructure(FPCdat_sim      = instance$FPC_data,
                                   reg_obj         = reg,
                                   eigenvalues_sim = instance$sim_settings$FPC_eigenvalues)
  
  # ----- 3.2) Larsson-Villani measure for evaluating the estimated FPC basis
  LV_FPCA <- eval_FPCstructure_LV(FPCdat_sim = instance$FPC_data,
                                  reg_obj    = reg)
  
  
  # 4) evaluation of the performance wrt. phase variation -------------------
  # ----- 4.1) Performance of the estimation of the domain lengths
  dat_trueDomainLengths <- instance$data %>%
    select(-index) %>% dplyr::rename(index = index_raw)
  
  MSE_domainLengths <- calcMSE_domainLengths(dat_trueDomainLengths, reg$Y)
  
  
  # ----- 4.2) Performance of the estimation of the warping functions
  MISE_warpings <- eval_warpingFunctions(dat_sim = instance$data,
                                         reg_obj = reg)
  
  
  # return results ----------------------------------------------------------
  return(list(convergence       = convergence,
              FPCs_nExtracted   = reg$fpca_obj$npc,
              FPCs_varExpl      = FPCs_varExpl,
              MISE_curves       = MISE_curves,
              MSE_globalMean    = MSE_globalMean,
              wMISE_FPCA                  = resList_FPC$wMISE_FPCA,
              wMISE_FPCA_onlyFirstFPC     = resList_FPC$wMISE_FPCA_onlyFirstFPC,
              wMISE_FPCA_onlyFirstTwoFPCs = resList_FPC$wMISE_FPCA_onlyFirstTwoFPCs,
              LV_FPCA           = LV_FPCA,
              MSE_domainLengths = MSE_domainLengths,
              MISE_warpings     = MISE_warpings,
              runtime           = runtime[[3]]))
}



#' Wrapper to call varEM with assumed completeness
#' 
#' Wrapper function to call varEM with Gaussian family and assumed completeness
#' in the simulation study. Written specifically for the use with the
#' \code{batchtools} package.
#' 
#' @inheritParams varEMWithIncompleteness_wrapper
#' 
#' @import registr
#' @export
#' 
#' @return List with performance criteria.
#' 
varEM_wrapper <- function(data, job, instance, max_iterations, n_cores, distribution_reg = NULL) {
  
  varEMWithIncompleteness_wrapper(data             = data,
                                  job              = job,
                                  instance         = instance,
                                  incompleteness   = FALSE,
                                  lambda_inc       = NULL,
                                  max_iterations   = max_iterations,
                                  n_cores          = n_cores,
                                  distribution_reg = distribution_reg)
}



#' Wrapper to call SRVF with assumed completeness
#' 
#' Wrapper function to call SRVF with assumed completeness in the simulation
#' study. Written specifically for the use with the \code{batchtools} package.
#' 
#' @inheritParams varEMWithIncompleteness_wrapper
#' 
#' @import registr
#' @export
#' 
#' @return List with performance criteria.
#' 
srvf_wrapper <- function(data, job, instance, max_iterations, n_cores) {
  
  # Note: Since the function can only handle complete curves on a regular grid,
  #       this function is only applicable to complete curve settings.
  #       Any observed grid is made regular by spline-based interpolation.
  
  # Note 2: A parallel call to fdasrvf::align_fPCA() lead to some errors
  #         on the server. Don't know why, but no errors showed up when only
  #         using one core. Accordingly, always use one core only
  n_cores <- 1
  
  npc <- instance$sim_settings$npc
  # set 'npc_criterion' to NULL if npc shouldn't be estimated based on the data
  npc_criterion <- NULL
  # npc_criterion <- .9
  
  
  # preparation: create regular grid ----------------------------------------
  regular_grid <- seq(min(instance$data$index), max(instance$data$index),
                      length.out = table(instance$data$id)[1])
  
  # linearly interpolate the randomly warped curves to obtain measurements on a
  # regular grid for align_fPCA.
  if (instance$sim_settings$distribution == "gamma") {
    family_mgcv <- mgcv::Tweedie(p = 2)
  } else { # gaussian
    family_mgcv <- gaussian()
  }
  dat_regular_list <- lapply(unique(instance$data$id), function(i) {
    dat_i     <- instance$data %>% filter(id == i)
    y_approx  <- approx(dat_i$index, y = dat_i$value, xout = regular_grid, method = "linear", rule = 2)
    y_regGrid <- y_approx$y
    return(data.frame(index     = regular_grid,
                      value     = y_regGrid,
                      id        = i,
                      row.names = NULL))
  })
  dat_regular <- dplyr::bind_rows(dat_regular_list)
  

  # estimation step ---------------------------------------------------------
  curve_matrix <- matrix(data = dat_regular$value, nrow = length(regular_grid),
                         byrow = FALSE)
  
  estimate_model <- function(n_FPCs) {
    align_fPCA(f        = curve_matrix,
               time     = regular_grid,
               num_comp = n_FPCs,
               MaxItr   = max_iterations,
               showplot = FALSE,
               parallel = FALSE) # only use one core; >1 cores lead to some errors on the server, don't know why
  }
  
  runtime <- system.time({
    # capture the messages since align_fPCA only prints the number of iterations to the console
    output <- capture.output({ reg <- estimate_model(npc) })
  })
  
  # choose the number of FPCs based on PVE.
  if (!is.null(npc_criterion)) {
    ev_varExpl <- reg$vfpca$latent / sum(reg$vfpca$latent)
    npc_pve    <- which(cumsum(ev_varExpl) >= npc_criterion)[1]
    
    if (npc_pve != npc) {
      npc     <- npc_pve
      runtime <- system.time({
        # capture the messages since align_fPCA only prints the number of iterations to the console
        output <- capture.output({ reg <- estimate_model(npc) })
      })
    }
  }

  # explained share of variance
  FPCs_varExpl <- reg$vfpca$latent[1:npc] / sum(reg$vfpca$latent)
  
  
  # 1) evaluation of the convergence ----------------------------------------
  n_iterations <- tail(output, 1)
  n_iterations <- as.numeric(substr(n_iterations, nchar(n_iterations), nchar(n_iterations)))
  
  convergence <- list(iterations = n_iterations)
  
  
  # 2) evaluation of the overall performance --------------------------------
  # --- 2.1) Performance of the estimation of the curve-specific mean structure
  #          (global mean + curve-specific FPC structure)
  MISE_curves <- eval_meanStructure(dat_sim          = instance$data,
                                    srvf_regularGrid = regular_grid,
                                    reg_obj          = reg)
  
  
  # --- 2.2) Performance of the estimation of the global mean
  #          (without the curve-specific FPC structure)
  MSE_globalMean <- eval_globalMean(dat_sim        = instance$data,
                                    dat_simRegular = dat_regular,
                                    distribution   = instance$sim_settings$distribution,
                                    reg_obj        = reg)
  
  
  # 3) evaluation of the performance wrt. amplitude variation ---------------
  # not applicable, since the approach doesn't yield orthonormal FPCs in the original space
  
  # ----- 3.1) weighted MISE for the estimation of the FPCs
  # resList_FPC <- eval_FPCstructure(FPCdat_sim       = instance$FPC_data,
  #                                  reg_obj          = reg,
  #                                  eigenvalues_sim  = instance$sim_settings$FPC_eigenvalues,
  #                                  srvf_regularGrid = regular_grid)
  
  # ----- 3.2) Larsson-Villani measure for evaluating the estimated FPC basis
  # LV_FPCA <- eval_FPCstructure_LV(FPCdat_sim       = instance$FPC_data,
  #                                 reg_obj          = reg,
  #                                 srvf_regularGrid = regular_grid)
    
  
  # 4) evaluation of the performance wrt. phase variation -------------------
  # ----- 4.1) Performance of the estimation of the domain lengths
  MSE_domainLengths <- 0 # align_fPCA is only applicable to complete curve settings
  
  
  # ----- 4.2) Performance of the estimation of the warping functions
  MISE_warpings <- eval_warpingFunctions(dat_sim          = instance$data,
                                         reg_obj          = reg,
                                         srvf_regularGrid = regular_grid)
  
  
  # return results ----------------------------------------------------------
  return(list(convergence       = convergence,
              FPCs_nExtracted   = npc,
              FPCs_varExpl      = FPCs_varExpl,
              MISE_curves       = MISE_curves,
              MSE_globalMean    = MSE_globalMean,
              wMISE_FPCA                  = NA, # resList_FPC$wMISE_FPCA,
              wMISE_FPCA_onlyFirstFPC     = NA, # resList_FPC$wMISE_FPCA_onlyFirstFPC,
              wMISE_FPCA_onlyFirstTwoFPCs = NA, # resList_FPC$wMISE_FPCA_onlyFirstTwoFPCs,
              LV_FPCA           = NA, # LV_FPCA,
              MSE_domainLengths = MSE_domainLengths,
              MISE_warpings     = MISE_warpings,
              runtime           = runtime[[3]]))
}
