
#' Wrapper for 'simulate_curves' for the simulation study
#' 
#' High-level wrapper function for \code{\link{simulate_curves}} for
#' simulating curves in the simulation study. Written specifically for the
#' use with the \code{batchtools} package.
#' 
#' @param data,job General arguments as requested by \code{batchtools}.
#' @param distribution One of \code{c("gaussian","gamma")}, specifying the
#' distribution of the data.
#' @param amplitude_rank One of \code{c("Rank 1","Rank 2-3","Rank 3-4")}.
#' @param incompleteness_strength One of \code{c("CC","WIC","SIC")},
#' specifying if complete curves (\code{"CC"}) should be simulated or curves
#' with weak incompleteness (\code{"WIC"}) or strong incompleteness
#' (\code{"SIC"}).
#' @param correlation_structure One of \code{c("ID","AP","AI")}, specifying
#' if amplitude, phase, and the level of incompleteness should be mutually
#' independent (\code{"ID"}), if amplitude and phase should be correlated
#' (\code{"AP"}) or if amplitude and the level of incompleteness should be
#' correlated (\code{"AI"}).
#' @param N Number of curves to simulate.
#' @param n_timeGrid Number of (regular) evaluation points over the complete
#' time grid.
#' 
#' @export
#' 
#' @return List with (i) the simulated data, (ii) a dataset containing the
#' underlying functional principal components, (iii) a list with the main
#' settings used for the data simulation.
#' 
simData_wrapper <- function(data, job, distribution, amplitude_rank,
                            incompleteness_strength, correlation_structure,
                            N, n_timeGrid) {
  
  # translate the simulation settings to some specifics
  if (amplitude_rank == "Rank 1") {
    npc            <- 1
    FPCA_structure <- list(n_FPCs = npc, eigenvalues = 1)
  } else if (amplitude_rank == "Rank 2-3") {
    npc            <- 3
    FPCA_structure <- list(n_FPCs = npc, eigenvalues = c(0.7,0.25,0.05))
  } else if (amplitude_rank == "Rank 3-4") {
    npc            <- 4
    FPCA_structure <- list(n_FPCs = npc, eigenvalues = c(0.4,0.3,0.2,0.1))
  }
  
  incompleteness      <- ifelse(incompleteness_strength == "CC", FALSE, TRUE)
  incompleteness_rate <- ifelse(incompleteness_strength == "WIC", 0.3, 0.6)
  
  corr_amplitude_phase          <- ifelse(correlation_structure == "AP", -1, 0)
  corr_amplitude_incompleteness <- ifelse(correlation_structure == "AI", 0.8, 0)
  
  
  # perform the simulation
  sim_list <- simulate_curves(N                             = N,
                              n_timeGrid                    = n_timeGrid,
                              distribution                  = distribution,
                              random_warping                = TRUE,
                              FPCA_structure                = FPCA_structure,
                              incompleteness                = incompleteness,
                              incompleteness_rate           = incompleteness_rate,
                              corr_amplitude_phase          = corr_amplitude_phase,
                              corr_amplitude_incompleteness = corr_amplitude_incompleteness)
  
  return(list(data          = sim_list$data,
              FPC_data      = sim_list$FPC_data,
              sim_settings  = list(distribution            = distribution,
                                   npc                     = npc,
                                   incompleteness_strength = incompleteness_strength,
                                   FPC_eigenvalues         = FPCA_structure$eigenvalues)))
}
