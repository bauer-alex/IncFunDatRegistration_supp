
library(dplyr)          # general data handling
library(microbenchmark) # benchmark functions

# data simulation
library(IncRegHelpers) # helper functions

# specify the directories where the package versions are installed
registrOld_dir <- "~/my_directory/registr_1.0.0/"
registrNew_dir <- "~/my_directory/registr_2.1.4/"
fdasrvf_dir    <- "~/my_directory/fdasrvf_1.9.4/"

# main settings
start_withSetting <- 1    # can be set to some number > 1 if the first few settings are already done
n_reps            <- 20   # number of replications for microbenchmark::microbenchmark()
max_iterations    <- 100  # max number of iterations of the joint approach
n_cores           <- 10   # only for fdasrvf and the new registr package version 

# simulation settings
sim_settings_gaussian <- data.frame(distribution   = "gaussian",
                                    incompleteness = "CC",
                                    N              = c(100, 100, 100, 1000, 1000, 1000, 2000, 3000),
                                    n_timeGrid     = c(50,  100, 150, 50,   100,  150,  50,   50))
sim_settings <- sim_settings_gaussian %>%
  mutate(distribution = "gamma") %>%
  dplyr::bind_rows(sim_settings_gaussian) %>%
  arrange(desc(distribution)) # perform the more expensive gamma computations at the end

# create a results dataset which gets filled throughout the code
results_dat <- data.frame(method         = character(),
                          distribution   = character(),
                          incompleteness = character(),
                          n_curves       = numeric(),
                          n_timeGrid     = numeric(),
                          n_cores        = numeric(),
                          runtime_min    = numeric(),
                          runtime_mean   = numeric(),
                          runtime_median = numeric(),
                          runtime_max    = numeric(),
                          iterations     = numeric(),
                          error_message  = character())



# main analyses -----------------------------------------------------------
if (start_withSetting > 1) {
  results_dat <- readRDS(paste0("../../results/simStudy_runtimeAnalysis_setting",start_withSetting - 1,".rds"))
}

# iterate over all simulation settings
for (i in start_withSetting:nrow(sim_settings)) {
  message("Run sim setting ", i, "/", nrow(sim_settings),
          " ('", sim_settings$distribution[i], "', '", sim_settings$incompleteness[i],
          "', N=", sim_settings$N[i], ", n_timeGrid=", sim_settings$n_timeGrid[i], ")")
  
  
  # simulate data -----------------------------------------------------------
  message("--- Simulating the data...")
  set.seed(07062021)
  sim_list <- simData_wrapper(distribution            = sim_settings$distribution[i],
                              amplitude_rank          = "Rank 2-3",
                              incompleteness_strength = sim_settings$incompleteness[i],
                              correlation_structure   = "ID",
                              N                       = sim_settings$N[i],
                              n_timeGrid              = sim_settings$n_timeGrid[i])
  dat <- sim_list$data
  
  # library(ggplot2)
  # ggplot(dat, aes(index, value, group = id)) + geom_line()
  
  
  
  # ensure that the relevant packages are not loaded in a newer version -----
  library(lme4)
  library(gamm4)
  library(registr)
  
  # unload the old package versions
  detach("package:registr", unload = TRUE)
  detach("package:gamm4",   unload = TRUE)
  detach("package:lme4",    unload = TRUE)
  
  
  
  # apply the methods -------------------------------------------------------
  ### call for the old registr version
  library("lme4",    lib.loc = registrOld_dir)
  library("gamm4",   lib.loc = registrOld_dir)
  library("registr", lib.loc = registrOld_dir)
  # packageVersion("registr"); packageVersion("gamm4"); packageVersion("lme4")
  
  #----- run registr with variationalEM FPCA
  if (sim_settings$distribution[i] == "gaussian") {
    
    message("--- Estimating 'registr old variationalEM complete'...")
    mb_old <- microbenchmark::microbenchmark({ try_old <- try({
      res_old <- register_fpca(Y = dat, family = "gaussian", npc = 3, max_iterations = max_iterations)
    }) }, times = n_reps)
    
    if (class(try_old) == "try-error")
      mb_old$time <- rep(0, n_reps)
    
    # add the results to the results dataset
    runtimes_dat <- t(as.data.frame(mb_old$time))
    row.names(runtimes_dat) <- NULL
    colnames(runtimes_dat) <- paste0("run",1:length(mb_old$time))
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(method = "registr old variationalEM complete",
                                               distribution   = sim_settings$distribution[i],
                                               incompleteness = sim_settings$incompleteness[i],
                                               n_curves       = sim_settings$N[i],
                                               n_timeGrid     = sim_settings$n_timeGrid[i],
                                               n_cores        = 1,
                                               runtime_min    = min(mb_old$time) / 10^9,
                                               runtime_mean   = mean(mb_old$time) / 10^9,
                                               runtime_median = median(mb_old$time) / 10^9,
                                               runtime_max    = max(mb_old$time) / 10^9,
                                               iterations     = res_old$iterations,
                                               error_message  = ifelse(class(try_old) == "try-error", as.character(try_old), "ran without error"),
                                               runtimes       = runtimes_dat))
    
  }
  
  
  
  ### call for the new registr version
  
  # unload the old package versions
  detach("package:registr", unload = TRUE)
  detach("package:gamm4",   unload = TRUE)
  detach("package:lme4",    unload = TRUE)
  
  library("lme4",    lib.loc = registrNew_dir)
  library("gamm4",   lib.loc = registrNew_dir)
  library("registr", lib.loc = registrNew_dir)
  # packageVersion("registr"); packageVersion("gamm4"); packageVersion("lme4")
  
  #----- run registr with variationalEM FPCA and without incompleteness
  if (sim_settings$distribution[i] == "gaussian") {
    message("--- Estimating 'registr new variationalEM complete'...")
    mb_new1 <- microbenchmark::microbenchmark({ try_new1 <- try({
      res_new1 <- register_fpca(Y = dat, family = "gaussian", fpca_type = "variationalEM",
                                npc = 3, cores = n_cores, max_iterations = max_iterations)
    }) }, times = n_reps)
    
    if (class(try_new1) == "try-error")
      mb_new1$time <- rep(0, n_reps)
    
    # add the results to the results dataset
    runtimes_dat <- t(as.data.frame(mb_new1$time))
    row.names(runtimes_dat) <- NULL
    colnames(runtimes_dat) <- paste0("run",1:length(mb_new1$time))
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(method = "registr new variationalEM complete",
                                               distribution   = sim_settings$distribution[i],
                                               incompleteness = sim_settings$incompleteness[i],
                                               n_curves       = sim_settings$N[i],
                                               n_timeGrid     = sim_settings$n_timeGrid[i],
                                               n_cores        = n_cores,
                                               runtime_min    = min(mb_new1$time) / 10^9,
                                               runtime_mean   = mean(mb_new1$time) / 10^9,
                                               runtime_median = median(mb_new1$time) / 10^9,
                                               runtime_max    = max(mb_new1$time) / 10^9,
                                               iterations     = res_new1$convergence$iterations,
                                               error_message  = ifelse(class(try_new1) == "try-error", as.character(try_new1), "ran without error"),
                                               runtimes       = runtimes_dat))
    
  }
  
  
  #----- run registr with variationalEM FPCA and incompleteness
  if (sim_settings$distribution[i] == "gaussian" & sim_settings$incompleteness[i] != "CC") {
    
    message("--- Estimating 'registr new variationalEM incomplete'...")
    mb_new2 <- microbenchmark::microbenchmark({ try_new2 <- try({
      res_new2 <- register_fpca(Y = dat, family = "gaussian", fpca_type = "variationalEM",
                                incompleteness = "trailing", lambda_inc = 1, npc = 3, cores = n_cores,
                                max_iterations = max_iterations)
    }) }, times = n_reps)
    
    if (class(try_new2) == "try-error")
      mb_new2$time <- rep(0, n_reps)
    
    # add the results to the results dataset
    runtimes_dat <- t(as.data.frame(mb_new2$time))
    row.names(runtimes_dat) <- NULL
    colnames(runtimes_dat) <- paste0("run",1:length(mb_new2$time))
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(method = "registr new variationalEM incomplete",
                                               distribution   = sim_settings$distribution[i],
                                               incompleteness = sim_settings$incompleteness[i],
                                               n_curves       = sim_settings$N[i],
                                               n_timeGrid     = sim_settings$n_timeGrid[i],
                                               n_cores        = n_cores,
                                               runtime_min    = min(mb_new2$time) / 10^9,
                                               runtime_mean   = mean(mb_new2$time) / 10^9,
                                               runtime_median = median(mb_new2$time) / 10^9,
                                               runtime_max    = max(mb_new2$time) / 10^9,
                                               iterations     = res_new2$convergence$iterations,
                                               error_message  = ifelse(class(try_new2) == "try-error", as.character(try_new2), "ran without error"),
                                               runtimes       = runtimes_dat))
  }
  
  
  #----- run registr with two-step FPCA and incompleteness
  message("--- Estimating 'registr new two-step incomplete'...")
  
  if (sim_settings$incompleteness[i] == "CC") {
    inc <- NULL
  } else {
    inc <- "trailing"
  }
  lambda_inc <- ifelse(sim_settings$distribution[i] == "gaussian", 1, 50)
  
  mb_new3 <- microbenchmark::microbenchmark({ try_new3 <- try({
    res_new3 <- register_fpca(Y = dat, family = sim_settings$distribution[i], fpca_type = "two-step",
                              incompleteness = inc, lambda_inc = lambda_inc, npc = 3, cores = n_cores,
                              max_iterations = max_iterations,
                              gradient = ifelse(sim_settings$distribution[i] == "gaussian", TRUE, FALSE)) # gradient is only available for gaussian
  }) }, times = n_reps)
  
  if (class(try_new3) == "try-error")
    mb_new3$time <- rep(0, n_reps)
  
  # add the results to the results dataset
  runtimes_dat <- t(as.data.frame(mb_new3$time))
  row.names(runtimes_dat) <- NULL
  colnames(runtimes_dat) <- paste0("run",1:length(mb_new3$time))
  results_dat <- dplyr::bind_rows(results_dat,
                                  data.frame(method = "registr new two-step incomplete",
                                             distribution   = sim_settings$distribution[i],
                                             incompleteness = sim_settings$incompleteness[i],
                                             n_curves       = sim_settings$N[i],
                                             n_timeGrid     = sim_settings$n_timeGrid[i],
                                             n_cores        = n_cores,
                                             runtime_min    = min(mb_new3$time) / 10^9,
                                             runtime_mean   = mean(mb_new3$time) / 10^9,
                                             runtime_median = median(mb_new3$time) / 10^9,
                                             runtime_max    = max(mb_new3$time) / 10^9,
                                             iterations     = res_new3$convergence$iterations,
                                             error_message  = ifelse(class(try_new3) == "try-error", as.character(try_new3), "ran without error"),
                                             runtimes       = runtimes_dat))
  
  
  
  ### call fdasrvf
  library("fdasrvf", lib.loc = fdasrvf_dir)
  
  if (sim_settings$distribution[i] == "gaussian" & sim_settings$incompleteness[i] == "CC") {
    
    # Note: Since the function can only handle complete curves on a regular grid,
    #       this function is only applicable to complete curve settings.
    #       Any observed grid is made regular by spline-based interpolation.
    
    # preparation: create regular grid
    regular_grid <- seq(min(dat$index), max(dat$index), length.out = table(dat$id)[1])
    
    # interpolate the randomly warped curves to obtain measurements on a regular
    # grid for align_fPCA.
    dat_regular_list <- lapply(unique(dat$id), function(i) {
      dat_i     <- dat %>% filter(id == i)
      model_i   <- gam(value ~ s(index, bs = "ps", k = 10, sp = 0), data = dat_i)
      y_regGrid <- predict.gam(model_i, newdata = data.frame(index = regular_grid))
      return(data.frame(index = regular_grid,
                        value = y_regGrid,
                        id    = i,
                        row.names = NULL))
    })
    dat_regular <- dplyr::bind_rows(dat_regular_list)
    
    curve_matrix <- matrix(data = dat_regular$value, nrow = length(regular_grid),
                           byrow = FALSE)
    
    message("--- Estimating 'fdasrvf'...")
    mb_srvf <- microbenchmark::microbenchmark({ try_srvf <- try({
      # capture the messages since align_fPCA only prints the number of iterations to the console
      output <- capture.output({
        res_srvf <- align_fPCA(f        = curve_matrix,
                               time     = regular_grid,
                               num_comp = 3,
                               MaxItr   = max_iterations,
                               showplot = FALSE,
                               parallel = TRUE,
                               cores    = n_cores)
      })
    }) }, times = n_reps)
    
    if (class(try_srvf) == "try-error") {
      mb_srvf$time <- rep(0, n_reps)
      n_iterations <- 0
    } else {
      iter_srvf <- tail(output, 1)
      iter_srvf <- as.numeric(substr(iter_srvf, nchar(iter_srvf), nchar(iter_srvf)))
    }
    
    
    # add the results to the results dataset
    runtimes_dat <- t(as.data.frame(mb_srvf$time))
    row.names(runtimes_dat) <- NULL
    colnames(runtimes_dat) <- paste0("run",1:length(mb_srvf$time))
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(method = "fdasrvf",
                                               distribution   = sim_settings$distribution[i],
                                               incompleteness = sim_settings$incompleteness[i],
                                               n_curves       = sim_settings$N[i],
                                               n_timeGrid     = sim_settings$n_timeGrid[i],
                                               n_cores        = n_cores,
                                               runtime_min    = min(mb_srvf$time) / 10^9,
                                               runtime_mean   = mean(mb_srvf$time) / 10^9,
                                               runtime_median = median(mb_srvf$time) / 10^9,
                                               runtime_max    = max(mb_srvf$time) / 10^9,
                                               iterations     = iter_srvf,
                                               error_message  = ifelse(class(try_srvf) == "try-error", as.character(try_srvf), "ran without error"),
                                               runtimes       = runtimes_dat))
  }
  
  
  # save the current temporary version of the results
  if (i > 1)
    file.remove(paste0("../../results/simStudy_runtimeAnalysis_setting",i-1,".rds"))
  saveRDS(results_dat, file = paste0("../../results/simStudy_runtimeAnalysis_setting",i,".rds"))
}


# save the results --------------------------------------------------------
# calculate the mean runtimes per iteration
results_dat <- results_dat %>%
  mutate(runtime_mean_perIteration   = runtime_mean / iterations,
         runtime_median_perIteration = runtime_median / iterations)

saveRDS(results_dat, file = "../../results/simStudy_runtimeAnalysis.rds")
