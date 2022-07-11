
library(dplyr)      # general data handling

library(batchtools) # framework for the simulation study

library(registr)    # our accompanying package
library(fdasrvf)    # SRVF approach

library(IncRegHelpers) # helper functions

### Note: If the simulation study should be run without fixing the number of FPCs,
#         but estimating it adaptively, the three rows in the helper file
#         'algorithm_wrappers.R' (under 'IncRegHelpers/R/') where 'npc_criterion'
#         is defined have to be uncommented manually.
sim_type <- "npcFixed" # either "npcFixed" or "npcSelection". Only for naming the results files.


# general settings
max_iterations      <- 100 # maximum number of joint iterations
n_cores             <- 20  # number of cores to use for parallelization
n_repls             <- 100 # number of replications of each simulation setting

# create combinations of all simulation settings
sim_design_settings    <- read.csv("simStudy_design_raw.csv", na.strings = "")
sim_design <- expand.grid(distribution            = sim_design_settings$distribution[!is.na(sim_design_settings$distribution)],
                          amplitude_rank          = sim_design_settings$amplitude_rank[!is.na(sim_design_settings$amplitude_rank)],
                          incompleteness_strength = sim_design_settings$incompleteness_strength[!is.na(sim_design_settings$incompleteness_strength)],
                          correlation_structure   = sim_design_settings$correlation_structure[!is.na(sim_design_settings$correlation_structure)],
                          stringsAsFactors        = FALSE) %>%
  filter(distribution == "t_strongerNoise")
sim_design$N <- 100         # number of curves
sim_design$n_timeGrid <- 50 # number of measurements per curve

# drop unrealistic combinations: complete curves with a correlation between amplitude and incompleteness
sim_design <- sim_design %>%
  filter(!(incompleteness_strength == "CC" & correlation_structure == "AI"))



# perform the simulations -------------------------------------------------
if (!checkmate::test_directory_exists("results_registry_t_strongerNoise")) {
  reg <- makeExperimentRegistry(file.dir = "results_registry_t_strongerNoise",
                                packages = c("registr"),
                                seed     = 2020)
}

reg <- loadRegistry("results_registry_t_strongerNoise", writeable = TRUE)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = n_cores)

addProblem(name = "joint_approach",
           fun  = simData_wrapper,
           seed = 2020)

addAlgorithm(name = "twoStepInc", fun = FGAMMWithIncompleteness_wrapper)
addAlgorithm(name = "twoStep",    fun = FGAMM_wrapper)
addAlgorithm(name = "registrInc", fun = varEMWithIncompleteness_wrapper)
addAlgorithm(name = "registr",    fun = varEM_wrapper)
addAlgorithm(name = "srvf",       fun = srvf_wrapper)

dat_algoSettings_general         <- data.frame(max_iterations  = max_iterations,
                                               n_cores         = n_cores)
dat_algoSettings_twoStep_IC      <- data.frame(lambda_inc      = .025,
                                               max_iterations  = max_iterations,
                                               n_cores         = n_cores)
# two-step and registr with incompleteness (runs only on the incomplete settings; for CC it would be identical to 'registr without incompleteness')
addExperiments(prob.designs = list(joint_approach = sim_design %>% filter(incompleteness_strength != "CC")),
               algo.designs = list(twoStepInc     = dat_algoSettings_twoStep_IC,
                                   registrInc     = dat_algoSettings_twoStep_IC),
               repls        = n_repls)
# two-step and registr without incompleteness (runs on all settings)
addExperiments(prob.designs = list(joint_approach = sim_design),
               algo.designs = list(twoStep        = dat_algoSettings_general,
                                   registr        = dat_algoSettings_general),
               repls        = n_repls)
# the SRVF approach only handles complete curves
addExperiments(prob.designs = list(joint_approach = sim_design %>% filter(incompleteness_strength == "CC")),
               algo.designs = list(srvf           = dat_algoSettings_general),
               repls        = n_repls)

ids_experiments <- findExperiments()
submitJobs(findNotDone(ids_experiments))
waitForJobs()



# gather the results of all simulations -----------------------------------
# look at result from one job
# loadResult(id = 2)

ids  <- findExperiments(prob.name = "joint_approach")
pars <- unwrap(getJobPars()) %>% as_tibble()
saveRDS(pars, file = paste0("../../results/simStudy_",sim_type,"_t_strongerNoise_pars.rds"))
res <- reduceResultsDataTable(ids=findDone(ids))
saveRDS(res, file = paste0("../../results/simStudy_",sim_type,"_t_strongerNoise.rds"))
