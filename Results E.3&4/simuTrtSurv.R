# ============================================================
# Main simulation driver for results in E.3 and E.4
#
# This script:
#   (1) loads all required packages and records package/R versions;
#   (2) sources the user-defined estimation functions;
#   (3) generates the simulated datasets and clustered Cox results;
#   (4) sets up parallel computing; and
#   (5) runs the proposed and comparator estimators across all
#       Monte Carlo simulation settings.
# ============================================================

### --- Packages and source functions preparation
# Clear unused objects from memory before running the simulation workflow.
gc()

# Install devtools if it is not already available.
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Please uncomment the following two lines to install CFsurvival and survSuperLearner packages from Github.
# devtools::install_github("tedwestling/CFsurvival")
# devtools::install_github("tedwestling/survSuperLearner")
library(CFsurvival)
library(survSuperLearner)
library(SuperLearner)
library(dplyr)
library(glmnet)
library(caret)
library(doParallel)
library(foreach)

# Record package and R versions to facilitate computational reproducibility.
pkg.names <- c("CFsurvival", "survSuperLearner", "SuperLearner", "dplyr", "glmnet", "caret", "doParallel", "foreach")
pkg.versions <- setNames(lapply(pkg.names, packageVersion), pkg.names)
r.version <- R.version.string
save(pkg.versions, r.version, file = "package_versions.RData")

# Load user-defined functions for influence-function-based survival estimation
# and the main federated survival estimation procedure.
source("EIFestimates.R")
source("FuseSurv.R")

### --- Run DGP.R to generate observed data in all simulation settings

# Generate and save all simulated datasets used below.
# Approximate runtime on the authors' computing environment is less than 3 minutes.
source("DGP.R")

### --- Run and save simulation results for CLCOX method first

# Run the clustered Cox benchmark separately.
# Approximate runtime on the authors' computing environment is less than 10 minutes.
source("ClusterCox.R")


### --- Set up parallel computing 

# Number of parallel workers. This value can be modified based on the
# available computing resources; the same value is also used as the batch size below.
n.workers <- 8L # as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# Create and register a PSOCK cluster for foreach parallelization.
cl <- parallel::makeCluster(n.workers, type = "PSOCK")
doParallel::registerDoParallel(cl)

# Export the current working directory to each worker.
working.directory <- getwd()
parallel::clusterExport(cl, varlist = "working.directory")

# Initialize each parallel worker with the same working directory,
# thread settings, packages, and user-defined source functions.
parallel::clusterEvalQ(cl, {
  setwd(working.directory)
  
  # Restrict lower-level numerical libraries to one thread per worker
  # to avoid nested parallelism and excessive CPU usage.
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  
  suppressPackageStartupMessages({
    library(CFsurvival)
    library(survSuperLearner)
    library(SuperLearner)
    library(dplyr)
    library(glmnet)
    library(caret)
  })
  
  source("EIFestimates.R")
  source("FuseSurv.R")
  NULL
})


### --- Function for conducting simulations for other methods
### --- Simulations are conducted in parallel within each setting

# ------------------------------------------------------------
# Run the survival estimators over a list of simulated datasets
#
# Inputs:
#   datlist    : list of independently generated Monte Carlo datasets
#   n.success  : desired number of successful simulation replicates
#   max.iter   : maximum number of datasets to attempt
#   save_file  : file name used to save the simulation results
#   save_path  : directory in which the result file is stored
#   start_iter : index of the first Monte Carlo replicate to attempt
#   batch.size : number of replicates submitted to parallel workers at once
#
# Output:
#   saves the successful simulation results, number of successes,
#   and corresponding replicate indices to an RData file;
#   invisibly returns the list of successful estimation results
# ------------------------------------------------------------
simulation.trtsurv <- function(datlist, n.success = 500L, max.iter = 600L, save_file = "sim_results.RData", save_path = ".", start_iter = 1L, batch.size = 5L) {
  
  # Create the output directory if needed and standardize control arguments.
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  
  n.success <- as.integer(n.success)
  max.iter <- min(as.integer(max.iter), length(datlist))
  start_iter <- max(1L, as.integer(start_iter))
  batch.size <- max(1L, as.integer(batch.size))
  
  # Initialize storage for successful simulation results and replicate indices.
  results <- list()
  success.iter <- integer(0)
  success.count <- 0L
  i <- start_iter
  
  # Continue until the requested number of successful replicates is obtained
  # or the maximum number of available iterations is reached.
  while (success.count < n.success && i <= max.iter) {
    batch.start <- Sys.time()
    batch.index <- seq.int(from = i, to = min(i + batch.size - 1L, max.iter))
    batch.data <- datlist[batch.index]
    
    # Run all replicates in the current batch in parallel.
    batch.results <- foreach(iter = batch.index, dat.i = batch.data, .inorder = TRUE) %dopar% {
      iteration.start <- Sys.time()
      
      # Use a replicate-specific seed for reproducibility.
      seed.i <- as.integer(iter * 11L)
      
      # Fit all estimators through FuseSurv(). Failed replicates are recorded
      # rather than terminating the complete simulation run.
      tryCatch({
        set.seed(seed.i)
        
        result.i <- FuseSurv(
          data = dat.i,
          covar.name = c("X1", "X2", "X3"),
          site.var = "site",
          trt.name = "A",
          time.var = "Y",
          event = "Delta",
          fit.times = 1:180,
          eval.times = c(30, 60, 90),
          prop.SL.library = c("SL.mean", "SL.glm"),
          event.SL.library = c("survSL.km", "survSL.coxph"),
          cens.SL.library = c("survSL.km", "survSL.coxph"),
          n.folds = 5,
          s = seed.i
        )
        
        elapsed <- as.numeric(difftime(Sys.time(), iteration.start, units = "secs"))
        list(success = TRUE, iteration = iter, result = result.i, elapsed = elapsed, error = NA_character_)
      }, error = function(e) {
        elapsed <- as.numeric(difftime(Sys.time(), iteration.start, units = "secs"))
        list(success = FALSE, iteration = iter, result = NULL, elapsed = elapsed, error = conditionMessage(e))
      })
    }
    
    # Retain successful replicates until the target number is reached
    # and report any failed iterations.
    for (output.i in batch.results) {
      if (isTRUE(output.i$success) && success.count < n.success) {
        success.count <- success.count + 1L
        results[[success.count]] <- output.i$result
        success.iter[success.count] <- output.i$iteration
        message(sprintf("Success %d at iteration %d: %.2f seconds", success.count, output.i$iteration, output.i$elapsed))
      } else if (!isTRUE(output.i$success)) {
        message(sprintf("Iteration %d failed: %s", output.i$iteration, output.i$error))
      }
    }
    
    # Report progress at the end of each parallel batch.
    batch.elapsed <- round(as.numeric(difftime(Sys.time(), batch.start, units = "mins")), 2)
    message(sprintf("Completed iterations %d-%d; %d/%d successful runs; %.2f minutes.", min(batch.index), max(batch.index), success.count, n.success, batch.elapsed))
    
    # Advance to the next batch and free temporary objects from memory.
    i <- max(batch.index) + 1L
    rm(batch.data, batch.results)
    gc(verbose = FALSE)
  }
  
  # Warn if the requested number of successful runs could not be completed.
  if (success.count < n.success) warning(sprintf("Only %d successful runs were completed before reaching max.iter or the end of datlist.", success.count))
  
  # Save successful results and the original Monte Carlo replicate indices.
  save_path_full <- file.path(save_path, save_file)
  save(results, success.count, success.iter, file = save_path_full)
  message(sprintf("Results saved to: %s", save_path_full))
  
  return(invisible(results))
}

### --- Run simulations of different settings and save results
### --- Main setting with nk = 300
# Load datasets with target-site n=300 and source-site n=300.
load("obsdata_s.Rdata")

# Run 500 successful Monte Carlo replicates for each distribution-shift scenario.
simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_s.RData", batch.size = n.workers)

# Remove the current datasets from memory before loading the next sample-size setting.
rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()

### --- Main setting with nk = 600 
# Load datasets with target-site n=300 and source-site n=600.
load("obsdata_l.Rdata")

simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_l.RData", batch.size = n.workers)

rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()

### --- Main setting with nk = 1000 
# Load datasets with target-site n=300 and source-site n=1000.
load("obsdata_l2.Rdata")

simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_l2.RData", batch.size = n.workers)

rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()

### --- Limited overlap setting
# Load the simulation datasets with poorer treatment overlap in the target site.
load("obsdata_limO.Rdata")

simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_limO.RData", batch.size = n.workers)

rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()
