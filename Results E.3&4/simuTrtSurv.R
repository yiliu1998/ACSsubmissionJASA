### --- Packages and source functions preparation
gc()
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Please uncomment the following two lines to install CFsurvival and survSuperLearner packages from Github
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

pkg.names <- c("CFsurvival", "survSuperLearner", "SuperLearner", "dplyr", "glmnet", "caret", "doParallel", "foreach")
pkg.versions <- setNames(lapply(pkg.names, packageVersion), pkg.names)
r.version <- R.version.string
save(pkg.versions, r.version, file = "package_versions.RData")

source("EIFestimates.R")
source("FuseSurv.R")

### --- Run DGP.R to generate observed data in all simulation settings
### --- This may take less than 3 minutes to finish
source("DGP.R")

### --- Run and save simulation results for CLCOX method first, which takes less than 10 minutes
source("ClusterCox.R")

### --- Set up parallel computing 
### --- You can modify the number of cores (also the batch size) below
n.workers <- 8L # as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
cl <- parallel::makeCluster(n.workers, type = "PSOCK")
doParallel::registerDoParallel(cl)

working.directory <- getwd()
parallel::clusterExport(cl, varlist = "working.directory")

parallel::clusterEvalQ(cl, {
  setwd(working.directory)
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
simulation.trtsurv <- function(datlist, n.success = 500L, max.iter = 600L, save_file = "sim_results.RData", save_path = ".", start_iter = 1L, batch.size = 5L) {
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  
  n.success <- as.integer(n.success)
  max.iter <- min(as.integer(max.iter), length(datlist))
  start_iter <- max(1L, as.integer(start_iter))
  batch.size <- max(1L, as.integer(batch.size))
  
  results <- list()
  success.iter <- integer(0)
  success.count <- 0L
  i <- start_iter
  
  while (success.count < n.success && i <= max.iter) {
    batch.start <- Sys.time()
    batch.index <- seq.int(from = i, to = min(i + batch.size - 1L, max.iter))
    batch.data <- datlist[batch.index]
    
    batch.results <- foreach(iter = batch.index, dat.i = batch.data, .inorder = TRUE) %dopar% {
      iteration.start <- Sys.time()
      seed.i <- as.integer(iter * 11L)
      
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
    
    batch.elapsed <- round(as.numeric(difftime(Sys.time(), batch.start, units = "mins")), 2)
    message(sprintf("Completed iterations %d-%d; %d/%d successful runs; %.2f minutes.", min(batch.index), max(batch.index), success.count, n.success, batch.elapsed))
    
    i <- max(batch.index) + 1L
    rm(batch.data, batch.results)
    gc(verbose = FALSE)
  }
  
  if (success.count < n.success) warning(sprintf("Only %d successful runs were completed before reaching max.iter or the end of datlist.", success.count))
  
  save_path_full <- file.path(save_path, save_file)
  save(results, success.count, success.iter, file = save_path_full)
  message(sprintf("Results saved to: %s", save_path_full))
  
  invisible(results)
}

### --- Run simulations of different settings and save results
### --- Main setting with nk = 300
load("obsdata_s.Rdata")

simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_s.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_s.RData", batch.size = n.workers)

rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()

### --- Main setting with nk = 600 
load("obsdata_l.Rdata")

simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_l.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_l.RData", batch.size = n.workers)

rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()

### --- Main setting with nk = 1000 
load("obsdata_l2.Rdata")

simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_l2.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_l2.RData", batch.size = n.workers)

rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()

### --- Limited overlap setting
load("obsdata_limO.Rdata")

simulation.trtsurv(datlist = dat.homo, save_file = "Res_homo_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffX, save_file = "Res_diffX_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffT, save_file = "Res_diffT_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffC, save_file = "Res_diffC_limO.RData", batch.size = n.workers)
simulation.trtsurv(datlist = dat.diffAll, save_file = "Res_diffAll_limO.RData", batch.size = n.workers)

rm(dat.homo, dat.diffX, dat.diffT, dat.diffC, dat.diffAll)
gc()