### --- packages and main function files
library(CFsurvival)
library(survSuperLearner)
library(SuperLearner)
library(dplyr)
library(glmnet)
library(caret)
library(foreach)
library(doParallel)

# Load the functions for treatment-specific survival estimation,
# extended causal contrasts, and influence-function calculations.
source("TrtSurvCurves.R")
source("Extends.R")
source("EIFestimates.R")

### --- generate observed data using DGP.R
source("DGP.R")

### --- number of cores used, also the batch size
### --- please customize the number
n.workers <- 10L # as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))


### --- run survival-curve estimation and extended contrasts in one pipeline

# Run the Section E.5 simulation pipeline.
# -----------------------------------------------------------------------------
# Inputs:
#   datlist    : list of simulated datasets
#   cl         : parallel cluster object
#   n.success  : target number of successful Monte Carlo replicates
#   max.iter   : maximum number of iterations attempted
#   save_file  : output RData file
#   save_path  : output directory
#   start_iter : first replicate to attempt
#   batch.size : number of replicates processed per parallel batch
#
# Output:
#   saves successful contrast-estimation results and replicate indices;
#   invisibly returns the list of successful results.
# -----------------------------------------------------------------------------
simulation.contrasts <- function(datlist, cl, n.success=500L, max.iter=600L,
                                 save_file="sim_contrasts.RData", save_path=".",
                                 start_iter=1L, batch.size=n.workers) {
  if (length(start_iter)!=1L || !is.numeric(start_iter) || start_iter<1L) stop("start_iter must be a positive scalar integer.")
  if (start_iter>length(datlist)) stop("start_iter exceeds length(datlist).")
  if (start_iter>max.iter) stop("start_iter exceeds max.iter.")
  
  dir.create(save_path, recursive=TRUE, showWarnings=FALSE)
  start_iter <- as.integer(start_iter); max.iter <- min(as.integer(max.iter), length(datlist))
  n.success <- as.integer(n.success); batch.size <- max(1L, as.integer(batch.size))
  results <- list(); success.iter <- integer(0)
  iter.list <- seq.int(start_iter, max.iter)
  batches <- split(iter.list, ceiling(seq_along(iter.list)/batch.size))
  
  for (batch.iter in batches) {
    if (length(results)>=n.success) break
    batch.data <- datlist[batch.iter]
    
    # Process one batch of Monte Carlo replicates in parallel.
    batch.results <- foreach(i=batch.iter, dat.i=batch.data, .inorder=TRUE) %dopar% {
      start <- Sys.time()
      stage <- "TrtSurvCurves"
      
      tryCatch({
        seed.i <- i*11L
        set.seed(seed.i)
        
        # First stage: construct cross-fitted survival and influence-function
        # quantities required for the causal contrasts.
        trt.fit <- TrtSurvCurves(
          data=dat.i,
          covar.name=c("X1","X2","X3"), site.var="site", tgt.name="0",
          trt.name="A", time.var="Y", event="Delta",
          fit.times  = 1:90,
          eval.times = 1:60,
          prop.SL.library  = c("SL.mean","SL.glm"),
          event.SL.library = c("survSL.km","survSL.coxph"),
          cens.SL.library  = c("survSL.km","survSL.coxph"),
          n.folds = 5,
          s = seed.i
        )
        
        trt.fit[c("df.TGT","df.FED","df.CCOD")] <- NULL
        
        # Second stage: estimate survival difference, survival ratio,
        # and RMST contrasts under TGT, CCOD, and FED.
        stage <- "FuseSurv_Extend"
        final.fit <- FuseSurv_Extend(
          site=trt.fit$site, eval.times=trt.fit$eval.times,
          IF.00=trt.fit$IF.00, IF.01=trt.fit$IF.01,
          S.00=trt.fit$S.00, S.01=trt.fit$S.01,
          Aug.00.mean=trt.fit$Aug.00.mean, Aug.01.mean=trt.fit$Aug.01.mean,
          Aug.R0.mean=trt.fit$Aug.R0.mean, Aug.R1.mean=trt.fit$Aug.R1.mean,
          Aug.R0.mean.sour=trt.fit$Aug.R0.mean.sour,
          Aug.R1.mean.sour=trt.fit$Aug.R1.mean.sour,
          IF.R0=trt.fit$IF.R0, IF.R1=trt.fit$IF.R1,
          IF.CCOD.0=trt.fit$IF.CCOD.0, IF.CCOD.1=trt.fit$IF.CCOD.1,
          ind.R1.ccod=trt.fit$ind.R1.ccod, s=seed.i
        )
        
        rm(trt.fit)
        list(success=TRUE, iter=i, result=final.fit,
             elapsed=round(as.numeric(difftime(Sys.time(), start, units="secs")), 2))
      }, error=function(e) {
        list(success=FALSE, iter=i, stage=stage, error=conditionMessage(e))
      })
    }
    
    # Retain successful replicates and report any failures.
    for (out in batch.results) {
      if (!out$success) {
        message(sprintf("Iteration %d failed in %s: %s", out$iter, out$stage, out$error))
      } else if (length(results)<n.success) {
        results[[length(results)+1L]] <- out$result
        success.iter <- c(success.iter, out$iter)
        message(sprintf("Success %d at iteration %d: %.2fs", length(results), out$iter, out$elapsed))
      }
    }
    
    rm(batch.data, batch.results)
    gc(verbose=FALSE)
  }
  
  success.count <- length(results)
  if (success.count<n.success) warning(sprintf("Only %d successful runs were completed.", success.count))
  
  # Save the successful contrast results and corresponding replicate indices.
  save_path_full <- file.path(save_path, save_file)
  save(results, success.count, success.iter, file=save_path_full)
  message(sprintf("Final contrast results saved to: %s", save_path_full))
  invisible(results)
}


### --- parallel backend

# Create a PSOCK cluster and initialize all workers with the same packages,
# working directory, and user-defined function files.
working.directory <- normalizePath(getwd())
master.libPaths <- .libPaths()
required.packages <- c("CFsurvival","survSuperLearner","SuperLearner","dplyr","glmnet","caret")
missing.packages <- required.packages[!vapply(required.packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing.packages)) stop("Missing packages: ", paste(missing.packages, collapse=", "))

cl <- parallel::makeCluster(n.workers, type="PSOCK")
doParallel::registerDoParallel(cl)

parallel::clusterExport(cl, c("working.directory","master.libPaths","required.packages"), envir=environment())
parallel::clusterEvalQ(cl, {
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
  setwd(working.directory); .libPaths(master.libPaths)
  invisible(lapply(required.packages, library, character.only=TRUE))
  source("TrtSurvCurves.R"); source("Extends.R"); source("EIFestimates.R")
  NULL
})


### --- final files for results on causal contrasts are generated

# Run all five simulation scenarios and save the resulting contrast estimates.
load("obsdata_s.Rdata")

simulation.contrasts(dat.homo, cl, save_file="homo_contrasts.RData")
gc(verbose=FALSE)

simulation.contrasts(dat.diffX, cl, save_file="diffX_contrasts.RData")
gc(verbose=FALSE)

simulation.contrasts(dat.diffT, cl, save_file="diffT_contrasts.RData")
gc(verbose=FALSE)

simulation.contrasts(dat.diffC, cl, save_file="diffC_contrasts.RData")
gc(verbose=FALSE)

simulation.contrasts(dat.diffAll, cl, save_file="diffAll_contrasts.RData")
gc(verbose=FALSE)

# Close the parallel backend after all scenarios are complete.
parallel::stopCluster(cl)
foreach::registerDoSEQ()
