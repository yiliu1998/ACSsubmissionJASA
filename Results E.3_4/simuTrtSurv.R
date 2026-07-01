gc()
library(CFsurvival)
library(survSuperLearner)
library(SuperLearner)
library(dplyr)
library(glmnet)
library(caret)
source("EIFestimates.R")
source("FuseSurv.R")

simulation.trtsurv <- function(datlist, 
                               n.success = 500, 
                               max.iter = 600,
                               save_file = "sim_results.RData",
                               save_path = ".",
                               start_iter = 1L) {    
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  results <- list()
  success.iter <- integer(0)
  success.count <- 0L
  
  i <- as.integer(start_iter)
  
  while (success.count < n.success && i <= max.iter && i <= length(datlist)) {
    
    start <- Sys.time()
    dat.i <- datlist[[i]]
    
    result.i <- tryCatch(
      FuseSurv(
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
        s = i * 11
      ),
      error = function(e) {
        message(sprintf("Iteration %d failed: %s", i, e$message))
        NULL
      }
    )
    
    if (!is.null(result.i)) {
      success.count <- success.count + 1L
      results[[success.count]] <- result.i
      success.iter[success.count] <- i
      
      elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2)
      message(sprintf("Success %d at iteration %d: %.2fs", 
                      success.count, i, elapsed))
    }
    
    i <- i + 1L
  }
  
  if (success.count < n.success) {
    warning(sprintf(
      "Only %d successful runs were completed before reaching max.iter or end of datlist.",
      success.count
    ))
  }
  
  save_path_full <- file.path(save_path, save_file)
  save(results, success.count, success.iter, file = save_path_full)
  message(sprintf("Results saved to: %s", save_path_full))
  
  invisible(results)
}

### --- Run DGP.R to generate observed data in all simulation settings
source("DGP.R")

### --- Run simulations of different settings, and save simulation results
load("obsdata_s.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_file="Res_homo_s.RData")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_file="Res_diffX_s.RData")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_file="Res_diffT_s.RData")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_file="Res_diffC_s.RData")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_file="Res_diffAll_s.RData")

load("obsdata_l.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_file="Res_homo_l.RData")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_file="Res_diffX_l.RData")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_file="Res_diffT_l.RData")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_file="Res_diffC_l.RData")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_file="Res_diffAll_l.RData")

load("obsdata_l2.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_file="Res_homo_l2.RData")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_file="Res_diffX_l2.RData")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_file="Res_diffT_l2.RData")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_file="Res_diffC_l2.RData")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_file="Res_diffAll_l2.RData")

load("obsdata_limO.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_file="Res_homo_limO.RData")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_file="Res_diffX_limO.RData")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_file="Res_diffT_limO.RData")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_file="Res_diffC_limO.RData")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_file="Res_diffAll_limO.RData")
