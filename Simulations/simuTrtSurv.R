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
                               save_every = 100,          
                               save_prefix = "sim_ckpt",  
                               save_path = ".",
                               start_iter = NULL) {    
  
  # ensure save directory exists
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  success.count <- 0L
  success.iter  <- integer(0)  
  
  # starting index
  if (is.null(start_iter)) {
    i <- 1L
  } else {
    if (!is.numeric(start_iter) || length(start_iter) != 1 || start_iter < 1)
      stop("start_iter must be a positive scalar integer.")
    if (start_iter > length(datlist))
      stop("start_iter exceeds length(datlist).")
    i <- as.integer(start_iter)
    if (i > max.iter) warning("start_iter > max.iter; loop will not run unless you increase max.iter.")
    message(sprintf("Starting fresh at iteration i=%d", i))
  }
  
  # main loop (stop if we hit datlist end as well)
  while (success.count < n.success && i <= max.iter && i <= length(datlist)) {
    start <- Sys.time()
    dat.i <- datlist[[i]]
    
    result.i <- tryCatch(
      FuseSurv(
        data = dat.i,
        covar.name=c("X1","X2","X3"), 
        site.var="site", 
        trt.name="A", 
        time.var="Y", 
        event="Delta", 
        fit.times=1:180, 
        eval.times=c(30,60,90),
        prop.SL.library=c("SL.mean", "SL.glm"), 
        event.SL.library=c("survSL.km", "survSL.coxph"), 
        cens.SL.library=c("survSL.km", "survSL.coxph"),
        n.folds=5,
        s = i*11
      ),
      error = function(e) {
        message(paste0("Iteration ", i, " failed: ", e$message))
        NULL
      }
    )
    
    if (!is.null(result.i)) {
      success.count <- success.count + 1L
      results[[success.count]] <- result.i
      length(success.iter) <- success.count
      success.iter[success.count] <- i
      
      print(paste0("Success ", success.count, " (iter ", i, 
                   "): ", round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2), "s"))
      
      # checkpoint every 'save_every' successes
      if (save_every > 0L && success.count %% save_every == 0L) {
        ckpt_file <- file.path(save_path, sprintf("%s_%03d.RData", save_prefix, success.count))
        save(results, success.count, success.iter, file = ckpt_file)
        message(sprintf("Checkpoint saved: %s", ckpt_file))
      }
    }
    i <- i + 1L
  }
  
  if (success.count < n.success) {
    warning(paste("Stopped after", success.count, 
                  "successful runs (max.iter or data exhausted)."))
  }
  invisible(results)
}

load("obsdata_s.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_prefix="Res_homo_s")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_prefix="Res_diffX_s")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_prefix="Res_diffT_s")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_prefix="Res_diffC_s")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_prefix="Res_diffAll_s")

load("obsdata_l.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_prefix="Res_homo_l")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_prefix="Res_diffX_l")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_prefix="Res_diffT_l")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_prefix="Res_diffC_l")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_prefix="Res_diffAll_l")

load("obsdata_l2.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_prefix="Res_homo_l2")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_prefix="Res_diffX_l2")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_prefix="Res_diffT_l2")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_prefix="Res_diffC_l2")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_prefix="Res_diffAll_l2")

load("obsdata_limO.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_prefix="Res_homo_limO")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_prefix="Res_diffX_limO")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_prefix="Res_diffT_limO")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_prefix="Res_diffC_limO")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_prefix="Res_diffAll_limO")

load("obsdata_imbal.Rdata")
r.imbal <- simulation.trtsurv(datlist=dat.imbal, save_prefix="Res_imbal")
r.moreK <- simulation.trtsurv(datlist=dat.moreK, save_prefix="Res_moreK")
