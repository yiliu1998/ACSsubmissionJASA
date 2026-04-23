library(CFsurvival)
library(survSuperLearner)
library(SuperLearner)
library(dplyr)
library(glmnet)
library(caret)
source("Extends.R")
source("EIFestimates.R")

simulation.extend <- function(datlist, 
                              n.success = 500, 
                              max.iter = 520,
                              save_every = 500,          
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
      FuseSurv_Extend(site = dat.i$site,
                      eval.times = dat.i$eval.times,
                      IF.00 = dat.i$IF.00, IF.01 = dat.i$IF.01,
                      S.00 = dat.i$S.00, S.01 = dat.i$S.01,
                      Aug.00.mean = dat.i$Aug.00.mean, Aug.01.mean = dat.i$Aug.01.mean,
                      Aug.R0.mean = dat.i$Aug.R0.mean, Aug.R1.mean = dat.i$Aug.R1.mean,
                      Aug.R0.mean.sour = dat.i$Aug.R0.mean.sour, 
                      Aug.R1.mean.sour = dat.i$Aug.R1.mean.sour,
                      IF.R0 = dat.i$IF.R0, IF.R1 = dat.i$IF.R1, 
                      IF.CCOD.0 = dat.i$IF.CCOD.0, IF.CCOD.1 = dat.i$IF.CCOD.1, 
                      ind.R1.ccod = dat.i$ind.R1.ccod,
                      s=4399),
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

load("homo_trtsurv.Rdata")
r.homo <- simulation.extend(datlist=results, save_prefix="homo")
load("diffX_trtsurv.Rdata")
r.diffX <- simulation.extend(datlist=results, save_prefix="diffX")
load("diffT_trtsurv.Rdata")
r.diffT <- simulation.extend(datlist=results, save_prefix="diffT")
load("diffC_trtsurv.Rdata")
r.diffC <- simulation.extend(datlist=results, save_prefix="diffC")
load("diffAll_trtsurv.Rdata")
r.diffAll <- simulation.extend(datlist=results, save_prefix="diffAll")

## Example
# data       = dat.diffT[[100]]
# covar.name = c("X1","X2","X3")
# site.var   = "site"
# trt.name   = "A"
# time.var   = "Y"
# event      = "Delta"
# fit.times  = 1:90
# eval.times = 1:45
# prop.SL.library  = c("SL.mean","SL.glm")
# event.SL.library = c("survSL.km","survSL.coxph")
# cens.SL.library  = c("survSL.km","survSL.coxph")
# n.folds = 5
# s = 1222
