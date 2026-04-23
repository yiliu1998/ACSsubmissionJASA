load("obsdata_main.Rdata")
library(CFsurvival)
library(survSuperLearner)
library(SuperLearner)
library(dplyr)
library(glmnet)
library(caret)
source("TrtSurvCurves.R")
source("EIFestimates.R")

simulation.trtsurv <- function(datlist, 
                               n.success = 500, 
                               max.iter = 520,
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
      TrtSurvCurves(
        data = dat.i,
        covar.name = c("X1","X2","X3"),
        site.var   = "site",
        tgt.name   = "0",
        trt.name   = "A",
        time.var   = "Y",
        event      = "Delta",
        fit.times  = 1:60,
        eval.times = 1:60,
        prop.SL.library  = c("SL.mean","SL.glm","SL.glm.interaction"),
        event.SL.library = c("survSL.km","survSL.coxph","survSL.weibreg"),
        cens.SL.library  = c("survSL.km","survSL.coxph","survSL.weibreg"),
        n.folds = 5,
        s = i*123
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

r.homo <- simulation.trtsurv(datlist=dat.homo, save_prefix="homo")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_prefix="diffX")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_prefix="diffT")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_prefix="diffC")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_prefix="diffAll")
