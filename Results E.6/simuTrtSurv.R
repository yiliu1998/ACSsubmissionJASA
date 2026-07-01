source("TrtSurvCurves.R")
source("EIFestimates.R")
simulation.trtsurv <- function(datlist, 
                               n.success = 500, 
                               max.iter = 520,
                               save_file = "sim_results.RData",
                               save_path = ".",
                               start_iter = 1L) {    
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  results <- list()
  success.count <- 0L
  success.iter  <- integer(0)
  
  i <- as.integer(start_iter)
  
  if (!is.numeric(start_iter) || length(start_iter) != 1 || start_iter < 1) {
    stop("start_iter must be a positive scalar integer.")
  }
  
  if (i > length(datlist)) {
    stop("start_iter exceeds length(datlist).")
  }
  
  if (i > max.iter) {
    warning("start_iter > max.iter; loop will not run unless you increase max.iter.")
  }
  
  while (success.count < n.success && i <= max.iter && i <= length(datlist)) {
    
    start <- Sys.time()
    dat.i <- datlist[[i]]
    
    result.i <- tryCatch(
      TrtSurvCurves(
        data = dat.i,
        covar.name = c("X1", "X2", "X3"),
        site.var   = "site",
        tgt.name   = "0",
        trt.name   = "A",
        time.var   = "Y",
        event      = "Delta",
        fit.times  = 1:60,
        eval.times = 1:60,
        prop.SL.library  = c("SL.mean", "SL.glm", "SL.glm.interaction"),
        event.SL.library = c("survSL.km", "survSL.coxph", "survSL.weibreg"),
        cens.SL.library  = c("survSL.km", "survSL.coxph", "survSL.weibreg"),
        n.folds = 5,
        s = i * 123
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
      "Only %d successful runs were completed before reaching max.iter or the end of datlist.",
      success.count
    ))
  }
  
  save_path_full <- file.path(save_path, save_file)
  save(results, success.count, success.iter, file = save_path_full)
  message(sprintf("Results saved to: %s", save_path_full))
  
  invisible(results)
}

### --- generate observed data using DGP.R file
### --- Results are already saved, if you wish to re-run, uncomment the following code
# source("DGP.R")
load("obsdata_main.Rdata")
r.homo <- simulation.trtsurv(datlist=dat.homo, save_file="homo_trtsurv.RData")
r.diffX <- simulation.trtsurv(datlist=dat.diffX, save_file="diffX_trtsurv.RData")
r.diffT <- simulation.trtsurv(datlist=dat.diffT, save_file="diffT_trtsurv.RData")
r.diffC <- simulation.trtsurv(datlist=dat.diffC, save_file="diffC_trtsurv.RData")
r.diffAll <- simulation.trtsurv(datlist=dat.diffAll, save_file="diffAll_trtsurv.RData")
