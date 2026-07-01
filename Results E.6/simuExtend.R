source("Extends.R")
source("EIFestimates.R")
simulation.extend <- function(datlist, 
                              n.success = 500, 
                              max.iter = 520,
                              save_file = "sim_extend_results.RData",
                              save_path = ".",
                              start_iter = 1L) {    
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  results <- list()
  success.count <- 0L
  success.iter  <- integer(0)
  
  if (!is.numeric(start_iter) || length(start_iter) != 1 || start_iter < 1) {
    stop("start_iter must be a positive scalar integer.")
  }
  
  if (start_iter > length(datlist)) {
    stop("start_iter exceeds length(datlist).")
  }
  
  i <- as.integer(start_iter)
  
  if (i > max.iter) {
    warning("start_iter > max.iter; loop will not run unless you increase max.iter.")
  }
  
  while (success.count < n.success && i <= max.iter && i <= length(datlist)) {
    
    start <- Sys.time()
    dat.i <- datlist[[i]]
    
    result.i <- tryCatch(
      FuseSurv_Extend(
        site = dat.i$site,
        eval.times = dat.i$eval.times,
        IF.00 = dat.i$IF.00, 
        IF.01 = dat.i$IF.01,
        S.00 = dat.i$S.00, 
        S.01 = dat.i$S.01,
        Aug.00.mean = dat.i$Aug.00.mean, 
        Aug.01.mean = dat.i$Aug.01.mean,
        Aug.R0.mean = dat.i$Aug.R0.mean, 
        Aug.R1.mean = dat.i$Aug.R1.mean,
        Aug.R0.mean.sour = dat.i$Aug.R0.mean.sour, 
        Aug.R1.mean.sour = dat.i$Aug.R1.mean.sour,
        IF.R0 = dat.i$IF.R0, 
        IF.R1 = dat.i$IF.R1, 
        IF.CCOD.0 = dat.i$IF.CCOD.0, 
        IF.CCOD.1 = dat.i$IF.CCOD.1, 
        ind.R1.ccod = dat.i$ind.R1.ccod,
        s = 4399
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

### --- Run "simuTrtSurv.R" to get the simulation results on treatment-specific survival curves first
### --- Results after running are already saved, if you wish to re-run the results, uncomment the code below
# source("simuTrtSurv.R")

load("homo_trtsurv.Rdata")
r.homo <- simulation.extend(datlist=results, save_file="homo.RData")
load("diffX_trtsurv.Rdata")
r.diffX <- simulation.extend(datlist=results, save_file="diffX.RData")
load("diffT_trtsurv.Rdata")
r.diffT <- simulation.extend(datlist=results, save_file="diffT.RData")
load("diffC_trtsurv.Rdata")
r.diffC <- simulation.extend(datlist=results, save_file="diffC.RData")
load("diffAll_trtsurv.Rdata")
r.diffAll <- simulation.extend(datlist=results, save_file="diffAll.RData")
