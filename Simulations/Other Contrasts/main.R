survFusion <- function(data, 
                       covar.name, 
                       site.var,
                       tgt.name,
                       trt.name, 
                       time.var, 
                       event, 
                       fit.times, 
                       eval.times,
                       prop.SL.library, 
                       event.SL.library, 
                       cens.SL.library,
                       n.folds=5, 
                       s=4399) {
  
  trt.surv.results <- TrtSurvCurves(data=data, 
                                    covar.name=covar.name, 
                                    site.var=site.var,
                                    tgt.name=tgt.name,
                                    trt.name=trt.name, 
                                    time.var=time.var, 
                                    event=event, 
                                    fit.times=fit.times, 
                                    eval.times=eval.times,
                                    prop.SL.library=prop.SL.library, 
                                    event.SL.library=event.SL.library, 
                                    cens.SL.library=cens.SL.library,
                                    n.folds=n.folds, 
                                    s=s)
  
  extend.results <- FuseSurv_Extend(eval.times=trt.surv.results$eval.times,
                                    site=trt.surv.results$site,
                                    IF.00=trt.surv.results$IF.00,
                                    IF.01=trt.surv.results$IF.01,
                                    S.00=trt.surv.results$S.00, 
                                    S.01=trt.surv.results$S.01,
                                    Aug.00.mean=trt.surv.results$Aug.00.mean, 
                                    Aug.01.mean=trt.surv.results$Aug.01.mean,
                                    Aug.R0.mean=trt.surv.results$Aug.R0.mean, 
                                    Aug.R1.mean=trt.surv.results$Aug.R1.mean,
                                    Aug.R0.mean.sour=trt.surv.results$Aug.R0.mean.sour, 
                                    Aug.R1.mean.sour=trt.surv.results$Aug.R1.mean.sour,
                                    IF.R0=trt.surv.results$IF.R0, 
                                    IF.R1=trt.surv.results$IF.R1, 
                                    IF.CCOD.0=trt.surv.results$IF.CCOD.0, 
                                    IF.CCOD.1=trt.surv.results$IF.CCOD.1,
                                    ind.R1.ccod=trt.surv.results$ind.R1.ccod,
                                    s=10*s)
  
  return(list(df.TGT=trt.surv.results$df.TGT, 
              df.CCOD=trt.surv.results$df.CCOD, 
              df.FED=trt.surv.results$df.FED,
              
              df.RD.TGT=extend.results$df.RD.TGT, 
              df.RD.CCOD=extend.results$df.RD.CCOD, 
              df.RD.FED=extend.results$df.RD.FED,
              
              df.SR.TGT=extend.results$df.SR.TGT, 
              df.SR.CCOD=extend.results$df.SR.CCOD, 
              df.SR.FED=extend.results$df.SR.FED,
              
              df.RMST.0.TGT=extend.results$df.RMST.0.TGT, 
              df.RMST.0.CCOD=extend.results$df.RMST.0.CCOD, 
              df.RMST.0.FED=extend.results$df.RMST.0.FED,
              
              df.RMST.1.TGT=extend.results$df.RMST.1.TGT, 
              df.RMST.1.CCOD=extend.results$df.RMST.1.CCOD, 
              df.RMST.1.FED=extend.results$df.RMST.1.FED,
              
              df.RMST.diff.TGT=extend.results$df.RMST.diff.TGT, 
              df.RMST.diff.CCOD=extend.results$df.RMST.diff.CCOD, 
              df.RMST.diff.FED=extend.results$df.RMST.diff.FED,
              
              eval.times=eval.times))
}

