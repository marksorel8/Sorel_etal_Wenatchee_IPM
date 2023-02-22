# made TMB model object with data, initial parameter values, vector of which parameters to treat as random effects, and map of which parameters to fix at initial values. Then optimize. 
initialize_model<-function(){
  if(file.exists(here("results","ipm_fit.Rdata"))){
    load(here("results","ipm_fit.Rdata"))
  }else{
    
  library(here)
  library(TMB)
  library(TMBhelper)
  
  
  setwd(here("src"))
  TMB::compile("Wen_spchk_IPM.cpp")
  dyn.load(dynlib("Wen_spchk_IPM"))
  

#intialize
  mod<-TMB::MakeADFun(data=input$dat_IPM,parameters = input_test2,
                      random=c(unique(input$rand_par_IPM)),
                      map=input$map[1:2],DLL ="Wen_spchk_IPM",silent = FALSE,
                      inner.control = list(maxit = 1000))
  


  gc()
  fit<-TMBhelper::fit_tmb(mod,newtonsteps = 1,getsd = TRUE,getReportCovariance = TRUE,
                          getJointPrecision = TRUE)
  

  
  obj<- list(mod=mod,fit=fit)
  save(obj,file=here("results","ipm_fit.Rdata"))

  }
  return(obj)
}