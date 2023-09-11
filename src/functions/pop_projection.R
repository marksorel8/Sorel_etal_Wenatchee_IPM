pop_projection<-function(){

  if(file.exists(here("results","sim_list_2_19.Rdata"))){
    load(here("results","sim_list_2_19.Rdata"))
  }else{
    
   here::i_am("data/broodstock_removal.csv")
  
  #load model object
  setwd(here("src"))
  TMB::compile("Wen_spchk_IPM.cpp")
  dyn.load(dynlib("Wen_spchk_IPM"))
 load(here("results","ipm_fit.Rdata"))
   mod<-obj$mod
   fit<-obj$fit
  rm(obj)
## grab data 
  dat_IPM<-mod$env$data
  
  #*************************************
  # model is fit to redds in years 1996-2019
  # projection is years 2020-2069
  #*************************************
  mod$env$data$proj_years<-50
  dat_IPM$proj_years<-50
  sim<-mod$simulate()  # make one simulaion

  # set hatchery management control rules parameters
  mod$env$data$bs_prop<-c(0.30,0.30) ##proportion of return to collect for broodstock
  mod$env$data$Hmax<-c(400,100) ## maximum number of hatchery fish to release onto spawning grounds
  mod$env$data$BS_Max<-c(74,64) ## broodstock size
  mod$env$data$Hslope<-c(1,1)  ##number of hatchery fish to release per natural origin fish forecasted to return
  
  
  #draw posterior samples
##function to sample from multivariate normal defined by mean and precision matrix
  rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
    set.seed( random_seed )
    z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L = Matrix::Cholesky(prec, super=TRUE)
    z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z = as.matrix(z)
    return(mu + z)
  }
  ##posterior samples
  list_of_draws<-t(rmvnorm_prec(mod$env$last.par.best,fit$SD$jointPrecision,200,1234))
  
  
  
  # Values used in simulation
  n_sims_per_par<-100  ## number of simulations per parameter
  total_years<-dat_IPM$last_t+mod$env$data$proj_years ## total years (retrospecitve and prospective)
  
  ## make arrays to hold simulations results
  sim_out<- #spawners
    sim_pHOS<- #pHOS
    sim_p_female<- #proportion female
    array(NA,dim=c(dim(sim$S_hat),dim(list_of_draws)[1]*n_sims_per_par))
  
  sim_NOB<- #natural origin broodstock collected
    array(NA,dim=c(dim(sim$broodstock_proj),dim(list_of_draws)[1]*n_sims_per_par))
  
  J_pred<- #juvenile emigrants
    J_LW<- #lower Wenatchee juveniles
    J_bon<- #Boneville Dam Juveniles
    array(NA,dim=c((total_years-3),dat_IPM$n_s,dat_IPM$n_l,dim(list_of_draws)[1]*n_sims_per_par))
  
  A_tum<-array(NA,dim=c((total_years-3),dat_IPM$n_s,dat_IPM$n_l,dat_IPM$n_ages,dim(list_of_draws)[1]*n_sims_per_par)) # adults at Tumwater juveniels

  
  # get design data for multi-state model so we can calculate survival rates based on different posterior samples
  here::i_am("src/Wen_spchk_IPM.cpp")
  setwd(here("Wenatchee-survival"))
  here::i_am("src/Wen_MSCJS_re_4_cont_cov.r")
  
  source(here("src","Wen_MSCJS_re_4.R"))
  
  #Just the years with mark recapture data
  mscjs_dat<- make_dat(mark_file_CH,sites=c( "LWe_J",
                                             "McN_J",
                                             #"JDD_J",
                                             "Bon_J",
                                             #"Est_J",
                                             "Bon_A","McN_A",
                                             #"PRa_A","RIs_A",
                                             "Tum_A"),cont_cov=c(),length_bin = 5,doy_bin = 5,inc_unk = FALSE,exc_unk=TRUE,start_year = 2018,end_year = 2018)
  
  
  
  dat<-mscjs_dat$Phi.design.dat %>% mutate(across(sum_flow:last_col(),function(x){x<-0})) 
  
  load(here("mscjs_fit_all_years_no_redd_new_1_6.Rdata"))
  
  design<-glmmTMB::glmmTMB(mscjs_fit_all_years_no_redd$Phi.design.glmmTMB$allForm$formula, data=dat,dispformula = ~0,doFit=FALSE,contrasts = list(LH="contr.sum",stream="contr.sum",age_class="contr.sum"))
  
  
  phi_design_dat<-design$data.tmb$X %>% as.matrix()
  
  rm(design)
  
  # indexes of parameters of interest
  p_fem_loc<-which(names(mod$env$last.par.best)=="mu_fem") # proportion females
   
  PSM_loc<-which(names(mod$env$last.par.best)=="mu_pss")   # prespawn mortality
  
  doy_loc<-which(names(mod$env$last.par.best)=="beta_DOY") # effect of emigration doy on fry survival from natal-reach emigration to the mouth of the Wenatchee River
  
  RRS_loc<-which(names(mod$env$last.par.best)=="logit_RRS") # relative reproductive success
  
  
  #names of the parameters to go into sensitivity analysis
  par_names<-c(
    ## juvenile production function shape parameters
    paste(rep(c("alpha","gamma","Jmax"),each=4),rep(c("fry","summer","fall","smolt"),times=3),sep="_"),
    ## survival from emigration to the mouth of the Wenatchee
    paste(c("fry","summer","fall","smolt"),"time1",sep="_"),
    ## downstream migration survival
    paste(c("sub","smolt"),"dowstream",sep="_"),
    ## smolt to adult returns
    paste(c("sub","smolt"),"SAR",sep="_"),
    ## upstream survival
    paste(c("three","four","five"),"Upstream",sep="_"),
    ## proportion adults returning at ages three and five
    paste(rep(c("sub","smolt"),each=2),c("percent_three","percent_five"),sep="_"),
    "PSM", ## prespawn mortality
    "P_female", ## average proportion female
    "logit_RRS") ## relative reproductive success
  
  #array to hold outputs of sensitivity analysis
  sim_pars_array<-array(NA,dim=c(length(par_names),3,(dim(list_of_draws)[1]*n_sims_per_par)),dimnames = list(par_names,c("Chiwawa","Nason","White"),1:(dim(list_of_draws)[1]*n_sims_per_par)))
  

  ## load simulated environmental covariate values
  here::i_am("data/broodstock_removal.csv")
  load(file=here("results","proj_arrays_8_10.Rdata"))
  
  gc()
  
#run simulations (this takes around an hour)
  start<-Sys.time()
  count<-1
  for ( i in 1:(dim(list_of_draws)[1])){ ## loop over posterior samples
    for(j in 1:n_sims_per_par){          ## loop over simulations per posterior samples
      
      # environmental covariates
      mod$env$data$X_proj <-( proj_arrays$x_proj_array[,,j] %>% as("dgTMatrix"))
      mod$env$data$X_phi_proj <- (proj_arrays$phi_proj_array[,,j] %>% as("dgTMatrix"))

      #conduct simulation with TMB model object using particular posterior samples
      sim<-mod$simulate(par=list_of_draws[i,])
      sim_out[,,count]<-sim$S_hat   #save natural-origin female spawners
      sim_NOB[,,count]<-sim$broodstock_proj # save projected broodstock collected
      
      
      # fill arrays with values of interes
      count_ts<-0
      count_ts2<-0
      count_tsl<-0
      for(t in 1:(total_years)){ ## loop over years
        for(s in 1:dat_IPM$n_s){ ## loop over natal streams
          if(t<(dat_IPM$first_t[s]+1)){next}
          
          count_ts<-count_ts+1
          sim_pHOS[t,s,count]<-sim$pHOS[count_ts]               #pHOS
          
       
            sim_p_female[t,s,count]<-sim$p_female[count_ts]       #pFemale
       
          
          if(t<=(total_years-3)){#end juvenile production calculation after fourth to last year
            for(l in 1:dat_IPM$n_l){
              count_tsl<-count_tsl+1
              J_pred[t,s,l,count]<-sim$J_pred[count_tsl] #Juvenile emigrants
              J_LW[t,s,l,count]  <-sim$J_LW[count_tsl]   #Juveniles at Lower Wenatchee
              J_bon[t,s,l,count] <-sim$JBon[count_tsl]   #Juveniles at Boneville
              
              for(a in 1:dat_IPM$n_ages){
                A_tum[t,s,l,a,count] <-sim$A_tum [count_tsl,a] #natural-origin adults at Tumwater
              }
              
            }
          }
        }
      } 
      
      
      #Juvenile production function parameters
      sim_pars_array[1:4,,count]<-sim$alpha
      sim_pars_array[5:8,,count]<-sim$gamma
      sim_pars_array[9:12,,count]<-sim$Jmax
      
      #Multistate model
      dat2<-dat %>% mutate(logit_phi=c(phi_design_dat %*% sim$beta_phi
      )) %>%  rbind(filter(.,time==1&LH=="summer") %>% mutate(logit_phi=logit_phi+dat_IPM$spr0_DOY*list_of_draws[i,doy_loc]) %>% mutate(LH="fry")) %>% mutate(phi=plogis(logit_phi))%>% mutate(LH=fct_relevel(LH,c("fry","summer","fall","smolt"))) %>% arrange(stream,LH)
      
      #time1 survival (emigration to Wenatchee mouth)
      sim_pars_array[13:16,,count]<-dat2 %>% filter(time==1) %>% select(LH,stream,phi) %>% pull(phi)
      #downstream survival
      sim_pars_array[17:18,,count]<-dat2 %>% filter(time%in%2:3 & LH %in%c("fall","smolt")) %>% select(LH,stream,phi,time) %>% group_by(stream,LH) %>% dplyr::summarize(phi=prod(phi),.groups ="keep") %>% pull(phi)
      #SAR
      sim_pars_array[19:20,,count]<-dat2 %>% filter(time%in%4 & LH %in%c("fall","smolt")) %>% select(LH,stream,phi) %>% pull(phi)
      #upstream
      sim_pars_array[21:23,,count]<-dat2 %>% filter(time%in%5:6 & LH %in%c("smolt"),stream=="Chiwawa") %>% select(stratum,phi,time) %>% group_by(stratum) %>% dplyr::summarize(phi=prod(phi)) %>% pull(phi)
      #percent three and five
      exp_alr_psi<-exp(sim$beta_psi[1:2])
      exp_alr_psi_sub<-exp(sim$beta_psi[1:2]+sim$beta_psi[3:4])
      ##subyearling emigrants
      sim_pars_array[24:25,,count]<-
        c(exp_alr_psi_sub)/(sum(exp_alr_psi_sub)+1)
      ##yearling emigrants
      sim_pars_array[26:27,,count]<-
        c(exp_alr_psi)/(sum(exp_alr_psi)+1)
      #PSM
      sim_pars_array[28,,count]<-plogis(list_of_draws[i,PSM_loc])
      #P_female
      sim_pars_array[29,,count]<-plogis(list_of_draws[i,p_fem_loc])
      #P_female
      sim_pars_array[30,,count]<-plogis(list_of_draws[i,RRS_loc])
      
      count<-count+1
    }
    # print(i)
    # gc()
  }
  end<-Sys.time()
  end-start

  #list of all of the simulation values
  sim_list<-list(sim_pars_array=sim_pars_array,
                 sim_out=sim_out,
                 sim_pHOS=sim_pHOS,
                 sim_p_female=sim_p_female,
                 A_tum=A_tum,
                 sim_NOB=sim_NOB)
  
  # save
  save(sim_list,file=here("results","sim_list_2_19.Rdata"))
  }


  return(sim_list)
}

