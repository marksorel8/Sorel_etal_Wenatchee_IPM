#function to make data, initial parameters, names of parameters to be integrated out of likelihood (rand_par), and "map" (i.e., fixed parameters) for model.

make_data<-function(){
   library(here)
  if(file.exists(here("results","inputs.Rdata"))){
    load(here("results","inputs.Rdata"))
  }else{
  
  
   #libraries

  library(TMB)
  library(TMBhelper)
  library(tidyverse)
  library(viridisLite)
  library(readxl)
  
  
  #---------
  #spawner-to-Juvenile model
  here::i_am("src/Wen_spchk_IPM.cpp")
  setwd(here("Wenatchee-screw-traps"))
  here::i_am("src/covariates.R")
  

  

  
  ##   Load data ##
  {
    #read redd data from Hatchery program annual report (Hillman et al 2020)
    redds<-read_csv(here("data","redd_counts.csv")) %>% pivot_longer(c("Chiwawa","Nason","White"),names_to="stream",values_to="redds") 
    
    
    ##emigrant abundance estimates (log means and standard deviations from screw trap model). See https://github.com/Quantitative-Conservation-Lab/Wenatchee-screw-traps
    ## estimates of e
    load(here("all_emigrants_estimates.rds"))
    
    #load stream flow covariate values (flow_covs)
    source(here("src","covariates.r"))
    
    # length of habitat surveyed for spawners in each tributary (based on motoring annual report, Hillman et al. 2020)
    trib_lengths<-c(32.5,15.4,16.1)
    names(trib_lengths)<-c("Chiwawa","Nason","White")
  }
  
  
  # source helper functions to make data object, make initial parameter values, and makes map (which fixes certain parameter),  all to feed to the TMB model "Stock_recruit_LVM"
  source(here("src","SR_helper_functions.R"))

  
  #--------------------------------------------------------------------- 
  setwd(here("src","TMB"))
  TMB::compile("Stock_recruit_LVM.cpp")
  dyn.load(dynlib("Stock_recruit_LVM"))
  # # 
  if(file.exists(here("fit_mod_result.Rdata"))){
    load(here("fit_mod_result.Rdata"))
  }else{
    fit_mod_result<-fit_mod_iter(rep(c(1,1,1,1),times=3),streams=0:2,LHs=1:4,n_f=1,no_rand_FF =0, fit_env =1, fit_attempts=5,additional_attempts = 0,rate=log(c(5,5,5,5,5)),fold_exclude=integer(0),rate_map=factor(c(NA,NA,NA,NA,NA)))
    # 

    
   save(fit_mod_result,file=here("fit_mod_result.Rdata"))
  }

  #make sure it converged
  fit_mod_result$BIC_vec
  
  dat<-fit_mod_result$dat
  
  juv_mat<-tibble(t=dat$BY,s=dat$s_i,l=dat$l_i,mean=dat$J_obs,sd=dat$J_obs_sd)
  
  #exclude yearling emigrant abundance data from first year and subyearlings from last year in each stream (because we don't have estimates for subs and yrlngs for that brood year so we aren't going ot model it in the IPM)
  range_years_stream<-tapply(dat$BY,dat$s_i,range)
  
  #adjust number of years
  dat$n_t<-diff(range_years_stream[[1]])-1
  
  # year by stream indices to include
  ind<-!((dat$s_i==0& dat$BY%in%range_years_stream[[1]])|
           (dat$s_i==1& dat$BY%in%range_years_stream[[2]])|
           (dat$s_i==2& dat$BY%in%range_years_stream[[3]]))
  
  #number of combinations of stream by life history by year to model.
  dat$n_slt<-sum(ind)
  
  #order that they will be used in model
  juv_mat2<-juv_mat[ind,] %>% arrange(t,s,l)
  
  dat$J_obs<-juv_mat2 %>% pull(mean)
  dat$J_obs_sd<-juv_mat2 %>% pull(sd)
  
  #first and last year for each stream
  dat$first_t<-tapply(dat$t_i[ind] %>% droplevels %>% as.numeric(), dat$s_i[ind],min)-1 
  dat$last_t<-max(dat$t_i[ind] %>% droplevels %>% as.numeric())+3
  #brood years to model including three years after the last juvenile observation for each stream
  BY_spawners<-lapply(tapply(dat$BY[ind],dat$s_i[ind],range),function(x)seq(from=(x[1]),to=(x[2]+3)))
  stream<-sapply(1:3,function(x)rep(c("Chiwawa","Nason","White")[x],times=length(BY_spawners[[x]])))
  redd_dat<-tibble(Year=BY_spawners %>% unlist(),stream=stream %>% unlist) %>% left_join(redds) %>% arrange(Year,stream)#%>% pull(redds)
  
  dat$log_S_obs<-redd_dat %>% pivot_wider(names_from=stream,values_from=redds) %>% select(2:4) %>% as.matrix() %>% log()
  dat$n_st<-sum(!is.na(dat$log_S_obs))
  
  # make design data for juvenile emigrant abundance process errors
 
  juv_mat2<-juv_mat2 %>% mutate(BY=t,Y0=BY+1,
                                stream=c("Chiwawa","Nason","White")[(s+1)],
                                LH=c("fry","summer","fall","smolt")[(l+1)])
  
  
  Cov_dat<-juv_mat2 %>% left_join(pivot_longer(flow_covs$winter_high,2:4,names_to = "stream", values_to="win_0") %>% mutate(BY=Year)%>% select(-Year)) %>% #winter is measured as fall year!
    left_join(pivot_longer(flow_covs$summer_low ,2:4,names_to = "stream", values_to="sum_0") %>% mutate(BY=Year-1) %>% select(-Year)) %>%
    left_join(pivot_longer(flow_covs$winter_high,2:4,names_to = "stream", values_to="win_1") %>% mutate(BY=(Year-1))%>% select(-Year))    %>%
    mutate(sum_0=ifelse(LH=="fry",NA,sum_0),win_1=ifelse(LH=="smolt",win_1,NA)) %>% 
    group_by(stream,LH) 
  
  
  
  
  Cov_dat_Z<-Cov_dat %>% 
    group_by(stream,LH) %>% dplyr::summarize(BY=(BY),win_0=c(scale(win_0)),sum_0=c(scale(sum_0)),win_1=c(scale(win_1)),.groups="keep") %>% mutate(across(win_0:win_1,replace_na,0)) %>% mutate(LH=factor(LH),LH=fct_relevel(LH,c("fry","summer","fall","smolt")))%>% arrange(BY,stream,LH)
  
  X2<-model.matrix(~-1+win_0:LH+sum_0:LH+win_1:LH,data=Cov_dat_Z) %>% as_tibble() %>% select(-c(`LHfry:sum_0`,`LHfall:win_1`:`LHfry:win_1`,`LHsummer:win_1`)) %>% as.matrix()
  
  dat$X<-as(X2,"dgTMatrix")
  dat$s_i<-juv_mat2$s
  dat$l_i<-juv_mat2$l
  dat$t_i<-juv_mat2$t
  
  # save map
  mod_map<-fit_mod_result$mod$env$map
           
  #--------------------------------------------------------------------- 
  
  ##intital parameters for IPM 
  params<-fit_mod_result$mod$env$parList(par=fit_mod_result$mod$env$last.par.best)
  # params<-make_params_func(dat)        # make initial parameter values
  params$eta<-params$eta[ind]
  params$Omega_xf <- params$Omega_xf[-c(1,23),,drop=FALSE] 
  random<-c("eta","Omega_xf","eps_alpha","eps_gamma","eps_Jmax","beta_gamma","beta_Jmax","Loadings_vec","beta_e")
  
  #end of juvenile production data and initial parameters
  #---------------------------------------
  #            Begining of survival and return age data and parameters   
  #----------------------------------------
  here::i_am("src/Wen_spchk_IPM.cpp")
  setwd(here("Wenatchee-survival"))
  here::i_am("src/Wen_MSCJS_re_4_cont_cov.r")
  
  
  source(here("src","Wen_MSCJS_re_4_cont_cov.R"))
  
  #Just the years with mark recapture data
  mscjs_dat_doy<- make_dat(mark_file_CH,sites=c( "LWe_J",
                                                 "McN_J",
                                                 #"JDD_J",
                                                 "Bon_J",
                                                 #"Est_J",
                                                 "Bon_A","McN_A",
                                                 #"PRa_A","RIs_A",
                                                 "Tum_A"),cont_cov=c("rel_DOY_bin"),length_bin = 5,doy_bin = 5,inc_unk = FALSE,exc_unk=TRUE)#,start_year = 1998,end_year = 2018)
  
  
  setwd(here("src"))
  TMB::compile("wen_mscjs_re_4.cpp")
  dyn.load("wen_mscjs_re_4")
  
  if(!file.exists(here("mscjs_fit_doy.rdata"))){
    mscjs_fit_doy<-fit_wen_mscjs(x=mscjs_dat,
                                 #*survival*
                                 phi_formula=par.index~
                                   
                                   ## to Lower Wenatchee
                                   -1+
                                   time1+
                                   
                                   time1:LH +
                                   time2:age_class+
                                   time3:age_class+
                                   time4:age_class +
                                   time1:stream+
                                   time1:LH:stream+

                                   time1:win_air:LHfall+
                                   time1:win_air:LHsummer +
                                  
                                   diag(0+time1+time1:LH|mig_year)+
                                   
                                   ## to McNary
                                   time2+
                                   time2:stream+
                                   time2:age_class:stream+
                                   diag(0+time2 +time2:age_class|mig_year)+
                                   
                                   ##to Bonneville
                                   time3 +
                                   
                                   time3:stream+
                                   time3:age_class:stream+
                                   diag(0+time3+time3:age_class|mig_year)+
                                   
                                   ## SAR to Bonneville adult
                                   time4 +
                                   
                                   time4:stream+
                                   time4:age_class:stream+
                                   time4:ersstWAcoast.sum:age_class +
                                   time4:cui.spr:age_class +

                                   diag(0+time4+time4:age_class|mig_year)+
                                   
                                   
                                   ## to McNary Adult
                                   time5 +
                                   time5:stratum +
                                   diag(0+time5|mig_year)+
                                   
                                   ## to Tumwater Adults
                                   time6 +
                                   time6:stratum+
                                   diag(0+time6|mig_year)
                                 ,
                                 
                                 #*detection*
                                 p_formula= par.index~
                                   
                                   ## Lower Wenatchee Trap
                                   -1+
                                   time2+
                                   time2:LH+
                                   time3:age_class+
                                   time4:age_class+
                                   
                                   time2:stream +
                                   time2:LH:stream+
                                   diag(0+time2+time2:LH|mig_year) +
                                   ## McNary Juveniles
                                   time3 +
                                   
                                   time3:stream +
                                   time3:age_class:stream+
                                   
                                   time3:McN_flow:age_class +
                                   time3:McN_spill:age_class +

                                   diag(0+time3+time3:age_class|mig_year) +
                                   
                                   # Bonneville Juvenile
                                   time4 +
                                   
                                   time4:stream +
                                   time4:age_class:stream+
                                   
                                   time4:Bon_flow:age_class+
                                   time4:Bon_spill:age_class+
                                   diag(0+time4+time4:age_class|mig_year) +
                                   
                                   
                                   # Bonn adult
                                   time5 +
                                   time5:stratum +
                                   diag(0+time5 |mig_year)+
                                   
                                   #McNary adult
                                   time6+
                                   time6:stratum+
                                   diag(0+time6 |mig_year)
                                 ,
                                 #maturation age
                                 psi_formula= par.index~
                                   -1+
                                   tostratum +
                                   tostratum:age_class+
                                   us(0+tostratum|mig_year)
                                 ,
                                 
                                 doFit = TRUE,silent=FALSE,
                                 sim_rand =0,REML=FALSE,pen=log(rep(5,5)), hypersd=1,map_hypers=c(FALSE,FALSE))
    
    save(mscjs_fit_DOY,file=here("mscjs_fit_DOY.rdata"))
  }else{
    load(here("mscjs_fit_DOY.rdata"))
  }
  
  
  #day of year effect on summer emigrant time 1 survival 
  ##location in variable vector
  loc<-which(names(mscjs_fit_DOY$fit$SD$value)=="beta_phi")[
    (mscjs_fit_DOY$dat_TMB$X_phi %>% colnames())=="time1:LHsummer:rel_DOY_bin"]
  
  DOY_effect_prior<-c(mu=mscjs_fit_DOY$fit$SD$value[loc],
                      sd=mscjs_fit_DOY$fit$SD$sd[loc])
  
  #Mean day of year of emigration of spring subyearlings scaled by the mean and sd of the summer subyearlings in the multistate model used to estimate the effect of DOY on survival
  # DOY_means_and_sds<-all_bio_data %>% as_tibble() %>% group_by(LH) %>% dplyr::summarize(bio_dat_mean=mean(DOY)) %>% left_join(mscjs_dat_doy$mean_DOY) %>% left_join(mscjs_dat_doy$sd_DOY)
  
  
  spring_sub_doy_Z<-
    (100.0414 #fry average DOY
     -
       222.8923# summer emigrant average DOY
     ) /
    20.03752 #summer emigrant doy standard deviation
  
  
  #Years for IPM
  source(here("src","Wen_MSCJS_re_4.R"))
  mscjs_dat<- make_dat(mark_file_CH,sites=c( "LWe_J",
                                             "McN_J",
                                             #"JDD_J",
                                             "Bon_J",
                                             #"Est_J",
                                             "Bon_A","McN_A",
                                             #"PRa_A","RIs_A",
                                             "Tum_A"),cont_cov=c(),length_bin = 5,doy_bin = 10,inc_unk = FALSE,exc_unk=TRUE,start_year = 1998,end_year = 2018)
  
  if(!file.exists(here("mscjs_fit_all_years_no_redd_new_1_6.rdata"))){ 
    
    
    mscjs_fit_all_years_no_redd<-fit_wen_mscjs(x=mscjs_dat,
                                               #*survival*
                                               phi_formula=par.index~
                                                 
                                                 ## to Lower Wenatchee
                                                 -1+
                                                 time1+
                                                 
                                                 time1:LH +
                                                 time2:age_class+
                                                 time3:age_class+
                                                 time4:age_class +
                                                 time1:stream+
                                                 time1:LH:stream+
                                                 
                                                 time1:win_air:LHfall+
                                                 time1:win_air:LHsummer +

                                                 
                                               diag(0+time1+time1:LH|mig_year)+
                                                 
                                                 ## to McNary
                                                 time2+
                                                 
                                                 time2:stream+
                                                 time2:age_class:stream+

                                                 diag(0+time2 +time2:age_class|mig_year)+
                                                 
                                                 ##to Bonneville
                                                 time3 +
                                                 
                                                 time3:stream+
                                                 time3:age_class:stream+

                                                 diag(0+time3+time3:age_class|mig_year)+
                                                 
                                                 ## SAR to Bonneville adult
                                                 time4 +
                                                 
                                                 time4:stream+
                                                 time4:age_class:stream+
                                                 
                                                 time4:ersstWAcoast.sum:age_class +
                                                 time4:cui.spr:age_class +
                                                 
                                                 diag(0+time4+time4:age_class|mig_year)+
                                                 
                                                 
                                                 ## to McNary Adult
                                                 time5 +
                                                 time5:stratum +
                                                 diag(0+time5|mig_year)+
                                                 
                                                 ## to Tumwater Adults
                                                 time6 +
                                                 time6:stratum+
                                                 diag(0+time6|mig_year)
                                               ,
                                               
                                               #*detection*
                                               p_formula= par.index~
                                                 
                                                 ## Lower Wenatchee Trap
                                                 -1+
                                                 time2+
                                                 time2:LH+
                                                 time3:age_class+
                                                 time4:age_class+
                                                 
                                                 time2:stream +
                                                 time2:LH:stream+

                                                 diag(0+time2+time2:LH|mig_year) +

                                                 ## McNary Juveniles
                                                 time3 +
                                                 
                                                 time3:stream +
                                                 time3:age_class:stream+
                                                 
                                                 time3:McN_flow:age_class +
                                                 time3:McN_spill:age_class +

                                                 diag(0+time3+time3:age_class|mig_year) +
                                                 
                                                 # Bonneville Juvenile
                                                 time4 +
                                                 
                                                 time4:stream +
                                                 time4:age_class:stream+
                                                 
                                                 time4:Bon_flow:age_class+
                                                 time4:Bon_spill:age_class+
                                                 diag(0+time4+time4:age_class|mig_year) +
                                                 
                                                 
                                                 # Bonn adult
                                                 time5 +
                                                 time5:stratum +
                                                 diag(0+time5 |mig_year)+
                                                 
                                                 #McNary adult
                                                 time6+
                                                 time6:stratum+
                                                 diag(0+time6 |mig_year)
                                               ,
                                               #maturation age
                                               psi_formula= par.index~
                                                 -1+
                                                 tostratum +
                                                 tostratum:age_class+

                                                 us(0+tostratum|mig_year)
                                               ,
                                               
                                               doFit = TRUE,silent=FALSE,
                                               sim_rand =0,REML=FALSE,pen=log(rep(5,5)), hypersd=1,map_hypers=c(FALSE,FALSE))
    
    
    
    save(mscjs_fit_all_years_no_redd,file=here("mscjs_fit_all_years_no_redd_new_1_6.rdata"))
  }else{
    load(here("mscjs_fit_all_years_no_redd_new_1_6.rdata"))
  }
  
  par_fits<-mscjs_fit_all_years_no_redd$mod$env$parList(par=mscjs_fit_all_years_no_redd$last_par_best)
 
  colnames_X_phi<-colnames(mscjs_fit_all_years_no_redd$dat_TMB$X_phi)
  
  
  mscjs_fit_all_years_no_redd$dat_TMB$X_phi<-mscjs_fit_all_years_no_redd$dat_TMB$X_phi[,]

  #---------
  
  ##other data
  ### (integer) number of rows of design matrix corresponding with summer subyearling time 1 survival
  sum_1<-mscjs_dat$Phi.design.dat %>% filter(LH=="summer"&time==1)
  range(sum_1$par.index) %>% diff()
  n_sum_0_1<-nrow(sum_1) 
  max_ind<-mscjs_dat$Phi.design.dat$par.index %>% max()
  ### (integer) index of first row of design matrix corresponding with summer subyearling time 1 survival
  ind_sum_0_1<-sum_1 %>% pull(par.index) %>% min()
  
  
  #indices for para
  stream<-c("Chiwawa","Nason","White")
  LH<-c("summer","fall","smolt")
  sea_Year_p<-as.factor(1998:2018)
  Time<-0:5
  stratum<-as.factor(1:3)
  
  
  #create indices for survival rates
  phi_ind<-expand.grid(Time=Time,stream=stream,LH=LH,sea_Year_p=sea_Year_p,stratum=stratum) %>% filter(!(Time<=3 & stratum!=1)) %>% left_join(mscjs_dat$Phi.design.dat %>% select("Time","stream","LH","sea_Year_p","stratum","par.index")) %>% 
    pivot_wider(names_from = c("Time","stratum"),values_from = "par.index") %>% 
    bind_rows(filter(.,LH=="summer") %>% mutate(LH="fry",`0_1`=seq(from=(max_ind+1),length=n_sum_0_1 ))) %>% left_join(x=juv_mat2 %>% mutate(sea_Year_p=as.factor(BY+2)),y=.)
  
  
  
  #create indices for return  age proportions
  psi_ind<-expand.grid(stream=stream,LH=LH,mig_year=sea_Year_p) %>% left_join( mscjs_dat$Psi.design.dat %>% filter(tostratum==2) %>% select("stream","LH","mig_year","par.index")) %>% 
    bind_rows(filter(.,LH=="summer") %>% mutate(LH="fry")) %>% left_join(x=juv_mat2 %>% mutate(mig_year=as.factor(BY+2)),y=.)
  
  
  
  # adult demographics and removals
  here::i_am("data/broodstock_removal.csv")
  redd_dat<-redd_dat %>% 
    left_join(read_xlsx(here("data","adult data.xlsx"),sheet="Carcass_origin") %>% rename(Year=Brood_year,stream=Stream)) %>% 
    left_join(read_xlsx(here("data","adult data.xlsx"),sheet="fish_per_redd") %>% rename(Year=Brood_year) %>% distinct()) %>% 
    left_join(read.csv(here("data","broodstock_removal.csv")) %>% rename(Year=Brood_year)) %>% mutate(across(n_collected:n_5,.fns=replace_na,0))
  
  #exclude first 5 years in each stream (init years) for data on carcass recoveries
  redd_dat2<-redd_dat %>% group_by(stream) %>% mutate(min_y=min(Year)+5) %>% filter(Year>=min_y) %>% ungroup()
  
  WH_carc<-redd_dat %>% select(n_wild:n_hatchery) %>% rowSums()  # number of wild and hatchery carcasses recovered
  H_carc<-redd_dat %>% pull(n_hatchery) # number of hatchery carcasses recovered
  W_carc<-redd_dat %>% pull(n_wild)
  p_F<-redd_dat %>% mutate(p_F=1/fish_per_redd) %>% pull(p_F)         # proportion of female spawners
  
  
  broodstock<-redd_dat %>%select(n_3:n_5) %>% as.matrix()  # natural origin fish removed for broodstock
  
  #fecundity of age 5 relative to age 5
  fecund_dat<-readxl::read_xlsx(here::here("data","fecundity.xlsx"),1)
  
  fecund_sum<-fecund_dat %>% group_by(origin) %>% dplyr::summarize(total_N_4 = sum(age4N),
                                                                   mean_f_4=sum(age4mean*age4N)/sum(age4N),
                                                                   sd_f_4=sum(age4SD*age4N)/sum(age4N),
                                                                   total_N_5 = sum(age5N),
                                                                   mean_f_5=sum(age5mean*age5N)/sum(age5N),
                                                                   sd_f_5=sum(age5SD*age5N)/sum(age5N))
  
  
  
  
  
  f_5_prior<-c(log(fecund_sum %>% filter(origin=="Wild") %>% summarise(mean_f_5/mean_f_4-1) %>% as.numeric()),.1)
  
  # ---------------------------------------------------------------------------------------
  #                      Projection data
  #---------------------------------------------------------------------------------------
  
  proj_years<-50
  n_slt_proj<-proj_years*dat$n_s*dat$n_l
  if(!file.exists(here("results","proj_arrays_8_10.Rdata"))){
    n_sim<-100
  }else{
    n_sim<-1
  }
  
  #simulated future covariate values from MARSS models. Generated in "scripts/Covariate simulation appendix.Rmd" file
  load(here("Wenatchee-survival","MARSS_sim_8_10_2022.Rdata"))
  
  #means and sds of tributary covariates
  center_scale_trib<-Cov_dat %>% group_by(stream) %>%  summarise(across(win_0:win_1,.fns=list(means=base::mean,sds=stats::sd),na.rm=T))%>% pivot_longer(cols = 2:7,names_to = c("var","val"),names_sep = 6  ) %>% pivot_wider(values_from=value,names_from = val) %>% mutate(var=substr(var,1,5)) %>% mutate(var=as.factor(var),var=fct_relevel(var,c("win_0","sum_0","win_1"))) %>% arrange(var)
  
  #names of tributary flow covariates
  var_names<-colnames(MARSS_sim$all_env %>% select(c("Chiwawa_win_flow":"White_sum_flow")))
  
  #row position in projection matrix (array really) of tributary covariates
  var_pos<-which(colnames(MARSS_sim$all_env) %in%var_names)-1
  
  #Z-scored projected tributary covariates
  x_proj_array<-array(NA,dim=c(n_slt_proj,ncol(dat$X),n_sim))
  
  
  for ( i in 1:n_sim){
    
    
    proj_i<-cbind(as.numeric(c((flow_covs$winter_high %>% filter(Year==max(juv_mat2$BY)) %>% select(Chiwawa:White)),
                               c(flow_covs$summer_low %>% filter(Year==max(juv_mat2$BY)+1) %>% select(Chiwawa:White)),
                               rep(NA,3))) , #include win_0 and sum_0 from the last observed year in the projection
                  MARSS_sim$transformed_sim [c(var_pos,var_pos[1:3]), #replicate winter for scaling by winter_1 mean and sd
                                             ,i]) %>% 
      
      
      `-`(center_scale_trib$means) %>% 
      `/`(center_scale_trib$sds) %>% 
      #arrange intro a longer format with stream column
      t() %>% as_tibble() %>% `colnames<-`(c(paste(var_names,"0",sep="_"),paste(var_names[1:3],"1",sep="_"))) %>%  
      #add a year column
      mutate(Year=seq(from=max(juv_mat2$BY)+1,length.out=nrow(.))) %>% 
      #make longer
      pivot_longer(cols=-Year,names_to = c("stream",".value"),names_pattern = "(.*)_(.*_.*_.)",) %>%
      rename(win_0=win_flow_0,sum_0=sum_flow_0,win_1=win_flow_1)
    
    
    #matrix for projections
    juv_mat_proj<-  expand.grid(BY=seq(from=max(juv_mat2$BY)+1,length=proj_years),stream=unique(juv_mat2$stream),LH=unique(juv_mat2$LH)) %>% arrange(BY,stream) %>% 
      
      left_join(proj_i %>% select(stream,Year,win_0:sum_0) %>% mutate(BY=Year-1) %>% select(-Year)) %>% # winter flow IN PROJECTION is year when winter ended
      left_join(proj_i %>% select(stream,Year,win_1) %>% mutate(BY=Year-2) %>% select(-Year))  # winter flow IN PROJECTION is year when winter ended
    
    
    
    
    X_proj<-model.matrix(~-1+win_0:LH+sum_0:LH+win_1:LH,data=juv_mat_proj) %>% as_tibble() %>% select(-c(`LHfry:sum_0`,`LHfall:win_1`:`LHfry:win_1`,`LHsummer:win_1`)) %>% as.matrix()
    
    x_proj_array[,,i]<-X_proj
  }
  
  
  dat$X_proj<-as(X_proj,"dgTMatrix")
  
  
  

  #create phi and psi design data with env cov proj
  #create design matrices and terms structures for projections
  #extract the individual design data for each parameter
  ##phi/s survival
  
  #indices for para
  stream<-c("Chiwawa","Nason","White")
  LH<-factor(x=c("fall","smolt","summer"),levels=c("fall","smolt","summer"),ordered=T)
  sea_Year_p<-as.factor(1998:2018)
  Time<-0:5
  stratum<-as.factor(1:3)
  
  #create phi design data
  Phi.design.dat.proj<-expand.grid(Time=Time,stream=stream,LH=LH,sea_Year_p=seq(2019,by=1,length.out=proj_years),stratum=stratum)%>% filter(!(Time<=3 & stratum!=1)) %>% arrange(desc(Time),desc(LH)) %>% mutate(par.index=((1:nrow(.))+max_ind+n_sum_0_1)) %>% mutate(time=Time+1,) %>% 
    
    mutate(time=as.factor(time),stratum=as.factor(stratum)) %>% 
    
    mutate(mig_year=as.factor(ifelse(Time<=(3),as.numeric(as.character(sea_Year_p)) ,as.numeric(as.character(sea_Year_p)) +as.numeric(as.character(stratum)))),
           mig_year_num=as.numeric(mig_year)) %>%   #add a column for the actual migration year, which is the seaward year for juveniles and the seaward + stratum for adults
    cbind(.,model.matrix(~time+stream+LH+stratum-1,data=.)) %>% 
    
    mutate(age_class=as.factor(ifelse(LH%in%c("fall","summer"),"sub",as.character(LH))), 
           LHsummer=as.numeric(LH=="summer"),LHfall=as.numeric(LH=="fall"),LHsmolt=as.numeric(LH=="smolt"),age_0=as.numeric(LH!="smolt")) %>% 
    #make length bin and release DOY numeric
    arrange(desc(Time),desc(LH)) %>% 
    mutate(redd_stan=0) 
  
  #format simulated environmental covariate values
  ##names of tributary flow covariates
  var_names<-colnames(MARSS_sim$all_env %>% select(c("win_air":"cui.spr")))
  
  ##row position in projection matrix (array really) of tributary covariates
  var_pos<-which(colnames(MARSS_sim$all_env) %in%var_names)-1
  
  ##Z-scored projected tributary covariates
  
  phi_proj_array<-array(NA,dim=c(nrow(Phi.design.dat.proj),mscjs_fit_all_years_no_redd$dat_TMB$X_phi %>% ncol(),n_sim))
  
  for (i in 1:n_sim){
    
    proj_i<-     MARSS_sim$transformed_sim [c(var_pos), ,i] %>% 
      
      #subbtract mean of covariate values for years used in fitting 
      `-`(sapply(var_names,function(x)pull(mscjs_dat$Phi.design.dat,x) %>% attr("scaled:center"))) %>% 
      #devide by sd of covariate values for years used in fitting
      `/`(sapply(var_names,function(x)pull(mscjs_dat$Phi.design.dat,x) %>% attr("scaled:scale"))) %>% 
      ####arrange intro a longer format with stream column
      t() %>% as_tibble() %>% `colnames<-`(var_names) %>%  
      ####add a year column
      mutate(Year=seq(from=max(juv_mat2$BY)+2,length.out=nrow(.))) 
    
    
    #add environmental covariates to design data  
    Phi.design.dat.proj_env  <-left_join(Phi.design.dat.proj,proj_i %>% rename(sea_Year_p=Year))          
    
    #create_design data
    Phi.design.glmmTMB.proj<-glmmTMB::glmmTMB(mscjs_fit_all_years_no_redd$Phi.design.glmmTMB$allForm$formula, data=Phi.design.dat.proj_env,dispformula = ~0,doFit=FALSE,contrasts = list(LH="contr.sum",stream="contr.sum",age_class="contr.sum"))
    
    
    
    phi_proj_array[,,i]<-Phi.design.glmmTMB.proj$data.tmb$X [,]
  }
  
  
  
  #create indices for survival rates
  ### (integer) number of rows of design matrix corresponding with summer subyearling time 1 survival
  sum_1_proj<-Phi.design.dat.proj %>% filter(LH=="summer"&Time==0)
  range(sum_1_proj$par.index) %>% diff()
  n_sum_0_1_proj<-nrow(sum_1_proj) 
  max_ind_proj<-max(Phi.design.dat.proj$par.index)
  ### (integer) index of first row of design matrix corresponding with summer subyearling time 1 survival
  ind_sum_0_1_proj<-sum_1_proj %>% pull(par.index) %>% min()
  
  #create indices for survival rates
  phi_ind_proj<-Phi.design.dat.proj %>% select(Time,stream,LH,sea_Year_p,stratum,par.index) %>% pivot_wider(names_from = c("Time","stratum"),values_from = "par.index") %>% 
    bind_rows(filter(.,LH=="summer") %>% mutate(LH="fry",`0_1`=seq(from=(max_ind_proj+1),length=n_sum_0_1_proj ))) %>% mutate(LH=as.factor(LH),LH=fct_relevel(LH,c("fry","summer","fall","smolt"))) %>% arrange(sea_Year_p,stream,LH) %>% select(`0_1`:`3_1`,`4_1`,`5_1`, `4_2`, `5_2`, `4_3`, `5_3`)
  
  
  
  #design data for psi
  Psi.design.dat.proj<-expand.grid(stream=stream,LH=LH,sea_Year_p=seq(2019,by=1,length.out=proj_years),tostratum=c(2,3)) %>% arrange(tostratum,sea_Year_p,stream,LH) %>%  mutate(par.index=((1:nrow(.))+mscjs_dat$Psi.design.dat %>% filter(tostratum==2) %>% pull(par.index) %>% max)) %>% rename(mig_year=sea_Year_p) %>% mutate( age_class=as.factor(ifelse(LH!="smolt","sub","yrlng")),tostratum=as.factor(tostratum))
  
  #create indices for return  age proportions
  psi_ind_proj<-Psi.design.dat.proj %>%
    bind_rows(filter(.,LH=="summer") %>% mutate(LH="fry")) %>% mutate(LH=as.factor(LH),LH=fct_relevel(LH,c("fry","summer","fall","smolt"))) %>% filter(tostratum==2) %>% arrange(mig_year,stream,LH) %>% mutate()
  
  
  Psi.design.glmmTMB.proj<-glmmTMB::glmmTMB(mscjs_fit_all_years_no_redd$Psi.design.glmmTMB$allForm$formula, data=Psi.design.dat.proj,dispformula = ~0,doFit=FALSE,contrasts = list(LH="contr.sum",stream="contr.sum",age_class="contr.sum"))
  
  
  proj_dat<-list(proj_years = proj_years,
                 n_slt_proj = n_slt_proj,
                 X_proj = X_proj,
                 phi_terms_proj = Phi.design.glmmTMB.proj$data.tmb$terms ,
                 psi_terms_proj = Psi.design.glmmTMB.proj$data.tmb$terms ,
                 X_phi_proj = Phi.design.glmmTMB.proj$data.tmb$X[,] %>% as("dgTMatrix"),
                 X_psi_proj = Psi.design.glmmTMB.proj$data.tmb$X %>% as("dgTMatrix"),
                 Z_phi_proj = Phi.design.glmmTMB.proj$data.tmb$Z,
                 Z_psi_proj = Psi.design.glmmTMB.proj$data.tmb$Z,
                 n_sum_0_1_proj = n_sum_0_1_proj ,
                 ind_sum_0_1_proj = ind_sum_0_1_proj)
  
  if(!file.exists(here("results","proj_arrays_8_10.Rdata"))){
    
    proj_arrays<-list(phi_proj_array=phi_proj_array,
                      x_proj_array=x_proj_array)
    save(proj_arrays,file=here("results","proj_arrays_8_10.Rdata"))
  }
  #---------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------
  # combine data into one list
  
  dat_IPM<-c(dat,
             mscjs_fit_all_years_no_redd$dat_TMB,
             nUS_OCC = mscjs_fit_all_years_no_redd$dat_TMB$n_OCC-mscjs_fit_all_years_no_redd$dat_TMB$nDS_OCC-1,
             proj_dat,
             list(
               stream_length=c(32.5,15.4,16.1),
               n_sum_0_1=n_sum_0_1,
               ind_sum_0_1=ind_sum_0_1,
               max_ind=max_ind,
               beta_DOY_prior= DOY_effect_prior , 
               spr0_DOY=spring_sub_doy_Z,
               phi_ind=rbind(phi_ind %>% select(`0_1`:`5_3`) %>% as.matrix(),
                             phi_ind_proj %>% as.matrix()),
               psi_ind=c(psi_ind %>% pull(par.index),
                         psi_ind_proj%>% pull(par.index)),
               WH_carc=WH_carc,
               H_carc=H_carc,
               W_carc=W_carc,
               broodstock=broodstock,
               n_ages=3,
               colnames_X_phi=colnames_X_phi,#for reference
               Hmax=c(200,200),       # maximum number of hatchery origin spawners on spawning ground
               Hslope=c(1,1),     # rate at which HO spawners decreases with increasing NO spawners
               stray_home=matrix(c(0.89,0.09,.02,
                                   0.01,.98,0.01),nrow=2,byrow = T), # proportions of fish from each hatchery program going to each trib
               bs_prop=c(.30,.20), #proportion of NO return taken for broodstock
               BS_Max=c(74,64),  #total number of fish needed for broodstock in each stream )
               mu_pss_prior_mean=qlogis(0.85),
               mu_pss_prior_sd=.15,
               logit_p_fem_obs =qlogis(p_F),
               p_fem_obs_sd=0.2,
               f_5_prior=f_5_prior
             ))
  dat_IPM$X_phi<-as(dat_IPM$X_phi,"dgTMatrix")
  dat_IPM$X_p<-as(dat_IPM$X_p,"dgTMatrix")
  dat_IPM$X_psi<-as(dat_IPM$X_psi,"dgTMatrix")
  dat_IPM$proj_years<-0
  dat_IPM$colnames_X_phi<-0
  dat_IPM$log_S_obs_CV_sd<-.25
  dat_IPM$log_S_obs_CV_mu<-log(.1)
  dat_IPM$do_tmbstan<-0
  dat_IPM$log_Jmax_prior_mean<-log(15000)
  dat_IPM$H_RRS_prior<-c(qlogis(.53),1)
  dat_IPM$p_fem_obs_CV<-0.25
  dat_IPM$mu_pss_prior_mean=qlogis(0.85)
  dat_IPM$mu_pss_prior_sd=.5
  dat_IPM$rate=3.5
  here::i_am("data/broodstock_removal.csv")

  #--------------------------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------------------------
  #    Parameters
  
  
  log_S_init<-dat_IPM$log_S_obs %>% apply(2,function(x)head(na.exclude(x),5)) %>% c() %>%  `+` (.25)
  

  logit_pHOS<-((dat_IPM$H_carc+1)/(dat_IPM$WH_carc+1)) %>% qlogis
  
  
  mu_pss <- qlogis(0.9)           # mean prespawn survival on logit scale
  log_pss_sd<-log(.1)       # log SD of random effects 
  eps_pss<-rep(.1,times=dat$n_st)
  eps_p_fem<-eps_pss
  
  par_IPM<-c(params,
             par_fits,
             list(
               log_S_init=log_S_init,
               mu_pss=mu_pss,
               log_pss_sd=log_pss_sd,
               mu_fem=0,
               eps_p_fem=eps_p_fem,
               log_p_fem_sd=log(.1),
               logit_pHOS=logit_pHOS+.5,
               logit_p_brood=rep(qlogis(.2),is.finite(broodstock[-(1:5),] %>% rowSums() %>% log) %>% sum),
               log_S_obs_CV=log(0.1),
               pen_new=log(rep(25,1)),
               beta_DOY= as.numeric(DOY_effect_prior[1])  ,
               b_phi_proj= c(Phi.design.glmmTMB.proj$parameters$b),
               b_psi_proj=c(Psi.design.glmmTMB.proj$parameters$b),
               mu_pHOS=0,
               log_sd_pHOS=0
             )
  )
  par_IPM$log_p_fem_sd<-log(.1)
  par_IPM$log_S_obs_CV<-log(.2)
  par_IPM$beta_phi_ints<-par_IPM$beta_phi_ints+.25
  par_IPM$mu_pHOS<-0
  par_IPM$log_sd_pHOS<-0
  par_IPM$log_f_5_fec<-log(.35)
  par_IPM$logit_RRS<-qlogis(.9)
  par_IPM$mu_pss<-qlogis(.9)
  par_IPM$eta<-par_IPM$eta+2
  par_IPM$log_S_init<-par_IPM$log_S_init+2
  par_IPM$theta_psi_cor<-par_IPM$theta_psi[3]
  par_IPM$theta_psi<-par_IPM$theta_psi[1:2]
  
  
  
  # vector of which paramaters to treat as random effects by integrating them out of marginal likelihoood with Laplace approximation
  rand_par_IPM<-c(random,"beta_phi_pen", "beta_p_pen" , "beta_psi_pen", "b_phi", "b_p" ,"b_psi","log_S_init","beta_DOY","mu_pss","logit_pHOS","eps_p_fem","logit_RRS","log_f_5_fec","theta_psi_cor")
    # vector of parameters to fix at initial values during fitting (projection parameters in this case)
  map_IPM<-c(list(b_phi_proj=rep(factor(NA),length(par_IPM$b_phi_proj)),b_psi_proj=rep(factor(NA),times=length(par_IPM$b_psi_proj))),mod_map)
  
 
  out<-list(dat_IPM=dat_IPM,
              par_IPM=par_IPM,
              rand_par_IPM=rand_par_IPM,
              map_IPM=map_IPM)
  save(out,file=here("results","inputs.Rdata"))

  }
  return(out)
}