# 4 helper functions for fitting spawner-to-juveniile models


##-------------------------------------------------------------------------------------------------------------------
##function to create data input for models. returns a list of data inputs to feed to TMB model "Stock_recruit_LVM". 
make_dat_func<-function(streams=c( 0:2), #streams to model 0 = chiwawa, 1 = nason, 2 = white
                        LHs=c(1:4),
                        n_f=1,        #life histories to model 1=fry, 2=sumemr, 3=fall, 4=smolt
                        no_rand_FF = 0,rate,fold_exclude){  # use hierarchical structure for functional form parameters if > 1 stream
  
  
  # concatenate log mean of observed emigrant distributions,
  # log standard devations of observed emigrant distributions,
  # and observed spawners for all streams, life histories, and years
  log_m_means<-log_m_sds<-numeric(0) #empty vectors
  sl_i<-st_i<-t_i<-s_i<-l_i<-BY<-redd_yr<-integer(0)
  S_obs<-matrix(nrow=0,ncol=4)
  X<-matrix(nrow = 0,ncol=24)         #design matrix
  
  
  #construction of vectors
  for (i in streams){ #loop through streams
    for(j in LHs){      #loop though life histories
      
      years<-sort(unique((all_emigrants_estimates[[i*2+1]]$dat$years)))#years (of observations)
      BY<-c(BY,(years-ifelse(j==4,2,1)))#brood years corresponding with observation years for given life history
      n_year<-length(years)#number of years
      
      # likelihood penalty mean and standard deviations for juvenile emigrant abundances
      if(j<4){
        means<-(all_emigrants_estimates[[i*2+1]]$boot$boot_log_means)[,j] #subyearlings
        sds<-(all_emigrants_estimates[[i*2+1]]$boot$boot_log_sds)[,j] #subyearlings
        
      }else{
        means<-all_emigrants_estimates[[i*2+2]]$boot$boot_log_means[,1] # yearlings
        sds<- all_emigrants_estimates[[i*2+2]]$boot$boot_log_sds[,1] # yearlings
      }
      # log means- concatenate current vector with all life histories
      log_m_means<-c(log_m_means,                                 #current vector
                     means)   
      # log standard deviations- concatenate current vector with all life histories
      log_m_sds<-c(log_m_sds,                                   #current vector
                   sds)   
      
      # observed redds
      redds_sub<-filter(redds,stream== c("Chiwawa","Nason","White")[(i+1)])
      redd_yr<-c(redd_yr,redds_sub$redds[match((years-ifelse(j==4,2,1)),redds_sub$Year)]/trib_lengths[(i+1)])
      
      # design matrix of environmental covariates
      X_s<-matrix(NA,nrow=n_year,ncol=24)
      #$ winter discharge brood year (incubation)
      X_s[(1:n_year),j]<-#(j+1)] <-##j+i*4] <-
        scale(flow_covs$winter_high[match((years-ifelse(j==4,2,1)), #Year of incubation
                                          flow_covs$winter_high$Year),(i+2)])
      
      #summer discharge brood year +1 (summer rearing)
      if(j>1){
        X_s[(1:n_year),(j+3)] <-##(j+3)]<-# (j+i*3)+11] <-
          scale(flow_covs$summer_low[match((years-ifelse(j==4,1,0)), #Year of incubation
                                           flow_covs$summer_low$Year),(i+2)])
      }
      
      #winter discharge brood year +1 (overwintering)
      if(j==4){
        X_s[(1:n_year),8] <- ##(i)+22] <-
          scale(flow_covs$winter_high[match((years-1), #Year of overwintering as parr
                                            flow_covs$winter_high$Year),(i+2)])
      }
      # add rows to design matrix
      X<-rbind(X,X_s)
      
      
      # stream by life history index 
      sl_i<- c(sl_i,
               rep(i*4+j,n_year))
      
      #stream by year index
      st_i<-c(st_i,
              (ifelse(j==4,1,2):(n_year+ifelse(j==4,0,1)))+(i*100))
      
      #brood year index
      t_i<-c(t_i,years-ifelse(j==4,2,1))
      
      #life history index (for reference/plotting only)
      l_i<- c(l_i,
              rep(j-1,n_year))
      
      #stream index (for reference/plotting only)
      s_i <- c(s_i,
               rep(i,n_year))
      
    }
  }
  
  
  
  
  #get rid of columns not used in envornmental cpovariate design matrix
  X<-X[,which(apply(X,2,function(x)sum(!is.na(x)))>0),drop=FALSE]
  #scale
  # X<-scale(X)
  #fill in NAs with zeros
  X[is.na(X)]<-0
  #check scaling
  apply(X,2,sd,na.rm=T)
  apply(X,2,mean,na.rm=T)
  
  # Cov_dat<-tibble(BY=BY,LH=factor(l_i,labels=c("fry","summer","fall","smolts")),stream=factor(s_i,labels=c("Chiwawa","Nason","White"))) %>% mutate(Y0=BY+1) %>% left_join(pivot_longer(flow_covs$winter_high,2:4,names_to = "stream", values_to="win_0") %>% rename("BY"="Year")) %>%
  #   left_join(pivot_longer(flow_covs$summer_low ,2:4,names_to = "stream", values_to="sum_0") %>% rename("Y0"="Year")) %>%
  #   left_join(pivot_longer(flow_covs$winter_high,2:4,names_to = "stream", values_to="win_1") %>% rename("Y0"="Year"))    %>%
  #   left_join(pHOS %>% rename("BY"="Brood_year","stream"="Stream") %>% select(-(n_carcasses:n_hatchery))) %>%
  #   group_by(stream,LH) %>% summarize(win_0=c(scale(win_0)),sum_0=c(scale(sum_0)),win_1=c(scale(win_1)),pHOS=c(scale(pHOS_weighted)))
  
  
  # X2<-model.matrix(BY~-1+win_0:LH:stream+sum_0:LH:stream+win_1:stream:LH+stream:LH:pHOS,data=Cov_dat) %>% as_tibble %>% select(-c(`LHfry:streamChiwawa:sum_0`,`LHfry:streamNason:sum_0`,`LHfry:streamWhite:sum_0`,`LHfry:streamChiwawa:win_1`:`LHfall:streamChiwawa:win_1`,`LHfry:streamNason:win_1`:`LHfall:streamNason:win_1`,`LHfry:streamWhite:win_1`:`LHfall:streamWhite:win_1`)) %>% as.matrix()
  
  # X2<-model.matrix(BY~-1+win_0:LH+sum_0:LH+win_1:LH,data=Cov_dat) %>% as_tibble() %>% select(-c(`LHfry:sum_0`,`LHfry:win_1`:`LHfall:win_1`)) %>% as.matrix()
  # +pHOS:LH:stream
  
  
  
  
  #construct data list
  dat<-list(n_sl = length(unique(sl_i)),           #n umber of distinct life history by stream combinations
            n_i = length(sl_i),                    # number of years x streams x life histories
            n_t = length(unique(t_i)),             # number of years
            n_l = length(LHs),             # number of life histories
            n_s = length(streams),             # number of life histories
            n_f = n_f,                             # number of latent variable factors
            mod= rep(1,length(unique(sl_i))),      # functional forms of process model
            J_obs = log_m_means,                   # Observed log recruits (posterior mean from juvenile model)
            J_obs_sd = log_m_sds,                  # Observed log recruits standard error (posterior sd from juvenile model)
            log_S_obs=log((redd_yr[!duplicated(st_i)] 
                           [order(st_i[!duplicated(st_i)])])), # log observed redds
            S_obs_CV = 0.1,                       # observed redd CV
            sl_i = factor(sl_i),                   # stream by life history index for observations
            st_i = factor(st_i),                   # stream by year index for observations
            t_i = factor(t_i),                     # year index of each observation
            X = as(X,"dgTMatrix"),                 # sparse design matrix of environmental covariates
            no_rand_FF = no_rand_FF,               # flag for whether to treat eps_alpha, eps_gamm, and eps_Jmax as random effects or fixed effects
            s_i = s_i,                             # stream index for reference (not used in model)
            l_i = l_i,                             # LHP index for reference (not used in model)
            BY=BY,                                 # Brood year for reference (not used in model)
            # beta_e_ind=as.numeric(as.factor(paste(c(rep(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."),times=3),
            #                                         rep(c("Sum. Sub","Fall Sub", "Spr. Yrl."),times=3),
            #                                         rep(c( "Spr. Yrl."),times=3),
            #                                         rep(c("Spr. Sub","Sum. Sub","Fall Sub", "Spr. Yrl."),times=3)),c(substr(colnames(X)[1:12],1,5),
            #                                                                                                          substr(colnames(X)[13:24],nchar(colnames(X)[13:24])-4,nchar(colnames(X)[13:24])))) %>% as.factor() %>% as.numeric()))-1
            rate=rate,
            rate_prior=c(0,25)
            # fold_exclude=fold_exclude
  )
  return(dat)
}   

##-------------------------------------------------------------------------------------------------------------------

#function to construct parameters list. returns a list of initial values for paramaters to feed to TMB model "Stock_recruit_DFA".
make_params_func<-function(dat){
  
  params<-list( beta_alpha=ifelse(rep(dat$no_rand_FF,dat$n_l),
                                  rep(0,dat$n_l),
                                  rnorm(dat$n_l,log(5),1)),   # log alpha interceptrs (by life-history)
                beta_Jmax= ifelse(rep(dat$no_rand_FF,dat$n_l),
                                  rep(0,dat$n_l),
                                  rnorm(dat$n_l,log(15000),1)), # log asymptotic maximum recruitment (Jmax) coefficients (by life-history) intercept
                beta_gamma= ifelse(rep(dat$no_rand_FF,dat$n_l),
                                   rep(0,dat$n_l),
                                   rnorm(dat$n_l,0,.2)),                    # log gamma coefficients (by life history)
                
                log_sigma_alpha=rep(0,dat$n_l),                    # log alpha random effect SD
                log_sigma_Jmax=rep(0,dat$n_l),                   # log Jmax random effect SD
                log_sigma_gamma=rep(0,dat$n_l),                    # log gamma random effect SD
                log_sigma_BJmax=rep(0,dat$n_l), 
                # log_loadings_sd=rep(0,dat$n_f*(dat$n_sl)-dat$n_f*(dat$n_f-1)/2),
                beta_e=rep(0,ncol(dat$X)),                     # environmental covariate coefficients for process error
                Loadings_vec=rnorm(dat$n_f*(dat$n_sl)-dat$n_f*(dat$n_f-1)/2,0,.1), # latent variable factor loadings
                log_sigma_eta=log(abs(rnorm(dat$n_sl,.5,.1))),# idiosyncratic process error log SDs
                
                #random effects
                eps_alpha=ifelse(rep(dat$no_rand_FF,dat$n_sl),
                                 rnorm(dat$n_sl,log(5),.1),
                                 rep(log(5),dat$n_sl)),                         # log alpha random effects
                eps_Jmax=ifelse(rep(dat$no_rand_FF,dat$n_sl),
                                rnorm(dat$n_sl,log(15000),.1),
                                rep(log(15000),dat$n_sl)),                        # log Jmax random effects
                eps_gamma= ifelse(rep(dat$no_rand_FF,dat$n_sl),
                                  rnorm(dat$n_sl,0,.1),
                                  rep(0,dat$n_sl)),                         # log gamma random effects
                log_S_hat = rnorm(length(dat$log_S_obs), dat$log_S_obs,.1), # latent spawners random effects
                Omega_xf=matrix(0,dat$n_t,dat$n_f),            #  latenct variable factors
                eta=rnorm(length(dat$J_obs),dat$J_obs,.1),#-mean(dat$J_obs),                         #  idisyncratic process error random effects
                # log_sd_beta_e=rep(log(1),times=ncol(dat$X)),
                # mu_beta_e=rep(0,times=ncol(dat$X)/3),
                # log_sd_beta_e=rep(0,times=ncol(dat$X)),
                log_sigma_Bgamma=rep(0,dat$n_l),
                log_sigma_beta_e=rep(-1,ncol(dat$X)),
                
                rate=dat$rate,
                log_sigma_loading=0
  )                       
  
  return(params)
}

##-------------------------------------------------------------------------------------------------------------------


#function to make map (controls which parameters are fit and which are held fixed and initials). Depending on functional form, some parameters need to be fixed (e.g., there are no gammas or Jmaxes in a straight line through the origin). Returns a list which can be included in the "map=" call of the `MakeADFun` function that initializes a TMB model "Stock_recruit_DFA" to fix paramaters. 
make_map<-function(mod_dat,fit_env=TRUE,rate_map=factor(NA)){
  
  # vectors of factors indicating all variables should be independently fit. 
  vec1<-factor(seq(1,mod_dat$n_l))    # beta_Jmax map
  vec1.1<-factor(seq(1,mod_dat$n_sl)) # eps_rm map
  vec2<-factor(seq(1,mod_dat$n_l))    # beta_gamma map
  vec2.1<-factor(seq(1,mod_dat$n_sl)) # eps_f map
  
  # change some levels to NA to fix those variables (depending on functional forms)
  ##  hyper means and random effect standard deviation
  for ( i in 1:mod_dat$n_l){
    if(mod_dat$mod[i]%in%c(3:4)){ vec1[i]<-factor(NA) # beta_Jmax map
    
    }
    if(mod_dat$mod[i]%in%c(4)){ vec2[i]<-factor(NA)  #beta_gamma map
    
    }
  } 
  ##random effects
  for ( i in 1:mod_dat$n_sl){
    if(mod_dat$mod[i]%in%c(3:4)){ vec1.1[i]<-factor(NA) #eps_rm map
    
    }
    if(mod_dat$mod[i]%in%c(2,4)){ vec2.1[i]<-factor(NA)  #eps_f map
    
    }
  }  
  
  #creat map list to return
  map=list(
    beta_alpha=factor(seq(1,mod_dat$n_l)),
    log_sigma_alpha=factor(seq(1,mod_dat$n_l)),
    
    eps_gamma=droplevels(vec2.1),
    beta_gamma=droplevels(vec2),
    log_sigma_gamma=droplevels(vec2),
    
    eps_Jmax=droplevels(vec1.1),
    beta_Jmax=droplevels(vec1),
    log_sigma_Jmax=droplevels(vec1),
    rate=rate_map
  )
  
  
  # if not sharing information about functional forms fix hyper parameters
  if(mod_dat$no_rand_FF){
    map$log_sigma_alpha<-map$log_sigma_Jmax<-map$log_sigma_gamma<-
      map$beta_alpha<-map$beta_gamma<-map$beta_Jmax<-factor(rep(NA,mod_dat$n_l))
  }
  
  #drop unused levels after replacing some values with NAs above
  map<-lapply(map,droplevels)
  
  #fix environmental covaraites at zero if fit_env= FALSE
  if(!fit_env){
    map$beta_e<-factor(rep(NA,ncol(mod_dat$X))) # fix environmental covariates at 0.0
  }
  # 
  # if(mod_dat$n_f==0){
  #   map$log_sigma_loading<-factor(NA) # fix environmental covariates at 0.0
  # }
  
  
  if(mod_dat$n_f>mod_dat$n_sl){stop("number of factors is > number of stream x LH combos")}
  return(map)
}

##-------------------------------------------------------------------------------------------------------------------

## function to attempt model fitting a specified number of times. Returns list of bets fit model (based on BIC), model and "fit", objects, and a vector of the BICs for all fitting attempts. Plus a "report" object and input data.
fit_mod_func<-function(mods, streams, LHs, n_f, no_rand_FF, fit_env, fit_attempts,rate,fold_exclude,rate_map){
  dat<-make_dat_func(streams,LHs,n_f,no_rand_FF,rate,fold_exclude) # make model data
  dat$mod<-mods                                  # specify functional forms
  mod_map<-make_map(mod_dat=dat,fit_env = fit_env,rate_map)                 # create map
  fit<-mod<-report<-NA # placeholders for converged model and "fit" objects
  BIC_vec<-rep(Inf,fit_attempts)
  for ( i in 1:fit_attempts){
    params<-make_params_func(dat)        # make initial parameter values
    if(no_rand_FF){
      mod_i<-TMB::MakeADFun(dat,params,random=c("log_S_hat","eta","Omega_xf"),DLL="Stock_recruit_LVM",map=mod_map,silent = TRUE)}else{
        mod_i<-TMB::MakeADFun(dat,params,random=c("log_S_hat","Omega_xf","eps_alpha","eps_gamma","eps_Jmax","beta_gamma","beta_Jmax","eta","rate","Loadings_vec","beta_e"),DLL="Stock_recruit_LVM",map=mod_map,silent = TRUE)#,,",,, "log_sigma"
      }
    fit_i<-NA # clear previous fit object 
    try(fit_i<-TMBhelper::fit_tmb(mod_i, newtonsteps = 1,getJointPrecision = TRUE)) # optimize
    BIC_mod<- NA #clear previous BIC
    try(BIC_mod<-(2*fit_i$objective+log(dat$n_t)*fit_i$number_of_coefficients[2])) # calculate  BIC
    try(if(min(BIC_vec,na.rm = TRUE)>BIC_mod){ # save model and "fit" object if BIC is lower than previous
      fit<-fit_i
      mod<-mod_i
      try(report<-mod$report())
    })
    try(BIC_vec[i]<-BIC_mod) # save BIC
  }
  if(is.na(mod))mod<-mod_i
  return(list(mod=mod,fit=fit,BIC_vec=BIC_vec,report=report,dat=dat))
}

#function to try a few more iterations if didn't converge in initial tries
fit_mod_iter<-function(mods, streams, LHs, n_f, no_rand_FF, fit_env, fit_attempts, additional_attempts,rate=c(.5,10,.1),fold_exclude,rate_map){
  x<-fit_mod_func(mods, streams, LHs, n_f, no_rand_FF, fit_env, fit_attempts,rate,fold_exclude,rate_map)
  i<-0
  while(i<=additional_attempts & (min(x$BIC_vec,na.rm=TRUE)==Inf)){
    x<-fit_mod_func(mods, streams, LHs, n_f, no_rand_FF, fit_env, fit_attempts=1,rate,fold_exclude)
    i<-i+1
  }
  
  return(x)
}



#function for simulation from a multivariate normal, taken fron J. Thorson's Fish Utils. Package
rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
  set.seed( random_seed )
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L = Matrix::Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}

