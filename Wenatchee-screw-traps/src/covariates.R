#This script compiles covariate values for a model of the transition from spawners to juvenile emigrants in the Chiwawa, Nason, and White Rivers. 

# covariates to be included

# *stream specific*
# - Weighted habitat area

# *year varying*
# - Max of daily flows Oct-March (Winter)
# - Min of daily flows June-Sept (Summer)

library(readr)
library(readxl)
library(here)


#load functions for downloading discharge data from internet
source(here("src","Discharge data funcs.R"))

#dowload Chiwawa data
if(file.exists(here("chiwDis.csv"))){
  load(here("chiwDis.csv"))}else{
    chiwDis<-Chiw_discharge_func()
    save(chiwDis,file=here("chiwDis.csv"))  
  }

#download Nason and White data
if(file.exists(here("nas_whi_dis.Rdata"))){
  load(here("nas_whi_dis.Rdata"))}else{
    nas_whi_dis<-Nason_White_Discharge_Func()
    save(nas_whi_dis,file=here("nas_whi_dis.Rdata"))  
  }





#functions to calculat winter maximum flow and summer minimum flow, code copied from M. Sheuerell. https://github.com/mdscheuerell/skagit_sthd/blob/master/analysis/App_1_Retrieve_covariates.pdf

#winter aximum (in year when winter begins)
winter_max_func<-function(dat_flow){
  ## autumn flows in year t
  flow_aut <- subset(dat_flow, (month>=10 & month<=12))
  ## spring flows in year t+1
  flow_spr <- subset(dat_flow,
                     (month>=1 & month<=2))
  ## change spr year index to match aut
  flow_spr[,"Year"] <- as.integer(flow_spr[,"Year"]) - 1
  ## combine flows indexed to winter start year & calculate max flow
  dat_flow_wtr <- aggregate(flow ~ Year, data = rbind(flow_aut,flow_spr),function(x)log (max(x)))
  dat_flow_wtr[,"flow"] <- round(dat_flow_wtr[,"flow"], 1) 
 
 return(dat_flow_wtr) 
}

#summer mean
summer_min_func<-function(dat_flow){
  ## autumn flows in year t
  flow_sum <- subset(dat_flow, (month>=6 & month<=9))
  ## calculate min flow
  dat_flow_sum <- aggregate(flow ~ Year, data = flow_sum, mean)

  return(dat_flow_sum) 
}

#compile flow data for all streams and both seasons
compile_flow_func<-function(){
##start with summer low flow
nad<-summer_min_func(as.data.frame(nas_whi_dis$Nason_Dis))
whd<-summer_min_func(as.data.frame(nas_whi_dis$White_Dis))
chd<-summer_min_func(as.data.frame(chiwDis))
summer_min_flow<-merge(chd,nad,by="Year",all=TRUE) %>% merge(whd,by="Year",all=TRUE)
colnames(summer_min_flow)[2:4]<-c("Chiwawa","Nason","White")

##winter max flow
nad<-winter_max_func(as.data.frame(nas_whi_dis$Nason_Dis))
whd<-winter_max_func(as.data.frame(nas_whi_dis$White_Dis))
chd<-winter_max_func(as.data.frame(chiwDis))
winter_max_flow<-merge(chd,nad,by="Year",all=TRUE) %>% merge(whd,by="Year",all=TRUE)
colnames(winter_max_flow)[2:4]<-c("Chiwawa","Nason","White")
return(list(summer_low=summer_min_flow,winter_high=winter_max_flow))
}

flow_covs<-compile_flow_func()


#----------------------------------------------------------------
# air temp
#https://w2.weather.gov/climate/xmacis.php?wfo=otx




#snowpack
