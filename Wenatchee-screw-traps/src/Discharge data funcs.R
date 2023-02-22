#function to load and download packages in necessary
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}#end of function 

#function to load Chiwawa discharge data
Chiw_discharge_func<-function(){

  
  pkgTest("dataRetrieval")
  require("dataRetrieval")
  pkgTest("here")
  require("here")
  Chiw.flowDV<- readNWISdv(siteNumber="12456500",
                           "00060", "1900-01-01", "2018-12-31",statCd="00003")
  #munging
  colnames(Chiw.flowDV)[4]<-"flow"
  Chiw.flowDV$date<-as.Date(Chiw.flowDV$Date,format="%Y-%m-%d")
  Chiw.flowDV$Year<-as.integer(format(Chiw.flowDV$date,"%Y"))
  Chiw.flowDV$week<-as.integer(format(Chiw.flowDV$date,"%W"))
  Chiw.flowDV$month<-as.integer(format(Chiw.flowDV$date,"%m"))
  
  Chiw.flowDV$day<-as.integer(format(Chiw.flowDV$date,"%d"))
  Chiw.flowDV$doy<-as.integer(format(Chiw.flowDV$date,"%j"))
  
  Chiw.flowDV<-Chiw.flowDV %>% filter(Year%in% (((pull(.,Year) %>% table())[(pull(.,Year) %>% table())>363]) %>% names()))
  
   return(Chiw.flowDV)
}


#####################################################


Nason_White_Discharge_Func<-function(){
  pkgTest("here")
  ##Nason Flow
  #downloaded from https://fortress.wa.gov/ecy/eap/flows/station.asp?sta=45J070
  
  #function to load discharge data: "loc" is the folder location with the files (one for each yar), and "breakYr" is the file where the format of the data changes (related to the number of rows of header). 
  loadDis<-function(loc){
    if(loc=="Nason"){
      site_num<-"45J070/45J070_"
      # "https://fortress.wa.gov/ecy/eap/flows/stafiles/"
      start_year<-2002
      break_year<-2004
    }else{
      site_num<-"45K090/45K090_"
      break_year<-2013
      start_year<-2003
    }
      
    out<-data.frame()
    for ( i in start_year:2019){
      if(i<=break_year) {test<-read_table2(paste0("https://apps.ecology.wa.gov/ContinuousFlowAndWQ/StationData/Prod/",site_num,i,"_DSG_DV.txt"),skip=ifelse(i==2002,250,4),n_max=366,col_names=c("date","time","flow","data.code")) 
      test<-test[,-2]}else{
        
        test<-read_table2(paste0("https://apps.ecology.wa.gov/ContinuousFlowAndWQ/StationData/Prod/",site_num,i,"_DSG_DV.txt"),skip=12,n_max=366,col_names=c("date","flow","data.code"))}
      
      out<-rbind(out,test)
    }
  
    #remove duplicate rows
    out<-out[!duplicated(out),]
    
    out$flow<-as.numeric(out$flow)
    out$date<-as.Date(out$date,format="%m/%d/%Y")
    out$Year<-as.integer(format(out$date,"%Y"))
    out$month<-as.integer(format(out$date,"%m"))
    out$doy<-as.integer(format(out$date,"%j"))
    
    return(out)
  }
  
  
  Nas<-loadDis(loc="Nason")
  
  
  ## White Flow
  #downloaded from https://fortress.wa.gov/ecy/eap/flows/station.asp?sta=45K090#block2
  
  White<-loadDis(loc="White")
  
  return(list(Nason_Dis=Nas,White_Dis=White))
  
}
