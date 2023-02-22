
#function to make a TS object out of catches
make_ts<-function(myDates,myCounts,main,plot=FALSE,plotDis=FALSE,dis_dates=NULL,dis_vals=NULL,range=NULL){
  
  ts_func<-function(my_dates,my_x){
    
    drop_leap_day<-function(x){
    x[!(format(x,"%m") == "02" & format(x, "%d") == "29")]
  }
  
  out<-ts(
    my_x[match(drop_leap_day(seq.Date(min(my_dates,na.rm=T),
                                          max(my_dates,na.rm=T),by=1)),my_dates)],start=as.numeric(c(format(min(my_dates,na.rm=T), "%Y"),format(min(my_dates,na.rm=T), "%j"))),frequency=365)

  return(out)
  }
  
 out2<-ts_func(myDates,myCounts)
 dis_ts<-NULL
 if(plotDis==TRUE){
 dis_ts<-ts_func(dis_dates,dis_vals)
 }
   
    if(isTRUE(plot)){
    plot_catch(out2,main,plotDis,dis_ts)
  }
  list(catch=out2,dis=dis_ts)

  
}



#function to plot a TS of catches
plot_catch<-function(catch,main,plot_dis=FALSE,dis_series=NULL,range=range){
  
  year<-start(catch)[1]
  year_ln<-end(catch)[1]-year
  par(mfrow=c(3,1),mar=c(3,5,1,2),oma=c(0,0,3,0))
  
  for ( i in 1:3){
    
    plot(catch,
         xlim=c(year,(year+floor(year_ln/3)+1)),ylab="",xlab="") 
    
    
    if (plot_dis==TRUE){
      par(new=T)
    plot(dis_series,
         xlim=c(year,(year+floor(year_ln/3)+1)),ylab="",xlab="",axes=F,col=rgb(.2,.1,.1,.4),lwd=.7) 
    }
    
    year<-year+floor(year_ln/3)+1

  }
  
  
  
  
  mtext("catch",2,-1.5,outer=T)
  mtext(main,3,0,outer=T,xpd=NA)
  
}


#function to plot efficiency data
plot_effic<-function(effic, dates, num_rel, main, cols, plot_Dis=FALSE, dis=NULL, dis_dat=NULL){
  dat_length<-diff(range(dates,na.rm = TRUE))
  year_ln<-diff(range(as.numeric(format(dates,form="%Y")),na.rm = TRUE))
  year<-min(as.numeric(format(dates,form="%Y")),na.rm = TRUE)
  par(mfrow=c(3,1),mar=c(3,5,1,2),oma=c(0,0,3,0))
  for ( i in 1:3){
    
    plot(dates,
         effic,cex=num_rel/100+.5,pch=19,col=cols,
         xlim=c(as.Date(paste0(year,"-01-01")),c(as.Date(paste0(year+floor(year_ln/3),"-12-30")))),ylab="") 
    
    if (plot_Dis==TRUE){
      par(new=T)
      plot(dis_dat,
           dis,type="l",
           xlim=c(as.Date(paste0(year,"-01-01")),c(as.Date(paste0(year+floor(year_ln/3),"-12-30")))),ylab="",xlab="",axes=F,col=rgb(.1,.1,.1,.5)) 
      
    }
    
    #points(supsmu(dates,effic,span=0.01),type="l",col="red")
    year<-year+floor(year_ln/3)+1
  }
  
  mtext("capture efficiency",2,-1.5,outer=T)
  mtext(main,3,0,outer=T)
  
}
