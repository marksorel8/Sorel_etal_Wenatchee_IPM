library(here)
library(tidyverse)
library(MARSS)

# read data
dat.z<-read.csv(here("data","dat_z.csv"),row.names = 1) %>% as.matrix()


# fit models

# set new control params
cntl.list <- list( maxit = 2000)


#------------------------------------------------------------
## log flow covariates

B<-c("identity","diagonal and equal","diagonal and unequal")
Q<-c("diagonal and equal","diagonal and unequal")
R<-c("diagonal and equal")
U<-c("zero","unequal")
Z <- matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), nrow = 6)


out_list<-list()
model.data <- data.frame(stringsAsFactors = FALSE)
# fit models & store results
for (b in B) {
  for (q in Q) {
    for(u in U){
      for (r in R){
      model <- list(R = r, B=b, Q=q, Z=Z, U = u)
      kemz <- MARSS::MARSS(dat.z[4:9,-(1:92)],
                           model = model, control = cntl.list)
      model.data <- rbind(
        model.data,
        data.frame(
          B = b,
          Q = q,
          U= u,
          R=r,
          logLik = kemz$logLik,
          K = kemz$num.params,
          AICc = kemz$AICc,
          stringsAsFactors = FALSE
        )
      )
      assign(paste("kemz", b,q,u, sep = "."), kemz)
      out_list[[paste("kemz", b,q,u, sep = ".")]]<-kemz
    } # end b loop
  } # end q loop
} # end of U loop
}#

model.data %>% arrange(AICc)
out_list$`kemz.diagonal and equal.diagonal and equal.zero`
#------------------------------------------------------------

## air temp 

B<-c("identity","diagonal and equal")
U<-c("zero","unequal")

out_list_temp<-list()
model.data_temp <- data.frame(stringsAsFactors = FALSE)
# fit models & store results
for (b in B) {
  for (u in U) {
    if(b!="identity"&u!="zero") next
    model <- list( B=b, U = u)
    kemz <- MARSS::MARSS(dat.z[1,-(1:59)],
                         model = model, control = cntl.list)
    model.data_temp <- rbind(
      model.data_temp,
      data.frame(
        B = b,
        U = u,
        logLik = kemz$logLik,
        K = kemz$num.params,
        AICc = kemz$AICc,
        stringsAsFactors = FALSE
      )
    )
    assign(paste("kemz", b,u, sep = "."), kemz)
    out_list_temp[[paste("kemz", b,u, sep = ".")]]<-kemz
  } # end b loop
} # end U loop

model.data_temp %>% arrange(AICc)

out_list_temp$`kemz.diagonal and equal.zero`
#------------------------------------------------------------
## upwelling  

B<-c("identity","diagonal and equal")
U<-c("zero","unequal")

out_list_upwelling <-list()
model.data_upwelling <- data.frame(stringsAsFactors = FALSE)
# fit models & store results
for (b in B) {
  for (u in U) {
    if(b!="identity"&u!="zero") next
    
    model <- list( B=b, U = u)
    kemz <- MARSS::MARSS(dat.z[3,-(1:47)],
                         model = model, control = cntl.list)
    model.data_upwelling <- rbind(
      model.data_upwelling,
      data.frame(
        B = b,
        U = u,
        logLik = kemz$logLik,
        K = kemz$num.params,
        AICc = kemz$AICc,
        stringsAsFactors = FALSE
      )
    )
    assign(paste("kemz", b,u, sep = "."), kemz)
    out_list_upwelling[[paste("kemz", b,u, sep = ".")]]<-kemz
  } # end b loop
} # end U loop


model.data_upwelling %>% arrange(AICc)

out_list_upwelling$kemz.identity.unequal
acf(c(out_list_upwelling$kemz.identity.unequal$states-
        out_list_upwelling$kemz.identity.unequal$ytT))


#------------------------------------------------------------
## Sea surface temperature 

B<-c("identity","diagonal and equal")
U<-c("zero","unequal")

out_list_SST <-list()
model.data_SST <- data.frame(stringsAsFactors = FALSE)
# fit models & store results
for (b in B) {
  for (u in U) {
    if(b!="identity"&u!="zero") next
    model <- list( B=b, U = u)
    kemz <- MARSS::MARSS(dat.z[2,],
                         model = model, control = cntl.list)
    model.data_SST <- rbind(
      model.data_SST,
      data.frame(
        B = b,
        U = u,
        logLik = kemz$logLik,
        K = kemz$num.params,
        AICc = kemz$AICc,
        stringsAsFactors = FALSE
      )
    )
    assign(paste("kemz", b,u, sep = "."), kemz)
    out_list_SST[[paste("kemz", b,u, sep = ".")]]<-kemz
  } # end b loop
} # end U loop


model.data_SST %>% arrange(AICc)
out_list_SST$kemz.identity.unequal

acf(c(out_list_SST$kemz.identity.unequal$states-
        out_list_SST$kemz.identity.unequal$ytT))



## obs model: y_t = alpha + beta * t + e_t
Z <- matrix(1)
A <- matrix("alpha")
D <- matrix("beta")
d <- matrix(seq(ncol(dat.z)), nrow = 1)
R <- matrix(0)

## autocorrelated errors: e_t = phi * e_{t-1} + w_t
B <- matrix("phi")
U <- matrix(0)
Q <- matrix("q")

model <- list(Z=Z,
              A=A,
              D=D,
              d=d,
              R=R,
              B=B,
              U=U,
              Q=Q
)
kemz <- MARSS::MARSS(dat.z[2,],
                     model = model, control = cntl.list)

