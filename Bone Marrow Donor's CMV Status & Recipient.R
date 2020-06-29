#STAT533 project R code
#Based on the bmt dataset in "KMsurv" package
#Hypothesis : Does bone marrow donor's cytomegalovirus(CMV) immune status (positive or nagative) affect the time to graft-verses-host disease (GVHD) of the recipient of the bone marrow from the donor?
#In the data set, there are two types of GVHD, which are:
#ta = acute GVHD
#tc = chronic GVHD
#These have different symptoms as well as general time in which they occur in the recipient's body.

##########################################################
setwd("/Users/haigelee/Hyejung/STAT533/")
install.packages("OIsurv")
library(KMsurv)
data(bmt)
attach(bmt)

######################################################################################

#Create non-parametric CGVHD data frame

##############creating NA estimators###########################################################################
NAestimator <- function (Z,delta) 
{ UZ<-unique(Z)      #Unique observed times;
N<-length(UZ)
UZ.order<-order(UZ)
UZ<-UZ[UZ.order]            #sort data;
N.A <- rep(0,N)             #cannot use NA (reserved by R) for N.A;
Y<-rep(0,N)
D<-rep(0,N)
D[1]<-sum(Z[delta==1]==UZ[1])
Y[1]<-sum(Z >= UZ[1])
N.A[1] <- D[1]/Y[1]            #this is for right continuous value
for (i in 2: N){
  D[i]<-sum(Z[delta==1]==UZ[i])
  Y[i]<-sum(Z >= UZ[i])
  N.A[i] <- N.A[i-1]+D[i]/Y[i]
}

# Calculate variance and sdandard error of N-A estimator;

sigma2.h<-rep(0,N)
for (i in 1: N){
  sigma2.h[i]<-sum( (UZ<=UZ[i])*(D/Y^2) )     
}
NA.var<-sigma2.h
NA.se<-sqrt(NA.var)

Full<-data.frame(UZ,D,Y,N.A,NA.var,NA.se)
Reduced<-subset(Full, (Full$D>0))
list(Full, Reduced)
}


options(width=60,length=200,digits=4)   # Set number of digits to be printed

NA.est1<-NAestimator(tc, dc)[[2]]
NA.est1
#Note that NA estimator is defined at time t=0. So we manually add one row to the NA.est.low for t=0.
test<- data.frame("UZ" = 0, "D" = 0, "Y" = 274, "N.A" = 0, "NA.var" = 0, "NA.se" = 0)
NA.est1<- rbind(test, NA.est1)
rownames(NA.est1) <- NULL
NA.est1
write.csv(NA.est1, "NA.CGVHD.csv", row.names = FALSE)
##################End. creating NA estimators#########################################
##################Creating KM estimators##############################################
KMestimator <- function (Z,delta) 
{ UZ<-unique(Z)      #Unique observed times;
N<-length(UZ)
UZ.order<-order(UZ)
UZ<-UZ[UZ.order]            #sort data;
KM <- rep(0,N)
Y<-rep(0,N)
D<-rep(0,N)
D[1]<-sum(Z[delta==1]==UZ[1])
Y[1]<-sum(Z >= UZ[1])
KM[1] <- 1-D[1]/Y[1]            #this is for right continuous value
for (i in 2: N){
  D[i]<-sum(Z[delta==1]==UZ[i])
  Y[i]<-sum(Z >= UZ[i])
  KM[i] <- KM[i-1]*(1-D[i]/Y[i])
}

# Calculate variance and sdandard error;
sigma2.s<-rep(0,N)
for (i in 1: N){
  ## sigma2.s[i]<-sum( (UZ<=UZ[i])*(D/(Y*(Y-D))) ) #old version
  sigma2.s[i]<-sum( (UZ[1:i]<=UZ[i])*(D[1:i]/(Y[1:i]*(Y[1:i]-D[1:i]))))   
  ##  Note the data is sorted by UZ;
  ##  Using this to avoid NaN for times smaller than the largest observation;
}
KM.var<-KM^2*sigma2.s
KM.se<-sqrt(KM.var)

Full<-data.frame(UZ,D,Y,KM,sigma2.s,KM.var,KM.se)
Reduced<-subset(Full, (Full$D>0))
list(Full, Reduced)
}

KM.all<-KMestimator(tc, dc)
write.csv(KM.all[[1]], "KM.CGVHD.not.censored.csv", row.names = FALSE) 
write.csv(KM.all[[2]], "KM.CGVHD.censored.csv", row.names = FALSE) 
###############################End. Creating KM estimators#################################


######################Creating H.hat data frame#####################################################
H.hat<- -log(KM.all$KM)
sigma2.h<-rep(0,length(H.hat))
for (i in 1: length(H.hat)){
  sigma2.h[i]<-sum( (KM.all$UZ<=KM.all$UZ[i])*(KM.all$D/KM.all$Y^2) )     
}
NA.var<-sigma2.h
NA.se<-sqrt(NA.var)
H.hat.all<-data.frame("UZ"=KM.all$UZ, "D"=KM.all$D, "Y"=KM.all$Y, H.hat, NA.var,NA.se)
#Note that Cum.hazard is defined at time t=0. So we manually add one row to the NA.est.low for t=0.
test<- data.frame("UZ" = 0, "D" = 0, "Y" = 274, "H.hat" = 0, "NA.var" = 0, "NA.se" = 0)
H.hat.all<- rbind(test, H.hat.all)
rownames(H.hat.all) <- NULL
write.csv(H.hat.all, "H.hat.CGVHD.csv", row.names = FALSE)  
######################End. Creating H.hat data frame##################################################
