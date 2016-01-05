## R version: 3.1.2
library(Bchron) #for calibration of 14C dates // version 4.1.1
library(zoo) #for plotting SPD outputs //version 1.7-12

setwd("~/github/jomonSPD")
source("./roca/src.R") #source functions 

##################
### Input Data ###
##################
sites=read.csv("./data/sites.csv")
c14dates=read.csv("./data/c14dates.csv")

c14dates.kanto=subset(c14dates,Region==1)
c14dates.aomori=subset(c14dates,Region==2)
c14dates.hokkaido=subset(c14dates,Region==3)

### Create Bins ###
binSize=200 ## Binsize set to 200yrs (cf Timpson et al 2015)
bins.kanto<-binPrep(sites=c14dates.kanto$SiteID,dates=c14dates.kanto$C14Age,h=binSize)
bins.aomori<-binPrep(sites=c14dates.aomori$SiteID,dates=c14dates.aomori$C14Age,h=binSize)
bins.hokkaido<-binPrep(sites=c14dates.hokkaido$SiteID,dates=c14dates.hokkaido$C14Age,h=binSize)



result005=roca(bins=bins.aomori,date=c14dates.aomori$C14Age,
    sd=c14dates.aomori$C14Error,marine=rep(FALSE,length(c14dates.aomori$C14Age)),
    DeltaR=rep(NA,length(c14dates.aomori$C14Age)),
    DeltaRsd=rep(NA,length(c14dates.aomori$C14Age)),
    yearRange=c(7000,3000),
    resolution=10,
    calCurves=c14dates.aomori$calCurve,
    nsim=5000,lag=5,rollrange=200)
result025=roca(bins=bins.aomori,date=c14dates.aomori$C14Age,
    sd=c14dates.aomori$C14Error,marine=rep(FALSE,length(c14dates.aomori$C14Age)),
    DeltaR=rep(NA,length(c14dates.aomori$C14Age)),
    DeltaRsd=rep(NA,length(c14dates.aomori$C14Age)),
    yearRange=c(7000,3000),
    resolution=10,
    calCurves=c14dates.aomori$calCurve,
    nsim=5000,lag=25,rollrange=200)
result050=roca(bins=bins.aomori,date=c14dates.aomori$C14Age,
    sd=c14dates.aomori$C14Error,marine=rep(FALSE,length(c14dates.aomori$C14Age)),
    DeltaR=rep(NA,length(c14dates.aomori$C14Age)),
    DeltaRsd=rep(NA,length(c14dates.aomori$C14Age)),
    yearRange=c(7000,3000),
    resolution=10,
    calCurves=c14dates.aomori$calCurve,
    nsim=5000,lag=50,rollrange=200)
result100=roca(bins=bins.aomori,date=c14dates.aomori$C14Age,
    sd=c14dates.aomori$C14Error,marine=rep(FALSE,length(c14dates.aomori$C14Age)),
    DeltaR=rep(NA,length(c14dates.aomori$C14Age)),
    DeltaRsd=rep(NA,length(c14dates.aomori$C14Age)),
    yearRange=c(7000,3000),
    resolution=10,
    calCurves=c14dates.aomori$calCurve,
    nsim=5000,lag=100,rollrange=200)


par(mfrow=c(2,2))
rocaplot(result005,main="lag=5yrs")
rocaplot(result025,main="lag=25yrs")
rocaplot(result050,main="lag=50yrs")
rocaplot(result100,main="lag=100yrs")


            
######### Codes ##########
