####################################
### Load Functions and Libraries ###
####################################

library(Bchron) #for calibration of 14C dates
library(zoo) #for plotting SPD outputs

source("./src.R")

##################
### Input Data ###
##################
sites=read.csv("./sites.csv")
c14dates=read.csv("./c14dates.csv")
save.image("./RDatas/intermediate1.RData")
##save.image("./RDatas/SI_26/intermediate1.RData") for DeltaC < -26 Threshold 


################################################################
### SPD Analysis Part 1: Exponential and Uniform Null Models ###
################################################################

### Kanto ###
set.seed(12345)
unif.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,sd=c14dates.kanto$C14Error,
    marine=rep(FALSE,nrow(c14dates.kanto)),DeltaR=rep(NA,nrow(c14dates.kanto)),DeltaRsd=rep(NA,nrow(c14dates.kanto)),
    yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,resolution=10,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,sd=c14dates.kanto$C14Error,
    marine=rep(FALSE,nrow(c14dates.kanto)),DeltaR=rep(NA,nrow(c14dates.kanto)),DeltaRsd=rep(NA,nrow(c14dates.kanto)),
    yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,resolution=10,edge=500,model="exponential",nsim=10000)

save.image("./part1.RData") 
### Aomori ###

set.seed(12345)
unif.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,sd=c14dates.aomori$C14Error,
    marine=rep(FALSE,nrow(c14dates.aomori)),DeltaR=rep(NA,nrow(c14dates.aomori)),DeltaRsd=rep(NA,nrow(c14dates.aomori)),
    yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,resolution=10,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,sd=c14dates.aomori$C14Error,
    marine=rep(FALSE,nrow(c14dates.aomori)),DeltaR=rep(NA,nrow(c14dates.aomori)),DeltaRsd=rep(NA,nrow(c14dates.aomori)),
    yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,resolution=10,edge=500,model="exponential",nsim=10000)

save.image("./part2.RData") 

### Hokkaido ###
set.seed(12345)
unif.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,sd=c14dates.hokkaido$C14Error,
    marine=rep(FALSE,nrow(c14dates.hokkaido)),DeltaR=rep(NA,nrow(c14dates.hokkaido)),DeltaRsd=rep(NA,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,resolution=10,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,sd=c14dates.hokkaido$C14Error,
    marine=rep(FALSE,nrow(c14dates.hokkaido)),DeltaR=rep(NA,nrow(c14dates.hokkaido)),DeltaRsd=rep(NA,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,resolution=10,edge=500,model="exponential",nsim=10000)

save.image("./part3.RData") 

load("./part1.RData") 
load("./part2.RData") 
load("./part3.RData") 
save.image("./RDatas/intermediate2.RData") 


#############################################
### SPD Analysis Part 2: Permutation Test ###
#############################################


##Hokkaido vs Aomori
HokAom<-subset(c14dates,Region==2|Region==3)
HokAom$Region=HokAom$Region-1
binsFullHokAom<-binPrep(sites=HokAom$SiteID,dates=HokAom$C14Age,h=binSize)

set.seed(12345)
HokAomResLong<-permutationTest(regions=HokAom$Region,bins=binsFullHokAom,
                               date=HokAom$C14Age,sd=HokAom$C14Error,
                               marine=rep(FALSE,nrow(HokAom)),
                               DeltaR=rep(NA,nrow(HokAom)),
                               DeltaRsd=rep(NA,nrow(HokAom)),
                               yearRange=c(7000,3000),resolution=10,
                               raw=FALSE,nsim=10000,calCurves=HokAom$calCurve)
set.seed(12345)
HokAomResShort<-permutationTest(regions=HokAom$Region,bins=binsFullHokAom,
                                date=HokAom$C14Age,sd=HokAom$C14Error,
                                marine=rep(FALSE,nrow(HokAom)),
                                DeltaR=rep(NA,nrow(HokAom)),
                                DeltaRsd=rep(NA,nrow(HokAom)),
                                yearRange=c(7000,4420),resolution=10,
                                raw=FALSE,nsim=10000,calCurves=HokAom$calCurve)

##Aomori vs Kanto
AomKanto<-subset(c14dates,Region==1|Region==2)
binsFullAomKanto<-binPrep(sites=AomKanto$SiteID,dates=AomKanto$C14Age,h=binSize)


set.seed(12345)
AomKantoResLong<-permutationTest(regions=AomKanto$Region,bins=binsFullAomKanto,
                                 date=AomKanto$C14Age,sd=AomKanto$C14Error,
                                 marine=rep(FALSE,nrow(AomKanto)),
                                 DeltaR=rep(NA,nrow(AomKanto)),
                                 DeltaRsd=rep(NA,nrow(AomKanto)),
                                 yearRange=c(7000,3000),resolution=10,
                                 raw=FALSE,nsim=10000,calCurves=AomKanto$calCurve)

set.seed(12345)
AomKantoResShort<-permutationTest(regions=AomKanto$Region,bins=binsFullAomKanto,
                                 date=AomKanto$C14Age,sd=AomKanto$C14Error,
                                 marine=rep(FALSE,nrow(AomKanto)),
                                 DeltaR=rep(NA,nrow(AomKanto)),
                                 DeltaRsd=rep(NA,nrow(AomKanto)),
                                 yearRange=c(7000,4420),resolution=10,
                                 raw=FALSE,nsim=10000,calCurves=AomKanto$calCurve)



##Hokkaido vs Kanto

HokkaidoKanto=subset(c14dates,Region==1|Region==3)
HokkaidoKanto$Region[which(HokkaidoKanto$Region==3)]=2

binsFullHokkaidoKanto<-binPrep(sites=HokkaidoKanto$SiteID,dates=HokkaidoKanto$C14Age,h=binSize)

set.seed(12345)
HokkaidoKantoResLong<-permutationTest(regions=HokkaidoKanto$Region,bins=binsFullHokkaidoKanto,
                                      date=HokkaidoKanto$C14Age,sd=HokkaidoKanto$C14Error,
                                      marine=rep(FALSE,nrow(HokkaidoKanto)),
                                      DeltaR=rep(NA,nrow(HokkaidoKanto)),
                                      DeltaRsd=rep(NA,nrow(HokkaidoKanto)),
                                      yearRange=c(7000,3000),resolution=10,
                                      raw=FALSE,nsim=10000,calCurves=HokkaidoKanto$calCurve)

set.seed(12345)
HokkaidoKantoResShort<-permutationTest(regions=HokkaidoKanto$Region,bins=binsFullHokkaidoKanto,
                                      date=HokkaidoKanto$C14Age,sd=HokkaidoKanto$C14Error,
                                      marine=rep(FALSE,nrow(HokkaidoKanto)),
                                      DeltaR=rep(NA,nrow(HokkaidoKanto)),
                                      DeltaRsd=rep(NA,nrow(HokkaidoKanto)),
                                      yearRange=c(7000,4420),resolution=10,
                                      raw=FALSE,nsim=10000,calCurves=HokkaidoKanto$calCurve)

save.image("./RDatas/final.RData") 
