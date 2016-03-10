####################################
### Load Functions and Libraries ###
####################################

library(Bchron) #for calibration of 14C dates
library(zoo) #for plotting SPD outputs

source("./src.R")

##################
### Input Data ###
##################
sites=read.csv("./data/sites.csv")
c14dates=read.csv("./data/c14dates.csv")

##c14dates=subset(c14dates,deltaC13<c(-26))
##sites=subset(sites,SiteID%in%unique(c14dates$SiteID))

#################
### Data Prep. ##
#################

### Subset dates ###

c14dates.kanto=subset(c14dates,Region==1)
c14dates.aomori=subset(c14dates,Region==2)
c14dates.hokkaido=subset(c14dates,Region==3)

### Create Bins ###
binSize=200 ## Binsize set to 200yrs (cf Timpson et al 2015)
bins.kanto<-binPrep(sites=c14dates.kanto$SiteID,dates=c14dates.kanto$C14Age,h=binSize)
bins.aomori<-binPrep(sites=c14dates.aomori$SiteID,dates=c14dates.aomori$C14Age,h=binSize)
bins.hokkaido<-binPrep(sites=c14dates.hokkaido$SiteID,dates=c14dates.hokkaido$C14Age,h=binSize)


save.image("./RDatas/intermediate1_24.RData")
##save.image("./RDatas/intermediate1_26.RData") # for DeltaC < -26 Threshold 

################################################################
### SPD Analysis Part 1: Exponential and Uniform Null Models ###
################################################################

### Kanto ###
set.seed(12345)
unif.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,error=c14dates.kanto$C14Error,
    DeltaR=rep(0,nrow(c14dates.kanto)),DeltaRsd=rep(0,nrow(c14dates.kanto)),yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,edge=500,model="uniform",nsim=10000)


set.seed(12345)
exp.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,error=c14dates.kanto$C14Error,
    DeltaR=rep(0,nrow(c14dates.kanto)),DeltaRsd=rep(0,nrow(c14dates.kanto)),
    yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part1.RData") 
### Aomori ###

set.seed(12345)
unif.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,error=c14dates.aomori$C14Error,
    DeltaR=rep(0,nrow(c14dates.aomori)),DeltaRsd=rep(0,nrow(c14dates.aomori)),
    yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,error=c14dates.aomori$C14Error,
    DeltaR=rep(0,nrow(c14dates.aomori)),DeltaRsd=rep(0,nrow(c14dates.aomori)),
    yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part2.RData") 

### Hokkaido ###
set.seed(12345)
unif.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,error=c14dates.hokkaido$C14Error,
    DeltaR=rep(0,nrow(c14dates.hokkaido)),DeltaRsd=rep(0,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,error=c14dates.hokkaido$C14Error,
    DeltaR=rep(0,nrow(c14dates.hokkaido)),DeltaRsd=rep(0,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,edge=500,model="exponential",nsim=10000)

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
                               date=HokAom$C14Age,error=HokAom$C14Error,
                               DeltaR=rep(0,nrow(HokAom)),
                               DeltaRsd=rep(0,nrow(HokAom)),
                               yearRange=c(7000,3000),
                               raw=FALSE,nsim=10000,calCurves=HokAom$calCurve)
set.seed(12345)
HokAomResShort<-permutationTest(regions=HokAom$Region,bins=binsFullHokAom,
                                date=HokAom$C14Age,error=HokAom$C14Error,
                                DeltaR=rep(0,nrow(HokAom)),
                                DeltaRsd=rep(0,nrow(HokAom)),
                                yearRange=c(7000,4420),
                                raw=FALSE,nsim=10000,calCurves=HokAom$calCurve)

##Aomori vs Kanto
AomKanto<-subset(c14dates,Region==1|Region==2)
binsFullAomKanto<-binPrep(sites=AomKanto$SiteID,dates=AomKanto$C14Age,h=binSize)


set.seed(12345)
AomKantoResLong<-permutationTest(regions=AomKanto$Region,bins=binsFullAomKanto,
                                 date=AomKanto$C14Age,error=AomKanto$C14Error,
                                 DeltaR=rep(0,nrow(AomKanto)),
                                 DeltaRsd=rep(0,nrow(AomKanto)),
                                 yearRange=c(7000,3000),
                                 raw=FALSE,nsim=10000,calCurves=AomKanto$calCurve)

set.seed(12345)
AomKantoResShort<-permutationTest(regions=AomKanto$Region,bins=binsFullAomKanto,
                                 date=AomKanto$C14Age,error=AomKanto$C14Error,
                                 DeltaR=rep(0,nrow(AomKanto)),
                                 DeltaRsd=rep(0,nrow(AomKanto)),
                                 yearRange=c(7000,4420),
                                 raw=FALSE,nsim=10000,calCurves=AomKanto$calCurve)



##Hokkaido vs Kanto

HokkaidoKanto=subset(c14dates,Region==1|Region==3)
HokkaidoKanto$Region[which(HokkaidoKanto$Region==3)]=2

binsFullHokkaidoKanto<-binPrep(sites=HokkaidoKanto$SiteID,dates=HokkaidoKanto$C14Age,h=binSize)

set.seed(12345)
HokkaidoKantoResLong<-permutationTest(regions=HokkaidoKanto$Region,bins=binsFullHokkaidoKanto,
                                      date=HokkaidoKanto$C14Age,error=HokkaidoKanto$C14Error,
                                      DeltaR=rep(0,nrow(HokkaidoKanto)),
                                      DeltaRsd=rep(0,nrow(HokkaidoKanto)),
                                      yearRange=c(7000,3000),
                                      raw=FALSE,nsim=10000,calCurves=HokkaidoKanto$calCurve)

set.seed(12345)
HokkaidoKantoResShort<-permutationTest(regions=HokkaidoKanto$Region,bins=binsFullHokkaidoKanto,
                                      date=HokkaidoKanto$C14Age,error=HokkaidoKanto$C14Error,
                                      DeltaR=rep(0,nrow(HokkaidoKanto)),
                                      DeltaRsd=rep(0,nrow(HokkaidoKanto)),
                                      yearRange=c(7000,4420),
                                      raw=FALSE,nsim=10000,calCurves=HokkaidoKanto$calCurve)


### Results p-value Matrix ###

resultPairwiseLong=matrix(NA,3,3)

row.names(resultPairwiseLong)=c("Kanto","Aomori","Hokkaido")
colnames(resultPairwiseLong)=c("Kanto","Aomori","Hokkaido")

resultPairwiseShort<-resultPairwiseLong


resultPairwiseLong[1,2]=AomKantoResLong$pValueList[2]
resultPairwiseLong[1,3]=HokkaidoKantoResLong$pValueList[1]

resultPairwiseLong[2,1]=AomKantoResLong$pValueList[1]
resultPairwiseLong[2,3]=HokAomResLong$pValueList[1]

resultPairwiseLong[3,1]=HokkaidoKantoResLong$pValueList[2]
resultPairwiseLong[3,2]=HokAomResLong$pValueList[2]


resultPairwiseShort[1,2]=AomKantoResShort$pValueList[2]
resultPairwiseShort[1,3]=HokkaidoKantoResShort$pValueList[1]

resultPairwiseShort[2,1]=AomKantoResShort$pValueList[1]
resultPairwiseShort[2,3]=HokAomResShort$pValueList[1]

resultPairwiseShort[3,1]=HokkaidoKantoResShort$pValueList[2]
resultPairwiseShort[3,2]=HokAomResShort$pValueList[2]


save.image("./RDatas/final.RData") 
