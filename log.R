####################################
### Load Functions and Libraries ###
####################################
## R version: 3.1.2
library(Bchron) #for calibration of 14C dates // version 4.1.1
library(rworldmap) #map making  //version 1.3-1
library(maptools) #map making  //version 0.8-36 
library(maps) #map making //version 2.3-9 
library(GISTools) #map making //version 0.7-4
library(mapdata) #map making //version 2.2-5
library(zoo) #for plotting SPD outputs //version 1.7-12

source("./src.R") #source functions 

##################
### Input Data ###
##################
sites=read.csv("./data/sites.csv")
c14dates=read.csv("./data/c14dates.csv")

## Subsetting using Delta13 C Threshold at -26
## c14dates=subset(c14dates,deltaC13<c(-26))
## sites=subset(sites,SiteID%in%unique(c14dates$SiteID))

#########################################
### Plot Site Distribution (Figure 1) ###
#########################################

detach("package:GISTools", unload=TRUE) #Unload GISTools temporarily#
basemap <- getMap(resolution = "high")
layout(matrix(c(1,2,2,2,2,2,2,2,2),3,3),height=c(1,0.8,0.8))
par(mar=c(0,0,0,0))
plot(basemap, col = "lightgrey", border = NA,xlim = c(114,149),ylim = c(24,50))
rect(xleft=135.2895, ybottom=35.16062, xright= 147, ytop=45)
box()   
plot(basemap, col = "lightgrey", border = NA,xlim = c(min(sites$Longitude), 142),ylim = c(min(sites$Latitude),45))
map.scale(x=142.0569, y=35.51394, ratio=FALSE, relwidth=0.25,cex=1.5)
library(GISTools)
north.arrow(xb=144.2259, yb=36.3012, len=0.25, lab="N",cex=2,col=1)
coordinates(sites)<-c("Longitude","Latitude")
latlong = "+init=epsg:4326"
map("japan","hokkaido",add=TRUE,col="white",lwd=1.5)
map("japan","aomori",add=TRUE,col="white",lwd=1.5)
map("japan","ibaraki",add=TRUE,col="white",lwd=1.5)
map("japan","chiba",add=TRUE,col="white",lwd=1.5)
map("japan","gunma",add=TRUE,col="white",lwd=1.5)
map("japan","tokyo",add=TRUE,col="white",lwd=1.5)
map("japan","saitama",add=TRUE,col="white",lwd=1.5)
map("japan","tokyo",add=TRUE,col="white",lwd=1.5)
map("japan","kanagawa",add=TRUE,col="white",lwd=1.5)
map("japan","tochigi",add=TRUE,col="white",lwd=1.5)
points(sites,pch=20,col=1)
text(x=138.9726,y=41.00709,"Aomori",cex=2.3)
text(x=140.0389,y=43.99951,"Hokkaido",cex=2.3)
text(x=141.4724,y=37.09934,"Kanto",cex=2.3)
lines(x=c(140.8214,141.6318),y=c(40.83573,41.15288))
lines(x=c(141.4925,142.0185),y=c(40.51684,40.64663))
text(x=142.628,y=41.24277,"Aomori City",cex=1.5)
text(x=143.0996,y=40.76567,"Hachinohe City",cex=1.5)

lines(x=c(140.4610,141.2926),y=c(36.28884,36.48591))
lines(x=c(139.7954,139.7954),y=c(36.76946,37.30036))
lines(x=c(139.0946,138.6011),y=c(36.73894,37.11231))
lines(x=c(138.9356,138.1129),y=c(35.97255,35.90670))
lines(x=c(140.3284,140.5558),y=c(35.42324,35.24928))
lines(x=c(139.0884,138.4872),y=c(35.80983,35.66287))
lines(x=c(139.0884,138.4872),y=c(35.80983,35.66287))
lines(x=c(139.0906,138.5351),y=c(35.27300,35.21326))
text(x=138.3618,y=37.30975,"Gunma",cex=1.5)
text(x=140.8808,y=35.09725,"Chiba",cex=1.5)
text(x=139.8198,y=37.47976,"Tochigi",cex=1.5)
text(x=141.5101,y=36.57375,"Ibaraki",cex=1.5)
text(x=137.4042,y=35.94444,"Saitama",cex=1.5)
text(x=138.0187,y=35.56576,"Tokyo",cex=1.5)
text(x=137.6419,y=35.16032,"Kanagawa",cex=1.5)
dev.print(device=pdf,"./figures/figure1.pdf") 


#############################
### Retreive Sample Sizes ###
#############################
nrow(c14dates)
##Total Sample Size (n=1442) Full Sample
##Total Sample Size (n=867) Reduced

table(c14dates$Region)
##Number of 14C Dates per Region (1=Kanto (n=276) ; 2=Aomori (n=259); 3=Hokkaido (n=332)) Reduced
##Number of 14C Dates per Region (1=Kanto (n=406) ; 2=Aomori (n=432); 3=Hokkaido (n=604)) Full Sample

c14dates.kanto=subset(c14dates,Region==1)
c14dates.aomori=subset(c14dates,Region==2)
c14dates.hokkaido=subset(c14dates,Region==3)

### Count Number of Sites ###
length(unique(c14dates.kanto$SiteID))
## n=41 Reduced
## n=47 Full Sample
length(unique(c14dates.aomori$SiteID))
## n=47 Reduced
## n=57 Full Sample
length(unique(c14dates.hokkaido$SiteID))
## n=71 Reduced
## n=82 Full Sample
     
### Create Bins ###
binSize=200 ## Binsize set to 200yrs (cf Timpson et al 2015)
bins.kanto<-binPrep(sites=c14dates.kanto$SiteID,dates=c14dates.kanto$C14Age,h=binSize)
bins.aomori<-binPrep(sites=c14dates.aomori$SiteID,dates=c14dates.aomori$C14Age,h=binSize)
bins.hokkaido<-binPrep(sites=c14dates.hokkaido$SiteID,dates=c14dates.hokkaido$C14Age,h=binSize)

### Retrieve Bin Counts ###
length(unique(bins.kanto))
##Kanto nb=75 Reduced
##Kanto nb=87 Full Sample
length(unique(bins.aomori))
##Aomori nb=89 Reduced
##Aomori nb=128 Full Sample
length(unique(bins.hokkaido))
##Hokkaido nb=137
##Hokkaido nb=187

################################################################
### SPD Analysis Part 1: Exponential and Uniform Null Models ###
################################################################


## N.B : These functions require considerable computing time (> c.3.5 hours each, with an Intel(R) Xeon(R) CPU E3-1280 V2 @ 3.60GHz, 32Gb RAM)

### Kanto ###
load("./intermediate1.RData")
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
load("./intermediate1.RData")
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
load("./intermediate1.RData")
set.seed(12345)
unif.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,sd=c14dates.hokkaido$C14Error,
    marine=rep(FALSE,nrow(c14dates.hokkaido)),DeltaR=rep(NA,nrow(c14dates.hokkaido)),DeltaRsd=rep(NA,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,resolution=10,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,sd=c14dates.hokkaido$C14Error,
    marine=rep(FALSE,nrow(c14dates.hokkaido)),DeltaR=rep(NA,nrow(c14dates.hokkaido)),DeltaRsd=rep(NA,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,resolution=10,edge=500,model="exponential",nsim=10000)
save.image("./part3.RData")



##### Plot RESULTS #######

load("~/github/jomonSPD/intermediate2.RData")
options(scipen=999) #Ensure p-values are displayed in non-scientific annotation 
par(mfrow=c(2,3))
plotSPDNull(unif.kanto)
title(paste("Kanto: uniform null model (p-value=",round(unif.kanto$pval,4),")",sep=""))
plotSPDNull(unif.aomori)
title(paste("Aomori: uniform null model (p-value=",round(unif.aomori$pval,4),")",sep=""))
plotSPDNull(unif.hokkaido)
title(paste("Hokkaido: uniform null model (p-value=",round(unif.hokkaido$pval,4),")",sep=""))

plotSPDNull(exp.kanto)
title(paste("Kanto: exponential null model (p-value<",round(exp.kanto$pval,4),")",sep=""))
legend("topleft",legend=c("SPD","RollingMean","CI","Positive Deviation","Negative Deviation"),col=c(1,1,"lightgrey","indianred","royalblue"),lty=c(1,2,1,1,1),lwd=c(1,2,5,5,5),cex=0.8)
plotSPDNull(exp.aomori)
title(paste("Aomori: exponential null model (p-value<",round(exp.aomori$pval,4),")",sep=""))
plotSPDNull(exp.hokkaido)
title(paste("Hokkaido: exponential null model (p-value<",round(exp.hokkaido$pval,4),")",sep=""))

dev.print(device=pdf,"./figure2.pdf") 

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



## Plot Matrix Example ##

layout(matrix(c(1:4,5,8:10,6,11:13,7,14:16),byrow=TRUE,4,4),width=c(0.15,1,1,1),height=c(0.15,1,1,1))
layout.show(n = 16)
par(mar=c(0,0,0,0))
plot(runif(1),axes=F,xlab="",ylab="",type="n")
plot(runif(1),axes=F,xlab="",ylab="",type="n")
mtext(side=1,"Kanto",line=-2)
plot(runif(1),axes=F,xlab="",ylab="",type="n")
mtext(side=1,"Aomori",line=-2)
plot(runif(1),axes=F,xlab="",ylab="",type="n")
mtext(side=1,"Hokkaido",line=-2)
plot(runif(1),axes=F,xlab="",ylab="",type="n")
mtext(side=2,"Kanto",line=-2)
plot(runif(1),axes=F,xlab="",ylab="",type="n")
mtext(side=2,"Aomori",line=-2)
plot(runif(1),axes=F,xlab="",ylab="",type="n")
mtext(side=2,"Hokkaido",line=-2)


##Kanto:
par(mar=c(1.5,1,1,1))
plot(runif(1),axes=F,xlab="",ylab="",type="n",xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,"NA",cex=2)
plotSPDSim(AomKantoResLong,index=2,main="",yMax=0.01)
plotSPDSim(HokkaidoKantoResLong,index=1,main="",yMax=0.01)

##Aomori
plotSPDSim(AomKantoResLong,index=1,main="",yMax=0.01)
plot(runif(1),axes=F,xlab="",ylab="",type="n",xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,"NA",cex=2)
plotSPDSim(HokAomResLong,index=1,main="",yMax=0.01)

##Hokkaido
plotSPDSim(HokkaidoKantoResLong,index=2,main="",yMax=0.01)
plotSPDSim(HokAomResLong,index=2,main="",yMax=0.01)
plot(runif(1),axes=F,xlab="",ylab="",type="n",xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,"NA",cex=2)
