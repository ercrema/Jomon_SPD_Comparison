## Load DATA ##
load("./RDatas/final.RData")


#######################################
## Figure 1: Map of Sample Locations ##
#######################################
library(rworldmap) 
library(maptools) 
library(maps) 
library(GISTools) 
detach("package:GISTools", unload=TRUE) #Unload GISTools momentarily#
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
points(subset(sites,Prefecture=="Aomori"),pch=20,col=2)
points(subset(sites,Prefecture=="Hokkaido"),pch=20,col=1)
points(subset(sites,Prefecture!="Hokkaido"&Prefecture!="Aomori"),pch=20,col=4)
legend(x=134,y=40.91,legend=c("Aomori","Hokkaido","Kanto"),col=c(2,1,4),pch=20,cex=1.7)
dev.print(device=pdf,"./figures/figure1.pdf") 


#####################################################
## Figure 2: SPD vs Exponential and Uniform Models ##
#####################################################
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

dev.print(device=pdf,"./figures/figure2.pdf") 

#####################################
## Figure 3: SPD Permutation Tests ##
#####################################


layout(matrix(c(1:4,5,8:10,6,11:13,7,14:16),byrow=TRUE,4,4),width=c(0.15,1,1,1),height=c(0.15,1,1,1))
##layout.show(n = 16)
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
dev.print(device=pdf,"./figures/figure3.pdf") 
