## Load DATA ##
##load("./RDatas/final_24.RData")
load("./RDatas/final_26.RData")

#######################################
## Figure 1: Map of Sample Locations ##
#######################################

postscript("./figures/figure1.eps", height = 6, width = 6,family = "Times", paper = "special", onefile = FALSE, horizontal = FALSE)

library(rworldmap) #map making  //version 1.3-1
library(maptools) #map making  //version 0.8-36 
library(maps) #map making //version 2.3-9 
library(GISTools) #map making //version 0.7-4
library(mapdata) #map making //version 2.2-5

detach("package:GISTools", unload=TRUE) #Unload GISTools temporarily#
basemap <- getMap(resolution = "high")
layout(matrix(c(1,2,2,2,2,2,2,2,2),3,3),height=c(1,0.8,0.8))
par(mar=c(0,0,0,0))
plot(basemap, col = "lightgrey", border = NA,xlim = c(114,149),ylim = c(24,50))
rect(xleft=135.2895, ybottom=35.16062, xright= 147, ytop=45)
box()
##par(new=T,bg="white")
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
dev.off()


#####################################################
## Figure 2: SPD vs Exponential and Uniform Models ##
#####################################################
library(Bchron) #for calibration of 14C dates // version 4.1.1
library(zoo) #for plotting SPD outputs //version 1.7-12
options(scipen=999) #Ensure p-values are displayed in non-scientific annotation 
par(mfrow=c(2,3),family="Times")
par(mar=c(4.5, 4, 1.5, 1))
plotSPDNull(unif.kanto,yMax=0.0009)
title("Kanto: uniform null model")
plotSPDNull(unif.aomori,yMax=0.0009)
title("Aomori: uniform null model")
plotSPDNull(unif.hokkaido,yMax=0.0009)
title("Hokkaido: uniform null model")

plotSPDNull(exp.kanto,yMax=0.0009)
title("Kanto: exponential null model")
plotSPDNull(exp.aomori,yMax=0.0009)
title("Aomori: exponential null model")
plotSPDNull(exp.hokkaido,yMax=0.0009)
title("Hokkaido: exponential null model")
legend("topleft",legend=c("SPD","RollingMean","CI","Positive Deviation","Negative Deviation"),col=c(1,1,"lightgrey","indianred","royalblue"),lty=c(1,1,1,1,1),lwd=c(0.5,2,5,5,5),cex=0.8,bg="white")
dev.print(device=pdf,useDingbats=FALSE,"~/github/jomonSPD/figures/figure2.pdf")
##dev.print(device=pdf,useDingbats=FALSE,"~/github/jomonSPD/figures/figure2_si.pdf")

#####################################
## Figure 3: SPD Permutation Tests ##
#####################################
library(Bchron) #for calibration of 14C dates // version 4.1.1
library(zoo) #for plotting SPD outputs //version 1.7-12



layout(matrix(c(1:4,5,8:10,6,11:13,7,14:16),byrow=TRUE,4,4),width=c(0.15,1,1,1),height=c(0.15,1,1,1))
##layout.show(n = 16)
par(mar=c(0,0,0,0),family="Times")
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
plotSPDSim(AomKantoResLong,index=2,main="",yMax=0.0009)
plotSPDSim(HokkaidoKantoResLong,index=1,main="",yMax=0.0009)

##Aomori
plotSPDSim(AomKantoResLong,index=1,main="",yMax=0.0009)
plot(runif(1),axes=F,xlab="",ylab="",type="n",xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,"NA",cex=2)
plotSPDSim(HokAomResLong,index=1,main="",yMax=0.0009)

##Hokkaido
plotSPDSim(HokkaidoKantoResLong,index=2,main="",yMax=0.0009)
plotSPDSim(HokAomResLong,index=2,main="",yMax=0.0009)
legend("topleft",legend=c("SPD","RollingMean","CI","Positive Deviation","Negative Deviation"),col=c(1,1,"lightgrey","indianred","royalblue"),lty=c(1,1,1,1,1),lwd=c(0.5,2,5,5,5),cex=0.8,bg="white")
plot(runif(1),axes=F,xlab="",ylab="",type="n",xlim=c(0,1),ylim=c(0,1))
text(0.5,0.5,"NA",cex=2)
dev.print(device=pdf,useDingbats=FALSE,"./figures/figure3.pdf")
##dev.print(device=pdf,useDingbats=FALSE,"~/github/jomonSPD/figures/figure3_si.pdf")
