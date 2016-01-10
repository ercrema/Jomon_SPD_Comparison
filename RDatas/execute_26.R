
load("intermediate2_26.RData")
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


save.image("final_26.RData")
