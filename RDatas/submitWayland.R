load("intermediate1_24.RData")

### Kanto ###
set.seed(12345)
unif.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,error=c14dates.kanto$C14Error,DeltaR=rep(0,nrow(c14dates.kanto)),DeltaRsd=rep(0,nrow(c14dates.kanto)),yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,edge=500,model="uniform",nsim=10000)


set.seed(12345)
exp.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,error=c14dates.kanto$C14Error,DeltaR=rep(0,nrow(c14dates.kanto)),DeltaRsd=rep(0,nrow(c14dates.kanto)),yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part1_24.RData")


### Aomori ###

load("intermediate1_24.RData")

set.seed(12345)
unif.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,error=c14dates.aomori$C14Error,DeltaR=rep(0,nrow(c14dates.aomori)),DeltaRsd=rep(0,nrow(c14dates.aomori)),yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,error=c14dates.aomori$C14Error,DeltaR=rep(0,nrow(c14dates.aomori)),DeltaRsd=rep(0,nrow(c14dates.aomori)),yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part2_24.RData") 



### Hokkaido ###

load("intermediate1_24.RData")

set.seed(12345)
unif.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,error=c14dates.hokkaido$C14Error,DeltaR=rep(0,nrow(c14dates.hokkaido)),DeltaRsd=rep(0,nrow(c14dates.hokkaido)),yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,error=c14dates.hokkaido$C14Error,DeltaR=rep(0,nrow(c14dates.hokkaido)),DeltaRsd=rep(0,nrow(c14dates.hokkaido)),yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part3_24.RData") 



load("./part1_24.RData") 
load("./part2_24.RData") 
load("./part3_24.RData") 



#################################################################


load("intermediate1_26.RData")

### Kanto ###
set.seed(12345)
unif.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,error=c14dates.kanto$C14Error,
    DeltaR=rep(0,nrow(c14dates.kanto)),DeltaRsd=rep(0,nrow(c14dates.kanto)),yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,edge=500,model="uniform",nsim=10000)


set.seed(12345)
exp.kanto=nullTest(bins=bins.kanto,date=c14dates.kanto$C14Age,error=c14dates.kanto$C14Error,
    DeltaR=rep(0,nrow(c14dates.kanto)),DeltaRsd=rep(0,nrow(c14dates.kanto)),
    yearRange=c(7000,3000),calCurves=c14dates.kanto$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part1_26.RData")


### Aomori ###

load("intermediate1_26.RData")

set.seed(12345)
unif.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,error=c14dates.aomori$C14Error,
    DeltaR=rep(0,nrow(c14dates.aomori)),DeltaRsd=rep(0,nrow(c14dates.aomori)),
    yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.aomori=nullTest(bins=bins.aomori,date=c14dates.aomori$C14Age,error=c14dates.aomori$C14Error,
    DeltaR=rep(0,nrow(c14dates.aomori)),DeltaRsd=rep(0,nrow(c14dates.aomori)),
    yearRange=c(7000,3000),calCurves=c14dates.aomori$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part2_26.RData") 



### Hokkaido ###

load("intermediate1_26.RData")

set.seed(12345)
unif.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,error=c14dates.hokkaido$C14Error,
    DeltaR=rep(0,nrow(c14dates.hokkaido)),DeltaRsd=rep(0,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,edge=500,model="uniform",nsim=10000)

set.seed(12345)
exp.hokkaido=nullTest(bins=bins.hokkaido,date=c14dates.hokkaido$C14Age,error=c14dates.hokkaido$C14Error,
    DeltaR=rep(0,nrow(c14dates.hokkaido)),DeltaRsd=rep(0,nrow(c14dates.hokkaido)),
    yearRange=c(7000,3000),calCurves=c14dates.hokkaido$calCurve,edge=500,model="exponential",nsim=10000)

save.image("./part3_26.RData") 





load("./part1_26.RData") 
load("./part2_26.RData") 
load("./part3_26.RData") 
