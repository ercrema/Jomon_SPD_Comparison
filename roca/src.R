
#################################################
## Binning Functions  for clustering 14C Dates ##
#################################################

### PARAMETERS:
## sites ... site IDs
## dates ... 14C dates
## h ... clustering range (in 14C yrs)

binPrep<-function(sites,dates,h=200)
    {
        clusters<-rep(NA,length(sites))
	for  (x in 1:length(unique(sites)))
            {
                index<-which(sites==unique(sites)[x])
                if (length(index)>1)
                    {
                        clusters[index]<-paste(unique(sites)[x],cutree(hclust(dist(dates[index])),h=h),sep="_")
                    }
                if (length(index)==1)
                    {
                        clusters[index]<-paste(unique(sites)[x],"1",sep="_")
                    }                                
            }
        return(clusters)
    }

#############################################################################
## Calibration Function (wrapper for BChron's  BchronCalibrate() function) ##
#############################################################################

### PARAMETERS:
## date ... 14C dates
## error ... 14C error
## calCurves ... calibration curve (see Bchron documentation)
## DeltaR ... DeltaR for marine curves
## DeltaRsd ... DeltaRsd for marine curves
## timeRange ... output time range

calibrate<-function(date, error, calCurves='intcal13', DeltaR=0 ,DeltaRsd=0,timeRange=c(10000,0)) 
{
    require(Bchron)
    date=date-DeltaR
    error=error+DeltaRsd
    tmp = BchronCalibrate(ages=date,ageSds=error,calCurves=calCurves,eps=0)
    calBP=rev(as.numeric(tmp[[1]][4][[1]]))
    prob=rev(as.numeric(tmp[[1]][[5]]))
    calBP.out=seq(50000,0,-1)
    prob.out=rep(0,length=length(calBP.out))
    index=which(calBP.out%in%calBP)
    prob.out[index]=prob
    res=cbind(calBP.out,prob.out)
    res=res[which(calBP.out<=timeRange[1]&calBP.out>=timeRange[2]),]
    return(res)
}

##############################
## "Uncalibration" Function ##
##############################

### PARAMETERS:
## dates ... calendar dates
## error ... obserbved error ranges
## calCurves ... calibration curve (see Bchron documentation)
## random ... if set to TRUE generate random 14C errors


uncalibrate<-function(dates,error,calCurves='intcal13',random=TRUE)
{
    require(Bchron) # Bchron v4.0 
    pathToCalCurves=system.file("data", package = "Bchron")
    calCurveFile = paste(pathToCalCurves, "/", calCurves,".txt.gz", sep = "")
    calcurve=as.matrix(read.table(calCurveFile))[,1:3]
    colnames(calcurve) <- c("CALBP", "C14BP", "Error")

    ## uncalibrate CAL BP dates, interpolating with approx
    dates <- data.frame(approx(calcurve, xout = dates))
    colnames(dates) <- c("CALBP", "C14BP")
    calcurve.error <- approx(calcurve[,c(1,3)], xout = dates$CALBP)$y
    dates$Error <- sqrt(error^2 + calcurve.error^2)
    if(random==TRUE){dates$C14.Age=round(rnorm(nrow(dates),mean=dates$C14BP,sd=dates$Error))}
    return(dates)
}

###################################
## Local Rate of Change Analysis ##
###################################

### PARAMETERS:
## bins ... output of the binPrep() function
## date ... 14C dates
## sd ... 14C error
## calCurves ... calibration curve (see Bchron documentation)
## resolution ... output resolution (in cal yrs)
## DeltaR ... DeltaR for marine curves
## DeltaRsd ... DeltaRsd for marine curves
## marine ... set to TRUE for marine dates
## yearRange ... output time range
## nsim ... number of bootstrap simulations
## lag ... lagged difference for computing local linear growth rate (in cal yrs)
## rollrange ... rolling mean interval (in cal yrs)

roca<-function(bins,date,sd,marine=FALSE,DeltaR=NA,DeltaRsd=NA,yearRange,
               resolution,calCurves,nsim=100,lag=10,rollrange=200)
{
    require(Bchron)
    require(zoo)
    ##Calibrate for each Date
    tmp=calibrate(date=0,sd=0,resolution=resolution,timeRange=yearRange) #retrieve size of the matrix
    individualDatesMatrix<-matrix(NA,nrow=nrow(tmp),ncol=length(date))

    print("Calibrating Individual Dates...")
    flush.console()

    pb <- txtProgressBar(min = 1, max = length(date), style=3)

    for (x in 1:length(date))
        {
            setTxtProgressBar(pb, x)
            individualDatesMatrix[,x]=calibrate(date=date[x],sd=sd[x],marine=marine[x],
                                     DeltaR=DeltaR[x],DeltaRsd=DeltaRsd[x],timeRange=yearRange,
                                     resolution=resolution,calCurves[x])[,2]
        }


    ##Aggregate by Bins
    binNames<-unique(bins)
    binnedMatrix<-matrix(NA,nrow=nrow(tmp),ncol=length(binNames))

    print("Binning...")
    flush.console()

    pb <- txtProgressBar(min = 1, max = length(binNames), style=3,title="Binning...")
    for (b in 1:length(binNames))
        {
            setTxtProgressBar(pb, b)
            index=which(bins==binNames[b])
            if (length(index)>1)
                {    
                    spd.tmp=apply(individualDatesMatrix[,index],1,sum)
                    binnedMatrix[,b]=spd.tmp/length(index)
                }
            else
                {
                    binnedMatrix[,b]=individualDatesMatrix[,index]
                }
        }
    close(pb)
    
    ##Create observed SPD
    obsSPD<-apply(binnedMatrix,1,sum)
    ##Reverse Order & Compute Rolled Mean SPD
    obsSPD<-rev(obsSPD)
    rollingMeanObsSPD=rollmean(obsSPD,k=rollrange/resolution,fill=NA)
    ##Compute Observed Lagged SPD
    laggedSPD=diff(rollingMeanObsSPD,lag)/(resolution*lag)


    ##Bootstrap Analysis##
    BOOTres=matrix(NA,nrow=nsim,ncol=length(laggedSPD))

    print("Bootstrapping...")
    flush.console()
    pb <- txtProgressBar(min = 1, max = nsim, style=3)
    for (s in 1:nsim)
        {
            setTxtProgressBar(pb, s)
            index=sample(1:ncol(binnedMatrix),replace=TRUE)
            tmp.simSPD<-apply(binnedMatrix[,index],1,sum)
            tmp.simSPD<-rev(tmp.simSPD)
            rollingMeanSimSPD=rollmean(tmp.simSPD,k=rollrange/resolution,fill=NA)
            BOOTres[s,]=diff(rollingMeanSimSPD,lag)/(resolution*lag)
        }
    close(pb)


    lo=apply(BOOTres,2,quantile,prob=0.025,na.rm=TRUE)
    hi=apply(BOOTres,2,quantile,prob=0.975,na.rm=TRUE)

    time=tmp[,1]
    timeSEQ=c(time+resolution/2*lag)[1:c(length(time)-lag)]

    obs=cbind(rev(timeSEQ),laggedSPD)
    envelope=cbind(lo,hi)
    return(list(obs=obs,envelope=envelope))
}


########################################
## Plot Local Rate of Change Analysis ##
########################################

### PARAMETERS:
## result ... output of roca()
## ...    ... other plot function parameters

rocaplot<-function(result,...)
{
    envelope=result$envelope
    obs=result$obs
    booms=which(envelope[,1]>0)
    busts=which(envelope[,2]<0)
    baseline=rep(0,nrow(obs))
    plot(obs[,1],obs[,2],xlim=c(max(obs[,1]),min(obs[,1])),ylim=c(min(envelope,na.rm=TRUE),max(envelope,na.rm=TRUE)),
         xlab="cal BP",ylab="Rate of Change",type="l",col=1,lwd=1,...)

    boomPlot=baseline
    boomPlot[booms]=obs[booms,2]
    bustPlot=baseline
    bustPlot[busts]=obs[busts,2]



    boomBlocks<-vector("list")
    counter=0
    state="off"
    for (x in 1:length(boomPlot))
        {
            if (boomPlot[x]>0&state=="off")
                {
                    counter=counter+1
                    boomBlocks=c(boomBlocks,vector("list",1))
                    boomBlocks[[counter]]=x
                    state="on"
                }
            if (state=="on")
                {
                    if (boomPlot[x]>0)
                        {
                            boomBlocks[[counter]]=c(boomBlocks[[counter]],x)
                        }
                    if (boomPlot[x]==0)
                        {
                            state="off"
                        }
                }
            
        }


    bustBlocks<-vector("list")
    counter=0
    state="off"
    for (x in 1:length(bustPlot))
        {
            if (bustPlot[x]<0&state=="off")
                {
                    counter=counter+1
                    bustBlocks=c(bustBlocks,vector("list",1))
                    bustBlocks[[counter]]=x
                    state="on"
                }
            if (state=="on")
                {
                    if (bustPlot[x]<0)
                        {
                            bustBlocks[[counter]]=c(bustBlocks[[counter]],x)
                        }
                    if (bustPlot[x]==0)
                        {
                            state="off"
                        }
                }
            
        }



    if (length(booms)>0) 
        {
            for (x in 1:length(boomBlocks))
                {
                    polygon(c(obs[boomBlocks[[x]],1],rev(obs[boomBlocks[[x]],1])),c(envelope[boomBlocks[[x]],1],rev(envelope[boomBlocks[[x]],2])),col=rgb(0.7,0,0,0.2),border=NA)
                }
        }
    
    if (length(busts)>0) #Fix this line
        {
            for (x in 1:length(bustBlocks))
                {
                    polygon(c(obs[bustBlocks[[x]],1],rev(obs[bustBlocks[[x]],1])),c(envelope[bustBlocks[[x]],1],rev(envelope[bustBlocks[[x]],2])),col=rgb(0,0,0.7,0.2),border=NA)
                }
        }

    lines(obs[,1],envelope[,1],lty=2)
    lines(obs[,1],envelope[,2],lty=2)
    abline(h=0,col=2,lty=3)

}
