
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


#####################################################################################################
## Local and Global MC-Hypothesis testing of SPDs with fitted Uniform or Expenonential null Models ##
#####################################################################################################


### PARAMETERS:
## bins ... output of the binPrep() function
## date ... 14C dates
## error ... 14C error
## calCurves ... calibration curve (see Bchron documentation)
## DeltaR ... DeltaR for marine curves
## DeltaRsd ... DeltaRsd for marine curves
## marine ... set to TRUE for marine dates
## yearRange ... output time range
## nsim ... number of simulations
## edge ... edge effect correction
## model ... option for type of null model
## raw ... if set to TRUE outputs individual simulations

nullTest<-function(bins,date,error,DeltaR=0,DeltaRsd=0,yearRange,calCurves,nsim=100,raw=FALSE,edge=500,model=c("uniform","exponential"))    
{
    require(Bchron) # Bchron v4.0 

    ##Calibrate for each Date
    tmp=calibrate(date=0,error=0,timeRange=yearRange) #retrieve size of the matrix
    individualDatesMatrix<-matrix(NA,nrow=nrow(tmp),ncol=length(date))

    print("Calibrating Individual Dates...")
    flush.console()

    pb <- txtProgressBar(min = 1, max = length(date), style=3)

    for (x in 1:length(date))
        {
            setTxtProgressBar(pb, x)
            individualDatesMatrix[,x]=calibrate(date=date[x],error=error[x],
                                     DeltaR=DeltaR[x],DeltaRsd=DeltaRsd[x],
                                     timeRange=yearRange,calCurves[x])[,2]
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
    finalSPD<-apply(binnedMatrix,1,sum)

    ##Normalise to 1
    finalSPD <- finalSPD/sum(finalSPD)

    min(finalSPD[finalSPD!=0])/10000    
    ##Fit Exponential Model ##

    plusoffset=min(finalSPD[finalSPD!=0])/10000 
    finalSPD=finalSPD+plusoffset #add positive jitter to avoid log(0)
    
    fit <- lm(log(finalSPD)~tmp[,1])
    time=seq(min(tmp[,1])-edge,max(tmp[,1])+edge,1)
    
    est <-  exp(fit$coefficients[1]) * exp(time*fit$coefficients[2]) #Null model density estimates
    pweights <- est/sum(est)       


    ##Simulation starts here...

    C14Interval=range(date)
    sim<-matrix(NA,nrow=length(finalSPD),ncol=nsim)
    print("Monte-Carlo test...")
    pb <- txtProgressBar(min = 1, max = nsim, style=3)
    flush.console()

    for (s in 1:nsim)
        {
            setTxtProgressBar(pb, s)
            ##simulate random dates in calendar date and uncalibrate to C14 dates (±edge yrs for edge effect)
            if (model=="uniform")
                {randomDates<-round(runif(length(unique(bins)),rev(yearRange)[1]-edge,rev(yearRange)[2]+edge))}
            if (model=="exponential")
                {randomDates<-round(sample(time,size=length(unique(bins)),prob=pweights))}        
            randomSDs<-sample(size=length(randomDates),error,replace=TRUE)
            simDates<-round(uncalibrate(randomDates,randomSDs,random=TRUE)[,2:3])
            randomDates<-simDates[,1]
            randomSDs<-simDates[,2]
            
            simDateMatrix=matrix(NA,nrow=nrow(tmp),ncol=length(randomDates))
            for (x in 1:length(randomDates))
                {
                    simDateMatrix[,x]=calibrate(date=randomDates[x],error=randomSDs[x],
                                     DeltaR=DeltaR[x],DeltaRsd=DeltaRsd[x],
                                     timeRange=yearRange,calCurves[x])[,2]
                }

            sim[,s]<-apply(simDateMatrix,1,sum)
            sim[,s]=sim[,s]/sum(sim[,s])
            sim[,s]=sim[,s]+plusoffset
        }

    ## Empirical 95% Intervals ##
    lo=apply(sim,1,quantile,prob=0.025)
    hi=apply(sim,1,quantile,prob=0.975)

    ## Z-score ##
    Zsim<-t(apply(sim,1,scale))
    zLo=apply(Zsim,1,quantile,prob=0.025,na.rm=TRUE)
    zHi=apply(Zsim,1,quantile,prob=0.975,na.rm=TRUE)

    ## z-score observed data ##
    Zscore_empirical <- (finalSPD - apply(sim, 1, mean))/apply(sim, 1, sd)
    busts=which(Zscore_empirical< zLo)
    booms=which(Zscore_empirical> zHi)

    busts2=which(finalSPD< lo)
    booms2=which(finalSPD> hi)


    ## compute global p-value ##
    observedStatistic=sum(c(zLo[busts] - Zscore_empirical[busts]),c(Zscore_empirical[booms]-zHi[booms]))

    expectedstatistic=abs(apply(Zsim,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=zLo))+
        apply(Zsim,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=zHi)

    pvalue <- 1 - c(length(expectedstatistic[expectedstatistic <= observedStatistic]))/c(length(expectedstatistic)+1)

    result=data.frame(calBP=tmp[,1],SPD=finalSPD,lo=lo,hi=hi)
    
    if(raw==FALSE)
        {
            return(list(result=result,pval=pvalue))
        }
    if(raw==TRUE)
        {
            return(list(result=result,sim=sim,pval=pvalue))
        }
}

##################################################
## Permutation Test for comparing multuple SPDs ##
##################################################


### PARAMETERS:
## regions ... value indicating membership to different sets 
## bins ... output of the binPrep() function
## date ... 14C dates
## error ... 14C error
## calCurves ... calibration curve (see Bchron documentation)
## DeltaR ... DeltaR for marine curves
## DeltaRsd ... DeltaRsd for marine curves
## marine ... set to TRUE for marine dates
## yearRange ... output time range
## nsim ... number of simulations
## raw ... if set to TRUE outputs individual simulations


permutationTest<-function(regions,bins,date,error,DeltaR=0,DeltaRsd=0,yearRange,nsim=1000,raw=TRUE,calCurves)
    {
        require(Bchron) # Bchron v4.0
        
        ##Execute Calibration for Each Date
        tmp=calibrate(date=0,error=0,timeRange=yearRange) #retrieve size of the matrix
        individualDatesMatrix<-matrix(NA,nrow=nrow(tmp),ncol=length(date))

        print("Calibrating Individual Dates...")
        flush.console()
        pb <- txtProgressBar(min = 1, max = length(date), style=3)

        for (x in 1:length(date))
            {
              setTxtProgressBar(pb, x)
              individualDatesMatrix[,x]=calibrate(date=date[x],error=error[x],
                                       DeltaR=DeltaR[x],DeltaRsd=DeltaRsd[x],
                                       timeRange=yearRange,calCurves=calCurves[x])[,2]
            }
        close(pb)


        ##Aggregate by Bins 
        binNames<-unique(bins)
        binnedMatrix<-matrix(NA,nrow=nrow(tmp),ncol=length(binNames))
        regionList<-numeric()
        
        print("Binning...")
        flush.console()
        
        pb <- txtProgressBar(min = 1, max = length(binNames), style=3)
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
                regionList[b]=regions[index][1]
            }
        close(pb)

        ##Calculate observed SPD

        observedSPD<-vector("list",length=length(unique(regionList)))
        names(observedSPD)<-unique(regionList)
        
        for (x in 1:length(unique(regionList)))
            {
                focus=unique(regionList)[x]
                index=which(regionList==focus)
                tmpSPD<-apply(binnedMatrix[,index],1,sum)
                tmpSPD=tmpSPD/sum(tmpSPD) #normalise to 1
                observedSPD[[x]]=data.frame(calBP=tmp[,1],SPD=tmpSPD)
            }
        
        ##Permutation Test
        simulatedSPD<-vector("list",length=length(unique(regionList)))
        for (x in 1:length(unique(regionList)))
            {
                simulatedSPD[[x]]=matrix(NA,nrow=nrow(tmp),ncol=nsim)
            }
        print("Permutation test...")
        flush.console()
        pb <- txtProgressBar(min = 1, max = nsim, style=3)
        for (s in 1:nsim)
            {
                setTxtProgressBar(pb, s)
                simRegionList=sample(regionList) #randomize Regions
                for (x in 1:length(unique(simRegionList)))
                    {
                        focus=unique(regionList)[x]
                        index=which(simRegionList==focus)
                        tmpSPD<-apply(binnedMatrix[,index],1,sum)
                        tmpSPD=tmpSPD/sum(tmpSPD) #normalise to 1
                        simulatedSPD[[x]][,s]=tmpSPD
                    }
            }
        close(pb)

        ##Retrieve Summary Stats

        simulatedCIlist<-vector("list",length=length(unique(regionList)))
        
        for (x in 1:length(unique(regionList)))
            {
                simulatedCIlist[[x]]=cbind(apply(simulatedSPD[[x]],1,quantile,prob=c(0.025)),
                                   apply(simulatedSPD[[x]],1,quantile,prob=c(0.975)))
            }


        ##Compute Global p-value:
        pValueList=numeric(length=length(simulatedSPD))
        for (x in 1:length(simulatedSPD))
            {
                ##Create Vector of Means
                zscoreMean=apply(simulatedSPD[[x]],1,mean)
                ##Create Vector of SDs
                zscoreSD=apply(simulatedSPD[[x]],1,sd)
                ##Z-Transform observed and simulated
                tmp.sim=t(apply(simulatedSPD[[x]],1,function(x){return((x - mean(x))/sd(x))}))
                tmp.obs=observedSPD[[x]]
                tmp.obs[,2]=(tmp.obs[,2]-zscoreMean)/zscoreSD
                ##Compute CI
                tmp.ci=t(apply(tmp.sim,1,quantile,prob=c(0.025,0.975)))

                expectedstatistic=abs(apply(tmp.sim,2,function(x,y){a=x-y;i=which(a<0);return(sum(a[i]))},y=tmp.ci[,1]))+
                    apply(tmp.sim,2,function(x,y){a=x-y;i=which(a>0);return(sum(a[i]))},y=tmp.ci[,2])

                lower=tmp.obs[,2]-tmp.ci[,1]
                indexLow=which(tmp.obs[,2]<tmp.ci[,1])
                higher=tmp.obs[,2]-tmp.ci[,2]
                indexHi=which(tmp.obs[,2]>tmp.ci[,2])
                observedStatistic= sum(abs(lower[indexLow]))+sum(higher[indexHi])

                pValueList[[x]]=1
                if (observedStatistic>0)
                    {    
                        pValueList[[x]] <- 1 - c(length(expectedstatistic[expectedstatistic <= observedStatistic]))/c(length(expectedstatistic)+1)
                    }
            }        
        if(raw==FALSE)
            {

                return(list(observed=observedSPD,envelope=simulatedCIlist,pValueList=pValueList))
            }
        if(raw==TRUE)
            {

                return(list(observed=observedSPD,envelope=simulatedCIlist,raw=simulatedSPD,pValueList=pValueList))
            }
        
    }


########################################
## Plot Functions: SPDs vs NULL Model ##
########################################

### PARAMETERS:
## data ... output of the nullTest() function

plotSPDNull<-function(data, ...)
    {
        require(zoo)
        
        obs=data$result[,1:2]
        resolution=1

        envelope=data$result[,3:4]
        yMax=max(envelope,obs[,2])
        
        booms=which(obs[,2]>envelope[,2])
        busts=which(obs[,2]<envelope[,1])
        baseline=rep(0,nrow(obs))
        plot(obs[,1],obs[,2],xlim=c(max(obs[,1]),min(obs[,1])),ylim=c(0,yMax),
             xlab="cal BP",ylab="SPD",type="l",col=1,lwd=0.5,...)

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
                        boomBlocks[[counter]]=vector("list",2)
                        boomBlocks[[counter]][[1]]=boomPlot[x]
                        boomBlocks[[counter]][[2]]=obs[x,1]
                        state="on"
                    }
                if (state=="on")
                    {
                        if (boomPlot[x]>0)
                            {
                                boomBlocks[[counter]][[1]]=c(boomBlocks[[counter]][[1]],boomPlot[x])
                                boomBlocks[[counter]][[2]]=c(boomBlocks[[counter]][[2]],obs[x,1])
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
                if (bustPlot[x]>0&state=="off")
                    {
                        counter=counter+1
                        bustBlocks=c(bustBlocks,vector("list",1))
                        bustBlocks[[counter]]=vector("list",2)
                        bustBlocks[[counter]][[1]]=bustPlot[x]
                        bustBlocks[[counter]][[2]]=obs[x,1]
                        state="on"
                    }
                if (state=="on")
                    {
                        if (bustPlot[x]>0)
                            {
                                bustBlocks[[counter]][[1]]=c(bustBlocks[[counter]][[1]],bustPlot[x])
                                bustBlocks[[counter]][[2]]=c(bustBlocks[[counter]][[2]],obs[x,1])
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
                        polygon(c(boomBlocks[[x]][[2]],rev(boomBlocks[[x]][[2]])),c(rep(+100,length(boomBlocks[[x]][[1]])),rep(-100,length(boomBlocks[[x]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
                    }
            }
        
        if (length(busts)>0) #Fix this line
            {
                for (x in 1:length(bustBlocks))
                    {
                        polygon(c(bustBlocks[[x]][[2]],rev(bustBlocks[[x]][[2]])),c(rep(+100,length(bustBlocks[[x]][[1]])),rep(-100,length(bustBlocks[[x]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
                    }
            }
        
        
        polygon(x=c(obs[,1],rev(obs[,1])),y=c(envelope[,1],rev(envelope[,2])),col=rgb(0,0,0,0.2),border=NA)
        spdSmooth<-rollmean(obs[,2],k=200/resolution,fill=NA)
        lines(obs[,1],spdSmooth,col=1,lwd=2,lty=2)
        axis(side=1,at=seq(max(obs[,1]),min(obs[,1]),-100),labels=NA,tck = -.01)
       
    }

###################################################
## Plot Functions: SPD Permutation test results  ##
########################################

### PARAMETERS:
## data ... output of the permutationTest() function
## index ... region index value for defining plot output
## yMax ... maximum value for the y-axis
plotSPDSim<-function(data,index,yMax=NA, ...)
    {
        require(zoo)
        obs=data$observed[[index]]
        resolution=1

        envelope=data$envelope[[index]]
        if (is.na(yMax))
            {yMax=max(as.numeric(envelope),obs[,2])}
        
        booms=which(obs[,2]>envelope[,2])
        busts=which(obs[,2]<envelope[,1])
        baseline=rep(0,nrow(obs))
        plot(obs[,1],obs[,2],xlim=c(max(obs[,1]),min(obs[,1])),ylim=c(0,yMax),
             xlab="cal BP",ylab="SPD",type="l",col=1,lwd=0.5,axes=FALSE,...)
        axis(side=1,padj=-1)
        axis(side=2,padj=1)
        box()
        
        boomPlot=baseline
        if (length(booms)>0){boomPlot[booms]=obs[booms,2]}
        bustPlot=baseline
        if (length(busts)>0){bustPlot[busts]=obs[busts,2]}                
        
        boomBlocks<-vector("list")
        counter=0
        state="off"
        for (x in 1:length(boomPlot))
            {
                if (boomPlot[x]>0&state=="off")
                    {
                        counter=counter+1
                        boomBlocks=c(boomBlocks,vector("list",1))
                        boomBlocks[[counter]]=vector("list",2)
                        boomBlocks[[counter]][[1]]=boomPlot[x]
                        boomBlocks[[counter]][[2]]=obs[x,1]
                        state="on"
                    }
                if (state=="on")
                    {
                        if (boomPlot[x]>0)
                            {
                                boomBlocks[[counter]][[1]]=c(boomBlocks[[counter]][[1]],boomPlot[x])
                                boomBlocks[[counter]][[2]]=c(boomBlocks[[counter]][[2]],obs[x,1])
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
                if (bustPlot[x]>0&state=="off")
                    {
                        counter=counter+1
                        bustBlocks=c(bustBlocks,vector("list",1))
                        bustBlocks[[counter]]=vector("list",2)
                        bustBlocks[[counter]][[1]]=bustPlot[x]
                        bustBlocks[[counter]][[2]]=obs[x,1]
                        state="on"
                    }
                if (state=="on")
                    {
                        if (bustPlot[x]>0)
                            {
                                bustBlocks[[counter]][[1]]=c(bustBlocks[[counter]][[1]],bustPlot[x])
                                bustBlocks[[counter]][[2]]=c(bustBlocks[[counter]][[2]],obs[x,1])
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
                        polygon(c(boomBlocks[[x]][[2]],rev(boomBlocks[[x]][[2]])),c(rep(+100,length(boomBlocks[[x]][[1]])),rep(-100,length(boomBlocks[[x]][[1]]))),col=rgb(0.7,0,0,0.2),border=NA)
                    }
            }

        if (length(busts)>0)
            {

                for (x in 1:length(bustBlocks))
                    {
                        polygon(c(bustBlocks[[x]][[2]],rev(bustBlocks[[x]][[2]])),c(rep(+100,length(bustBlocks[[x]][[1]])),rep(-100,length(bustBlocks[[x]][[1]]))),col=rgb(0,0,0.7,0.2),border=NA)
                    }
            }
        
        polygon(x=c(obs[,1],rev(obs[,1])),y=c(envelope[,1],rev(envelope[,2])),col=rgb(0,0,0,0.2),border=NA)
        spdSmooth<-rollmean(obs[,2],k=200/resolution,fill=NA)
        lines(obs[,1],spdSmooth,col=1,lwd=2,lty=2)
        axis(side=1,at=seq(max(obs[,1]),min(obs[,1]),-100),labels=NA,tck = -.01)  
    }
