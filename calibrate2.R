calibrate2<-function (date, error, calibrationCurve, ids = NULL, positions = NULL, 
    pathToCalCurves = system.file("data", package = "Bchron"), 
    dfs = 100,timeRange,resolution=1,DeltaR=0,DeltaRsd=0) 
{
    date=date-DeltaR
    error=error+DeltaRsd
    
    calCurve = calBP = c14BP = calSd = ageGrid = mu = tau1 = list()
    calCurveFile = paste(pathToCalCurves, "/", calibrationCurve, 
        ".txt.gz", sep = "")
    if (!file.exists(calCurveFile)) 
        stop(paste("Calibration curve file", calCurveFile, 
                   "not found"))
    calCurve = as.matrix(read.table(calCurveFile))
    calBP = calCurve[, 1]
    c14BP = calCurve[, 2]
    calSd = calCurve[, 3]
    ageGrid = seq(min(calBP), max(calBP),by = resolution)
    mu = approx(calBP, c14BP, xout = ageGrid)$y
    tau1 = approx(calBP, calSd, xout = ageGrid)$y

    if (date > max(ageGrid) | date < min(ageGrid)) 
        stop(paste("Date", ids[i], "outside of calibration range"))
    tau = error^2 + tau1
    currAgeGrid = ageGrid
    dens = dt((date - mu)/sqrt(tau),df = dfs)
    dens = dens/sum(dens)

    out = cbind(currAgeGrid, dens)
    if (any(!is.na(timeRange)))
        {
            yearRange=seq(timeRange[1],timeRange[2],-resolution)
            yearRangeIndex=which(currAgeGrid%in%yearRange)
            out = cbind(currAgeGrid[yearRangeIndex], dens[yearRangeIndex])
        }    
    return(out)
}
