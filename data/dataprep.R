#### Data Preparation Code (Not to be made public in the final version) ####

source("~/github/SPDC14/src/calibrate.R")
source("~/github/SPDC14/src/binFuns.R")
source('~/github/SPDC14/src/permutationTest.R')
source('~/github/SPDC14/src/plotFuns.R')
source('~/github/SPDC14/src/uniformTest.R')
source('~/github/SPDC14/src/spd.R')
source('~/github/SPDC14/src/nullTest.R')
source('~/github/SPDC14/src/uncalibrate.R')

### Load Data ###
library("RSQLite")
con=dbConnect(SQLite(),dbname="~/Dropbox/nicoss/C14Database/JomonC14.sqlite")
alltables = dbListTables(con)
C14Samples = dbGetQuery( con,'SELECT * FROM C14Samples' )
Sites = dbGetQuery( con,'SELECT * FROM Sites' )

### General Parameters ###
deltaC13Cutoff=-24 
temporalScopeData = c(7500,2500)  #in C14Age
temporalScopeAnalysisLong = c(7000,3000) #in calYr 
temporalScopeAnalysisShort = c(7000,4420) #in calYr (Early~Late Jomon period; Kobayashi's Scheme)
binSize = 200 #years



### Clean Data ###
KantoSitesIndex<-Sites$SiteID[which(Sites$Prefecture%in%c("Ibaragi","Gunma","Saitama","Kanagawa","Chiba","Tokyo","Tochigi"))]
AomoriSitesIndex<-Sites$SiteID[which(Sites$Prefecture=="Aomori")]
HokkaidoSitesIndex<-Sites$SiteID[which(Sites$Prefecture=="Hokkaido")]

## Subset to Region and Chronological Range
dataSet=subset(C14Samples,SiteID%in%c(KantoSitesIndex,AomoriSitesIndex,HokkaidoSitesIndex)&C14Age>=temporalScopeData[2]&C14Age<=temporalScopeData[1]) #subset to 

##Use corrected C14 dates when available##
for (x in 1:nrow(dataSet))
    {
        if(!is.na(dataSet$C14AgeCorrected[x]))
            {
                dataSet$C14Age[x]=dataSet$C14AgeCorrected[x]
                dataSet$C14Error[x]=dataSet$C14ErrorCorrected[x]
            }

    }

##Remove Marine##
dataSet$marine=FALSE
dataSet$marine[which(dataSet$Material=="Shell"
                     |dataSet$Material=="Fishbone")]=TRUE
dataSet=subset(dataSet,marine==FALSE)
dataSet$DeltaR=NA
dataSet$DeltaRsd=NA

## Define Regions ##
dataSet$Region=NA
dataSet$Region[which(dataSet$SiteID%in%KantoSitesIndex)]=1
dataSet$Region[which(dataSet$SiteID%in%AomoriSitesIndex)]=2
dataSet$Region[which(dataSet$SiteID%in%HokkaidoSitesIndex)]=3

## Define Calibration Curve ##
dataSet$calCurve='intcal13'

## Define delta13C subset (using 26 threshold) ##
dataSet = dataSet[which(dataSet$deltaC13<c(deltaC13Cutoff)),]


## Extract only columns that are needed ##

dataSet=data.frame(C14ID=dataSet$C14ID,
    SiteID=dataSet$SiteID,
    C14Age=dataSet$C14Age,
    C14Error=dataSet$C14Error,
    LabCode=dataSet$LabCode,
    deltaC13=dataSet$deltaC13,
    deltaC13Error=dataSet$deltaC13Error,
    Material=dataSet$Material,
    calCurve=dataSet$calCurve,
    Region=dataSet$Region,
    Ref=dataSet$Ref)

Sites=Sites[which(Sites$SiteID%in%dataSet$SiteID),]

Sites=data.frame(SiteID=Sites$SiteID,
    SiteName=Sites$SiteName,
    Latitude=Sites$Latitude,
    Longitude=Sites$Longitude,
    Prefecture=Sites$Prefecture)


C14Dates=dataSet
C14Dates$Material=as.character(C14Dates$Material)
C14Dates$Material[which(C14Dates$Material=="Carbon")]="Charcoal"
C14Dates$Material[which(C14Dates$Material=="Pottery Residue")]="Pottery Food Residue"
C14Dates$Material[which(C14Dates$Material=="Carbonised Material")]="Charred Remain"
C14Dates$Material[which(C14Dates$Material=="Seed")]="Seed/Nut"
C14Dates$Material[which(C14Dates$Material=="Wood")]="Wood"
C14Dates$Material[which(C14Dates$Material=="Lacquer")]="Lacquer"
C14Dates$Material[which(C14Dates$Material=="Plant")]="Organic"
C14Dates$Material[which(C14Dates$Material=="Carbonised Seed")]="Seed/Nut"
C14Dates$Material[which(C14Dates$Material=="Organic")]="Organic"
C14Dates$Material[which(C14Dates$Material=="Carbonised Wood")]="Wood"
C14Dates$Material[which(C14Dates$Material=="Charcoal")]="Charcoal"
C14Dates$Material[which(C14Dates$Material=="Others")]="Others"
C14Dates$Material[which(C14Dates$Material=="Nutshell")]="Seed/Nut"


##Remove samples without Labcode
C14Dates=C14Dates[-which(is.na(C14Dates$LabCode)),]
##C14Dates[which(duplicated(C14Dates$LabCode)),]

##Identifying Instances of Duplicated Labcodes
anyDuplicated(C14Dates$LabCode)
which(duplicated(C14Dates$LabCode))


##Sites with no locations:
site2remove=Sites[which(is.na(Sites$Latitude)),]$SiteID
#   SiteID  SiteName Latitude Longitude Prefecture
#20     29  Kenyoshi       NA        NA     Aomori
#80    161 Shimonone       NA        NA   Kanagawa

## Remove these Sites and C14Dates linked to them
C14Dates=subset(C14Dates,!SiteID%in%site2remove)
Sites=subset(Sites,!SiteID%in%site2remove)


write.csv(C14Dates[,-11],"~/github/jomonSPD/data/c14dates.csv",row.names=FALSE)
write.csv(Sites,"~/github/jomonSPD/data/sites.csv",row.names=FALSE)
