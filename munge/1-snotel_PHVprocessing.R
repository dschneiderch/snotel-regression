library('ProjectTemplate')
library(R.matlab)
library(matlab)
library(spatial.tools)
setwd('~/GoogleDrive/snotel-regression_project')
#load.project()

# Read snotel, phv data from matlab output.
data=readMat('data/snotel_phv_matrices.mat')
stationid=readMat('data/StationID.mat')
stationid=unlist(unlist(stationid))
stationid=gsub(' ','',stationid)
#
phvrec=data[[1]]
#Assign appropriate column names for phv
colnames(phvrec)=c('Lat','Long','Elev','Eastness','Northness','Slope','RegionalSlope','RegionalEastness','RegionalNorthness','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')
phvrec=as.data.frame(phvrec)
phvrec$Station_ID=stationid
### --- Subset phvrec based on stations found in snotellocs
# snotellocs was created selecting snotel station inside the recon domain in geographic coords and then further subset based on snotel file existence and span of record. phvrec$Station_ID is from matlab code that subsets stations in geographic coords using ~100km buffer around recon domain and then snotel file existence and span of record. they currently overlap with 237 stations.
ind=which(phvrec$Station_ID %in% snotellocs$Station_ID)
phvrec=phvrec[ind,]


### Get UCO domain terrain variables
### LOAD UCO TERRAIN VARIABLES
uco_phv=readMat('data/uco_variables_MASTER.mat')
names(uco_phv)=c('Lat','Long','RegionalSlope','RegionalAspect','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','Waz','NWd2ocean','NWaz','SWd2ocean','SWaz','Aspect','Slope','Elev')
### Stack each variable
Lat=raster(uco_phv$Lat,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Long=raster(uco_phv$Long,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
RegionalSlope=raster(uco_phv$RegionalSlope,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
RegionalEastness=raster(uco_phv$RegionalAspect,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
RegionalNorthness=sin(RegionalSlope)*RegionalEastness
FtprtW=raster(uco_phv$FtprtW,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Wbdiff=raster(uco_phv$Wbdiff,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
NWbdiff=raster(uco_phv$NWbdiff,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
SWbdiff=raster(uco_phv$SWbdif,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Wd2ocean=raster(uco_phv$Wd2ocean,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
NWd2ocean=raster(uco_phv$NWd2ocean,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
SWd2ocean=raster(uco_phv$SWd2ocean,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Eastness=raster(uco_phv$Aspect,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Slope=raster(uco_phv$Slope,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Northness=Eastness*sin(Slope)
Elev=raster(uco_phv$Elev,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
ucophv.stack=stack(list(Lat,Long,RegionalSlope,RegionalEastness,RegionalNorthness,FtprtW,Wbdiff,NWbdiff,SWbdiff,Wd2ocean,NWd2ocean,SWd2ocean,Eastness,Northness,Slope,Elev))
names(ucophv.stack)=c('Lat','Long','RegionalSlope','RegionalEastness','RegionalNorthness','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean','Eastness','Northness','Slope','Elev')
projection(ucophv.stack) <- '+proj=longlat +datum=NAD83'

#ucophv.stack=projectRaster(ucophv.stack,crs=CRS('+init=epsg:5070'))

phvfn='data/ucophv_variables_stack.grd'
writeRaster(ucophv.stack,filename=phvfn,overwrite=T)
ucophv.stack=stack(phvfn)
#### Create dataframe of stack for modeling later
ucophv=as.data.frame(getValues(ucophv.stack))

#### --- with PHV variables of enlarged domain, could get snotel pixel attributes just by extracting from stack. currentyl PHV only extends for RCN domain. but no recon outside recon domain
#phvrec=extract(ucophv.stack,snotellocs,df=T)

save(file='cache/ucophv.RData',list=c('phvrec','ucophv','ucophv.stack'))


