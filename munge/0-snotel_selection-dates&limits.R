library('ProjectTemplate')
setwd('~/Documents/R/snotel-regression_project')
load.project()


## data must have historical/<STATE>
## set dates
cutoff_date=as.POSIXct('1999-10-01',tz='MST')
current_wy=as.POSIXct('2012-10-01',tz='MST')
## snotel format date
cutoff_recorddate=strftime(cutoff_date,'%m%d%y')
current_wateryear= strftime(current_wy,'%m%d%y')

## geographic extents
Wlimit=-112.25;
Elimit=-104.125;
Nlimit=43.75;
Slimit=33;
boxpoly=Polygon(matrix(c(Wlimit,Elimit,Elimit,Wlimit,Wlimit,Slimit,Slimit,Nlimit,Nlimit,Slimit),nrow=5),hole=F)
boxpolylist=Polygons(list(boxpoly),'uppercolorado')
geobox=SpatialPolygons(list(boxpolylist),proj4string=CRS('+proj=longlat +datum=NAD83'))
#usgsbox=spTransform(geobox,CRS('+init=epsg:5070'))

#increase box by 100km -- gaah don't bother because you don't have recon data outside
#usgsbox_over=gBuffer(usgsbox,width=100000,capStyle='SQUARE')
#geobox_over=gBuffer(geobox,width=0.00416666667/500*100000)

## get snotel coordinates
snotellist=read.csv('data/SNOTEL_MASTER.csv',sep=',',header=T,stringsAsFactors=F,strip.white=T)
coordinates(snotellist) <- ~Longitude+Latitude
crs(snotellist) <- crs(geobox)
#snotellist=spTransform(snotellist,CRS('+init=epsg:5070'))

## compare snotel locations to uppercolorado polygon
snotel2keep=which(!is.na(over(snotellist,geobox)))
snotellocs=snotellist[snotel2keep,]
#snotellocs=spTransform(snotellocs,CRS('+init=epsg:5070'))

## get data from snotels in snotellocs
states=c('arizona','colorado','new_mexico','wyoming','idaho','utah')
filedir=file.path('data/historical',states)
fn=list.files(filedir,full.names=T)
temp1=matrix(unlist(strsplit(fn,'/',fixed=T)),ncol=4,byrow=T)
temp2=matrix(unlist(strsplit(temp1[,4],'_',fixed=T)),ncol=2,byrow=T)
fn=fn[order(temp2[,1])]
fnshort=list.files(filedir)
#subset snotel files based on geographic location from snotellocs
stationid=laply(strsplit(fnshort,'_',fixed=T),'[',1)
fn_ind=which(toupper(stationid) %in% snotellocs$Station_ID)
stalist=fn[fn_ind]#
stationid=stationid[fn_ind]
dat=llply(stalist,read.table,sep='\t',fill=T,strip.white=T)
#
##**** these have incompatible first lines (usually a # in name)
##FUTURE WORK: use textscan to parse first line
## [1] "data/historical/colorado/05m03s_all.txt"
## [1] "data/historical/colorado/07m35s_all.txt"
## [1] "data/historical/new_mexico/05n11s_all.txt"
## [1] "data/historical/new_mexico/06p10s_all.txt"
## [1] "data/historical/utah/09j01s_all.txt"
## [1] "data/historical/utah/10j10s_all.txt"
## [1] "data/historical/utah/10j12s_all.txt"
## [1] "data/historical/utah/10k02s_all.txt"
## [1] "data/historical/utah/11j01s_all.txt"
## [1] "data/historical/utah/11j02s_all.txt"
## [1] "data/historical/utah/11j52s_all.txt"
## [1] "data/historical/utah/11k15s_all.txt"
## [1] "data/historical/utah/11k21s_all.txt"
## [1] "data/historical/utah/11k22s_all.txt"
## [1] "data/historical/utah/11m03s_all.txt"

# *** make sure the snotel record has swe. remove any radiation measurements
dat2=list()
stdlabels=c('pill','prec','tmax','tmin','tavg','prcp')
j=1
for(i in seq(1,length(dat))){
    x=dat[[i]]
    temp=x[1,]
    labind=which(temp %in% stdlabels)
    if(length(labind)==6){
        x=x[-1,c(1,labind)]
        colnames(x) <- c('snoteldate','swe','apcp','tmax','tmin','tavg','precip')
        x=as.data.frame(lapply(x,as.numeric))
        dat2[[j]]=x
        j=j+1
    } else {
        print('WARNING: check the file')
        print(stalist[i])
        print(x[1,])
    }
}
#

# ** function to check for full record of interest
check_snotelrec <- function(x){
# check dates of record
    indstart=which(x$snoteldate==as.numeric(cutoff_recorddate))
    indend=which(x$snoteldate==as.numeric(current_wateryear)-7100)#last day of previous water year
    if(length(indstart)==1 & length(indend)==1){
        x=x[indstart:indend,]
        return(x)
    }
}


# ** create dataframe of snotel stations only with full temporal record
snotelrec=data.frame(matrix(NA,ncol=8))
colnames(snotelrec) <- c('snoteldate','swe','apcp','tmax','tmin','tavg','precip','Station_ID')
for(i in seq(1,length(dat2))){
    df=check_snotelrec(dat2[[i]])
    if(is.data.frame(df)){
        df$Station_ID=toupper(stationid[i])
        snotelrec=rbind(snotelrec,df)
    }
}
snotelrec=snotelrec[-1,]
snotelrec$swe=snotelrec$swe*2.54/100#convert to meters
snotelrec$apcp=snotelrec$apcp*0.0254
snotelrec$precip=snotelrec$precip*0.0254
snotelrec$date=as.POSIXct(strptime(sprintf(fmt='%06d',snotelrec$snoteldate),format='%m%d%y',tz='MST'))

#finally subset snotellocs based on snotel that had records for our dates
ind=which(snotellocs$Station_ID %in% unique(snotelrec$Station_ID) )
snotellocs=snotellocs[ind,]


save(list=c('snotelrec', 'snotellocs','current_wy','cutoff_date','usgsbox','geobox'),file='cache/snoteldata.RData')
