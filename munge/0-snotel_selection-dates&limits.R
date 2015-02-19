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
     #adjust Temp by 1 day. NRCS reports Temp from day i-1 as day i
     x$tmin=c(NA,x$tmin[2:nrow(x)])
     x$tmax=c(NA,x$tmax[2:nrow(x)])
     x$tavg=c(NA,x$tavg[2:nrow(x)])
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
#
snotelrec=snotelrec[order(snotelrec$Station_ID,snotelrec$date),]
snotelrec$yr=as.numeric(strftime(snotelrec$date,'%Y'))
snotelrec$wy=water_yr(snotelrec$date)
snotelrec$mth=as.numeric(strftime(snotelrec$date,'%m'))
snotelrec$dy=as.numeric(strftime(snotelrec$date,'%d'))

tmp=snotelrec

snotelrec=tmp
## output snotel files we are interested in with qa and qc
snotelrec=ddply(snotelrec,.inform=F,.(Station_ID),function(dF){
     #
     #save raw snotelrec for plotting later
     dFraw=dF
     dFraw$version='raw'
     
     ##---------- temperature
     #QA temperature
     tempind=with(dF,
                  which(tmin < -40 | tmax < -40 | tavg < -20 | tmin > 40 | tmax > 40 | tavg > 20))
     dF$tmin[tempind]=NA
     dF$tmax[tempind]=NA
     dF$tavg[tempind]=NA
     #QC temperature
     dF$tmin=fillMissing(dF$tmin,max.fill=14) #I probably don't need a complete series
     dF$tmax=fillMissing(dF$tmax,max.fill=14)
     dF$tavg=fillMissing(dF$tavg,max.fill=14)
     
     ## ----------- swe and precip
     #set summer SWE NAs to 0
     summerNA=which(dF$mth>=7 & dF$mth<11 & is.na(dF[,'swe']))
     dF[summerNA,'swe']=0
     
     #check for and convert neg swe to 0
     negs=which(dF$swe<0)
     if(length(negs)>0)  dF[negs,'swe']=0
     
     #make continuous swe timeseries 
     dFqc1=ddply(dF,.(wy),fixSWE,.inform=T)#does some minor precip QA, then uses precip to fill swe. should never see "end of ts not filled if summerNA sets Sept NAs to 0
     
     #adjust precip to fix some of the undercatch compared to the pillow.
     precipmeth=2
     dFqc1=precip_adjust(dFqc1,precipmeth)#script transcribed from J.Barsugli's matlab script. adjusts precip for undercatch based on pillow. several methods available
     dFqc1$version='qc1'
     #write new record to file
     dirname=paste0('data/snotelQAQC')
     fn=file.path(dirname,'qc',paste0(unique(dFqc1$Station_ID),'-qc1_precipmeth',precipmeth,'.txt'))
     write.table(dFqc1[,c('snoteldate','swe','apcp','tmax','tmin','tavg','precip')],fn,sep='\t',row.names=F,quote=F)
          
     ##CAVEAT: precip and aprecip are used to fill SWE. precip_adjust then increases apcp and precip based on max(dswe,precip). since there are no NA's in swe, rerunning fixSWE on dFqc1 wouldn't do anything. rerun fixSWE with raw swe but updated apcp and precip? -- there are so few days where SWE =NA. plus with precipmeth=1 swe needs to increase for precip to be adjusted so it might not change swe much to rerun fixSWE...? maybe more with precipmeth=2
     dFqc2=dF
     dFqc2$precip=dFqc1$precip
     dFqc2$apcp=dFqc1$apcp
     dFqc2=ddply(dFqc2,.(wy),fixSWE,.inform=T)
     dFqc2$version='qc2'
     #write new record to file
     dirname=paste0('data/snotelQAQC')
     fn=file.path(dirname,'qc',paste0(unique(dFqc2$Station_ID),'-qc2_precipmeth',precipmeth,'.txt'))
     write.table(dFqc1[,c('snoteldate','swe','apcp','tmax','tmin','tavg','precip')],fn,sep='\t',row.names=F,quote=F)
     
     #combine raw and qc snotelrec and plot together.
     dFall=rbind(dFraw,dFqc1,dFqc2)
     dFplot=melt(dFall[,c('date','swe','apcp','version')],.(date,version))
     g=ggplot(dFplot)+
          geom_line(aes(x=date,y=value,colour=version,linetype=variable),alpha=0.5)+
          labs(title=paste(unique(dF$Station_ID),'\nprecip method',precipmeth))+
          theme_minimal()
     ggsave(plot=g,filename=paste0(dirname,'/plots/snotelrec_',unique(dF$Station_ID),'_precipmeth',precipmeth,'.pdf'),height=3,width=6)
     ## ggplot is giving a warning that 23 rows contain missing values but snotelrec doesn't have any missing values so something in dFplot after melt?
     return(dFqc2)
})

#finally subset snotellocs based on snotel that had records for our dates
ind=which(snotellocs$Station_ID %in% unique(snotelrec$Station_ID) )
snotellocs=snotellocs[ind,]

save(list=c('snotelrec', 'snotellocs','current_wy','cutoff_date','usgsbox','geobox'),file='cache/snoteldata.RData')
