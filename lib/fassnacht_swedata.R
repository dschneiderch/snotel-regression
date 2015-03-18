fassnacht_swedata=function(snotelrecon,qc,precipmeth){
          
snotel1993=read.table(paste0('data/fassnacht/snotelrec_1993_',qc,'_precipmeth',precipmeth,'.txt'),header=T)
snotel1998=read.table(paste0('data/fassnacht/snotelrec_1998_',qc,'_precipmeth',precipmeth,'.txt'),header=T)
snotel1999=read.table(paste0('data/fassnacht/snotelrec_1999_',qc,'_precipmeth',precipmeth,'.txt'),header=T)

snotelrec=rbind(snotel1993,snotel1998,snotel1999)
snotelrec$date=as.POSIXct(snotelrec$date,tz='MST')
# create temporally continuous snotelrecon record ----------------
swedata=merge(snotelrecon[,c('date','Station_ID','recon')],snotelrec[,c('swe','date','Station_ID')],all=T)
swedata$yrdoy=as.numeric(strftime(swedata$date,'%Y%j'))
swedata$dy=as.numeric(strftime(swedata$date,'%d'))
swedata$mth=as.numeric(strftime(swedata$date,'%m'))
swedata$yr=as.numeric(strftime(swedata$date,'%Y'))
swedata$mnth=factor(strftime(swedata$date,'%b'),levels=unique(strftime(swedata$date[order(swedata$date)],'%b')))
#
# scale the x variables -----------------------------------------
phv_scatt=scale(phvrec[,-ncol(phvrec)])
phvrec.sc=data.frame(phv_scatt,Station_ID=phvrec$Station_ID)
phvrec.sc=phvrec.sc[order(phvrec$Station_ID),]
#
# merge phv and snotel data
swedata=merge(swedata,phvrec.sc,by='Station_ID',all=T)
colnames(swedata) <- gsub(pattern='swe',replacement='snotel',x=colnames(swedata))
#

# ## get snotel coordinates for this dataset. even if all years don't have the same stations, there is more subsetting in doDOYfit 
# snotellist=read.csv('data/SNOTEL_MASTER.csv',sep=',',header=T,stringsAsFactors=F,strip.white=T)
# coordinates(snotellist) <- ~Longitude+Latitude
# crs(snotellist) <- CRS('+proj=longlat +datum=WGS84'))
# #snotellist=spTransform(snotellist,CRS('+init=epsg:5070'))
# ind=which(snotellocs$Station_ID %in% unique(swedata$Station_ID) )
# snotellocs=snotellocs[ind,]

# return(list(swedata,snotellocs))

return(swedata)
}
