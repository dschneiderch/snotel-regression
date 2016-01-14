setwd('~/Documents/snotel-regression_project')
library('ProjectTemplate')
load.project()

fordensource='umd_forden'

## ----- calculate daily average snotel scaled by fsca for relative calcs below.
    snotellocs.sub=snotellocs[snotellocs$Station_ID %in% snotelrec$Station_ID,]
    snotellocs.sub=snotellocs.sub[order(snotellocs.sub$Station_ID),]

## ----- need to get simulated ddates from the crossvalidation files.
recon.version='v3.1'
covrange='idp1'
snotelscale='scale'
fscaMatch='wofsca'
cost='r2'
style='real-time'
#
fnpath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/')
fnpattern=paste0(style,'_recondate_selection_',fscaMatch,'_',cost)
fns=list.files(path=fnpath,pattern=fnpattern,full.names=T)
skill=ldply(fns,read.table,sep='\t',header=T,stringsAsFactors=F)
skill$date=as.POSIXct(skill$date,tz='MST')
skill$phvrcn_recondate=as.POSIXct(skill$phvrcn_recondate,tz='MST')
skill$recon_costdate=as.POSIXct(skill$recon_costdate,tz='MST')
#
swe_doi=snotelrec[snotelrec$date %in% skill$date,c('date','swe','Station_ID')]
colnames(swe_doi)=c('date','snotel','Station_ID')
 
dateind=which(strftime(swe_doi$date,'%m')<3)

# mdate=unique(swe_doi[,'date'])[1]
fsca_doi=ldply(unique(swe_doi[,'date']),function(mdate) {
    print(mdate)
 if(as.numeric(strftime(mdate,'%m'))<3) {
         sca=get_sca(mdate,'fsca')
         scaraster=raster(matrix(sca,nrow=2580,byrow=T),xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
      projection(scaraster)='+proj=longlat +datum=WGS84'
      snotelsca=extract(scaraster, snotellocs.sub,cellnum=F)#get values at snotel pixel locations. named sca so we can use its length later
      snotelfsca=data.frame(date=mdate,fsca=snotelsca,Station_ID=snotellocs.sub$Station_ID)
} else {
yr=as.numeric(strftime(mdate,'%Y'))
rdfn=paste0('snotelrecon',yr)
if(!exists(rdfn)){
    load(paste0('data/recon_',recon.version,'/snotelfscadata_',yr,'_',recon.version,'.RData'))
    snotelfsca=get(rdfn)
    snotelfsca$date=as.POSIXct(snotelfsca$date,tz='MST')
    colnames(snotelfsca)[ncol(snotelfsca)]='fsca'
    snotelfsca=snotelfsca[snotelfsca$date==mdate,c('date','fsca','Station_ID')]
}
}
return(snotelfsca)
})

fsca_doi$yr=as.numeric(strftime(fsca_doi$date,'%Y'))
doi2=ddply(fsca_doi,.(yr),function(dF){
    yr=unique(dF$yr)

if(!exists('forden'))   forden=stack(paste0('data/forestdensity/',fordensource,'.nc'))
      if(fordensource=='umd_forden'){        
          layerind=yr-2000+1
        } else { 
          layerind=1
      }
snotelforden=extract(forden[[layerind]],snotellocs.sub)
snotelforden[snotelforden==1]=0.999# avoid divide by zero     
#--
dF[dF$fsca==0,'fsca']=.1
dF[,'fsca']=dF[,'fsca']/(1-snotelforden)#if 1-forden is > sca from satellite, then resultant sca is >1. different products don't always match
dF[dF$fsca>1,'fsca']=1##if clouded, assume 100%sca
return(dF)
})

data_oi=merge(doi2,swe_doi)
data_oi=arrange(data_oi,date,Station_ID)
data_oi$swe=data_oi$snotel*data_oi$fsca
saveRDS(data_oi,'data/station_swe_scaled.rds')

dailyavg_scaled=ddply(data_oi,.(date),.inform=T,function(dF){
    summarise(dF, avgsnotel=mean(swe))
})
# write.table(x=dailyavg,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_dailysnotelswe.txt'),quote=F,row.names=F)
saveRDS(dailyavg_scaled,'data/dailyavg_snotel_scaled.rds')