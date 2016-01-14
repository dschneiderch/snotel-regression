library('ProjectTemplate')
setwd('/Volumes/Dominik/Documents/snotel-regression_project')
load.project()
library(gtable)
library(xtable)
library(xlsx)

cost='r2'
style='real-time'
recon.version='v3.1'
covrange='idp1'
snotelscale='scale'#'noscale','scale'
fscaMatch='wofsca'#'wofsca', 'fsca'
fordensource='umd_forden'
dateflag='B'
#get daily average snotel swe that was scaled by fsca 
dailyavg_scaled=readRDS('data/dailyavg_snotel_scaled.rds')

### for comparison to fassnacht
davg=readRDS('data/station_swe_scaled.rds')
daily=dcast(davg,date+yr~.,value.var='snotel',fun.agg=mean)
colnames(daily)=c('date','yr','domainavg')
dcast(daily,yr~.,value.var='domainavg',fun.agg=max)

qc='qc1'
precipmeth='1'
snotel1993=read.table(paste0('data/fassnacht/snotelrec_1993_',qc,'_precipmeth',precipmeth,'.txt'),header=T)
snotel1998=read.table(paste0('data/fassnacht/snotelrec_1998_',qc,'_precipmeth',precipmeth,'.txt'),header=T)
snotel1999=read.table(paste0('data/fassnacht/snotelrec_1999_',qc,'_precipmeth',precipmeth,'.txt'),header=T)
snotelrec=rbind(snotel1993,snotel1998,snotel1999)
snotelrec$date=as.POSIXct(snotelrec$date,tz='MST')
swedata=snotelrec
swedata$yr=as.numeric(strftime(swedata$date,'%Y'))
str(swedata)
dailyf=dcast(swedata[,c('swe','date','yr')],date+yr~.,value.var='swe',fun.agg=mean)
colnames(dailyf)=c('date','yr','domainavg')
dcast(dailyf,yr~.,value.var='domainavg',fun.agg=max)

#### ------


dFall=data.frame()
for(cost in c('rmse','r2')){
    for(style in 'real-time'){
        for(covrange in c('idp1')){
            for(snotelscale in c('scale','noscale')){
                for(fscaMatch in c('wofsca','fsca')){

print(cost)
print(style)
print(covrange)
print(snotelscale)
print(fscaMatch)
print('---------')
# graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',fscaMatch)

#skill=read.table(paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/fassnacht/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T)
fnpath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag)
fnpattern=paste0(style,'_recondate_selection_',fscaMatch,'_',cost)
fns=list.files(path=fnpath,pattern=fnpattern,full.names=T)
skill=ldply(fns,read.table,sep='\t',header=T,stringsAsFactors=F)
skill$date=as.POSIXct(skill$date,tz='MST')
skill$phvrcn_recondate=as.POSIXct(skill$phvrcn_recondate,tz='MST')
skill$recon_costdate=as.POSIXct(skill$recon_costdate,tz='MST')
##
skill$fulldiff=skill$skill_phvrcnfull-skill$skill_phvfull
skill$diff=skill$skill_phvrcn-skill$skill_phv
skill.m=melt(skill[,c('date','yr','skill_phv','skill_phvrcn','skill_phvfull','skill_phvrcnfull','fulldiff','diff','yrdoy')],.(date,yr,yrdoy))
skill.m$set=ifelse(grepl('skill',as.character(skill.m$variable)),'data','diff')
skill.m$set=as.factor(skill.m$set)
skill.m$modeltype=ifelse(grepl('phvrcn',as.character(skill.m$variable)),'phvrcn','phv')
skill.m$modeltype=ifelse(grepl('diff',as.character(skill.m$variable)),'diff',skill.m$modeltype)
skill.m$modeltype=as.factor(skill.m$modeltype)
skill.m$mnthdy=strftime(skill.m$date,'%b%d')
skill.m$mnth=factor(strftime(skill.m$date,'%b'),levels=c('Jan','Feb','Mar','Apr','May','Jun'))
skill.m$doy=as.numeric(strftime(skill.m$date,'%j'))


## --- set up labels based on skill metric and compute relative rmse if applicable
if(cost=='rmse'){
    
    costshort='rmse'
    costlong=expression('RMSE')
    diffcostlong=expression(Delta * RMSE)

if(snotelscale=='noscale'){
dailyavg=ddply(snotelrec[snotelrec$date %in% skill$date,],.(date),function(x){
    summarise(x,
        avgsnotel=mean(swe)
        )})
} else if(snotelscale=='scale'){
    dailyavg=dailyavg_scaled
}

    dind=dailyavg$date %in% skill.m$date    
    skill.m=cbind(skill.m,avgsnotel=dailyavg$avgsnotel[dind])
    skill.m$Rrmse=skill.m$value/skill.m$avgsnotel
    summrel=ddply(skill.m,.(variable,yr),summarise,
        avg=mean(Rrmse,na.rm=T),
        med=median(Rrmse,na.rm=T),
        min=min(Rrmse,na.rm=T),
        max=max(Rrmse,na.rm=T))
# write.table(x=summrel,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',costshort,'.txt'),quote=F,row.names=F)

} else if(cost=='r2'){
    # skill.m=subset(skill.m,value<0.5)
    costshort='r2'
    costlong=expression(r^2)
    diffcostlong=expression(Delta * r^2)   
}

## --------------
### r2 or rmse
summ=ddply(skill.m,.(variable,yr), function(x){
    summarise(x,
    avg=mean(value,na.rm=T),
    med=median(value,na.rm=T),
    min=min(value,na.rm=T),
    max=max(value,na.rm=T),
    sd=sd(value,na.rm=T),
    cv=sd/avg)
})
# write.table(x=summ,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',cost,'.txt'),quote=F,row.names=F)
dcast(summ,yr~variable,value.var='avg')

### relative rmse
if(cost=='rmse'){
summ=ddply(skill.m,.(variable,yr), function(x){
    summarise(x,
    avg=mean(Rrmse,na.rm=T),
    med=median(Rrmse,na.rm=T),
    min=min(Rrmse,na.rm=T),
    max=max(Rrmse,na.rm=T),
    sd=sd(Rrmse,na.rm=T),
    cv=sd/avg)
})
tmp=dcast(summ,yr~variable,value.var='avg')
}
mean(tmp$diff)

###----


dFallstat=data.frame()
## stat refers to summarise columns above
for(stat in c('min','max','med','avg')){
    dF=dcast(summ,variable~yr,value.var=stat)
    dF$stat=stat
    dF$cost=cost
    dF$covrange=covrange
    dF$fscaMatch=fscaMatch
    dF$snotelscale=snotelscale
    if(cost=='rmse'){
        dFrel=dcast(summrel,variable~yr,value.var=stat)
        dFrel$stat=stat
        dFrel$cost='Relrmse'
        dFrel$covrange=covrange
        dFrel$fscaMatch=fscaMatch
        dFrel$snotelscale=snotelscale        
        dF=rbind(dF,dFrel)
    }
    colnames(dF)[1]='model'
    dFallstat=rbind(dFallstat,dF)
#
    dFmdl=subset(dFallstat,!grepl('diff',as.character(model)))
    dFm=melt(dFmdl,.(cost,covrange,fscaMatch,snotelscale,stat,model))
#
}
dFall=rbind(dFall,dFm)
                }
            }
        }
    }
}

str(dFall)

#### Compare to fassnacht. 
#1993 bigyear, rmse=155
#1998 avgyr, rmse=126
#1999 lowyr, rmse=124
tmp=subset(dFall,variable==2010 | variable==2011 | variable==2012)
tmp2=subset(tmp,stat=='avg')
tmp3=subset(tmp2,cost=='rmse')
format(dcast(tmp3,variable+model~fscaMatch+snotelscale,value.var='value'),digits=2)



## -- median values
# this is indicative of overall performance to determine with parameterization is best
fn='reports/xval_reports/median.xlsx'
dFm.stat=subset(dFall,stat=='med')
dFmed=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale,model),function(x){
    summarise(x,
            minval=min(value),
            maxval=max(value),
            avgval=mean(value))
})
format(dFmed,digits=1,scientific=F)

write.xlsx(dFmed,sheetName='all_values',row.names=F,file=fn)
wb=loadWorkbook(file=fn)

dFmed.mdl=ddply(dFmed,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
        minmodel=model[which.min(avgval)],#lowest average median worst model config, highest average median best model config
        maxmodel=model[which.max(avgval)])
})
dFmed.mdl

dFmed.vals=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
            # minval=min(value),#gives min median value for all the years.
            # maxval=max(value),
            avgval=mean(value))
})
format(dFmed.vals,digits=1,scientific=F)
dFmed.all=merge(dFmed.mdl,dFmed.vals,by=c('cost','covrange','fscaMatch','snotelscale'))

cs=createSheet(wb,sheetName='best_values')
addDataFrame(dFmed.all,sheet=cs,row.names=F)
saveWorkbook(wb,fn)



## -- average values
# this is indicative of overall performance to determine with parameterization is best
fn='reports/xval_reports/mean.xlsx'
dFm.stat=subset(dFall,stat=='avg')
dFmed=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale,model),function(x){
    summarise(x,
            minval=min(value),
            maxval=max(value),
            avgval=mean(value))
})
format(dFmed,digits=1,scientific=F)

write.xlsx(dFmed,sheetName='all_values',row.names=F,file=fn)
wb=loadWorkbook(file=fn)

dFmed.mdl=ddply(dFmed,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
        minmodel=model[which.min(avgval)],#lowest average median worst model config, highest average median best model config
        maxmodel=model[which.max(avgval)])
})
dFmed.mdl

dFmed.vals=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
            # minval=min(value),#gives min median value for all the years.
            # maxval=max(value),
            avgval=mean(value))
})
format(dFmed.vals,digits=1,scientific=F)
dFmed.all=merge(dFmed.mdl,dFmed.vals,by=c('cost','covrange','fscaMatch','snotelscale'))

cs=createSheet(wb,sheetName='best_values')
addDataFrame(dFmed.all,sheet=cs,row.names=F)
saveWorkbook(wb,fn)



## -- min values
#only care about min of min
fn='reports/xval_reports/min.xlsx'
dFm.stat=subset(dFall,stat=='min')
dFmin=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale,model),function(x){
    summarise(x,
            val=min(value))
})
format(dFmin,digits=1,scientific=F)
write.xlsx(dFmin,sheetName='all_values',row.names=F,file=fn)
wb=loadWorkbook(file=fn)

dFmin.mdl=ddply(dFmin,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
        minmodel=model[which.min(val)])
})
dFmin.mdl

dFmin.vals=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
        minval=min(value))
})
format(dFmin.vals,digits=1,scientific=F)
dFmin.all=merge(dFmin.mdl,dFmin.vals,by=c('cost','covrange','fscaMatch','snotelscale'))

cs=createSheet(wb,sheetName='best_values')
addDataFrame(dFmin.all,sheet=cs,row.names=F)
saveWorkbook(wb,fn)


##initialiize max values
fn='reports/xval_reports/max.xlsx'
dFm.stat=subset(dFall,stat=='max')
dFmax=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale,model),function(x){
    summarise(x,
            val=max(value))
})
format(dFmax,digits=1,scientific=F)

write.xlsx(dFmax,sheetName='all_values',row.names=F,file=fn)
wb=loadWorkbook(file=fn)

dFmax.mdl=ddply(dFmax,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
        maxmodel=model[which.max(val)])
})
dFmax.mdl

dFmax.vals=ddply(dFm.stat,.(cost,covrange,fscaMatch,snotelscale),function(x){
    summarise(x,
        maxval=max(value))
})
dFmax.vals
dFmax.all=merge(dFmax.mdl,dFmax.vals,by=c('cost','covrange','fscaMatch','snotelscale'))

cs=createSheet(wb,sheetName='best_values')
addDataFrame(dFmax.all,sheet=cs,row.names=F)
saveWorkbook(wb,fn)



# dte=as.POSIXct('2011-05-01',tz='MST')
# startday=as.POSIXct('2011-03-01',tz='MST')
# startdoy=as.numeric(strftime(startday,'%j'))
# seq(dte,dte+7*24*3600, 3600*24)
# datpath=paste0('/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/recon/v4.0_senlat.nldas.usace.modscag.umdforden.fscamasking/snow_output/2011/')
# system(paste0('mkdir -p ',datpath,'/swe_tif'))
# var='fsca'
# if(var=='swe') fid=file('/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/recon/v4.0_senlat.nldas.usace.modscag.umdforden.fscamasking/snow_output/2011/swe.dat','rb')
# if(var=='fsca') fid=file('/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/recon/v4.0_senlat.nldas.usace.modscag.umdforden.fscamasking/snow_output/2011/fsca.dat','rb')
# for(day in seq(dte,dte+7*24*3600, 3600*24)){
#     # day=1301641200
#     day= as.POSIXct(day,origin=as.POSIXct('1970-01-01 00:00:00'),tz='MST')
#     doy=as.numeric(strftime(day,'%j'))-startdoy+1
#     swe=readgrads(doy,fid)    
#     writeRaster(swe,paste0(datpath,var,'_tif/',strftime(day,'%d%b%Y'),'.tif'),overwrite=T)
# }

# library(rasterVis)
# fid=file('/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/recon/v4.0_senlat.nldas.usace.modscag.umdforden.fscamasking/snow_output/2011/swe.dat','rb')
# swe=readgrads(71,fid)
# swe2=sampleRegular(swe,2e6,asRaster=T)
# swe.df=as.data.frame(swe2,xy=T)
# swe.df$swecut=cut(swe.df$layer,breaks=c(0,.1,.5,1,2,12))
# dev.set(2)
# ggplot(swe.df)+geom_raster(aes(x,y,fill=swecut))+scale_color_brewer(palette='PuBu')
# fid2=file('/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/recon/v4.0_senlat.nldas.usace.modscag.umdforden.fscamasking/snow_output/2011/fsca.dat','rb')
# fsca=readgrads(70,fid2)
# fsca2=sampleRegular(fsca,2e6,asRaster=T)
# fsca.df=as.data.frame(fsca2,xy=T)
# fsca.df$fscacut=cut(fsca.df$layer,breaks=seq(0,1,.1))
# dev.set(3)
# ggplot(fsca.df)+geom_raster(aes(x,y,fill=fscacut))+scale_color_brewer(palette='BuYlGn')