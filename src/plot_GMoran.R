library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(dplyr)

cost='rmse'
costshort='globalI'
style='real-time'#reanalysis'
recon.version='v3.2'
covrange='300km'

## **
# you'll need to run src/check_spatialcorrelation.R first
## **
skill=read.table(paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/',style,'_moran_info_for_recondate_selection_',cost,'.txt'),sep='\t',header=T)
skill$date=as.POSIXct(skill$date,tz='MST')
skill$yr=strftime(skill$date,'%Y')
skill$recondate=as.POSIXct(skill$recondate,tz='MST')
# str(skill)
##

# set plot title info about GMoran calc
plttitle2=paste0('\nIDW p=0.5, binary neighborhood') # from calc_GMoran.R

#set alpha for pvalue cutoff
aval=0.05

colid=colnames(skill)[grep('globalI',colnames(skill))]
skill.I=melt(skill[,c('yrdoy','date','yr',colid)],.(date,yr,yrdoy))
colnames(skill.I)=gsub('variable','globalI',colnames(skill.I))
colnames(skill.I)=gsub('value','I',colnames(skill.I))
skill.I$model='phv'
skill.I$model=ifelse(grepl('phvrcn',skill.I$globalI),'phvrcn','phv')
skill.I$blend='noblend'
skill.I$blend=ifelse(grepl('full',skill.I$globalI),'blended','noblend')
#
colid=colnames(skill)[grep('pvalue',colnames(skill))]
skill.p=melt(skill[,c('yrdoy','date','yr',colid)],.(date,yr,yrdoy))
colnames(skill.p)=gsub('variable','pvalue',colnames(skill.p))
colnames(skill.p)=gsub('^value','p',colnames(skill.p))
skill.p$model='phv'
skill.p$model=ifelse(grepl('phvrcn',skill.p$pvalue),'phvrcn','phv')
skill.p$blend='noblend'
skill.p$blend=ifelse(grepl('full',skill.p$pvalue),'blended','noblend')
#
skill.m=merge(skill.I,skill.p,by=c('date','yr','yrdoy','model','blend'))
skill.m$sig=paste0('<',aval)
skill.m$sig=ifelse(skill.m$p<aval,skill.m$sig,paste0('>=',aval))
skill.m=melt(skill.m,.(date,yr,yrdoy,model,blend,sig,globalI),measure.var='I')

# str(skill.m)
summ_yr=ddply(skill.m,.(globalI,yr,model,blend,sig), function(x){
    summarise(x,
    avg=mean(value,na.rm=T),
    med=median(value,na.rm=T),
    min=min(value,na.rm=T),
    max=max(value,na.rm=T))
})
# str(summ_yr)
# write.table(x=summ,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',cost,'.txt'),quote=F,row.names=F)
#

indsig=which(summ_yr$sig==paste0('<',aval))
#calc difference in mean Moran I by model, by blend
ddply(.data=summ_yr[indsig,],.(blend,sig),summarise,avg=mean(avg))
ddply(.data=summ_yr[indsig,],.(model,sig),summarise,avg=median(avg))
ddply(.data=summ_yr[indsig,],.(blend,model,sig),summarise,avg=mean(avg))

#yearly averages
summ_avg=dcast(summ_yr[indsig,],blend+model+sig~yr,value.var='avg',mean,na.rm=T)
summ_avg

skill.m$mnthdy=strftime(skill.m$date,'%b%d')
skill.m$mnth=factor(strftime(skill.m$date,'%b'),levels=c('Jan','Feb','Mar','Apr','May'))
skill.m$doy=strftime(skill.m$date,'%j')

g=ggplot(skill.m)+
    geom_point(aes(x=date,y=value,colour=globalI,shape=sig))+
    scale_shape_manual(values=c(16,1))+
    labs(x='Date',y='Global Moran\'s I',title=paste0('Global Moran I Comparison\nPrediction Residuals',plttitle2))+
    theme_minimal()+
    theme(strip.text=element_text(face='bold',size=14),
        panel.background=element_rect(colour='grey80'))+
    facet_wrap(~yr,scale='free_x')
show(g)

ggsave(plot=g,filename=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/scatterplot_globalMoranI_resid_p<',aval,'_',cost,'_',style,'_byyear_blended&unblended.pdf'),width=12,height=6)