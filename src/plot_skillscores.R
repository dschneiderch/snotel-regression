library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(dplyr)

cost='rmse'
style='reanalysis'
recon.version='v3.2'
skill=read.table(paste0('diagnostics/rswe_',recon.version,'/fassnacht/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T)
# skill=read.table(paste0('diagnostics/rswe_',recon.version,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T)
skill$date=as.POSIXct(skill$date,tz='MST')
skill$phvrcn_recondate=as.POSIXct(skill$phvrcn_recondate,tz='MST')
skill$recon_costdate=as.POSIXct(skill$recon_costdate,tz='MST')
# str(skill)
##
skill.m=melt(skill[,c('date','yr','skill_phv','skill_phvrcn','skill_phvfull','skill_phvrcnfull','yrdoy')],.(date,yr,yrdoy))

# str(skill.m)
summ=ddply(skill.m,.(variable), function(x){
    summarise(x,
    avg=mean(value,na.rm=T),
    med=median(value,na.rm=T),
    min=min(value,na.rm=T),
    max=max(value,na.rm=T))
})
# write.table(x=summ,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',cost,'.txt'),quote=F,row.names=F)
summ
dF=dcast(summ,variable~yr,value.var='avg')
dF[3,]-dF[4,]

dailyavg=snotelrec[snotelrec$date %in% skill$date,] %>%
group_by(date) %>%
select(swe) %>%
summarise(  
    avgsnotel=mean(swe) )
# write.table(x=dailyavg,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_dailysnotelswe.txt'),quote=F,row.names=F)

if(cost=='rmse'){
	
    costshort='rmse'
    costlong='RMSE'

    skill.m=cbind(skill.m,avgsnotel=dailyavg$avgsnotel)
    skill.m$Rrmse=skill.m$value/skill.m$avgsnotel

    summrel=ddply(skill.m,.(variable,yr),summarise,
        avg=mean(Rrmse,na.rm=T),
    med=median(Rrmse,na.rm=T),
    min=min(Rrmse,na.rm=T),
    max=max(Rrmse,na.rm=T))
# write.table(x=summrel,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',costshort,'.txt'),quote=F,row.names=F)
 
} else if(cost=='r2'){
    costshort='r2'
    costlong=bquote(r^2)
    
}


skill.m$mnthdy=strftime(skill.m$date,'%b%d')
skill.m$mnth=factor(strftime(skill.m$date,'%b'),levels=c('Jan','Feb','Mar','Apr','May'))
skill.m$doy=strftime(skill.m$date,'%j')

#plot rmse and r2
if(cost=='rmse' | cost=='r2'){

ggplot(skill.m)+
    geom_boxplot(aes(x=as.factor(yr),y=value,colour=variable))+
    labs(x='Year',y=costlong)+
    scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
    theme_minimal()+
    theme(axis.line=element_line(colour='grey10'))
ggsave(filename=paste0('graphs/rswe_',recon.version,'/boxplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=10,height=4.5)

###

ggplot(skill.m)+
     geom_freqpoly(aes(x=value,colour=variable),alpha=0.75,size=2,binwidth=0.01)+
     labs(x=costlong)+
     scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
     theme_minimal()+
     theme(axis.line=element_line(colour='grey10'))
ggsave(filename=paste0('graphs/rswe_',recon.version,'/freqpoly_',costshort,'_',style,'_blended&unblended.pdf'),width=10,height=4.5)

###

ggplot(skill.m)+      
     geom_point(aes(x=as.numeric(doy),y=value,colour=variable),alpha=.75)+    
     facet_grid(.~yr)+
     labs(y=costlong)+
     scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
     scale_x_continuous('day of year',limits=c(0,160),labels=c(1,seq(30,160,30)),breaks=c(1,seq(30,160,30)))+
     #scale_y_continuous(limits=c(0,.5))+
     guides(colour=guide_legend(override.aes=list(alpha=1)))+
     theme_bw()+
     theme(axis.line=element_line(colour='grey10'),
           strip.background=element_rect(fill='white',colour='white'),
           strip.text=element_text(face='bold',size='12'),
           axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=.5))

ggsave(filename=paste0('graphs/rswe_',recon.version,'/scatterplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=11,height=3)
  }

#plot relative rmse
if(cost=='rmse'){

    costshort='RelRMSE'
    costlong='Relative RMSE'

ggplot(skill.m)+
    geom_boxplot(aes(x=as.factor(yr),y=Rrmse,colour=variable))+
    labs(x='Year',y=costlong)+
    scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
    theme_minimal()+
    theme(axis.line=element_line(colour='grey10'))
ggsave(filename=paste0('graphs/rswe_',recon.version,'/boxplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=10,height=4.5)

###

ggplot(skill.m)+
     geom_freqpoly(aes(x=Rrmse,colour=variable),alpha=0.75,size=2,binwidth=0.01)+
     labs(x=costlong)+
     scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
     theme_minimal()+
     theme(axis.line=element_line(colour='grey10'))
ggsave(filename=paste0('graphs/rswe_',recon.version,'/freqpoly_',costshort,'_',style,'_blended&unblended.pdf'),width=10,height=4.5)

###

ggplot(skill.m)+      
     geom_point(aes(x=as.numeric(doy),y=Rrmse,colour=variable),alpha=.75)+    
     facet_grid(.~yr)+
     labs(y=costlong)+
     scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
     scale_x_continuous('day of year',limits=c(0,160),labels=c(1,seq(30,160,30)),breaks=c(1,seq(30,160,30)))+
     #scale_y_continuous(limits=c(0,.5))+
     guides(colour=guide_legend(override.aes=list(alpha=1)))+
     theme_bw()+
     theme(axis.line=element_line(colour='grey10'),
           strip.background=element_rect(fill='white',colour='white'),
           strip.text=element_text(face='bold',size='12'),
           axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=.5))

ggsave(filename=paste0('graphs/rswe_',recon.version,'/scatterplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=11,height=3)
}
