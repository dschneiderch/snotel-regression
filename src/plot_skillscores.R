library('ProjectTemplate')
load.project()

cost='rmse'
style='real-time'
recon.version='v3.2'
skill=read.table(paste0('diagnostics/rswe_',recon.version,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T)
skill$date=as.POSIXct(skill$date,tz='MST')
skill$phvrcn_recondate=as.POSIXct(skill$phvrcn_recondate,tz='MST')
skill$recon_costdate=as.POSIXct(skill$recon_costdate,tz='MST')
str(skill)
##
if(cost=='rmse'){
	costshort='RelRMSE'
	costlong='Relative RMSE'
library(dplyr)
dailyavg=snotelrec[snotelrec$date %in% skill$date,] %>%
group_by(date) %>%
select(swe) %>%
summarise(	avgsnotel=mean(swe)	)

skill.m=melt(skill[,c('date','yr','skill_phv','skill_phvrcn','skill_phvfull','skill_phvrcnfull','yrdoy')],.(date,yr,yrdoy))
skill.m=cbind(skill.m,dailyavg)
skill.m$Rrmse=with(skill.m,value/avgsnotel)
} else if(cost=='r2'){
    costshort='r2'
    costlong=bquote(r^2)
    skill.m=melt(skill[,c('date','yr','skill_phv','skill_phvrcn','skill_phvfull','skill_phvrcnfull','yrdoy')],.(date,yr,yrdoy))
}

ggplot(skill.m)+
    geom_boxplot(aes(x=as.factor(yr),y=Rrmse,colour=variable))+
    labs(x='year',y=costlong)+
    scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
    theme_minimal()+
    theme(axis.line=element_line(colour='grey10'))
ggsave(filename=paste0('graphs/rswe_',recon.version,'/boxplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=10,height=4.5)

###

skill.m$mnthdy=strftime(skill.m$date,'%b%d')
skill.m$mnth=factor(strftime(skill.m$date,'%b'),levels=c('Jan','Feb','Mar','Apr','May'))
ggplot(skill.m)+
     geom_freqpoly(aes(x=Rrmse,colour=variable),alpha=0.5,size=2,binwidth=0.01)+
     labs(x=costlong)+
     scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
     theme_minimal()+
     theme(axis.line=element_line(colour='grey10'))
ggsave(filename=paste0('graphs/rswe_',recon.version,'/freqpoly_',cost,'_',style,'_blended&unblended.pdf'),width=10,height=4.5)

###
skill.m$doy=strftime(skill.m$date,'%j')
ggplot(skill.m)+      
     geom_point(aes(x=as.numeric(doy),y=Rrmse,colour=variable),alpha=.5)+    
     facet_grid(.~yr)+
     labs(y=costlong)+
     scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
     scale_x_continuous('day of year',limits=c(30,160),breaks=seq(30,160,30),labels=seq(30,160,30))+
     scale_y_continuous(limits=c(0,1))+
     theme_bw()+
     theme(axis.line=element_line(colour='grey10'),
           strip.background=element_rect(fill='white',colour='white'),
           axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
ggsave(filename=paste0('graphs/rswe_',recon.version,'/scatterplot_',cost,'_',style,'_byyear_blended&unblended.pdf'),width=9,height=3)
 
