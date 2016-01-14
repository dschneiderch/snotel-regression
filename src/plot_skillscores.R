library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project()
library(gtable)
library(xtable)
#plotting functions
source('lib/skill_score_plot_functions.R')
#get daily average snotel swe that was scaled by fsca 
dailyavg_scaled=readRDS('data/dailyavg_snotel_scaled.rds')

cost='r2'
fscaMatch='wofsca'
snotelscale='scale'

style='real-time'
recon.version='v3.1'
covrange='idp1'
dateflag='B'

for(cost in c('rmse')){
for(snotelscale in c('scale')){
for(fscaMatch in c('wofsca')){

graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/',fscaMatch,'/',style)

#skill=read.table(paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/fassnacht/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T)
fnpath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag)
fnpattern=paste0(style,'_recondate_selection_',fscaMatch,'_',cost)
fns=list.files(path=fnpath,pattern=fnpattern,full.names=T)
skill=ldply(fns,read.table,sep='\t',header=T,stringsAsFactors=F)
skill$date=as.POSIXct(skill$date,tz='MST')
skill$phvrcn_recondate=as.POSIXct(skill$phvrcn_recondate,tz='MST')
skill$recon_costdate=as.POSIXct(skill$recon_costdate,tz='MST')
# str(skill)
##
skill$fulldiff=skill$skill_phvrcnfull-skill$skill_phvfull
skill$diff=skill$skill_phvrcn-skill$skill_phv
skill.m=melt(skill[,c('date','yr','skill_phv','skill_phvrcn','skill_phvfull','skill_phvrcnfull','fulldiff','diff','yrdoy')],.(date,yr,yrdoy))
skill.m$set=ifelse(grepl('skill',as.character(skill.m$variable)),'data','diff')
skill.m$set=as.factor(skill.m$set)
skill.m$modeltype=ifelse(grepl('phvrcn',as.character(skill.m$variable)),'phvrcn','phv')
skill.m$modeltype=ifelse(grepl('diff',as.character(skill.m$variable)),'diff',skill.m$modeltype)
skill.m$modeltype=as.factor(skill.m$modeltype)
skill.m$residblend=ifelse(grepl('full',as.character(skill.m$variable)),'blended','unblended')
skill.m$residblend=factor(skill.m$residblend,levels=c('unblended','blended'))
skill.m$mnthdy=strftime(skill.m$date,'%b%d')
skill.m$mnth=factor(strftime(skill.m$date,'%b'),levels=c('Jan','Feb','Mar','Apr','May','Jun'))
skill.m$doy=as.numeric(strftime(skill.m$date,'%j'))

## -- calculate daily average snotel for relative calcs below.
if(snotelscale=='noscale'){
dailyavg=ddply(snotelrec[snotelrec$date %in% skill$date,],.(date),function(x){
    summarise(x,
        avgsnotel=mean(swe)
        )})
} else if(snotelscale=='scale'){
    dailyavg=dailyavg_scaled
}
# write.table(x=dailyavg,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_dailysnotelswe.txt'),quote=F,row.names=F)

## --- set up labels based on skill metric and compute relative rmse if applicable
if(cost=='rmse'){
    
    costshort='rmse'
    costlong=expression('RMSE')
    diffcostlong=expression(Delta * RMSE)
    
    skill.m=cbind(skill.m,avgsnotel=dailyavg$avgsnotel)
    skill.m$Rrmse=skill.m$value/skill.m$avgsnotel

    summrel=ddply(skill.m,.(variable,yr),summarise,
        avg=mean(Rrmse,na.rm=T),
        med=median(Rrmse,na.rm=T),
        min=min(Rrmse,na.rm=T),
        max=max(Rrmse,na.rm=T))
# write.table(x=summrel,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',costshort,'.txt'),quote=F,row.names=F)

# ddply(skill.m,.(variable,yr),function(x){
#     summarise(x,
#     min=doy[which.min(Rrmse)],
#     max=doy[which.max(Rrmse)]
#     )
# })
} else if(cost=='r2'){
    costshort='r2'
    costlong=expression(r^2)
    diffcostlong=expression(Delta * r^2)
    # bquote(r^2)
    
}

}}}

## --------------
if(cost=='r2'){ 
summ=ddply(skill.m,.(variable,yr), function(x){
    summarise(x,
    avg=mean(value,na.rm=T),
    med=median(value,na.rm=T),
    min=min(value,na.rm=T),
    max=max(value,na.rm=T),
    sd=sd(value,na.rm=T),
    cv=sd/avg)
})
} else if(cost=='rmse'){
summ=ddply(skill.m,.(variable,yr), function(x){
    summarise(x,
    avg=mean(Rrmse,na.rm=T),
    med=median(Rrmse,na.rm=T),
    min=min(Rrmse,na.rm=T),
    mindoy=doy[which.min(Rrmse)],
    max=max(Rrmse,na.rm=T),
    sd=sd(Rrmse,na.rm=T),
    cv=sd/avg)
})
}

head(summ)
ddply(summ,.(variable),function(dF){
    summarise(dF,
        avg=mean(avg,na.rm=T),
        min=min(min),
        min(mindoy),
        max(mindoy)
        )
})
range(summ$mindoy)
summ2=subset(summ,!grepl('full',as.character(variable)))
mean(dcast(summ2,yr~variable,value.var='sd')$skill_phv)
mean(dcast(summ2,yr~variable,value.var='sd')$skill_phvrcn)
dcast(summ,variable~.,value.var='mindoy',fun.agg=mean)
dcast(summ,variable~.,value.var='mindoy',fun.agg=sd)
dcast(summ2,variable~.,value.var='sd',fun.agg=mean)

str(summ)

# write.table(x=summ,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',cost,'.txt'),quote=F,row.names=F)

monavg=dcast(skill.m,mnth~variable,value.var='Rrmse',fun.agg=mean,na.rm=T)

# stat='med'
# dF=dcast(summ,variable~yr,value.var=stat,margins=T,fun.aggregate=mean)
# pdf(paste0(graphbase,'/',stat,'_',costshort,'_',style,'_byyear_blended&unblended&diff.pdf'),
#     width=12,
#     height=3)
# grid.arrange(tableGrob(format(dF,digits=2)))
# dev.off()

# if(cost=='rmse'){
# dFrel=dcast(summrel,variable~yr,value.var=stat,margins=T,fun.aggregate=mean)
# pdf(paste0(graphbase,'/',stat,'_RelRMSE_',style,'_byyear_blended&unblended&diff.pdf'),
#     width=12,
#     height=3)
# grid.arrange(tableGrob(format(dFrel,digits=2)))
# dev.off()
#  }

## -- find min and max doy 
summdoy=ddply(skill.m,.(variable,yr), function(x){
    summarise(x,
    min=doy[which.min(value)],
    max=doy[which.max(value)]
)})
dcast(summdoy,yr~variable,value.var='min')
dcast(summdoy,yr~variable,value.var='max')
sd(dcast(summdoy,yr~variable,value.var='max')$skill_phvrcn)
sd(dcast(summdoy,yr~variable,value.var='max')$skill_phv)

head(summdoy)
temp=subset(summdoy,variable=='skill_phvrcn')
sd(temp$min)
sd(temp$max)

summdoy$model=ifelse(grepl('diff',summdoy$modeltype),'modeldiff','model')
ddply(summdoy,.(model),function(x){
    summarise(x,
        minrange=range(min),
        maxrange=range(max)
        )
})
# ddply(summdoy,.(modeltype),function(x){
#     summarise(x,
#         avgmin=mean(min),
#         avgmax=mean(max),
#         sdmin=sd(min),
#         sdmax=sd(max)
#         )
# })





## ----- Visual inspection where rRMSE drastically changes.
# test=skill.m[!is.na(skill.m$Rrmse),]
# test=subset(test,modeltype=='phvrcn' & variable=='skill_phvrcnfull')
# qplot(test$doy,test$Rrmse,colour=as.factor(test$yr))+geom_line()+scale_colour_brewer(palette='Paired')

# spldf=ddply(test,.(yr),function(x){ 
#     doy=x$doy
#     splfit1=predict(smooth.spline(x$doy,x$Rrmse),deriv=1)$y
#     splfit2=predict(smooth.spline(x$doy,x$Rrmse),deriv=2)$y
#     # splwin=rollapply(splfit,7,median)
#     data.frame(doy,splfit1,splfit2)
# })

# rollapply <- function(x, n, f, ...) {
#   out <- rep(NA, length(x))

#   offset <- trunc(n / 2)
#   for (i in (offset + 1):(length(x) - n + offset + 1)) {
#     out[i] <- f(x[(i - offset):(i + offset)], ...)
#   }
#   out
# }
# dev.new()
# ggplot(spldf,aes(x=doy,y=splfit1,colour=as.factor(yr)))+
# geom_point(size=2)+
# geom_line(size=1,alpha=.8)+
# geom_hline(aes(yintercept=0.01))+
# scale_colour_brewer(palette='Paired')+
# scale_x_continuous(breaks=seq(1,150,2),labels=seq(1,150,2))+
# coord_cartesian(ylim=c(-0.01,0.05),xlim=c(45,150))+
# theme(panel.background=element_rect(fill='grey10'),
#     panel.grid=element_line(colour='grey50'))
      
# ## Days of the year that rRMSE turns up. Visual inspection.
# spikedoy=data.frame(yr=seq(2001,2012),doy=c(92,90,90,68,91,91,75,108,112,91,65,78))#some initial guesses
# output=ddply(test,.(yr),function(x){
#     doi=spikedoy[spikedoy$yr==unique(x$yr),'doy']
#     dind=which.min(abs(x$doy/doi-1))#find the closest date to the guess that was made from the graph
#     summarise(x[1:dind,],
#         avg=mean(Rrmse,na.rm=T),
#         sd=sd(Rrmse,na.rm=T),
#         med=median(Rrmse,na.rm=T),
#         min=min(Rrmse,na.rm=T),
#         doy=doy[dind])
# })
# spikedoy=output[,c('yr','doy','avg','min')]
# ## -------





# ## ---- plot rmse and r2 boxplots
# if(cost=='rmse'){
#   mg=myboxplot(skill.m,'value')+coord_cartesian(ylim=c(0,0.5))
#   #scale_y_continuous(lim=c(0,5))
# } else if( cost=='r2'){
#   mg=myboxplot(skill.m,'value')
#   # scale_y_continuous(lim=c(0,1))
# }
# ggsave(plot=mg,filename=paste0(graphbase,'/boxplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=6.6,height=3)
    
###
# ggplot(skill.m)+
#      geom_freqpoly(aes(x=value,colour=variable),alpha=0.75,size=2,binwidth=0.01)+
#      labs(x=costlong)+
#      scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
#      theme_minimal()+
#      theme(axis.line=element_line(colour='grey10'))
# ggsave(filename=paste0('graphs/rswe_',recon.version,'/freqpoly_',costshort,'_',style,'_blended&unblended.pdf'),width=10,height=4.5)


## plot monthly skill scores as barplot
skill.m.sub=subset(skill.m,residblend=='unblended' & set=='data' & modeltype!='diff')
monthskill=dcast(skill.m.sub,modeltype+residblend+set~mnth,value.var='value',fun.agg=mean,na.rm=T)
monthskillplot=melt(monthskill,.(modeltype,residblend,set))
monthskillplot$modeltype=ifelse(monthskillplot$modeltype=='phv','PHV-baseline','PHV-RCN')
if(cost=='r2'){
ggmnth=monthlybarplot(monthskillplot,'value')+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
        scale_fill_manual('Model',values=c('blue','red'))+
        annotate('text',label='a',x=.8,y=0.95)
} else if(cost=='rmse'){
ggmnth=monthlybarplot(monthskillplot,'value')+
        scale_y_continuous(expand=c(0,0))+
        annotate('text',label='b',x=.8,y=0.19)
}
ggsave(plot=ggmnth,paste0(graphbase,'/barplot_',costshort,'_',style,'_monthly_unblended.pdf'),width=6.6,height=3)
  
if(cost=='rmse'){
    costshort='RelRMSE'
    costlong='Relative RMSE'
    diffcostlong=expression(Delta * 'Rel RMSE')

monthskill=dcast(skill.m.sub,modeltype+residblend+set~mnth,value.var='Rrmse',fun.agg=mean,na.rm=T)
monthskillplot=melt(monthskill,.(modeltype,residblend,set))
monthskillplot$modeltype=ifelse(monthskillplot$modeltype=='phv','PHV-baseline','PHV-RCN')
ggmnth=monthlybarplot(monthskillplot,'value')+
        scale_y_continuous(expand=c(0,0))+
        coord_cartesian(ylim=c(0,2))+
        annotate('text',label='b',x=.8,y=1.9)+##.95*2=1.9 -- use y value from r2 plot * max of ylim to position text label.
        scale_fill_manual('Model' ,values=c('blue','red'),guide=F)

ggsave(plot=ggmnth,paste0(graphbase,'/barplot_',costshort,'_',style,'_monthly_unblended.pdf'),width=6.6,height=3)
}



###
## --scatter plot of skill scores
#### create a facet version so we can steal the legend.
lp=multifacetplot(skill.m,'value',shapeguide='F') 

## create 1st plot for combined figure
skill.mplot=subset(skill.m,modeltype!='diff' & residblend=='unblended')
if(cost=='rmse'){
  maxrmse=round(max(skill.mplot$value,na.rm=T),1)
  annotedF=data.frame(yr=2001,label='b',x=60,y=0.9*maxrmse)
  dap=dataplot(skill.mplot,'value')+
       geom_text(data=annotedF,aes(x,y,label=label),size=4,face='bold')+
       coord_cartesian(ylim=c(0,maxrmse))
   # +   scale_y_continuous(limits=c(0,5),breaks=seq(0,5,1),labels=seq(0,5,1))
   } else if(cost=='r2'){
  annotedF=data.frame(yr=2001,label='a',x=60,y=0.9)
  dap=dataplot(skill.mplot,'value') +
        geom_text(data=annotedF,aes(x,y,label=label),size=4,face='bold')+
        coord_cartesian(ylim=c(0,1))
    # +  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1),labels=seq(0,1,.1))
   }

## create 2nd plot for combined figure
skill.mplot=subset(skill.m,modeltype=='diff' & residblend=='unblended')
if(cost=='rmse'){
dfp=diffplot(skill.mplot,'value')
# +    scale_y_continuous(limits=c(0,5),breaks=seq(0,5,1),labels=seq(0,5,1))
} else if(cost=='r2'){
dfp=diffplot(skill.mplot,'value') +  scale_y_continuous(limits=c(-.5,.5),breaks=seq(-.5,.5,.5),labels=seq(-.5,.5,.5))
}

## combine plots and save
gA <- ggplot_gtable(ggplot_build(dap))
    gB <- ggplot_gtable(ggplot_build(dfp))
if(cost=='r2'){
    maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
    gA$widths[2:3] <- as.list(maxWidth)
    gB$widths[2:3] <- as.list(maxWidth)
    saveRDS(maxWidth,paste0(graphbase,'/scatterplot_grobwidth.rds'))
} else if(cost=='rmse'){
    maxWidth=readRDS(paste0(graphbase,'/scatterplot_grobwidth.rds'))
    gA$widths[2:3] <- as.list(maxWidth)
    gB$widths[2:3] <- as.list(maxWidth)
}

pdf(paste0(graphbase,'/scatterplot_',costshort,'_',style,'_byyear_unblended&diff.pdf'),
    width=6.6,
    height=3.5)
# grid.arrange(gA,gB,heights=c(0.65,0.35))
if(cost=='r2'){ #was trying to make r2 figure without the legend
grid.arrange(
        arrangeGrob(gA, gB,nrow=2,heights=c(.8,.2)),
            glegend(lp),nrow=2,heights=c(.65,.2))
} else if(cost=='rmse'){
grid.arrange(
    arrangeGrob(gA, gB,nrow=2,heights=c(.75,.25)),
        glegend(lp), nrow=2,heights=c(.65,.2))    
}
 dev.off()
     
 
#plot relative rmse
if(cost=='rmse'){

    costshort='RelRMSE'
    costlong='Relative RMSE'
    diffcostlong=expression(Delta * 'Rel RMSE')

    mbp=myboxplot(skill.m,'Rrmse')+coord_cartesian(ylim=c(0,2))
    ggsave(plot=mbp,filename=paste0(graphbase,'/boxplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=6.6,height=3)

###

# ggplot(skill.m)+
#      geom_freqpoly(aes(x=Rrmse,colour=variable),alpha=0.75,size=2,binwidth=0.01)+
#      labs(x=costlong)+
#      scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
#      theme_minimal()+
#      theme(axis.line=element_line(colour='grey10'))
# ggsave(filename=paste0('graphs/rswe_',recon.version,'/freqpoly_',costshort,'_',style,'_blended&unblended.pdf'),width=10,height=4.5)

###

skill.mplot=subset(skill.m,modeltype!='diff' & residblend=='unblended')
annotedF=data.frame(yr=2001,label='b',x=60,y=0.9*2)#0.9*2 should give the same relative height as 0.9*1 in r2 plot. 2 should match ylim max
dap2=dataplot(skill.mplot,'Rrmse')+
        geom_text(data=annotedF,aes(x,y,label=label),size=4,face='bold')+
        coord_cartesian(ylim=c(0,2))

skill.mplot=subset(skill.m,modeltype=='diff' & residblend=='unblended')
#str(skill.mplot)
range(skill.mplot$Rrmse,na.rm=T)
dfp2=diffplot(skill.mplot,'Rrmse')#+coord_cartesian(ylim=c(-0.2,0))

gA <- ggplot_gtable(ggplot_build(dap2))
gB <- ggplot_gtable(ggplot_build(dfp2))
maxWidth=readRDS(paste0(graphbase,'/scatterplot_grobwidth.rds'))
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)

pdf(paste0(graphbase,'/scatterplot_',costshort,'_',style,'_byyear_unblended&diff.pdf'),
    # res=300,
    width=6.6,
    height=3)
# grid.arrange(gA,gB,heights=c(0.65,0.35))
grid.arrange(
    arrangeGrob(gA, gB,nrow=2,heights=c(.75,.25)),
        glegend(lp), nrow=2,heights=c(.65,.2))
dev.off()
   


  

# ggplot(skill.m)+      
#      geom_point(aes(x=as.numeric(doy),y=Rrmse,colour=variable),alpha=.75,size=1)+    
#      facet_grid(set~yr,scale='free',space='free')+
#      labs(y='')+
#      scale_colour_brewer(name='',palette='Set1',labels=c('PHV (unblended)','PHV+RCN (unblended)','PHV (blended)','PHV+RCN (blended)','Blended Skill Difference'))+
#      scale_x_continuous('day of year',limits=c(1,160),labels=c(1,seq(30,160,30)),breaks=c(1,seq(30,160,30)))+
#      scale_y_continuous(breaks=seq(0,.7,.1),labels=seq(0,.7,.1))+
#      guides(colour=guide_legend(ncol=3, byrow=F,override.aes=list(alpha=1,size=2)))+
#      theme_bw()+
#      theme(axis.line=element_line(colour='grey10'),
#            strip.background=element_rect(fill=NA,colour=NA),
#            strip.text=element_blank(),#element_text(face='bold',size='12'),
#            axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=.5),
#            axis.text.y=element_text(size=6,angle=0),
#            legend.key=element_rect(colour='white'),
#            legend.position = 'bottom', 
#            legend.box = "horizontal")
# 
#      
#      scale_colour_discrete(name='Model',labels=c('PHV Regression (unblended)','PHV+RCN Regression (unblended)','PHV Regression (blended)','PHV+RCN Regression (blended)'))+
#      scale_x_continuous('day of year',limits=c(0,160),labels=c(1,seq(30,160,30)),breaks=c(1,seq(30,160,30)))+
#      #scale_y_continuous(limits=c(0,.5))+
#      guides(colour=guide_legend(override.aes=list(alpha=1)))+
#      theme_bw()+
#      theme(axis.line=element_line(colour='grey10'),
#            strip.background=element_rect(fill='white',colour='white'),
#            strip.text=element_text(face='bold',size='12'),
#            axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=.5))

# ggsave(filename=paste0(graphbase,'/scatterplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=11,height=3)
}
