library('ProjectTemplate')
# setwd('~/Documents/snotel-regression_project')
reload.project()
library(cowplot)
library(grid)
# library(gridExtra)
# library(gtable)
library(dplyr)
#plotting functions
source('lib/skill_score_plot_functions.R')
#get daily average snotel swe that was scaled by fsca 
dailyavg_scaled=readRDS('data/dailyavg_snotel_scaled.rds')

cost='rmse'
fscaMatch='wofsca'
snotelscale='scale'
style='real-time'
recon.version='v3.1'
covrange='idp1'
dateflag='B'

# for(cost in c('rmse')){
# for(snotelscale in c('scale')){
# for(fscaMatch in c('wofsca')){

graphbase=paste0('graphs/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/',fscaMatch,'/',style)

plotdata=setup_plotdata(cost,recon.version,covrange,snotelscale,dateflag,fscaMatch,style)
costshort=plotdata[[1]]
costlong=plotdata[[2]]
diffcostlong=plotdata[[3]]
skill.m=plotdata[[4]]

## summarise skill.m --------
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
}
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

##  find min and max doy ----
# summdoy=ddply(skill.m,.(variable,yr), function(x){
#     summarise(x,
#     min=doy[which.min(value)],
#     max=doy[which.max(value)]
# )})
# dcast(summdoy,yr~variable,value.var='min')
# dcast(summdoy,yr~variable,value.var='max')
# sd(dcast(summdoy,yr~variable,value.var='max')$skill_phvrcn)
# sd(dcast(summdoy,yr~variable,value.var='max')$skill_phv)
# 
# head(summdoy)
# temp=subset(summdoy,variable=='skill_phvrcn')
# sd(temp$min)
# sd(temp$max)
# 
# summdoy$model=ifelse(grepl('diff',summdoy$modeltype),'modeldiff','model')
# ddply(summdoy,.(model),function(x){
#     summarise(x,
#         minrange=range(min),
#         maxrange=range(max)
#         )
# })
# ddply(summdoy,.(modeltype),function(x){
#     summarise(x,
#         avgmin=mean(min),
#         avgmax=mean(max),
#         sdmin=sd(min),
#         sdmax=sd(max)
#         )
# })



## Visual inspection where rRMSE drastically changes. ----
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


## ---- plot rmse and r2 boxplots ----
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


## plot monthly skill scores as barplot ----
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
# ggsave(plot=ggmnth,paste0(graphbase,'/barplot_',costshort,'_',style,'_monthly_unblended.pdf'),width=6.6,height=3)

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
     
     # ggsave(plot=ggmnth,paste0(graphbase,'/barplot_',costshort,'_',style,'_monthly_unblended.pdf'),width=6.6,height=3)
}



## scatter plot of skill scores ----
## ----
## create 1st plot for combined figure ----
skill.mplot=subset(skill.m,modeltype!='diff' & residblend=='unblended')
yrbreak=2007
if(cost=='rmse'){
     maxrmse=round(max(skill.mplot$value,na.rm=T),1)
     annotedF=data.frame(yr=2001,label='b',x=30,y=0.8*maxrmse)
     dap1=filter(skill.mplot,yr<yrbreak) %>%
          dataplot(.,'value')+
          geom_text(data=annotedF,aes(x,y,label=label),size=6,fontface='bold')+
          coord_cartesian(ylim=c(0,maxrmse))+
          scale_y_continuous(expand=c(0,0))
     dap2=dplyr::filter(skill.mplot,yr>=yrbreak) %>%
          dataplot(.,'value')+
          scale_y_continuous(expand=c(0,0))
} else if(cost=='r2'){
     annotedF=data.frame(yr=2001,label='a',x=30,y=0.8)
     dap1=dplyr::filter(skill.mplot,yr<yrbreak) %>%
          dataplot(.,'value')+
          geom_text(data=annotedF,aes(x,y,label=label),size=6,fontface='bold')+
          scale_y_continuous(expand=c(0,0))+
          coord_cartesian(ylim=c(0,1))
     dap2=dplyr::filter(skill.mplot,yr>=yrbreak) %>%
          dataplot(.,'value')+
          scale_y_continuous(expand=c(0,0))
     # dap=dataplot(skill.mplot,'value') +
     #       geom_text(data=annotedF,aes(x,y,label=label),size=4,fontface='bold')+
     #       coord_cartesian(ylim=c(0,1))
     #   # +  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.1),labels=seq(0,1,.1))
}

## create 2nd plot for combined figure ----
skill.mplot=subset(skill.m,modeltype=='diff' & residblend=='unblended')
if(cost=='rmse'){
     dfp1=skill.mplot %>% 
          filter(yr<yrbreak) %>%
          diffplot(.,'value') +
          scale_y_continuous(expand=c(0,0))
     dfp2=skill.mplot %>% 
          filter(yr>=yrbreak) %>%
          diffplot(.,'value')  +
          scale_y_continuous(expand=c(0,0))
     # +    scale_y_continuous(limits=c(0,5),breaks=seq(0,5,1),labels=seq(0,5,1))
} else if(cost=='r2'){
     dfp1=skill.mplot %>% 
          filter(yr<yrbreak) %>%
          diffplot(.,'value') + 
          scale_y_continuous(limits=c(-.5,.5),breaks=seq(-.5,.5,.5),labels=seq(-.5,.5,.5),expand=c(0,0))
     dfp2=skill.mplot %>% 
          filter(yr>=yrbreak) %>%
          diffplot(.,'value') + 
          scale_y_continuous(limits=c(-.5,.5),breaks=seq(-.5,.5,.5),labels=seq(-.5,.5,.5),expand=c(0,0))
}

## create a facet version so we can steal the legend. ----
lp=multifacetplot(skill.m,'value',shapeguide='F')
grobs <- ggplotGrob(lp)$grobs
legend_b <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

## combine plots and save ----
#make sure all plots have same width
if(cost=='r2'){
     gb=ggplotGrob(dap1)
     gb2=ggplotGrob(dfp1)
     gb$widths=grid:::unit.list(gb$widths)
     gb2$widths=grid:::unit.list(gb2$widths)
     maxWidth=as.list(unit.pmax(gb$widths,gb2$widths))
     saveRDS(maxWidth,paste0(graphbase,'/scatterplot_grobwidth.rds'))
} else {
     maxWidth=readRDS(paste0(graphbase,'/scatterplot_grobwidth.rds'))
}

dap1=unify_grobwidth(dap1,maxWidth)
dap2=unify_grobwidth(dap2,maxWidth)
dfp1=unify_grobwidth(dfp1,maxWidth)
dfp2=unify_grobwidth(dfp2,maxWidth)
plotA=plot_grid(dap1,dfp1,dap2,dfp2,ncol=1,rel_heights = c(1.5,1,1.5,1))
plotA=plot_grid(plotA,legend_b,ncol=1,rel_heights=c(1.1,.1))
# plotA
save_plot(base_height=3,
          base_aspect_ratio = 2,
          nrow=2,
          plot=plotA,
          filename=paste0(graphbase,'/scatterplot_',costshort,'_',style,'_byyear_unblended&diff.pdf'))


# plot relative rmse ----
if(cost=='rmse'){
     
     costshort='RelRMSE'
     costlong='RelRMSE'
     diffcostlong=expression(Delta * 'RelRMSE')
     
     #### boxplots rrmse ----
     mbp=myboxplot(skill.m,'Rrmse')+coord_cartesian(ylim=c(0,2))
     # ggsave(plot=mbp,filename=paste0(graphbase,'/boxplot_',costshort,'_',style,'_byyear_blended&unblended.pdf'),width=6.6,height=3)
     
     #### scatter plots Rrmse ----
     skill.mplot=subset(skill.m,modeltype!='diff' & residblend=='unblended')
     annotedF=data.frame(yr=2001,label='b',x=30,y=0.8*2)#0.9*2 should give the same relative height as 0.9*1 in r2 plot. 2 should match ylim max
     dap1=filter(skill.mplot,yr<yrbreak) %>%
          dataplot(.,'Rrmse')+
          geom_text(data=annotedF,aes(x,y,label=label),size=6,fontface='bold')+
          coord_cartesian(ylim=c(0,2))+
          scale_y_continuous(expand=c(0,0))
     dap2=filter(skill.mplot,yr>=yrbreak) %>%
          dataplot(.,'Rrmse')+
          coord_cartesian(ylim=c(0,2))+
          scale_y_continuous(expand=c(0,0))
     # +   scale_y_continuous(limits=c(0,5),breaks=seq(0,5,1),labels=seq(0,5,1))
     
     ## create 2nd plot for combined figure
     skill.mplot=subset(skill.m,modeltype=='diff' & residblend=='unblended')
     dfp1=skill.mplot %>% 
          filter(yr<yrbreak) %>%
          diffplot(.,'Rrmse')+
          scale_y_continuous(expand=c(0,0))
     dfp2=skill.mplot %>% 
          filter(yr>=yrbreak) %>%
          diffplot(.,'Rrmse')  +
          scale_y_continuous(expand=c(0,0))
     
     
     ## combine plots and save
     dap1=unify_grobwidth(dap1,maxWidth)
     dap2=unify_grobwidth(dap2,maxWidth)
     dfp1=unify_grobwidth(dfp1,maxWidth)
     dfp2=unify_grobwidth(dfp2,maxWidth)
     plotA=plot_grid(dap1,dfp1,dap2,dfp2,ncol=1,rel_heights = c(1.5,1,1.5,1))
     plotA=plot_grid(plotA,legend_b,ncol=1,rel_heights=c(1.1,.1))
     # plotA
     # ggsave(plotA,filename=paste0(graphbase,'/scatterplot_test2.pdf'),width=6.5)
     save_plot(base_height=3,
               base_aspect_ratio = 2,
               nrow=2,
               plot=plotA,
               filename=paste0(graphbase,'/scatterplot_',costshort,'_',style,'_byyear_unblended&diff.pdf'))
     
     
     
     ## combine plots and save
     # pdf(paste0(graphbase,'/scatterplot_',costshort,'_',style,'_byyear_unblended&diff-part1.pdf'), width=6.6, height=3, onefile = FALSE)
     # arrange_2x1_shared_legend(dap1,dfp1,NULL)
     # dev.off()
     
     # pdf(paste0(graphbase,'/scatterplot_',costshort,'_',style,'_byyear_unblended&diff-part2.pdf'), width=6.6, height=3, onefile = FALSE)
     # arrange_2x1_shared_legend(dap2,dfp2,lp)
     # dev.off()
}


