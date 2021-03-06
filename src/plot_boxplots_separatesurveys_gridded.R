# setwd('~/Documents/snotel-regression_project')
# setwd('~/Documents/snotel-regression_project')

library(raster)
library(reshape2)
library(plyr)
library(ggplot2)
library(doMC)
library(RColorBrewer)
library(gridExtra)
library(dplyr)

rswe='v3.1'
covrange='idp1'
fscaMatch='wofsca'
scalesnotel='scale'
cost='r2'
dateflag='surveyvalidation'
blending='unblended'
style='real-time'
gridding=''#_snodasgrid'
postscaled=''#_postscaled'

swe=readRDS(file=paste0('data/spatialvalidation/',blending,'/',dateflag,'/',style,'/surveygriddedswe',gridding,'_rswe',rswe,'_',covrange,'_snotel',scalesnotel,postscaled,'_',cost,'_',fscaMatch,'.rds'))
swe_all=swe
cellcount=ddply(swe_all,.(model.type),function(dF) plyr::count(dF$cell))
freqdf=dcast(cellcount,model.type~x,value.var='freq')
ind=which(freqdf[1,2:ncol(freqdf)]<140)#about 278 30m cells in an 500m pixel #rows 1,2,3 will be the same
cell2rm=as.numeric(colnames(freqdf)[-1][ind])
ind2=which(freqdf[nrow(freqdf),2:ncol(freqdf)]<560)#about 1111 30m cells in 1km pixel #last row is snodas
cell2rm2=as.numeric(colnames(freqdf)[-1][ind2])
swe=subset(swe_all,(model.type!='SNODAS' & !(cell %in% cell2rm)) | (model.type=='SNODAS' & !(cell %in% cell2rm2)))
swe$swe.model[swe$swe.model<0.001]=0

## define cost function
docost=function(y,yhat,cost){
if(cost=='mae') { #statistics to minimize
          costfun<-function(y, yhat) mean(abs(yhat-y),na.rm=T)
} else if(cost=='rmse') {
          costfun<-function(y,yhat) sqrt(mean((yhat-y)^2,na.rm=T))
} else if (cost=='relrmse'){
          costfun<-function(y,yhat) sqrt(mean((yhat-y)^2,na.rm=T))/mean(y,na.rm=T)
} else if (cost=='bias'){
          costfun<-function(y,yhat) mean(yhat-y,na.rm=T)
} else if (cost=='pctbias'){
              costfun<-function(y,yhat) mean(yhat-y,na.rm=T)/mean(y,na.rm=T)*100
          }

     costfun(y,yhat)
}

cellavg=ddply(swe,.(model.type,date,site,yr,mth,cell),function(dF){
    summarise(dF,
            swe.obs.avg=mean(swe.obs,na.rm=T),#the cell number in swe represents the model cell number so there can be muliple obs cell values (30m)
                                        # na.rm=F will average only for model pixels that are completely covered by observations, otherwise there can be overlap
              swe.model=unique(swe.model))# there should only be one value for each cell
    })


## average rmse of each model for each survey ----
# sweerror=ddply(cellavg,.(model.type,site,date,mth,yr),function(dF) with(dF,docost(swe.obs.avg,swe.model,'rmse')))
# names(sweerror)[length(sweerror)]='objcost'#c('source','site','date','mth','yr','objcost')
sweerror=cellavg %>%
     group_by(model.type,site,date,mth,yr) %>%
     summarise(objcost=docost(swe.obs.avg,swe.model,'rmse'),
               swe.obs.avg=mean(swe.obs.avg,na.rm=TRUE))

filter(sweerror,grepl('PHV',model.type)) %>%
     ggplot(.,aes(x=swe.obs.avg,y=objcost/swe.obs.avg,colour=model.type)) +
          geom_point()+
          geom_smooth(method='glm',formula=y~x, method.args = list(family = gaussian(link='identity')))


## combine modeled and observed swe for boxplots ----
obs=swe[,c('site','date','yr','mth','swe.obs')]
obs$source='OBS'
mdl=swe[,c('site','date','yr','mth','swe.model','model.type')]
names(obs) <- names(mdl) <- c('site','date','yr','mth','swe','source')
sweboxplot=rbind(obs,mdl)

## average swe of each source for each survey ----
sweavg=ddply(sweboxplot,.(site,date,mth,yr,source),function(dF) summarise(dF, avg=mean(swe,na.rm=T)))

## meromy bias ----
meromybias=read.csv('data/spatialvalidation/snow_survey_table.csv',stringsAsFactor=F)
colnames(meromybias)=c('site','date','survey.swe.avg.m','snotel.swe.m','bias.pct')
meromybias$date=as.POSIXct(strptime(as.character(meromybias$date),'%m/%d/%Y'))
meromybias=arrange(meromybias,date,site)


## subset sweboxplot for phv and phvrcn ----
sweboxplot=sweboxplot %>%
     filter(source!='RCN', source!='SNODAS') %>%
     mutate(surveytype=ifelse(site=='Green Lakes Valley','alpine','forested'))
sweboxplot_max=ddply(sweboxplot,.(date,site,source),summarise,maxswe=max(swe,na.rm=T))
sweboxplot_max=dcast(sweboxplot_max,date+site~.,value.var='maxswe',fun.agg=max)
colnames(sweboxplot_max)=c('date','site','maxvalue')
sweboxplot_max

sweerror=sweerror %>%
     filter(grepl('PHV',model.type)) %>%
     mutate(surveytype=ifelse(site=='Green Lakes Valley','alpine','forested')) %>%
     group_by(surveytype,model.type) %>%
     summarise(sweavg=mean(swe.obs.avg),
               erroravg=mean(objcost))
ggplot(sweboxplot,aes(x=source,y=swe))+
     geom_boxplot(outlier.colour = NA)+
     geom_point(data=sweavg,aes(x=model.type,y=sweavg),colour='blue')
     geom_point(data=sweerror,aes(x=model.type,y=erroravg),colour='red')+
     facet_wrap(~surveytype)

sweavg=subset(sweavg,source!='RCN' & source!='SNODAS')
sweerror=subset(sweerror,source!='RCN' & source!='SNODAS')
sweerror=arrange(sweerror,source,date,site)
sweerrorbest=reshape2::dcast(sweerror,date+site~source,value.var='objcost')
sweerrorbest$best=ifelse(sweerrorbest$PHV>sweerrorbest$PHVRCN,'PHVRCN','PHV')
sweerrorbest=merge(sweerrorbest,sweboxplot_max,by=c('date','site'))
sweerrorbest
sweerror$snotelbias=meromybias$bias.pct
# sweerror
# dcast(sweerror,source~.,value.var='objcost',fun.agg=mean,na.rm=T)


## boxplots ----
sweboxplot2=sweboxplot %>%
     mutate(surveytype=ifelse(site=='Green Lakes Valley','alpine','forested'))
sweavg2=sweavg %>%
     mutate(surveytype=ifelse(site=='Green Lakes Valley','alpine','forested'))
sweerror2=sweerror %>%
     mutate(surveytype=ifelse(site=='Green Lakes Valley','alpine','forested'))
sweerrorbest2=sweerrorbest %>%
     mutate(surveytype=ifelse(site=='Green Lakes Valley','alpine','forested'))

# nwt=subset(sweavg,site=='Niwot')

st='forested'
makeBoxplot=function(st,sweboxplot2,sweavg2,sweerror2,swerrorbest2){
     if(st=='alpine') mypal='#FFB90F'
     if(st=='forested') {
          mypal=c(brewer.pal(9,'Set1'),'black')
          mypal[5]='#EE7600'
          mypal=mypal[-6]
}
     sweboxplot=filter_(sweboxplot2,~surveytype==st)
     sweavg=filter_(sweavg2,~surveytype==st)
     sweerror=filter_(sweerror2,~surveytype==st)
     sweerrorbest=filter_(sweerrorbest2,~surveytype==st)
g=ggplot()+
     geom_boxplot(data=sweboxplot,aes(x=source,y=swe),outlier.shape=NA)+
     geom_point(data=sweavg,aes(x=source,y=avg,colour=site),size=3)+
     geom_point(data=sweerror,aes(x=source,y=objcost,colour=site),shape=18,size=3)+
     # facet_grid(mth~yr,scales='free_y')+
     facet_wrap(~date+site,scales='free_y',ncol=5)+
     geom_hline(yintercept = 0)+
     scale_x_discrete(labels=c('OBS','PHV\nbaseline','PHV\nRCN'))+
     scale_y_continuous(expand=c(0,0))+
     scale_colour_manual(values=mypal)+
     geom_text(data=sweerrorbest,aes(x=best,y=(maxvalue*0.9),label='*'),size=6)+
     #     geom_text(data=meromybias,aes(x=1,y=0.05,label=bias.pct),size=4)+
     labs(x='',y='SWE [m]')+#,title=paste(rswe,covrange,scalesnotel,fscaMatch))+
     guides(colour=guide_legend('Site',ncol=2,override.aes=list(shape=19,size=5)))+
     theme_minimal(base_size=14)+
     expand_limits(y=0)+
     expand_limits(y=0)+
     theme(legend.position=c(1,-.01),
           legend.justification=c(1,0),
           legend.text=element_text(size=8),
           strip.text=element_text(face='bold'),
           axis.line=element_line(size=.5),
           axis.text.x=element_text(angle=45,hjust=1,vjust=1),
           panel.border=element_rect(fill=NA,colour='grey50'))

return(g)
}
rbind_gtable_max <- function(...){

     gtl <- list(...)
     stopifnot(all(sapply(gtl, is.gtable)))
     bind2 <- function (x, y)
     {
          stopifnot(ncol(x) == ncol(y))
          if (nrow(x) == 0)
               return(y)
          if (nrow(y) == 0)
               return(x)
          y$layout$t <- y$layout$t + nrow(x)
          y$layout$b <- y$layout$b + nrow(x)
          x$layout <- rbind(x$layout, y$layout)
          x$heights <- gtable:::insert.unit(x$heights, y$heights)
          x$rownames <- c(x$rownames, y$rownames)
          x$widths <- grid::unit.pmax(x$widths, y$widths)
          x$grobs <- append(x$grobs, y$grobs)
          x
     }
     Reduce(bind2, gtl)
}

st='forested'
gf=makeBoxplot(st,sweboxplot2,sweavg2,sweerror2,swerrorbest2)
st='alpine'
ga=makeBoxplot(st,sweboxplot2,sweavg2,sweerror2,swerrorbest2)

gfgrob=ggplotGrob(gf)#+theme(legend.position='none'))
gagrob=ggplotGrob(ga+theme(legend.position='none'))

# g=rbind(gagrob,gfgrob,size='first')
# g$heights
# g$widths=unit.pmax(gagrob$widths,gfgrob$widths)
# grid.newpage()
# grid.draw(g)

pdf('test4.pdf',height=350/25.4,width=250/25.4)
grid.draw(rbind_gtable_max(gagrob, gfgrob))
dev.off()

# ggsave('test.pdf',width=12,height=10)
# maxwidth <- grid::unit.pmax(gfgrob$widths[2:5],gagrob$widths[2:5])
# gfgrob$widths[2:5]=as.list(maxwidth)
# gagrob$widths[2:5]=as.list(maxwidth)
# grid.arrange(gagrob,gfgrob,nrow=2)
#
# ggsave(plot=g,paste0('test_',st,'.pdf'),width=12,height=10)

# paste0('graphs/rswe_',rswe,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/boxplot_',blending,'models_vs_griddedobs_facet_datesite_',cost,'-rcnselect_noRCN_TEST.pdf'),width=18.6,height=10.1)



## the rest -----
glvswe=subset(cellavg,site=='Green Lakes Valley')
# glvswe
objf='rmse'
objcost=ddply(glvswe,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overallrmse=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overallrmse)[ncol(overallrmse)]=objf
overallrmse

# dcast(objcost,date+site+yr+mth~model.type,value.var='cost')
cost.df=dcast(subset(objcost,model.type=='PHV' | model.type=='PHVRCN'),model.type+site+date~.,value.var='cost')
colnames(cost.df)=c('model.type','site','date','avg_cost')
cost.df
bestcost=ddply(cost.df,.(site,date),function(dF){
  as.character(dF$model.type[which.min(dF$avg_cost)])
})
# bestcost
# tmpdf=cost.df[c(1,8),]
#   as.character(tmpdf$model.type[which.min(tmpdf$avg_cost)])
count(bestcost$V1)

# dcast(objcost,date+site+yr+mth~model.type,value.var='cost')
cost.df=dcast(subset(objcost,yr>2003),model.type+site+date~.,value.var='cost')
colnames(cost.df)=c('model.type','site','date','avg_cost')
cost.df
bestcost=ddply(cost.df,.(site,date),function(dF){
  as.character(dF$model.type[which.min(dF$avg_cost)])
})
bestcost
count(bestcost$V1)

glvswe=subset(cellavg,site=='Green Lakes Valley')
objf='relrmse'
objcost=ddply(glvswe,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overallrmse=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overallrmse)[ncol(overallrmse)]=objf
overallrmse

objf='pctbias'
glvswe=subset(cellavg,site=='Green Lakes Valley')
objcost=ddply(glvswe,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overall=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overall)[ncol(overall)]=objf
overall

forestswe=subset(cellavg,site!='Green Lakes Valley')
objf='rmse'
objcost=ddply(forestswe,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overallrmse=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overallrmse)[ncol(overallrmse)]=objf
overallrmse

forestswe=subset(cellavg,site!='Green Lakes Valley')
# forestswe=subset(cellavg,site=='Niwot' | site=='Green Lakes Valley' | site=='Dry Lake' | site=='Joe Wright' | site=='South Brush Creek')
objf='relrmse'
objcost=ddply(forestswe,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overallrmse=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overallrmse)[ncol(overallrmse)]=objf
overallrmse

cost.df=dcast(subset(objcost,model.type=='PHV' | model.type=='PHVRCN'),model.type+site+date~.,value.var='cost')
colnames(cost.df)=c('model.type','site','date','avg_cost')
bestcost=ddply(cost.df,.(site,date),function(dF){
  as.character(dF$model.type[which.min(dF$avg_cost)])
})
# bestcost
count(bestcost$V1)

objf='pctbias'
forestswe=subset(cellavg,site!='Green Lakes Valley')
forestswe=subset(cellavg,site=='Niwot' | site=='Green Lakes Valley' | site=='Dry Lake' | site=='Joe Wright' | site=='South Brush Creek')
objcost=ddply(forestswe,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overall=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overall)[ncol(overall)]=objf
overall

objf='rmse'
objcost=ddply(cellavg,.(model.type,date,site,yr,mth),function(dF){
    cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
    data.frame(cost)
})
# objcost
overallrmse=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overallrmse)[ncol(overallrmse)]=objf
overallrmse
# write.table(format(overallrmse,digits=1,trim=F,nsmall=2),paste0('graphs/rswe_',rswe,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/boxplot_overallrmse_',blending,'models_vs_griddedobs_facet-datesite_',cost,'-rcnselect.txt'),sep='\t',row.names=F,quote=F)

cost.df=dcast(subset(objcost,model.type=='PHV' | model.type=='PHVRCN'),model.type+site+date~.,value.var='cost')
colnames(cost.df)=c('model.type','site','date','avg_cost')
# cost.df
bestcost=ddply(cost.df,.(site,date),function(dF){
  as.character(dF$model.type[which.min(dF$avg_cost)])
})
# bestcost
count(bestcost$V1)

cost.df=dcast(subset(objcost,yr>2003),model.type+site+date~.,value.var='cost')
colnames(cost.df)=c('model.type','site','date','avg_cost')
#cost.df
bestcost=ddply(cost.df,.(site,date),function(dF){
  as.character(dF$model.type[which.min(dF$avg_cost)])
})
# bestcost
count(bestcost$V1)

objf='relrmse'
objcost=ddply(cellavg,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overallrmse=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overallrmse)[ncol(overallrmse)]=objf
overallrmse
#write.table(format(overallrmse,digits=1,trim=F,nsmall=2),paste0('graphs/rswe_',rswe,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/boxplot_overallrmse_',blending,'models_vs_griddedobs_facet-datesite_',cost,'-rcnselect.txt'),sep='\t',row.names=F,quote=F)

objf='bias'
objcost=ddply(cellavg,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
overall=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
# dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overall)[ncol(overall)]=objf
overall


objf='pctbias'
objcost=ddply(cellavg,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
dcast(objcost,date+site~model.type,value.var='cost')

overall=dcast(objcost,model.type~.,fun.aggregate=mean,na.rm=T,value.var='cost')
colnames(overall)[ncol(overall)]=objf
overall

dcast(objcost,yr+mth~model.type,value.var='cost',fun.agg=mean,na.rm=T)
# phv=c(-1,32,33,3)
# phvrcn=c(-7,25,24,-5)
# snotel=c(17,65,-3,23)
# cor(phv,snotel,method='spearman')
# cor(phvrcn,snotel,method='spearman')
dcast(objcost,site+date~model.type,value.var='cost')

#leave out survey station for mean bias because survey station wasn't used in regression
weeklysurveyavg=ddply(objcost,.(model.type,yr,mth),function(dF){
#print(dF)
#     str(dF)
    dFo=data.frame()
    for(s in unique(dF$site)){
  #      print(s)
        sid=which(dF$site==s)
 #       str(sid)
       # str(dF$cost)
        avg=mean(dF$cost[-sid],na.rm=T)
#          print(avg)
        dFo=rbind(dFo,data.frame(site=s,date=dF$date[sid],avg))
        }
#    print(dFo)
    return(dFo)
})
dcast(weeklysurveyavg,site+date~model.type,value.var='avg')
# weekavg=dcast(weeklysurveyavg,site+date+yr+mth~model.type,value.var='avg',fun.agg=mean,na.rm=T)
# weekavg$pattern=c(rep('a',15),'a1','a','m','ac','m','a','m','ac','m','a','m','a','m','ac')
# mean(weekavg$PHV-weekavg$PHVRCN,na.rm=T)
# data.frame(mean(weekavg$PHV,na.rm=T),mean(weekavg$PHVRCN,na.rm=T))
# weekavg

dcast(weeklysurveyavg,yr+mth~model.type,value.var='avg',fun.agg=mean,na.rm=T)

ddply(weekavg,.(yr,pattern),function(dF){
    summarise(dF,
              avg_phv=mean(PHV,na.rm=T),
              avg_phvrcn=mean(PHVRCN,na.rm=T),
              avg_rcn=mean(RCN,na.rm=T),
              avg_snodas=mean(SNODAS,na.rm=T)
              )})

objf='pctbias'
objcost=ddply(cellavg,.(model.type,date,site,yr,mth),function(dF){
  cost=docost(dF$swe.obs.avg,dF$swe.model,objf)
  data.frame(cost)
})
# print(objcost)
swebias=dcast(objcost,site+date~model.type,value.var='cost')
swebias
summarise(swebias,count(PHV>0),count(PHVRCN>0))
summarise(swebias,mean(PHV),mean(PHVRCN))

# average rmse of each model for each survey
sweerror=ddply(cellavg,.(model.type,site,date,mth,yr),function(dF) with(dF,docost(swe.obs.avg,swe.model,'rmse')))
names(sweerror)=c('source','site','date','mth','yr','objcost')
# combine modeled and observed swe for boxplots
obs=swe[,c('site','date','yr','mth','swe.obs')]
obs$source='OBS'
mdl=swe[,c('site','date','yr','mth','swe.model','model.type')]
names(obs) <- names(mdl) <- c('site','date','yr','mth','swe','source')
sweboxplot=rbind(obs,mdl)
# calculate average swe of each source for each survey
sweavg=ddply(sweboxplot,.(site,date,mth,yr,source),function(dF) summarise(dF, avg=mean(swe,na.rm=T)))

tmp=dcast(sweavg,site+date~source,value.var='avg')
docost(tmp$OBS,tmp$PHV,'rmse')
docost(tmp$OBS,tmp$PHVRCN,'rmse')
docost(tmp$OBS,tmp$SNODAS,'rmse')

meromybias=read.csv('data/spatialvalidation/snow_survey_table.csv',stringsAsFactor=F)
colnames(meromybias)=c('site','date','survey.swe.avg.m','snotel.swe.m','bias.pct')
meromybias$date=as.POSIXct(strptime(as.character(meromybias$date),'%m/%d/%Y'))
meromybias=arrange(meromybias,date,site)
data.frame(med=median(meromybias$bias.pct,na.rm=T),avg=mean(meromybias$bias.pct,na.rm=T))

mb.avg=dcast(meromybias,site~.,value.var='bias.pct',fun.agg=sd,na.rm=T)
colnames(mb.avg)=c('site','bias.pct.sd')
mb.avg
mean(mb.avg$bias.pct.sd,na.rm=T)

sweboxplot=subset(sweboxplot,source!='RCN' & source!='SNODAS')
sweboxplot_max=ddply(sweboxplot,.(date,site,source),summarise,maxswe=max(swe,na.rm=T))
sweboxplot_max=dcast(sweboxplot_max,date+site~.,value.var='maxswe',fun.agg=max)
colnames(sweboxplot_max)=c('date','site','maxvalue')
sweboxplot_max

sweavg=subset(sweavg,source!='RCN' & source!='SNODAS')
sweerror=subset(sweerror,source!='RCN' & source!='SNODAS')
sweerror=arrange(sweerror,source,date,site)
sweerrorbest=dcast(sweerror,date+site~source,value.var='objcost')
sweerrorbest$best=ifelse(sweerrorbest$PHV>sweerrorbest$PHVRCN,'PHVRCN','PHV')
sweerrorbest=merge(sweerrorbest,sweboxplot_max,by=c('date','site'))
sweerrorbest
sweerror$snotelbias=meromybias$bias.pct
# sweerror
# dcast(sweerror,source~.,value.var='objcost',fun.agg=mean,na.rm=T)

sweboxplot=sweboxplot %>% mutate(surveytype=ifelse(site=='Green Lakes Valley','alpine','forested'))


head(sweboxplot)

# options(repr.plot.width = 18.6, repr.plot.height = 10.1, jupyter.plot_mimetypes = 'image/png')
mypal=c(brewer.pal(9,'Set1'),'#000000')
mypal[5]='#EE7600'
mypal[6]='#FFB90F'
nwt=subset(sweavg,site=='Niwot')
g=ggplot()+
    geom_boxplot(data=sweboxplot,aes(x=source,y=swe),outlier.shape=NA)+
    geom_point(data=sweavg,aes(x=source,y=avg,colour=site),size=3)+
    geom_point(data=sweerror,aes(x=source,y=objcost,colour=site),shape=18,size=3)+
    geom_hline(yintercept = 0)+
    scale_x_discrete(labels=c('OBS','PHV-baseline','PHV-RCN'))+
    scale_y_continuous(expand=c(0,0))+
    scale_colour_manual(values=mypal)+
    geom_text(data=sweerrorbest,aes(x=best,y=(maxvalue-0.05),label='*'))+
#     geom_text(data=meromybias,aes(x=1,y=0.05,label=bias.pct),size=4)+
    labs(x='',y='SWE [m]')+#,title=paste(rswe,covrange,scalesnotel,fscaMatch))+
    guides(colour=guide_legend('Site',ncol=2,override.aes=list(shape=19,size=5)))+
    theme_minimal()+
    expand_limits(y=0)+
    facet_wrap(~date+site,scales='free_y')+
    expand_limits(y=0)+
    theme(legend.position=c(1,-.01),
          legend.justification=c(1,0),
          legend.text=element_text(size=8),
          strip.text=element_text(face='bold'),
          axis.line=element_line(size=.5),
          panel.border=element_rect(fill=NA,colour='grey50'))
# show(g)

ggsave(plot=g,paste0('graphs/rswe_',rswe,'/covrange',covrange,'/snotel',scalesnotel,'/',dateflag,'/',fscaMatch,'/',style,'/boxplot_',blending,'models_vs_griddedobs_facet_datesite_',cost,'-rcnselect_noRCN_TEST.pdf'),width=18.6,height=10.1)

##plot error vs bias
options(repr.plot.width = 18.6, repr.plot.height = 10.1, jupyter.plot_mimetypes = 'image/png')
mypal=c(brewer.pal(9,'Set1'),'#000000')
ggplot(sweerror)+
    geom_point(aes(snotelbias,objcost,colour=site),size=5)+
    scale_colour_manual(values=mypal,drop=TRUE)+
    coord_fixed(ratio=1000,ylim=c(0,0.3))+
    facet_wrap(~source)+
    theme_bw(base_size=24)


