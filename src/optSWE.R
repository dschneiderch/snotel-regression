setwd('~/GoogleDrive/snotel-regression_project')
optdata=read.table('diagnostics/rswe_v3.2/covrange300km/fullpreds/xval/real-time_snotel_xval_bestpreds_r2.txt',header=T,stringsAsFactors=F)
optdata$date=as.POSIXct(optdata$date,tz='MST')
optdata$recondate=as.POSIXct(optdata$recondate,tz='MST')

cost='rmse'
if(cost=='cor') { #statistics to maximize
          costfun<-function(y,yhat) cor(y,yhat)
          flag='max'
     } else if(cost=='r2') {     
          costfun<-function(y,yhat) 1-sum((y-yhat)^2)/sum((y-mean(y,na.rm=T))^2)
          flag='max'
     } else if(cost=='mae') { #statistics to minimize
          costfun<-function(y, yhat) mean(abs(yhat-y),na.rm=T)
          flag='min'
     } else if(cost=='rmse') {
          costfun<-function(y,yhat) sqrt(mean((yhat-y)^2,na.rm=T))
          flag='min'
     }

library(plyr)
library(ggplot2)
optskill=ddply(optdata,.(yrdoy),function(dF){
	skill=costfun(dF$snotel,dF$reconopt)
	relskill=skill/mean(dF$snotel)
	dte=unique(dF$date)
	yr=strftime(dte,'%Y')
	data.frame(date=dte,yr,skill,relskill)
})


ggplot(optskill)+
geom_point(aes(x=date,y=relskill))+
coord_cartesian(ylim=c(0,2))+
facet_grid(~yr,scale='free_x')

ggsave('graphs/rswe_v3.2/covrange300km/scatterplot_Relrmse_real-time_byyear_optrecon.pdf')