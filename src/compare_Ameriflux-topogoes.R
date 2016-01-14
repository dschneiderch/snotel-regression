library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project()

valdata=read.csv('data/metvalidation/Ameriflux/Level2-with_gaps/Niwot_Glees.csv',stringsAsFactor=F)
valdata$date=as.POSIXct(strptime(valdata$hour,'%Y%j%H'))
valdata$siteid=ifelse(valdata$Site=='Niwot_Ridge','US-NR1',
	ifelse(valdata$Site=='GLEES','US-GLE',NA))

modelpath='/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/recon/v3.1_senlat.nldas.usace.modscag.umdforden/met_output/'
yrs=seq(2000,2012)
snodismet=ldply(yrs,function(yr){
	fn=file.path(modelpath,paste0('topogoes_eblocs_',yr,'.txt'))
	dat=read.table(fn,header=T,sep='\t')
	dat2=subset(dat,siteid=='US-NR1' | siteid=='US-GLE')
	return(dat2)
})
snodismet$date=as.POSIXct(strptime(snodismet$date,'%Y-%m-%d %H:%M',tz='MST'))

obsdF=valdata[,c('siteid','date','Rg')]
obsdF$type='obs'

mdldF=snodismet[,c('siteid','date','value')]
colnames(mdldF)=c('siteid','date','Rg')
mdldF$type='mdl'
globalR=rbind(obsdF,mdldF)
globalR$Rg=as.numeric(globalR$Rg)
globalR$yr=as.numeric(strftime(globalR$date,'%Y'))
globalR$mth=as.numeric(strftime(globalR$date,'%m'))

globalR.plot=subset(globalR, siteid=='US-GLE' & Rg!=0 & mth<7)
ggplot(globalR.plot)+
	geom_point(aes(x=date,y=Rg,colour=type))+
	# scale_x_continuous(limits=c(min(mdldF$date),max(mdldF$date)))+
	facet_wrap(~yr,scales='free_x')

gRp=dcast(globalR.plot,date+yr~type,mean,value.var='Rg',na.rm=T)
dev.new()
ggplot(gRp)+
	geom_point(aes(x=obs,y=mdl))+
	geom_abline()+
	facet_wrap(~yr)