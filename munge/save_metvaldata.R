library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')

sites=c('GLEES','Niwot_Ridge','Valles_Caldera_Mixed_Conifer')#c('Chimney_Park','GLEES','Niwot_Ridge','Valles_Caldera_Mixed_Conifer','Valles_Caldera_Ponderosa_Pine')
pn='data/metvalidation/Ameriflux/Level2-with_gaps'

metdata=ldply(sites,function(sitename){
sitepn=file.path(pn,sitename)
if(sitename=='Niwot_Ridge'){
fn=file.path(sitepn,paste0('AMF_USNR1_',seq(2000,2012),'_L2_WG_V008.csv'))
} else if(sitename=='GLEES'){
fn=file.path(sitepn,paste0('AMF_USGLE_',seq(2004,2012),'_L2_WG_V005.csv'))
} else if(sitename=='Valles_Caldera_Mixed_Conifer'){
fn=file.path(sitepn,paste0('AMF_USVcm_',seq(2007,2010),'_L2_WG_V003.csv'))
}
csvdata=ldply(fn,function(f){
print(sitename)
all_content = readLines(f)
skip_second = all_content[-c(1:17,19:20)]
csvdata=read.csv(textConnection(skip_second), header = TRUE, stringsAsFactors = FALSE)
	# read.csv(f,skip=17,header=T,stringsAsFactor=F)
})
csvdata$Site=sitename
return(csvdata)
})
metdata[metdata==-6999]=NA
metdata[metdata==-9999]=NA
str(metdata)
metdata$date=as.POSIXct(strptime(with(metdata,paste0(YEAR,sprintf('%03d',DOY),sprintf('%04d',HRMIN))),'%Y%j%H%M'),tz='MST')
metdata$hour=strftime(metdata$date,'%Y%j%H')
metdata=metdata[,c('date','hour','Site','H','LE','Rg','RgOut','Rgl','RglOut','RH','TA')]
metdata.m=melt(metdata,id=c('date','hour','Site'))
metdata2=ddply(metdata.m,.(hour,variable,Site),function(x){
	summarise(x,value_avg=mean(value,na.rm=T))
})
metavg=dcast(metdata2,hour+Site~variable)	

write.csv(metavg,'data/metvalidation/Ameriflux/Level2-with_gaps/Niwot_Glees_VMC.csv',row.names=F, quote=F)

