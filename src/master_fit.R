library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
#setwd('~/Documents/snotel-regression_project')
library(doMC)
load.project()

# select recon version ----------------------------------------------------
recon.version='v3.2'
ryrs=seq(2000,2012)
snotelrecon=get_snotelrecondata(recon.version,ryrs)

# generate swedata dataframe
swedata=generate_swedata(snotelrecon)
#subset snotel locations to those we have.
snotellocs=snotellocs[snotellocs$Station_ID %in% unique(swedata$Station_ID),]

# create opt lookup table
optlist=opt_LUT()
optLUT=optlist[[1]]
stascale=optlist[[2]]

#generate real-time optimized data
optdata=ddply(snotelrec,.(wy,Station_ID),function(dF){
  cddavg=stascale[stascale$Station_ID==unique(dF$Station_ID),'cddavg']
  cddstd=stascale[stascale$Station_ID==unique(dF$Station_ID),'cddstd']
  apcpavg=stascale[stascale$Station_ID==unique(dF$Station_ID),'apcpavg']
  apcpstd=stascale[stascale$Station_ID==unique(dF$Station_ID),'apcpstd']
  # degday=ifelse(dF$tavg<0,0,dF$tavg)
  degday=dF$tavg
  degday=ifelse(is.na(degday),0,degday)
  cumdegday=(cumsum(degday)-cddavg)/cddstd
  apcp=(dF$apcp-apcpavg)/apcpstd
  date=as.POSIXct(strftime(strptime(sprintf('%06d',dF$snoteldate),'%m%d%y'),'%Y-%m-%d'),tz='MST')
  yr=strftime(date,'%Y')
  doy=strftime(date,'%j')
  return(data.frame(date,yr,cumdegday,apcp,doy))
})

# # setup parallel processing -----------
# getnodes <- function() {
#        f <- Sys.getenv('PBS_NODEFILE')
#        x <- if (nzchar(f)) readLines(f) else rep('localhost', 2)
#        as.data.frame(table(x), stringsAsFactors=FALSE)
#       }
#   nodes <- getnodes()
#   cl <- makeSOCKcluster(nodes$x)

# print(cl)
# print(nodes)

# registerDoSNOW(cl)

# setcores <- function(cl, nodes) {
#    cores=nodes$Freq
#    f <- function(cores) assign('allocated.cores', cores, envir=.GlobalEnv)
#    clusterApply(cl, nodes$Freq, f)
#  }
# setcores(cl, nodes)
allocated.cores=2

# get modelling data --------------------------------------------------------
mdldata=swedata[swedata$mth<6,]#if changing mth<6, need to change dates2model too
mdldata$recondate=mdldata$date

# select dates to model ----------
# Option A will only model dates for dates selected from modscag images.
# Option B will model selected dates from modscag for months prior to March and then 1st, 8th, 15th, 22nd of March, April, May
dateselect=dates2model(opt='B')
ind=mdldata$date %in% dateselect
doidata=arrange(mdldata[ind,],date,Station_ID)#!!!important

# select recon dates to evaluate model with
#recondata is used to iterate through recondates. 
#these dates can be different than olddata/mdldata if wanted but must include at least 1 date before each mdl date if used in real-time mode
recondata=subset(mdldata, (dy==1 | dy==15) ) 

# run model -------------------
cost='r2'#cor, r2, mae, rmse
style='real-time'#real-time'#'real-time','reanalysis'
spatialblend='blend'#flag to do geostatistical blending or not (prediction stage only; always does in the CV stage). blending takes long time for the entire domain..
output='surface'#points' #'surface' #just predict at snotel pixels #for 'points' spatialblend must also be 'blend'
covrange='300km'

print(cost)
print(style)
print(spatialblend)
print(output)
print(covrange)


if(style=='real-time'){
     doidata=subset(doidata,yr>2000)
}
# dte=strptime('20010524','%Y%m%d','MST')
# ldte=strptime('20020202','%Y%m%d','MST')
# doidata=subset(doidata,date>=dte)
# doidata=subset(doidata,date<ldte)
#find which recondate gives best estimate.
#output dF of GLobal Moran I and objective functions for best model estimates for each yrdoy. do NOT parallelize
doidata=subset(doidata,yr==2002 & (mth==3))
doylist=dlply(doidata,.(yrdoy),doDOYfit,cost,style,.parallel=F, .drop=F,.inform=F)
cleandF=function(dF){
    mutate(dF,
            yrdoy=attr(doylist,'split_labels')$yrdoy,
            date=as.POSIXct(strptime(yrdoy,'%Y%j','MST')),
            yr=strftime(date,'%Y'))
}
which_recon_date=cleandF(ldply(doylist,'[[',1))
xvalpreds=ldply(doylist,'[[',2,.inform=F)
optweights=ldply(doylist,'[[',3)

write.table(which_recon_date,paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',row.names=F,quote=F)     
write.table(xvalpreds,paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/fullpreds/xval/',style,'_snotel_xval_bestpreds_',cost,'.txt'),sep='\t',row.names=F,quote=F)
write.table(optweights,paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/',style,'_optmodel_',cost,'.txt'),sep='\t',row.names=F,quote=F)




#quit(save='no')
