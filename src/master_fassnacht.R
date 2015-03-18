library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
setwd('~/Documents/snotel-regression_project')
library(doMC)
load.project()

recon.version='v3.2'
ryrs=2000#seq(2000,2012) ## need at least 1 year to keep things simpler
snotelrecon=get_snotelrecondata(recon.version,ryrs)

# generate swedata dataframe
swedata=fassnacht_swedata(snotelrecon,'qc1',precipmeth=1)
#subset snotel locations to those we have.
snotellocs=snotellocs[snotellocs$Station_ID %in% unique(swedata$Station_ID),]

# setup parallel processing -----------
parfun=function(){
  setup_parallel=function(){
    getnodes <- function() {
       f <- Sys.getenv('PBS_NODEFILE')
       x <- if (nzchar(f)) readLines(f) else rep('localhost', 2)
       as.data.frame(table(x), stringsAsFactors=FALSE)
      }
  nodes <- getnodes()
  cl <- makeSOCKcluster(nodes$x)

print(cl)
print(nodes)

registerDoSNOW(cl)

setcores <- function(cl, nodes) {
   cores=nodes$Freq
   f <- function(cores) assign('allocated.cores', cores, envir=.GlobalEnv)
   clusterApply(cl, nodes$Freq, f)
 }
setcores(cl, nodes)
}
}

# get modelling data --------------------------------------------------------
mdldata=swedata[swedata$mth<6,]#if changing mth<6, need to change dates2model too
mdldata$recondate=mdldata$date

# select dates to model ----------
# Option A will only model dates for dates selected from modscag images.
# Option B will model selected dates from modscag for months prior to March and then 1st, 8th, 15th, 22nd of March, April, May
dateselect=dates2model(opt='F')
ind=mdldata$date %in% dateselect
doidata=mdldata[ind,]

# select recon dates to evaluate model with
#recondata is used to iterate through recondates. 
#these dates can be different than olddata/mdldata if wanted but must include at least 1 date before each mdl date if used in real-time mode
recondata=subset(mdldata, (dy==1) )

# run model -------------------
cost='rmse'#cor, r2, mae, rmse
style='reanalysis'#real-time'#'real-time','reanalysis'  #reanalysis only. should have at least 1 year in ryrs above
spatialblend='blend'#flag to do geostatistical blending or not (prediction stage only; always does in the CV stage). blending takes long time for the entire domain..
output='points'#points' #'surface' #just predict at snotel pixels #for 'points' spatialblend must also be 'blend'

  # dte=strptime('20010524','%Y%m%d','MST')
  # ldte=strptime('20020202','%Y%m%d','MST')
  # doidata=subset(doidata,date>=dte)
  # doidata=subset(doidata,date<ldte)
#find which recondate gives best estimate.
#output dF of GLobal Moran I and objective functions for best model estimates for each yrdoy. do NOT parallelize
doylist=dlply(doidata,.(yrdoy),doDOYfit,cost,style,.parallel=F, .drop=F,.inform=F)
cleandF=function(dF){
     mutate(dF,
            yrdoy=attr(doylist,'split_labels')$yrdoy,
            date=as.POSIXct(strptime(yrdoy,'%Y%j','MST')),
            yr=strftime(date,'%Y'))
}
which_recon_date=cleandF(ldply(doylist,'[[',1))
moran.df=cleandF(ldply(doylist,'[[',2))
xvalpreds=ldply(doylist,'[[',3,.inform=F)

write.table(which_recon_date,paste0('diagnostics/rswe_',recon.version,'/fassnacht/',style,'_recondate_selection_',cost,'.txt'),sep='\t',row.names=F,quote=F)     
write.table(moran.df,paste0('diagnostics/rswe_',recon.version,'/fassnacht/',style,'_moran_info_for_recondate_selection_',cost,'.txt'),sep='\t',row.names=F,quote=F)     
write.table(xvalpreds,paste0('diagnostics/rswe_',recon.version,'/fassnacht/',style,'_snotel_xval_bestpreds_',cost,'.txt'),sep='\t',row.names=F,quote=F)
