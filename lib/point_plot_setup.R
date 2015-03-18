point_plot_setup=function(recon.version,style,cost,covrange){
# # Map domain and snotel stations
# #*** For Upper CO River Basin domain
upleft=c(43.75,-112.25)
lowright=c(33,-104.125)
lon.low=-112.25
lon.high=-104.125
lat.low=33
lat.high=43.75
lim=c(lon.low, lon.high, lat.low, lat.high)

## Get state boundaries
assign('states',map_data('state'),envir=.GlobalEnv)

## Get shaded relief for background
ucorelief=raster('data/shaded.relief.grey/shadedrelief_uco.tif')
ucorelief.agg=aggregate(ucorelief,2)
ucorelief.df=data.frame(rasterToPoints(ucorelief.agg))
names(ucorelief.df)=c('long','lat','ter')
assign('ucorelief.df',ucorelief.df,envir=.GlobalEnv)

snotelxval=read.table(paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/fullpreds/xval/',style,'_','snotel_xval_bestpreds_',cost,'.txt'),header=T)
snotelxval=snotelxval[!is.na(snotelxval$Station_ID),]
snotelxval$date=as.POSIXct(snotelxval$date,tz='MST')
snotelxval$recondate=as.POSIXct(snotelxval$recondate,tz='MST')
snotelxval=arrange(snotelxval,date,Station_ID)

snotellocs.df=as.data.frame(snotellocs)[snotellocs$Station_ID %in% unique(snotelxval$Station_ID),]

snotelxval=ddply(snotelxval,.(Station_ID),function(dF){ 
	sid=unique(dF$Station_ID)
	dF$lat=snotellocs.df[snotellocs.df$Station_ID==sid,'y']
	dF$long=snotellocs.df[snotellocs.df$Station_ID==sid,'x']
	return(dF)
})


swediffmap=mutate(snotelxval,
	swediff.phv=phv.pred-snotel,
	swediff.phvrcn=phvrcn.pred-snotel,
	swediff.reconrt=reconrt-snotel,
	swediff.phvfull=phv.fullpred-snotel,
	swediff.phvrcnfull=phvrcn.fullpred-snotel,
	swediff.reconrtfull=reconrt.fullpred-snotel,
	swediffpct.phv=swediff.phv/snotel,
	swediffpct.phvrcn=swediff.phvrcn/snotel,
	swediffpct.reconrt=swediff.reconrt/snotel,
	swediffpct.phvfull=swediff.phvfull/snotel,
	swediffpct.phvrcnfull=swediff.phvrcnfull/snotel,
	swediffpct.reconrtfull=swediff.reconrtfull/snotel)
swediffmap$yr=strftime(swediffmap$date,'%Y')
swediffmap$mnth=strftime(swediffmap$date,'%b')
swediffmap$mnthdy=strftime(swediffmap$date,'%b%d')


## **** this other stuff was for the origional snotelxval  output. originally output every recondate-xval combo but its very big and slow. limiting to the best makes more sense i think.
# ## Add snotel swe measurements to snotel prediction dF
# snotelrec=arrange(snotelrec,Station_ID,date)
# snotelfullpreds=read.table(paste0('diagnostics/rswe_',recon.version,'/fullpreds/xval/',style,'_','snotel_xval_preds_',cost,'.txt'),header=T)
# snotelfullpreds$date=as.POSIXct(strptime(snotelfullpreds$yrdoy,'%Y%j',tz='MST'))
# snotelfullpreds$recondate=as.POSIXct(snotelfullpreds$recondate,tz='MST')
# snotelfullpreds=arrange(snotelfullpreds,Station_ID,date)
# snotelfullpreds=ddply( snotelfullpreds,.(date),function(dF){
# 	ind=which(with(snotelrec,date==unique(dF$date)))
# 	 dF$snotel=snotelrec[ind,'swe']
# 	 return(dF)
# })

# which_recon_date=read.table(paste0('diagnostics/rswe_',recon.version,'/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T,stringsAsFactors=F)
# which_recon_date$date=as.POSIXct(strptime(which_recon_date$date,'%Y-%m-%d',tz='MST'))
# which_recon_date$phvrcn_recondate=as.POSIXct(strptime(which_recon_date$phvrcn_recondate,'%Y-%m-%d',tz='MST'))
# which_recon_date$recon_costdate=as.POSIXct(strptime(which_recon_date$recon_costdate,'%Y-%m-%d',tz='MST'))
# # str(which_recon_date)

# test=ddply(snotelfullpreds,.(yrdoy),function(dF){
#   wr=subset(which_recon_date,date==unique(dF$date))
#   ind=which(dF$recondate %in% wr$phvrcn_recondate)
#   newdF=dF[ind,]
#   return(newdF)
# })
# test=arrange(test,date,Station_ID)

# # fnms=list.files(paste0('diagnostics/rswe_',recon.version,'/fullpreds/snotellocs'),pattern='.txt$',full.name=T)
# # snotelfullpreds=ldply(basename(fnms),function(f) read.table(f,header=T))
# # snotelfullpreds$date=as.POSIXct(strptime(snotelfullpreds$yrdoy,'%Y%j',tz='MST'))
# # snotelfullpreds=arrange(snotelfullpreds,Station_ID,date)
# # snotelfullpreds=ddply( snotelfullpreds,.(date),function(dF){
# # 	ind=which(with(snotelrec,date==unique(dF$date)))
# # 	 dF$snotel=snotelrec[ind,'swe']
# # 	 return(dF)
# # })

# snotellocs.df=arrange(as.data.frame(snotellocs),Station_ID)
# swediffmap=mutate(test,
# 	swediff.phv=phv.pred-snotel,
# 	swediff.phvrcn=phvrcn.pred-snotel,
# 	swediff.recon=recon-snotel,
# 	swediff.phvfull=phv.fullpred-snotel,
# 	swediff.phvrcnfull=phvrcn.fullpred-snotel,
# 	swediff.reconfull=recon.fullpred-snotel,
# 	swediffpct.phv=swediff.phv/snotel,
# 	swediffpct.phvrcn=swediff.phvrcn/snotel,
# 	swediffpct.recon=swediff.recon/snotel,
# 	swediffpct.phvfull=swediff.phvfull/snotel,
# 	swediffpct.phvrcnfull=swediff.phvrcnfull/snotel,
# 	swediffpct.reconfull=swediff.reconfull/snotel,
# 	lat=rep(snotellocs.df$y,length.out=nrow(test)),
# 	long=rep(snotellocs.df$x,length.out=nrow(test)))
# swediffmap$yr=strftime(swediffmap$date,'%Y')
# swediffmap$mnth=strftime(swediffmap$date,'%b')
# swediffmap$mnthdy=strftime(swediffmap$date,'%b%d')


return(swediffmap)
}