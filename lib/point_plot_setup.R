point_plot_setup=function(recon.version){
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


## Add snotel swe measurements to snotel prediction dF
snotelrec=arrange(snotelrec,Station_ID,date)
fnms=list.files(paste0('diagnostics/rswe_',recon.version,'/fullpreds/snotellocs'),pattern='.txt$',full.name=T)
snotelfullpreds=ldply(fnms,function(f) read.table(f,header=T))
snotelfullpreds$date=as.POSIXct(strptime(snotelfullpreds$yrdoy,'%Y%j',tz='MST'))
snotelfullpreds=arrange(snotelfullpreds,Station_ID,date)
snotelfullpreds=ddply( snotelfullpreds,.(date),function(dF){
	ind=which(with(snotelrec,date==unique(dF$date)))
	 dF$snotel=snotelrec[ind,'swe']
	 return(dF)
})

snotellocs.df=arrange(as.data.frame(snotellocs),Station_ID)
swediffmap=mutate(snotelfullpreds,
	swediff.phv=phv.fullpred-snotel,
	swediff.phvrcn=phvrcn.fullpred-snotel,
	swediff.recon=recon-snotel,
	swediffpct.phv=swediff.phv/snotel,
	swediffpct.phvrcn=swediff.phvrcn/snotel,
	swediffpct.recon=swediff.recon/snotel,
	lat=rep(snotellocs.df$y,length.out=nrow(snotelfullpreds)),
	long=rep(snotellocs.df$x,length.out=nrow(snotelfullpreds)))
swediffmap$yr=strftime(swediffmap$date,'%Y')
swediffmap$mnth=strftime(swediffmap$date,'%b')
swediffmap$mnthdy=strftime(swediffmap$date,'%b%d')
return(swediffmap)
}