blend_recon=function(rdtsdF,newdatalocs.usgs,newdatalocs.agg.usgs,spatialblend,recon.version,covrange){

print(paste('**',unique(rdtsdF$yr),'**'))	
basepath=paste0('diagnostics/rswe_',recon.version,'/fullpreds')

get_reconraster=function(dte){
    yr=strftime(dte,'%Y')
    rfn=paste0('data/recon_',recon.version,'/recondata_',yr,'_',recon.version,'.nc')
    ncstack=stack(rfn)
    stckdate=strftime(dte,'X%Y%m%d')
   	rind=grep(stckdate,names(ncstack))
   	rvals=tryCatch({getValues(ncstack[[rind]])},error=function(x) {#this error should never happen
   	print(paste0('the date ',dte,' was not available in the selected sca image'))})
	return(rvals)
}

#load water mask
	watermask=getValues(raster('data/cs_NHD_MOD44_water_mask.tif'))  
#setup netcdf dimensions
    dim1=ncdim_def('Long','degree',seq(-112.25,-104.125,0.00416666667))
    dim2=ncdim_def('Lat','degree',seq(43.75,33,-0.00416666667))
    dim3=ncdim_def('time','yrdoy',unlim=T,vals=as.numeric(unique(rdtsdF$yrdoy)))
    var=ncvar_def('swe','meters',dim=list(dim1,dim2,dim3),missval=-99,longname='snow water equivalent',compression=6)
#create netcdf file
    rcn.fn=file.path(bpth,'netcdf',paste0('recon_',yr,'_',spatialblend,'.nc'))
    rcn.nc=nc_create(rcn.fn,var)
#add projection info as global attribute
    ncatt_put(rcn.nc,0,'proj4string','+proj=longlat +datum=WGS84')

#blend each date and save into netcdf file
d_ply(rdtsdF,.(date),function(rdatedF){
	dte=unique(rdatedF$date)
	yr=unique(rdtsdF$yr)
	print(dte)

#get recon surface
	rvals=get_reconraster(dte)
	#get swe data of snotels, includes recon values. calculate residual
	rcndata=subset(swedata,date==dte)
	reconresid=with(rcndata, recon-snotel)

#fit kriging model to reconresid
	locsmat=as.matrix(as.data.frame(snotellocs.usgs)[,c('x','y')])
    rlambda=summary(Krig(locsmat,reconresid,theta=covrange,cov.function="wendland.cov"))$lambda
   	rresidkr=mKrig(locsmat,reconresid,theta=covrange,lambda=rlambda,cov.function="wendland.cov")

#predict geostat model at larger resolution and then interpolate bilinearly.
	xpred=as.data.frame(newdatalocs.agg.usgs)[,c('x','y')]
    grid.list=list(x=newdatapreds[,'x'],y=newdatapreds[,'y'])
    predagg=predict(rresidkr,x=xpred)
    obj=list(x=xpred[,'x'],y=xpred[,'y'],z=predagg)
    rcn.fullpred=rvals-interp.surface.grid(obj,grid.list)

#use order of yrdoy in the dataframes to add data by layer          
	tind=unique(rdatedF$yrdoy)         
    val=matrix((rcn.fullpred*watermask),nrow=2580,byrow=T)
    ncvar_put(rcn.nc,var,vals=val,start=c(1,1,which(tind==unique(rdtsdF$yrdoy))),count=c(-1,-1,1))
})
     nc_close(rcn.nc)
}