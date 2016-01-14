setwd('/Volumes/Dominik/Documents/snotel-regression_project')
library('ProjectTemplate')
load.project()


recon.version='v3.1'
covrange='idp1'
cost='r2'
product='phvrcn'
fscaMatch='fsca'
snotelscale='noscale'
residblending='unblended'
if(residblending=='unblended') {
		resid='' 
	} else {
		resid='full'
	}
config=''#bic-v2.5-removed_NWbdiff'#'bic-v2.2-nosnoteltransform'
dateflag='B'
style='reanalysis'#'real-time'
dte=as.POSIXct('1900-04-01',tz='MST')

watermask=raster('data/cs_NHD_MOD44_water_mask.tif')
pathout='/Volumes/shareProjects/WSC/data/swe_regression'
pathin='/Volumes/Dominik/Documents/snotel-regression_project/diagnostics'
parampath=paste0('rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag,'/fullpreds/',cost,'/netcdf/',style,'/',config,'/')
fordenstack=stack('data/forestdensity/umd_forden.nc')
names(fordenstack)=seq(2000,2012)
#
Getfsca=function(dte){
if( as.numeric(strftime(dte,'%m')) >= 3){
     rfn=paste0('data/recon_',recon.version,'/reconfscadata_',yr,'_',recon.version,'.nc')
     ncstack=stack(rfn)
     stckdate=strftime(dte,'X%Y%m%d')
    rind=grep(stckdate,names(ncstack))
    vals=ncstack[[rind]]
    } else {
     print('using modscag fsca for snotel pixel sca. may not be completely cloudfree')
     sfn=file.path('data','selectdates','modscag',paste0('fsca',yr,'.nc'))
     ncstack=stack(sfn)
     stckdate=strftime(dte,'X%Y%j')
    rind=grep(stckdate,names(ncstack))
    vals=tryCatch({ncstack[[rind]]},error=function(x) {#this error should never happen
    print(paste0('the date ',dte,' was not available in the selected sca image'))})
    vals=shift(vals,x=0.00416666667/2,y=-0.00416666667/2) ## raw modscag images were incorrectly written with 1/2 pixel shift
    valind=vals[vals<=100]
    vals[valind]=vals[valind]/100
  }
  return(vals)
}

yr=2001
registerDoMC(3)
for(yr in seq(2001,2012),.parallel=T,function(yr){
	print(yr)

	predfile=file.path(pathin,parampath,paste0(resid,'preds-',product,'_',yr,'_blend-',fscaMatch,'.nc'))
	swestack=stack(predfile)
	datevec=as.POSIXct(strptime(names(swestack),'X%Y%m%d'),tz='MST')
	
	datei=1
	for(datei in 1:length(datevec)){
		dte=datevec[datei]
		fsca=Getfsca(dte)
		forden=fordenstack[[grep(yr,names(fordenstack))]]
		forden[forden==1]=.999
		fsca_corr=fsca
		fsca_corr[fsca_corr>1]=NA
		fsca_corr=fsca_corr/(1-forden)
		fsca_corr[fsca_corr>1]=1
		sweraster=swestack[[datei]]
		sweraster=sweraster*fsca_corr
		sweraster[fsca>=200]=fsca[fsca>=200]
		sweraster=mask(sweraster,watermask)

	}
	layerid=paste0('X',yr-2011+12)
	layernum=grep(layerid,names(fordenstack))
	forden=fordenstack[[layernum]]
	
	
	
		
	
	s=stack(predfile)
	band=grep(paste0(yr,mth,dy),names(s))
	r=s[[band]]
	r[r==253]=NA
	if(snotelscale=='scale'){
		r=r*fsca
	}
	#
	r.ext=as.data.frame(zonal(r,huc4raster,'mean'))#swe avg is meters
	colnames(r.ext)=c('zone','sweavg')
	r.ext$vol=r.ext$sweavg*basinsize$count*500*500/10^9#cu km
	if(units=='metric'){
		basin_stats=merge(basin_static,r.ext,by='zone')
		swedF=data.frame(date=as.POSIXct(paste0(yr,'-',mth,'-',dy),tz='MST'),yr=yr,model=product,basin_stats)
	}
	swetmp=rbind(swetmp,swedF)
}
return(swetmp)
})