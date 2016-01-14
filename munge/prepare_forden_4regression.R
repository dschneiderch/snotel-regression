setwd('~/Documents/snotel-regression_project')
library('ProjectTemplate')
load.project()

pn='/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/model_code/input_static/UpperColoradoRiver'
fid=file(file.path(pn,'forden_nlcd2006.dat'),'rb')
nlcd=readgrads(1,fid)
writeRaster(nlcd,'data/forestdensity/nlcd_forden.nc')

forden=stack()
yr=2000
fns=list.files(pn,full.names=T,pattern=glob2rx('forden*umd.dat$'))
for(x in fns){
	fid=file(x,'rb')
	r=readgrads(1,fid)
	names(r)=yr
	forden=addLayer(forden,r)
	yr=yr+1
}
forden
writeRaster(forden,'data/forestdensity/umd_forden.nc')
