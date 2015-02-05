library(raster)
library(plyr)

recon.version='v3.2'
setwd(paste0('~/GoogleDrive/snotel-regression_project/data/recon_',recon.version))


yrs=as.list(seq(2000,2011))
l_ply(yrs,function(yr){
print(yr)
	load(paste0('recondata_',yr,'_',recon.version,'.RData'))
	rm(list=paste0('recon',yr,'df'))
	save(file=paste0('recondata_',yr,'_',recon.version,'.RData'),
         list=ls(pattern=glob2rx(pattern='recon2*')))
rm(list=ls(pattern=glob2rx(pattern='recon2*')))
}    ,.progress=T )