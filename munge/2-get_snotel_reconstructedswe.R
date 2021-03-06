library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project()

#: reads in swe.dat from recon output from selected version, saves usgs  projected stack and snotel values to recondata_v#.RData
## ----  Which Recon version? -----

for(varname in c('swe','fsca')){
# varname='fsca'
yrs=2000:2012
recon.version='v3.2'
reconbase='/Volumes/hydroProjects/SWE/Rockies/SWE_SNODIS/recon'

## ------

# set up out directory ------
pathout=paste0('data/recon_',recon.version)
if(is.na(file.info(pathout)$isdir))  dir.create(pathout,recursive=T)

get_recon_rasterstack <- function(yr,filename){
     tr=file(filename,'rb')
     s=stack()
     nx=1950#columns
     ny=2580#rows
     nt=184#days
     dsize=4

     dt=strptime(paste0(yr,'-Mar-01'),format='%Y-%B-%d',tz='MST')
     for(i in seq(1,nt)){#don't make this an apply function. list concatenation messes structure of output so extract doens't work. could stack(s)? no speed gain
          data=readBin(tr,what='numeric',size=dsize,n=nx*ny)#grid file grid by grid
          mat=matrix(data,nrow=ny,byrow=T)#convert to matrix with nrow
          mat=mat[nrow(mat):1,]#flip in y-direction to make right side map
          r=raster(mat)#rasters use column major cellnums instead row major, i.e. cell nums count up left to right, top to bottom instead of normal matrix way of top to bottom, left to right.
          extent(r) <- c(-112.25,-104.125,33,43.75)
          projection(r) <- '+proj=longlat +datum=WGS84'
          names(r) <- as.character(strftime(dt+24*3600*(i-1),format='%Y%m%d'))
          s=addLayer(s,r)
     }
close(tr)
return(s)
}

get_snotel_reconswe <- function(x){
     snotelrecon=as.data.frame(extract(x,snotellocs,sp=T))
     snotelrecon=melt(snotelrecon,id.vars=c(names(snotellocs),coordnames(snotellocs)),value.name='recon',variable.name='date')
     snotelrecon$date=strptime(as.character(snotelrecon$date),format='X%Y%m%d',tz='MST')
     return(snotelrecon)
}

fn=list.files(reconbase,full.name=T)
fileind=which(file.info(fn)$isdir)
recondir=fn[fileind][grep(recon.version,x=fn[fileind])]
#
fn=list.files(recondir,recursive=T,full.names=T,pattern=paste0(varname,'.dat$'))
#
## yrs=strsplit(fn,'/',fixed=T)
## pathlen=unlist(lapply(yrs,length))[1]
## yrs=as.numeric(unlist(lapply(yrs,'[',pathlen-1)) )
nx=0.0041666666667
ny=0.0041666666667
xdim=-112.25+nx/2+nx*(seq(1,1950)-1)
ydim=43.75-ny/2-ny*(seq(1,2580)-1)
dim1=ncdim_def('Long','degree',xdim)
dim2=ncdim_def('Lat','degree',ydim)

for(yr in yrs){
     print(yr)
     ifn=fn[grep(glob2rx(paste0('*/',yr,'/*')),fn)]
     recon.stack=get_recon_rasterstack(yr, ifn)
     #recon.stack=projectRaster(recon.stack,crs=CRS('+init=epsg:5070'),method='ngb')
     snotelrecon=get_snotel_reconswe(recon.stack)
     assign(paste0('snotelrecon',yr),snotelrecon)
     #
#      registerDoMC(3)#bug in aaply. parallel fails :(
#      mm=getValues(recon.stack)
#      vals=aaply(mm,2,.parallel=F,.inform=F,function(m){
#           mat=matrix(m,nrow=2580,byrow=T)#fill byrow because raster pkg
#           mat=mat[nrow(mat):1,]#invert to match netcdf dims
#           return(mat)
#      })
#      valwrite=aperm(vals,c(3,2,1))#permute to lat, long, time
valwrite=getValues(recon.stack)#this is easy. make sure that y dim is backwards for netcdf
    
    if(varname=='swe'){
#netcdf write path
fnnc=file.path(pathout,paste0('recondata_',yr,'_',recon.version,'.nc'))
#setup netcdf time dimension. different every year (20000301...20000831,20010301...20010831,etc)
dim3=ncdim_def('time','yrdoy',unlim=T,vals=as.numeric(gsub('X','',names(recon.stack))))
var=ncvar_def('swe','meters',dim=list(dim1,dim2,dim3),missval=-99,longname='reconstructed swe',compression=6)
#write netcdf files
     ncnew=nc_create(fnnc,var)
     ncvar_put(ncnew, var, valwrite)
     ncatt_put(ncnew,0,'projection#proj4text',projection(recon.stack))#'+proj=longlat +datum=NAD83')
     ncatt_put(ncnew,0,'dates','dates are defined in dim3 (time). should be daily Mar01 - Aug31')
     nc_close(ncnew)
     
#save snotelrecon data
     save(file=paste0(pathout,'/snotelrecondata_',yr,'_',recon.version,'.RData'),
          list=c(ls(pattern='snotelrecon'), ls(pattern=glob2rx(pattern='*df$'))))
     rm(list=c(ls(pattern='snotelrecon'), ls(pattern=glob2rx(pattern='*df$'))))

} else if(varname=='fsca'){
#netcdf write path
fnnc=file.path(pathout,paste0('reconfscadata_',yr,'_',recon.version,'.nc'))
#setup netcdf time dimension. different every year (20000301...20000831,20010301...20010831,etc)
dim3=ncdim_def('time','yrdoy',unlim=T,vals=as.numeric(gsub('X','',names(recon.stack))))
var=ncvar_def('fsca','/100',dim=list(dim1,dim2,dim3),missval=-99,longname='fsca from reconstructed swe',compression=6)
#write netcdf files
     ncnew=nc_create(fnnc,var)
     ncvar_put(ncnew, var, valwrite)
     ncatt_put(ncnew,0,'projection#proj4text',projection(recon.stack))#'+proj=longlat +datum=NAD83')
     ncatt_put(ncnew,0,'dates','dates are defined in dim3 (time). should be daily Mar01 - Aug31')
     nc_close(ncnew)
     
#save snotelrecon data
     save(file=paste0(pathout,'/snotelfscadata_',yr,'_',recon.version,'.RData'),
          list=c(ls(pattern='snotelrecon'), ls(pattern=glob2rx(pattern='*df$'))))
     rm(list=c(ls(pattern='snotelrecon'), ls(pattern=glob2rx(pattern='*df$'))))
}
}
}