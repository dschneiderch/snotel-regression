library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()

#: reads in swe.dat from recon output from selected version, saves usgs  projected stack and snotel values to recondata_v#.RData
## ----  Which Recon version? -----

yrs=seq(2000,2011)
recon.version='v3.1'
reconbase='/Volumes/PIA/snodis_output'

## ------

# set up out directory ------
pathout=paste0('data/recon_',recon.version)
if(is.na(file.info(pathout)$isdir))  dir.create(pathout,recursive=T)

get_recon_rasterstack <- function(filename){
    tr=file(filename,'rb')
    s=stack()
    nx=1950#columns
    ny=2580#rows
    nt=184#days
    dsize=4
    yr=unlist(strsplit(filename,'/',fixed=T))
    yr=yr[length(yr)-1]
    dt=strptime(paste0(yr,'-Mar-01'),format='%Y-%B-%d',tz='MST')
    for(i in seq(1,nt)){
        data=readBin(tr,what='numeric',size=dsize,n=nx*ny)#grid file grid by grid
        mat=matrix(data,nrow=ny,byrow=T)#convert to matrix with nrow
        mat=mat[nrow(mat):1,]#flip in y-direction to make right side map
        r=raster(mat)#rasters use column major cellnums instead row major, i.e. cell nums count up left to right, top to bottom instead of normal matrix way of top to bottom, left to right.
        extent(r) <- c(-112.25,-104.125,33,43.75)
        projection(r) <- '+proj=longlat +datum=NAD83'
        names(r) <- as.character(strftime(dt+24*3600*(i-1),format='%Y%m%d'))
        s=addLayer(s,r)
    }
    close(tr)
    return(s)
}

get_snotel_reconswe <- function(x){
    snotelrecon=as.data.frame(extract(x,snotellocs,sp=T))[,-c(1,2)]
    snotelrecon=melt(snotelrecon,id.vars=names(snotellocs),value.name='recon',variable.name='date')
    snotelrecon$date=strptime(as.character(snotelrecon$date),format='X%Y%m%d',tz='MST')
    return(snotelrecon)
}
    
fn=list.files(reconbase,full.name=T)
fileind=which(file.info(fn)$isdir)
recondir=fn[fileind][grep(recon.version,x=fn[fileind])]
#
fn=list.files(recondir,recursive=T,full.names=T,pattern='swe.dat')
#
## yrs=strsplit(fn,'/',fixed=T)
## pathlen=unlist(lapply(yrs,length))[1]
## yrs=as.numeric(unlist(lapply(yrs,'[',pathlen-1)) )
for(iyr in 1:length(yrs)){
    print(yrs[iyr])
    recon.stack=get_recon_rasterstack(fn[iyr])
    #recon.stack=projectRaster(recon.stack,crs=CRS('+init=epsg:5070'),method='ngb')
    snotelrecon=get_snotel_reconswe(recon.stack)
#
    yr=as.numeric(yrs[iyr])
    assign(paste0('snotelrecon',yr),snotelrecon)
    assign(paste0('recon',yr),recon.stack)
    #assign(paste0('recon',yr,'df'),getValues(recon.stack))

    save(file=paste0(pathout,'/recondata_',yr,'_',recon.version,'.RData'),
         list=ls(pattern=glob2rx(pattern='recon2*')))
    rm(list=ls(pattern=glob2rx(pattern='recon2*')))
    
    save(file=paste0(pathout,'/snotelrecondata_',yr,'_',recon.version,'.RData'),
         list=c(ls(pattern='snotelrecon'), ls(pattern=glob2rx(pattern='*df$'))))
    rm(list=c(ls(pattern='snotelrecon'), ls(pattern=glob2rx(pattern='*df$'))))
}
