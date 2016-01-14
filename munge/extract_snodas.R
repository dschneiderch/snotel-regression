library(ncdf4)
library(raster)

source('lib/read_surveydates.R')
fn='data/spatialvalidation/snow_survey_dates.txt'
dtes=read_surveydates(fn)

dts=as.Date(dtes$date,tz='MST')
dts.snodas=as.numeric(dts)-as.numeric(as.Date('1601-01-01'))

fns=list.files('data/spatialvalidation/snodas/orig_data',pattern='sub.nc',full.names=T)

snodas.stack=llply(fns,function(x){
nid=nc_open(x)
xcoord=nid$dim[[1]]$vals
ycoord=nid$dim[[2]]$vals
varname='Snow_Water_Equivalent_with_State_Model_variable_type'
tind=which(floor(nid$dim$time$vals) %in% dts.snodas)
print(tind)
date.snodas=nid$dim$time$vals[tind]
date=date.snodas+as.Date('1601-01-01')
nc_close(nid)
s=raster::stack(x)
if(length(tind)){
print('??')
# snodas=ncvar_get(nid,varname,start=c(1,1,tind),count=c(length(xcoord),length(ycoord),1))
r=s[[tind]]
names(r)=date
print(r)
return(r)
} else {
return(NULL)
}
})
swe.snodas=stack(unlist(snodas.stack))
swe.snodas.ucrb=crop(swe.snodas,extent(c(-112.25,-104.125,33,43.75)))

x2=xmin(swe.snodas.ucrb)+res(swe.snodas.ucrb)[1]/2+res(swe.snodas.ucrb)[1]*(seq(1,ncol(swe.snodas.ucrb))-1)
y2=ymax(swe.snodas.ucrb)-res(swe.snodas.ucrb)[2]/2-res(swe.snodas.ucrb)[2]*(seq(1,nrow(swe.snodas.ucrb))-1)
ncdimdates=as.numeric(gsub('.','',gsub('X','',names(swe.snodas)),fixed=T))

dim1=ncdim_def('long','degree',x2)
dim2=ncdim_def('lat','degree',y2)
dim3=ncdim_def('time','ansi date',ncdimdates)
var=ncvar_def('swe','meters',dim=list(dim1,dim2,dim3),missval=-9999)
ncnew=nc_create('data/spatialvalidation/snodas/snodas.survey.nc',var)
ncvar_put(ncnew,var,getValues(swe.snodas.ucrb))
ncatt_put(ncnew,0,'projection#proj4text',projection(swe.snodas.ucrb))
ncatt_put(ncnew,0,'dates',names(swe.snodas.ucrb))
ncatt_put(ncnew,0,'dates_misc','snodas unavailable before 2003')
nc_close(ncnew)

