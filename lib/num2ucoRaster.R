num2ucoRaster=function(v){
     rpred=raster(matrix(v,nrow=2580,byrow=T))
     projection(rpred)=CRS(projection(ucophv.stack))
     extent(rpred)=extent(ucophv.stack)
     return(rpred)
}