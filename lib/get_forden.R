get_forden=function(fordensource,yr){
#
      forden=stack(paste0('data/forestdensity/',fordensource,'.nc'))
      if(fordensource=='umd_forden'){
          names(forden)=seq(2000,2012)
          if(yr>2012) yr=2012
          layerind=grep(as.character(yr),names(forden))
        } else {
          layerind=1
      }
      return(forden[[layerind]])
  }