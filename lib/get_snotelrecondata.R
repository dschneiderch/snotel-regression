get_snotelrecondata<- function(rversion,yrs){
    tmp=new.env()
    for(fn in paste0('data/recon_',rversion,'/snotelrecondata_',yrs,'_',rversion,'.RData')){
        load(fn,envir=tmp)
    }
     snotelrecon=do.call(rbind,lapply(ls(pattern='snotelrecon2',envir=tmp),get,envir=tmp))
     return(snotelrecon)
}
