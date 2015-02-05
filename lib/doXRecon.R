#takes static model and fits all iterations of dyn model to get cross validated predictions.

doXRecon <- function(doydF,static_mdl,style){
     print(paste0('model date: ',doydF$date[1]))
     rdates=unique(recondata[!is.na(recondata$recon),'recondate'])
     if(style=='real-time') {
          myr=unique(doydF$yr)#get year of model sim
          ryr=strftime(rdates,'%Y')#vector of year of recon dates
          rdates=rdates[ryr<myr]#subset recondates to feed to fitting process so only previous years are available.
     }
     if(length(rdates!=0)){
          ldply(as.list(rdates),CVwrapper,doydF,static_mdl, .parallel=F)
     }
}