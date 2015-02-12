#this is the workhorse fitting function. with ddply, subset master dataframe for a single doy, fit static model, then iterate through all dynamic models FOR THAT DOY, and produce crossvalidated estimates. Then I'm fitting a spatial model to the residuals to create a composite prediction. then use some cost function to figure out which dynamic model is best and output a dataframe 'which_recon_date'.
doDOYfit<-function(doydF,cost,style) {  ## from ddply
     #set up empty dataframes to output below if there are no results.
     emptydf=function(){
          colvec=c('yrdoy','date','yr','phvrcn_recondate','recon_costdate','skill_phv','skill_phvrcn','skill_recon','skill_phvfull','skill_phvrcnfull','skill_reconfull')
          c1=data.frame(matrix(NA,ncol=length(colvec)))
          colnames(c1)=colvec
#           #
colvec=c('yrdoy','date','recondate','phv.globalI','phvrcn.globalI','phv.pvalue','phvrcn.pvalue','phv.full.globalI','phvrcn.full.globalI','phv.full.pvalue','phvrcn.full.pvalue')
          c2=data.frame(matrix(NA,ncol=length(colvec)))
          colnames(c2)=colvec
          return(list(c1,c2))
     }
     #     
#      withOptions <- function(optionList, expr) {
#           oldOpts <- options(optionList)
#           on.exit(options(oldOpts))
#           expr # lazily evaluate
#      }
     
     print(paste0('model date: ',doydF$date[1]))
#      static_mdl=withOptions(list(warn=2),doPHVfit(doydF))
     static_mdl<-doPHVfit(doydF)
     #      print(static_mdl)
     if(length(static_mdl)!=1){
          cvpreds=doXRecon(doydF,static_mdl,style)#give all relevant dates to do xrecon. realtime only uses dates from previous years.
          #print(head(cvpreds))
          snotellocs.usgs=spTransform(snotellocs,CRS('+init=epsg:5070'))
          if(!empty(cvpreds)){
               fullpreds=ddply(cvpreds,.(recondate),doSPATIAL,snotellocs.usgs,'xval',.inform=T,.parallel=F)
               which_recon_date=doCOST(fullpreds,cost,recon.version)#
               morandf=plot_localmoran(fullpreds[fullpreds$recondate==which_recon_date$phvrcn_recondate,],snotellocs.usgs,recon.version)
               return(list(which_recon_date,morandf))    
          }
     } else {
          print('--phv model didn\'t converge')
          return(emptydf())
#           return(list(data.frame(),data.frame()))
     }
}

