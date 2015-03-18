#this is the workhorse fitting function. with ddply, subset master dataframe for a single doy, fit static model, then iterate through all dynamic models FOR THAT DOY, and produce crossvalidated estimates. Then I'm fitting a spatial model to the residuals to create a composite prediction. then use some cost function to figure out which dynamic model is best and output a dataframe 'which_recon_date'.
doDOYfit<-function(doydF,cost,style) {  ## from ddply
     #set up empty dataframes to output below if there are no results.
     emptydf=function(yrdoy,dte){
          colvec=c('yrdoy','date','yr','phvrcn_recondate','recon_costdate','skill_phv','skill_phvrcn','skill_recon','skill_phvfull','skill_phvrcnfull','skill_reconfull')
          c1=data.frame(matrix(NA,ncol=length(colvec)))
          colnames(c1)=colvec
          c1$yrdoy=yrdoy
          c1$date=dte
#           #
          colvec=c('date','yrdoy','snotel','recondate','reconrt','phvresid','phvrcnresid','reconrtresid','phv.pred','phvrcn.pred','phv.fullpred','phvrcn.fullpred','reconrt.fullpred','Station_ID')
          c3=data.frame(matrix(NA,ncol=length(colvec)))
          colnames(c3)=colvec 
          c3$yrdoy=yrdoy
          c3$date=dte

          colvec=c('date','ryr','Station_ID','historical_apcp','historical_cumdegday','apcp','cumdegday','dist','weight','mth','dy','recondate','recon','reconopt')
          c4=data.frame(matrix(NA,ncol=length(colvec)))
          colnames(c4)=colvec 
          c4$date=dte

          return(list(c1,c3,c4))
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
          #
          snotellocs=snotellocs[snotellocs$Station_ID %in% doydF$Station_ID,]#this is incase we do another year such as fassnacht years. keep in mind that only years availble 2000-2012 will be included regardless if there wre others.
          snotellocs.usgs=spTransform(snotellocs,CRS('+init=epsg:5070'))
          if(!empty(cvpreds)){
          
               registerDoMC(allocated.cores)
               # registerDoMC(2)
               fullpreds=ddply(cvpreds,.(recondate),doSPATIAL,snotellocs.usgs,covrange,'xval',.inform=F,.parallel=F)
               which_recon_date=doCOST(fullpreds,cost,recon.version)#
               bestind=which(fullpreds$recondate %in% which_recon_date$phvrcn_recondate)
               fullpreds=fullpreds[bestind,]
               
               optlist=doOptRecon(doydF,static_mdl,style)
               optpreds=optlist[[1]]
               wmat=optlist[[2]]

               fullpreds$phvrcn.opt=optpreds$phvrcn.predictions
               fullpreds$reconopt=optpreds$reconrt

               return(list(which_recon_date,fullpreds,wmat))    
          }
     } else {
          print('--phv model didn\'t converge')
          return(emptydf(unique(doydF$yrdoy),unique(doydF$date)))
#           return(list(data.frame(),data.frame()))
     }
}

