#this is the workhorse fitting function. with ddply, subset master dataframe for a single doy, fit static model, then iterate through all dynamic models FOR THAT DOY, and produce crossvalidated estimates. Then I'm fitting a spatial model to the residuals to create a composite prediction. then use some cost function to figure out which dynamic model is best and output a dataframe 'which_recon_date'.
doDOYfit<-function(doydF,cost,style) {  ## from ddply
     static_mdl<-doPHVfit(doydF) 
     cvpreds=doXRecon(doydF,static_mdl,style)#give all relevant dates to do xrecon. realtime only usees dates from previous years.
     #print(head(cvpreds))
          locs=spTransform(snotellocs,CRS('+init=epsg:5070'))
          fullpreds=ddply(cvpreds,.(recondate),doSPATIAL,locs,'xval',.inform=T)
          which_recon_date=doCOST(fullpreds,cost,recon.version)#
      
          return(which_recon_date)
     #return(fullpreds)
}
