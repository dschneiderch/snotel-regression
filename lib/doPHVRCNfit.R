doPHVRCNfit <- function(xrdate,doydF,static_mdl){
     #get appropriate recon data
     doydF$recon=mdldata[mdldata$recondate==xrdate,'recon']
     doydF$recondate=mdldata[mdldata$recondate==xrdate,'recondate']
     doydF$recon=scale(doydF$recon)
     
     # fix snotel record so there are no zeros
     snotelmdl=doydF$snotel
     snotelmdl=snotelmdl+runif(1,0,0.00001) # doesn't converge usually if there are 0s
     #   snotelmdl=log(snotelmdl)
     ## snotelmdl[is.infinite(snotelmdl)] <- NA
     doydF$snotel=snotelmdl
     #
          
     if((!is.character(static_mdl) && !is.na(static_mdl))) { #check no error messages in previous steps
          
          ## PHV and RECON
          newformula=paste(formula(static_mdl)[2], formula(static_mdl)[1], formula(static_mdl)[3],' + recon')
      #    print(newformula)
          mdl.wrecon=tryCatch({
               glm(newformula,data=doydF,family=gaussian(link='log'))},
               error=function(e) {
                    print(paste0('returning empty phvrcn model ', doydF$date[1], ', rswe: ',doydF$recondate))
               })
          
          ### Check significance of regression models
          mdl.anova=anova(static_mdl,mdl.wrecon,test='F')
          #
          if(mdl.anova$P[2]<0.05){
               rcn.sig.phv=1
          } else {
               rcn.sig.phv=0
          }
     } else {
          mdl.stepaic.phv=NA
          mdl.wrecon=NA
          rcn.sig.phv=NA
     }
     #
     return(list(mdl.wrecon,rcn.sig.phv,doydF))
}
