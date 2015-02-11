doPHVRCNfit <- function(xrdate,doydF,static_mdl){
     oldopt=options()
#      print(options('warn'))
     on.exit(options(oldopt))
     options(warn=2)
     #get appropriate recon data
     doydF$recon=mdldata[mdldata$recondate==xrdate,'recon']
     doydF$recondate=mdldata[mdldata$recondate==xrdate,'recondate']
     doydF$recon=scale(doydF$recon)
     
     # fix snotel record so there are no zeros
     snotelmdl=doydF$snotel
     snotelmdl=snotelmdl+runif(1,0,0.0001) # doesn't converge usually if there are 0s
     #   snotelmdl=log(snotelmdl)
     ## snotelmdl[is.infinite(snotelmdl)] <- NA
     doydF$snotel=snotelmdl
     #     
     ## PHV and RECON
     newformula=paste(formula(static_mdl)[2], formula(static_mdl)[1], formula(static_mdl)[3],' + recon')
     #    print(newformula)
     mdl.wrecon=try(
          glm(newformula,data=doydF,family=gaussian(link='log')),
          TRUE)
     failed=inherits(mdl.wrecon,'try-error')
     if(!failed || mdl.wrecon$converged) {
          #           print(mdl.wrecon)
          mdl.anova=tryCatch({
               anova(static_mdl,mdl.wrecon,test='F')},error=function(e) mdl.anova=NA)
     } else {
          mdl.wrecon=NA
     }
     
     if(length(mdl.anova)>1){
          ### Check significance of regression models
          rcn.sig.phv=tryCatch({
               if(mdl.anova$P[2]<0.05){
                    rcn.sig.phv=1
               } else {
                    rcn.sig.phv=0}
          }, error=function(e){ 
               print('----anova did not converge')
               return(rcn.sig.phv=NA)})
     } else { 
          rcn.sig.phv=NA
     }
          #
     return(list(mdl.wrecon,rcn.sig.phv,doydF))
}
