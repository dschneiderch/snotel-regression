library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(MASS)
library(doMC)

# select recon version ----------------------------------------------------
rversion='v3.2'
ryrs=seq(2000,2011)
snotelrecon=get_recondata(rversion,ryrs)

# generate swedata dataframe
swedata=generate_swedata(snotelrecon)

# SETUP MODEL DATA --------------------------------------------------------
mdldata=swedata[swedata$mth<6,]
mdldata$recondate=mdldata$date
mdldata$recondate[mdldata$mth<3] = as.POSIXct(strptime(paste0(mdldata$yr[mdldata$mth<3],'0301'),'%Y%m%d',tz='MST'))
#
mdlrecon <- function(x){
    # print(x$yr[1])
     x$recon=ifelse(x$recondate==x$date,
                    x$recon,
                    x$recon[x$mth==3 & x$dy==1])
     return(x)
}
mdldata=ddply(mdldata[mdldata$yr<2012,],.(yr,Station_ID),mdlrecon,.inform=T)

#dF=subset(mdldata,yrdoy==2001151)
# DEFINE MODEL FITTING FUNCTION -------------------------------------------
doFIT <- function(xrdate,dF, rdata){
     print(dF$yrdoy[1])   
     print(xrdate)
#    
# fix snotel record so there are no zeros
     snotelmdl=dF$snotel
     snotelmdl=snotelmdl+runif(1,0,0.00001)
     #   snotelmdl=log(snotelmdl)
     ## snotelmdl[is.infinite(snotelmdl)] <- NA
     dF$snotel=snotelmdl
     #
# because all dates before march 1 use mar1 recon, there are multiple sets of records with recondate==xrdate. sort by date then station and then take the first 237 since it is a repeating vector.
     ord=order(rdata$date,rdata$Station_ID)
     rdata2=rdata[ord,]
     rmdl=rdata2$recon[rdata2$recondate==xrdate]
     rmdl=rmdl[1:length(unique(rdata2$Station_ID))]
     ## rmdl[rmdl==0]=runif(1,0,0.00001)
     ## rmdl=log(rmdl)
     rmdl=scale(rmdl)
     rmdl[is.infinite(rmdl)] <- NA
     dF$recon=rmdl
     #
     ## print(head(dF))
     
     mdl=tryCatch({
          ### Base model
          ## start.param=fitdistr(swe2model$snotel,'gamma')
          ## lmmodel=lm(log(snotel) ~ Lat+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+Wd2ocean+1,data=dF)
          # startvals=c(log(mean(dF$snotel+runif(1,0,0.00001))),rep(0,11))
          #
          glm(snotel ~ Lat+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+Wd2ocean+1,family=gaussian(link='log'),data=dF)},#,start=startvals)},
          error= function(e) {
               print(paste0('returning empty phv model ', dF$date[1]))})
     if(!is.character(mdl)) {
          #
          ### Reduce with stepaic
          ## PHV only
          mdl.stepaic.phv=stepAIC(mdl, scope=list(upper= ~ Lat + Elev + Eastness + Northness + Slope + RegionalSlope + RegionalEastness + RegionalNorthness + FtprtW + Wbdiff + NWbdiff + Wd2ocean + 1, lower=~1), direction='both',trace=F,k=log(nrow(dF)))
          #   
          ## PHV and RECON
          newformula=paste(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv$formula[3],'+ recon',sep=' ')
          print(newformula)
          mdl.wrecon=tryCatch({
               glm(newformula,data=dF,family=gaussian(link='log'))},
               error=function(e) {
                    print(paste0('returning empty phvrcn model ', dF$date[1], ', rswe: ',dF$recondate))})
          
          ### Check significance of regression models
          mdl.anova=anova(mdl.stepaic.phv,mdl.wrecon,test='F')
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
     ## ### Predict surface from regression
     ## #PHV only surface
     ##     phv=predict(mdl.stepaic.phv,type='response',newdata=dF)
     ## #PHV and RECON surface
     ##     phvrcn=predict(mdl.wrecon,type='response',newdata=dF)     
     ## #
     return(list(mdl.stepaic.phv,mdl.wrecon,rcn.sig.phv))
}

# run model ---------------------------------------------------------------

registerDoMC(3)
olddata=subset(mdldata,(mth==1 | mth==2 | mth==3 | mth == 4 | mth==5) & (dy==1) & yr>2000)
rundata=olddata
doXRecon <- function(datasub,rdata){
     llply(as.list(unique(rdata$recondate)), doFIT, datasub, rdata)
}    
test=dlply(rundata,.(yrdoy),.fun=doXRecon,rundata[,c('Station_ID','date','recondate','recon')],.parallel=T,.drop=F,.inform=T, .paropts=list(.export=c('rundata','doFIT'),.packages(all.available = TRUE)))
#
mdl.stepaic.phv=lapply(test,lapply,'[[',1)  
mdl.wrecon=lapply(test,lapply,'[[',2)
rcn.sig.phv=lapply(test,lapply,'[[',3)

# save models -------------------------------------------------------------
pn='data/latest.model'
fullpath=file.path(pn,paste0('rswe_',rversion))
fn='mdl_1st_of_month.stepaic.gauss_log.eachday.snotel.rswescaled.xrecon.Ftest.RData'
if(is.na(file.info(fullpath)$isdir)){
     dir.create(fullpath,recursive=T)
}
save(file=file.path(fullpath,fn),
     list=c('mdl.stepaic.phv','mdl.wrecon','rcn.sig.phv','originaldata','swedata'))


## ---------------------   end of ddply stuff

# 
# 
# #fit models
# for(yr in 2001:2012){
#     stdate=strptime(paste0((yr),'-01-01'),tz='MST',format='%Y-%m-%d')
#     enddate=strptime(paste0(yr,'-05-15'),tz='MST',format='%Y-%m-%d')
#     workdate=stdate
#     j=1
#     mdl.stepaic.phv[[yr-2000]]=replicate(numdays,list)
#     mdl.wrecon[[yr-2000]]=replicate(numdays,list)
#     while(workdate <= enddate) {#
#         print(paste0('iteration - ',strftime(workdate,'%B-%d-%Y')))
# 
#         wdate=as.numeric(strftime(workdate,'%Y%j'))
#         swedatasub=subset(swedata,yrdoy==wdate)
# 
#         ## snotel_scatt=scale(swedatasub$snotel)
#         ## dev.set(3)
#         ## hist(snotel_scatt,main='scaled snotel')
#         swedatasub$snotel[swedatasub$snotel==0]=0.001
#         snotel_log=log(swedatasub$snotel)
#         snotel_log[is.infinite(snotel_log)] <- NA
#         ## dev.set(4)
#         ## hist(snotel_log,main='log snotel')
#         ## dev.set(5)
#         ## snotel_log_sc=(snotel_log-mean(snotel_log,na.rm=T))/sd(snotel_log,na.rm=T)
#         ## hist(snotel_log_sc,main='scaled log  snotel')
#         
#         snotelraw=swedatasub$snotel#save for later
#         swedatasub$snotel=snotel_log
# 
#         ## z2keep=which(snotelraw>0)#thinkning recon/fsca should be the definer here
#         z2keep=1:length(snotelraw)
#         swedata2model=swedatasub[z2keep,]
# 
#         
#         ### Base model           
#             mdl=glm(snotel ~ Lat+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+Wd2ocean+1,family=gaussian,data=swedata2model)
# 
#         mdl.wrecon[[yr-2000]][[j]]=replicate(numryrs,list)
#         mdl.stepaic.phv[[yr-2000]][[j]]=replicate(numryrs,list)
#         for(ryr in min(unique(year(snotelrecon$date))):max(unique(year(snotelrecon$date)))){
#             rdate=strptime(paste0(yr,'0301'),format='%Y%m%d',tz='MST') #do use yr - used to check if before march 1
# 
#             if(workdate<rdate){
#                 rswemax=as.numeric(strftime(strptime(paste0(ryr,'0301'),format='%Y%m%d',tz='MST'),'%Y%j'))
#                 rswe=subset(swedata[,'recon'],swedata$yrdoy==rswemax)
#             } else {
#                 rswedate=as.numeric(strftime(strptime(paste0(ryr,strftime(workdate,'%m%d')),format='%Y%m%d',tz='MST'),'%Y%j'))
#                 rswe=subset(swedata[,'recon'],swedata$yrdoy==rswedate)
#             }
# 
#             rswe[rswe==0]=0.0001 #model won't fit to NA (from logtransform below) so make some small non-zero number since we don't have fsca data for jan1-feb28
#             rswe.log=log(rswe)
#             rswe.log[is.infinite(rswe.log)] <- NA
#             rswe.log=scale(rswe.log)
#             swedata2model$recon=rswe.log[z2keep]
#                  
# 
# ### Reduce with stepaic
#             ## PHV only
#             mdl.stepaic.phv[[yr-2000]][[j]][[ryr-1999]]=stepAIC(mdl, scope=list(upper= ~ Lat + Elev + Eastness + Northness + Slope + RegionalSlope + RegionalEastness + RegionalNorthness + FtprtW + Wbdiff + NWbdiff + Wd2ocean + 1, lower=~1), direction='both',trace=F,k=log(nrow(swedata2model)))
#             
#             oldformula=paste0(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv[[yr-2000]][[j]][[ryr-1999]]$formula[3])
#             newformula=paste0(oldformula,' + recon')
#             print(newformula)
#             print(paste('ryr = ',ryr))
#             
# ### Reduce with stepaic
#             ## PHV and RECON
#             mdl.wrecon[[yr-2000]][[j]][[ryr-1999]]=glm(newformula,data=swedata2model,family=gaussian)
#             mdl.anova=anova(mdl.stepaic.phv[[yr-2000]][[j]][[ryr-1999]],mdl.wrecon[[yr-2000]][[j]][[ryr-1999]],test='F')
#             
# ### Predict surface from regression 
#             if(mdl.anova$P[2]<0.05){
#                 rcn.sig.phv[yr-2000,j,ryr-1999]=1
#                 print('significant? y')
#             } else {
#                 rcn.sig.phv[yr-2000,j,ryr-1999]=0
#                 print('significant? n')
#             }
# 
#         }
#         j=j+1
#         workdate=workdate+days(1)
#     }
# }
# 
# pn='data/latest.model'
# fullpath=file.path(pn,paste0('rswe_',rversion))
# if(file.info(fullpath)$isdir){
#     save(file=file.path(fullpath,'mdl.stepaic.gauss.eachday.snotellog.rswelog.xrecon.Ftest.RData'),list=c('mdl.stepaic.phv','mdl.wrecon','rcn.sig.phv','originaldata','swedata'))
# } else {
#     dir.create(fullpath,recursive=T)
#     save(file=file.path(fullpath,'mdl.stepaic.gauss.eachday.snotellog.rswelog.xrecon.Ftest.RData'),list=c('mdl.stepaic.phv','mdl.wrecon','rcn.sig.phv','originaldata','swedata'))
# }
