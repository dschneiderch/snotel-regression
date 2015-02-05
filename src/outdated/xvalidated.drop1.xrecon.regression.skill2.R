library('ProjectTemplate')
#setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(MASS)
library(ipred)
library(doMC)
library(doSNOW)

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
mdldata=ddply(mdldata[mdldata$yr<2012,],.(yr,Station_ID),mdlrecon,.inform=F)

#dF=subset(mdldata,yrdoy==2001091)
#xrdate=strptime(20010401,'%Y%m%d')
#rdata=mdldata[,c('Station_ID','date','recondate','recon')]
# DEFINE MODEL FITTING FUNCTION -------------------------------------------
doXVAL <- function(xrdate,dF, rdata){
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
               print(paste0('returning empty phv model ', dF$date[1]))
               }
          )
     if(!is.character(mdl)) {
          #
          ### Reduce with stepaic
          ## PHV only
          mdl.stepaic.phv=stepAIC(mdl, scope=list(upper= ~ Lat + Elev + Eastness + Northness + Slope + RegionalSlope + RegionalEastness + RegionalNorthness + FtprtW + Wbdiff + NWbdiff + Wd2ocean + 1, lower=~1), direction='both',trace=F,k=log(nrow(dF)))
          #   
          ## PHV and RECON
          newformula=paste(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv$formula[3],'+ recon',sep=' ')
          #    print(newformula)
          
          mdl.wrecon=tryCatch({
               glm(newformula,data=dF,family=gaussian(link='log'))},
               error=function(e) {
                    print(paste0('returning empty phvrcn model ', dF$date[1], ', rswe: ',dF$recondate))
               }
          )
          
          ### Check significance of regression models
          mdl.anova=anova(mdl.stepaic.phv,mdl.wrecon,test='F')
          #
          if(mdl.anova$P[2]<0.05){
               rcn.sig.phv=1
          } else {
               rcn.sig.phv=0
          }
          
          ### cross validation
          myglm=function(formula,data) {
               glm(formula,data,family=gaussian(link='log'))
          }
          mypredict.glm <- function(object, newdata){
               predict(object, newdata,type='response')
          }
          #PHV only surface
          phv=errorest(formula=formula(mdl.stepaic.phv),data=dF, model=myglm, predict=mypredict.glm, estimator='cv',est.para=control.errorest(k=nrow(dF), predictions=T))
          #PHV and RECON surface
          phvrcn=errorest(formula=formula(mdl.wrecon),data=dF, model=myglm, predict=mypredict.glm, estimator='cv',est.para=control.errorest(k=nrow(dF), predictions=T))
          
     } else {
          phv=rep(NA,nrow(dF))
          phvrcn=rep(NA,nrow(dF))
          rcn.sig.phv=rep(NA,nrow(dF))
     }
     #
     
     return(data.frame(recondate=xrdate,Station_ID=dF$Station_ID,snotel=dF$snotel,phv$predictions,phvrcn$predictions,recon=dF$recon,sig=rcn.sig.phv))
}

# run model ---------------------------------------------------------------

registerDoMC(3)
olddata=subset(mdldata,(mth==1 | mth==2 |  mth==3 | mth == 4 | mth==5) & (dy==1) & yr>2000)
rundata=olddata
doXRecon <- function(datasub,rdata){
     ldply(as.list(unique(rdata$recondate)), doXVAL, datasub, rdata)
}    
cvpreds=ddply(rundata,
              .(yrdoy),
              .fun=doXRecon,rundata[,c('Station_ID','date','recondate','recon')],
              .parallel=T,
              .drop=F,
              .inform=F, 
              .paropts=list(.export=c('rundata','doXVAL'),
                            .packages=c('ipred','MASS')
              )
)



# save models -------------------------------------------------------------
pn='diagnostics'
fullpath=file.path(pn,paste0('rswe_',rversion))
if(is.na(file.info(fullpath)$isdir)) dir.create(fullpath, recursive=T)
#
fn='cvdrop1.skillscores_1st_of_month.stepaic.gauss_log.eachday.snotel.rswescaled.xrecon.RData'
if(is.na(file.info(fullpath)$isdir)){
     dir.create(fullpath,recursive=T)
}
save(file=file.path(fullpath,fn),
     list='cvpreds')

