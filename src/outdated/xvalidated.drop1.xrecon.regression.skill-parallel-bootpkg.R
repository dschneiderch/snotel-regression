library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(spdep)
library(MASS)
library(foreach)
library(doSNOW)
library(doParallel)
library(doMC)
library(boot)

# --- Load recon version
rversion='v3.1'
pn='data/latest.model'
fullpath=file.path(pn,paste0('rswe_',rversion))
load(file.path(fullpath,'mdl.stepaic.gauss.eachday.snotellog.rswelog.xrecon.Ftest.RData'))
# contains models and modeldata
# --- 

## --- SETUP PARALLEL PROCESSING
getnodes <- function() {
  f <- Sys.getenv('PBS_NODEFILE')
  x <- if (nzchar(f)) readLines(f) else rep('localhost', 3)
  as.data.frame(table(x), stringsAsFactors=FALSE)
}

nodes <- getnodes()
cl <- makeCluster(nodes$x, outfile='')
registerDoParallel(cl)

registerDoParallel(cl=1,cores=3)


#difference between swedata and original data is phvrec was scaled in swedata
swedata$mnth=factor(strftime(swedata$date,'%b'),levels=unique(strftime(swedata$date[order(swedata$date)],'%b')))
originaldata$mnth=swedata$mnth


## INITIALIZE CROSS-VALIDATION
mdlyrs=seq(2001,2002)
nummdlyrs=length(mdlyrs)
mdldoys=3#will always start on jan 1.
numstations=length(unique(swedata$Station_ID))
numryrs=2#always starts from 2000. 
#
pr.phvrcn=NULL
pr.phv=NULL
SSphv=NULL
SSphvrcn=NULL
SSrecon=NULL
swediff.phv=NULL
swediff.phvrcn=NULL
swediffpct.phv=NULL
swediffpct.phvrcn=NULL
#
swediff.phvmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
swediff.phvrcnmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
swediffpct.phvmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
swediffpct.phvrcnmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
#
r2.phvmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs))
r2.phvrcnmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs))
r2.reconmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs))
#
valvec=as.numeric(paste0(
    rep(mdlyrs,each=mdldoys),
    rep(sprintf('%03d',seq(1,mdldoys)),length=(diff(range(mdlyrs))+1)*mdldoys)))

## comb <- function(x, ...) {
##   lapply(seq_along(x),
##     function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
## }

oper=foreach(yj=valvec, .packages=c('doMC','MASS')) %dopar% {
    #,.combine='comb',.multicombine=T,.init=list(list(),list(),list())) %dopar% {
    print(paste('iteration: ',yj,sep=''))
#       #select day being modeled.
    swe2model.all=swedata[swedata$yrdoy==yj,]
    swe2model.all$snotel[swe2model.all$snotel==0]=runif(1,0,0.00001)
                                        #
    for(ryr in seq(2000,(2000+numryrs-1))){
        print(paste0('recon year: ',ryr))
        ryj=as.numeric(paste0(ryr,substring(yj,5,7)))

        if(strptime(ryj,'%Y%j') < strptime(paste0(ryr,'0301'),'%Y%m%d',tz='MST')){
            rswemax=as.numeric(strftime(strptime(paste0(ryr,'0301'),format='%Y%m%d',tz='MST'),'%Y%j'))
            rmdl=swedata$recon[swedata$yrdoy==rswemax]
            rmdl[rmdl==0]=runif(1,0,0.00001)
            rmdl=scale(log(rmdl))
            rmdl[is.infinite(rmdl)] <- NA
            swe2model.all$recon=rmdl
        } else {              
            rmdl=swedata$recon[swedata$yrdoy==ryj]
            rmdl[rmdl==0]=runif(1,0,0.00001)
            rmdl=scale(log(rmdl))
            rmdl[is.infinite(rmdl)] <- NA
            swe2model.all$recon=rmdl
        }
  #
        ### Base model
            ## start.param=fitdistr(swe2model$snotel,'gamma')
            mdl=glm(snotel ~ Lat+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+Wd2ocean+1,family=gaussian,data=swe2model.all)
#                  
### Reduce with stepaic
            ## PHV only
            mdl.stepaic.phv=stepAIC(mdl, scope=list(upper= ~ Lat + Elev + Eastness + Northness + Slope + RegionalSlope + RegionalEastness + RegionalNorthness + FtprtW + Wbdiff + NWbdiff + Wd2ocean + 1, lower=~1), direction='both',trace=F,k=log(nrow(swe2model.all)))
   #                 
### Reduce with stepaic
            ## PHV and RECON
        newformula=paste(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv$formula[3],'+ recon',sep=' ')
                                        #print(newformula)
        mdl.wrecon=glm(newformula,data=swe2model.all,family=gaussian)

difffunc <- function(y,yhat){
    return(mean(yhat-y,na.rm=T))
}
pctdifffunc <- function(y,yhat){
    return(mean((yhat-y)^2,na.rm=T))
}
spearmanfunc <- function(y,yhat){
return(1- (6*sum((order(y)-order(yhat))^2))/ (length(y)*(length(y)^2-1)))
    }

    return(cor(y,yhat,method='spearman'))
}
        
        cv.glm(swe2model.all,mdl.wrecon,cost=pctdifffunc,K=nrow(swe2model.all))$delta
        cv.glm(swe2model.all,mdl.wrecon,cost=difffunc,K=nrow(swe2model.all))$delta
        cv.glm(swe2model.all,mdl.wrecon,cost=spearmanfunc,K=nrow(swe2model.all))$delta
        
        for(dropiter in seq(1,numstations)) {
            z2keep=seq(1,numstations,1)[-dropiter]
            swe2model=swe2model.all[z2keep,]
            swe2model$snotel=log(swe2model$snotel)
            swe2model$snotel[is.infinite(swe2model$snotel)] <- NA
                                        #
### Base model
            ## start.param=fitdistr(swe2model$snotel,'gamma')
            mdl=glm(snotel ~ Lat+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+Wd2ocean+1,family=gaussian,data=swe2model)
#                  
### Reduce with stepaic
            ## PHV only
            mdl.stepaic.phv=stepAIC(mdl, scope=list(upper= ~ Lat + Elev + Eastness + Northness + Slope + RegionalSlope + RegionalEastness + RegionalNorthness + FtprtW + Wbdiff + NWbdiff + Wd2ocean + 1, lower=~1), direction='both',trace=F,k=log(nrow(swe2model)))
   #                 
### Reduce with stepaic
            ## PHV and RECON
            newformula=paste(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv$formula[3],'+ recon',sep=' ')
                                        #print(newformula)
            mdl.wrecon=glm(newformula,data=swe2model,family=gaussian)
                                        #                
### Predict surface from regression
#PHV only surface
            phv=predict(mdl.stepaic.phv,type='response',newdata=swe2model.all[-z2keep,])
#PHV and RECON surface
            phvrcn=predict(mdl.wrecon,type='response',newdata=swe2model.all[-z2keep,])           
#
                                        # transform data back from log space
            pr.phv[dropiter]=exp(phv)
            pr.phvrcn[dropiter]=exp(phvrcn)
        }
                                        #           
### Compute skill
        yobs=swe2model.all$snotel
        yrecon=as.numeric(swedata$recon[swedata$yrdoy==yj])
        SSrecon = cor(yrecon,yobs)^2
  #          
        pr.phv[yrecon==0]=0
        SSphv = cor(pr.phv,yobs)^2#1 - SSE/SSNull
        swediff.phv = pr.phv-yobs
        swediffpct.phv = (pr.phv-yobs)/yobs
        #swediffpct.phv[is.na(swediffpct.phv)]=0
        swediffpct.phv[is.infinite(swediffpct.phv)]=1000
   #         
        pr.phvrcn[yrecon==0]=0
        SSphvrcn = cor(pr.phvrcn,yobs)^2 #1 - SSE/SSNull
        swediff.phvrcn = pr.phvrcn-yobs
        swediffpct.phvrcn = (pr.phvrcn-yobs)/yobs
        #swediffpct.phvrcn[is.na(swediffpct.phvrcn)]=0
        swediffpct.phvrcn[is.infinite(swediffpct.phvrcn)]=1000
                                        #           
    ## if(length(which(is.na(SSrecon)))!=0) print(paste('NA values - recon ',mth,yr,'- recon',ryr,sep=''))
 #
        dayind=as.numeric(substring(yj,5,7))
        yrind=as.numeric(substring(yj,1,4))-2000
        #
        r2.phvmodel[dayind, yrind,ryr-1999]=SSphv
        r2.phvrcnmodel[dayind, yrind,ryr-1999]=SSphvrcn
        r2.reconmodel[dayind, yrind,ryr-1999]=SSrecon
        swediff.phvmodel[dayind, yrind,ryr-1999,]=swediff.phv
        swediff.phvrcnmodel[dayind, yrind,ryr-1999,]=swediff.phvrcn
        swediffpct.phvmodel[dayind, yrind,ryr-1999,]=swediffpct.phv
        swediffpct.phvrcnmodel[dayind, yrind,ryr-1999,]=swediffpct.phvrcn
                                        #          
    }
    list(r2.phvmodel,r2.phvrcnmodel,r2.reconmodel,
         swediff.phvmodel,swediff.phvrcnmodel,
         swediffpct.phvmodel,swediffpct.phvrcnmodel)
} 

split_dopar <- function(x,arnum){
# x is output from dopar without any combining
# arnum is the array number you want to extract. corresponds with location in list from each iteration
    tmp1=lapply(x,'[[',arnum)
    ar=do.call(pmax,args=c(tmp1,na.rm=T))
    return(ar)
}
r2.phvmodel=split_dopar(oper,1)
r2.phvrcnmodel=split_dopar(oper,2)
r2.reconmodel=split_dopar(oper,3)
swediff.phvmodel=split_dopar(oper,4)
swediff.phvrcnmodel=split_dopar(oper,5)
swediffpct.phvmodel=split_dopar(oper,6)
swediffpct.phvrcnmodel=split_dopar(oper,7)
 


## comb <- function(x, ...) {
##   lapply(seq_along(x),
##     function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
## }

registerDoParallel(2)
## In a much smaller loop, find swediff for recon vs snotel
swediff.reconmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
swediffpct.reconmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
oper=foreach(yj = valvec, .packages=c('doMC')) %dopar% {
#,.combine='comb',.multicombine=T,.init=list(list(),list()))
        print(paste('iteration: ',yj,sep=''))
#       
        swe2model.all=swedata[swedata$yrdoy==yj,]
        swe2model.all$snotel[swe2model.all$snotel==0]=runif(1,0,0.00001)
        yobs=swe2model.all$snotel
        #
        registerDoMC(2)
      foreach(ryr = seq(2000,(2000+numryrs-1))) %dopar% {
            print(paste('ryr = ',ryr))
            ryj=as.numeric(paste0(ryr,substring(yj,5,7)))
            if(strptime(ryj,'%Y%j') < strptime(paste0(ryr,'0301'),'%Y%m%d',tz='MST')){
                rswemax=as.numeric(strftime(strptime(paste0(ryr,'0301'),format='%Y%m%d',tz='MST'),'%Y%j'))
                yrecon=swedata$recon[swedata$yrdoy==rswemax]
            } else {              
                yrecon=swedata$recon[swedata$yrdoy==ryj]
            }
#                  
#            print(paste('yobs: ',which(is.na(yobs)),sep=''))
 #           print(paste('yrecon: ',which(is.na(yrecon)),sep=''))
                                        #
            dayind=as.numeric(substring(yj,5,7))
            yrind=as.numeric(substring(yj,1,4))-2000
            #
            swediff.reconmodel[dayind,yrind,ryr-1999,]=yrecon-yobs
            swediffpct.reconmodel[dayind,yrind,ryr-1999,]=(yrecon-yobs)/yobs
            ## swediffpct.reconmodel[icount,ryr-1999,is.na(swediffpct.reconmodel[icount,ryr-1999,])]=0
            swediffpct.reconmodel[dayind,yrind,ryr-1999,is.infinite(swediffpct.reconmodel[dayind,yrind,ryr-1999,])]=1000
#
            list(swediff.reconmodel,swediffpct.reconmodel)                                   #
        }
    }

split_dopar <- function(x,arnum){
# x is output from dopar without any combining
# arnum is the array number you want to extract. corresponds with location in list from each iteration
    tmp1=lapply(unlist(x,recursive=F),'[[',arnum)
    ar=do.call(pmax,args=c(tmp1,na.rm=T))
    return(ar)
}
swediff.reconmodel1=split_dopar(oper,1)
swediffpct.reconmodel2=split_dopar(oper,2)
#swediff.reconmodel=Reduce(function(x,y) ifelse(!is.na(x),x,y),ar1)
#swediffpct.reconmodel=Reduce(function(x,y) ifelse(!is.na(x),x,y),ar2)
#or can do do.call(pmax,c(ar1,na.rm=T))

pn='diagnostics'
fullpath=file.path(pn,paste0('rswe_',rversion))
save(file=file.path(fullpath,'cvdrop1.skillscores.xrecon.regression.gauss.snotellog.reconlog.RData'),
     list=c('r2.phvmodel','r2.phvrcnmodel','r2.reconmodel',
         'swediff.phvmodel','swediff.phvrcnmodel','swediff.reconmodel',
         'swediffpct.reconmodel','swediffpct.phvmodel','swediffpct.phvrcnmodel'))
