library('ProjectTemplate')
setwd('~/GoogleDrive/snotel-regression_project')
load.project()
library(spdep)
library(MASS)

# --- Load recon version
rversion='v3.1'
pn='data/latest.model'
fullpath=file.path(pn,paste0('rswe_',rversion))
load(file.path(fullpath,'mdl.stepaic.gauss.eachday.snotellog.rswelog.xrecon.Ftest.RData'))
# contains models and modeldata
# --- 

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
icount=1
for(yj in valvec){
    print(paste('iteration: ',yj,sep=''))
#       #select day being modeled.
    swe2model.all=swedata[swedata$yrdoy==yj,]
    swe2model.all$snotel[swe2model.all$snotel==0]=runif(1,0,0.00001)
#
    for (ryr in seq(2000,(2000+numryrs-1))){
        print(paste('ryr = ',ryr))
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
        for(dropiter in seq(1,numstations)){
            z2keep=seq(1,numstations,1)[-dropiter]
                                        # x.scaled=scale(swe2model.all[z2keep,c('recon','Lat','Long','Elev','Eastness','Northness','Slope','RegionalSlope','RegionalEastness','RegionalNorthness','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')])
                                        #  swe2model=swe2model.all[z2keep,c('Station_ID','date','snotel','yrdoy','dy','mth','yr')]
            swe2model=swe2model.all[z2keep,]
            swe2model$snotel=log(swe2model$snotel)
            swe2model$snotel[is.infinite(swe2model$snotel)] <- NA
                                        #    swe2model=data.frame(swe2model,x.scaled)
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
} 


#
## In a much smaller loop, find swediff for recon vs snotel
swediff.reconmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
swediffpct.reconmodel=array(NA,dim=c(mdldoys,nummdlyrs,numryrs,numstations))
for(yj in valvec){
        print(paste('iteration: ',yj,sep=''))
#       
        swe2model.all=swedata[swedata$yrdoy==yj,]
        swe2model.all$snotel[swe2model.all$snotel==0]=runif(1,0,0.00001)
        yobs=swe2model.all$snotel
        #
        for (ryr in seq(2000,(2000+numryrs-1))){
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
        }
    }

pn='diagnostics'
fullpath=file.path(pn,paste0('rswe_',rversion))
save(file=file.path(fullpath,'cvdrop1.skillscores.xrecon.regression.gauss.snotellog.reconlog.RData',
     list=c('r2.phvmodel','r2.phvrcnmodel','r2.reconmodel',
         'swediff.phvmodel','swediff.phvrcnmodel','swediff.reconmodel',
         'swediffpct.reconmodel','swediffpct.phvmodel','swediffpct.phvrcnmodel'))



## # ****** Reassign names so we can use existing code.
## r2phvmed=r2.phvmodel
## r2phvrcnmed=r2.phvrcnmodel
## r2reconmed=r2.reconmodel
## #
## swediff.phv.avg=swediff.phvmodel
## swediff.phvrcn.avg=swediff.phvrcnmodel
## swediffpct.phv.avg=swediffpct.phvmodel
## swediffpct.phvrcn.avg=swediffpct.phvrcnmodel

## ##### Dataframe of skill scores by yr, mth, and ryr for boxplot
## r2df=data.frame(yrdoy=rep(valvec,times=numryrs))
## r2df$reconyear=rep(seq(2000,(2000+numryrs-1)),each=length(unique(r2df$yrdoy)))
## #
## r2df=cbind(r2df,r2phv=as.vector(r2phvmed),r2phvrcn=as.vector(r2phvrcnmed),r2recon=as.vector(r2reconmed))
## r2df$date=as.POSIXct(strptime(r2df$yrdoy,'%Y%j',tz='MST'))
## r2df$mnth=factor(strftime(r2df$date,'%b'),levels=levels(swedata$mnth))


## #
## #Find the mean r2 for phv and recon. Each recon year has numstations iterations for each timestep. the mean is used to represent r2 skill for a given timestep since phv and recon modesl don't change with the rcnyear used.
## r2dfmelt=melt(r2df,id=c('date','mnth','yrdoy','reconyear'))
## r2avg=dcast(r2dfmelt,yrdoy~variable,fun.aggregate=mean,na.rm=T)

## ggplot(r2dfmelt)+
##     geom_boxplot(aes(x=date,y=value,colour=variable))




## r2.max=ddply(.data=r2df,.(yrdoy),summarize,
##     reconyear=which.max(r2phvrcn)+1999,
##     maxr2=max(r2phvrcn,na.rm=T))

## r2.max.recon=ddply(.data=r2df,.(yrdoy),summarize,
##     reconyear=which.max(r2recon)+1999,
##     maxr2=max(r2recon,na.rm=T))
        
## names(r2.max)=c('mth','yr','modelnum','maxr2')
## names(r2.max.recon)=c('mth','yr','modelnum','maxr2')
## r2.max$reconyear=r2.max$modelnum+1999
## r2.max.recon$reconyear=r2.max.recon$modelnum+1999
## r2.max=r2.max[order(r2.max$yr),] #places it in order mth, yr same as the stack. assumes factor mth is ordered correctly.
## r2.max.recon=r2.max.recon[order(r2.max.recon$yr),] #places it in order mth, yr same as the stack. assumes factor mth is ordered correctly.
## rcnyr.phvrcn[yrd-1999,mthi]=subset(r2.max,mth==mthd & yr==yrd)$reconyear
## rcnyr.recon[yrd-1999,mthi]=subset(r2.max.recon,mth==mthd & yr==yrd)$reconyear
## mthi=mthi+1



## bestphvrcn=data.frame(rcnyr.phvrcn)
## rownames(bestphvrcn) <- 2000:2011
## colnames(bestphvrcn) <- c('Mar','Apr','May')
## write.table(x=bestphvrcn,file='../swe.validation/best_recon.phvrcn.drop1.txt',sep=' ')

## bestrecon=data.frame(rcnyr.recon)
## rownames(bestrecon) <- 2000:2011
## colnames(bestrecon) <- c('Mar','Apr','May')
## write.table(x=bestphvrcn,file='../swe.validation/best_recon.recon.drop1.txt',sep=' ')




## ## save.image('cvdrop1.xrecon.regression.gauss.snotelunscaled.skill.image.RData')
## load('cvdrorp1.xrecon.regression.gauss.snotelunscaled.skill.image.RData')






## ##############-------------- MAPS ---------------------------

## require(lubridate)
## states=map_data('state')
## #### Map of Avg SWE Difference in drop10 prediction
## # and percent swe difference
## swediffmap=data.frame()
## for(yri in 2000:2011){
##     for(mthi in 3:5){
##         mths=sprintf('0%d',mthi)
##         mths=strftime(ymd(paste(yri,mths,'01',sep=''),tz='MST'),'%b')
##         for(ryri in 2000:2011){
##             df=data.frame(mth=mths,
##                 yr=yri,
##                 reconyear=ryri,
##                 swediff.phv=swediff.phv.avg[yri-1999,mthi-2,ryri-1999,],
##                 swediff.phvrcn=swediff.phvrcn.avg[yri-1999,mthi-2,ryri-1999,],
##                 swediff.recon=swediff.reconmodel[yri-1999,mthi-2,ryri-1999,],
##                 swediffpct.phv=swediffpct.phv.avg[yri-1999,mthi-2,ryri-1999,],
##                 swediffpct.phvrcn=swediffpct.phvrcn.avg[yri-1999,mthi-2,ryri-1999,],
##                 swediffpct.recon=swediffpct.reconmodel[yri-1999,mthi-2,ryri-1999,],
##                 lat=snotel.locs$lat,
##                 long=snotel.locs$long)
##             swediffmap=rbind(swediffmap,df)
##         }
##     }
## }
## #
## #

## ## Get shaded relief for background
## relief=raster('shaded.relief.grey/GRAY_HR_SR.tif')
## e=extent(x=c(-112.5,-104.25),y=c(32,44))
## ucorelief=crop(relief,e)
## ucorelief.agg=aggregate(ucorelief,2)
## ucorelief.df=data.frame(rasterToPoints(ucorelief.agg))
## names(ucorelief.df)=c('long','lat','ter')
## #
## ## Get info about best rcnyr for phvrcn model from drop1 r2
## bestmodel=read.table('~/Documents/R/swe.validation/best_recon.phvrcn.drop1.txt')




## ##--------------- PERCENT DIFFERENCES IN SWE INSTEAD OF SWE DIFF (METERS)
## dev.set(2)
## ## Map regression phv % swe differences with observed facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## swediffmapplot$swediffpct.phv=swediffmapplot$swediffpct.phv*100
## swediffmapplot$cuts=cut(swediffmapplot$swediffpct.phv,breaks=seq(-200,200,length.out=17))
## #
## swediffmapplot$extreme=NA
## swediffmapplot$extreme[swediffmapplot$swediffpct.phv > 200] ='> 200%'
## swediffmapplot$extreme[swediffmapplot$swediffpct.phv < -200] = '< -200%'
## #
## extremesz=2
## szmin=2
## szstep=1
## szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
## #
## gphv=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,swediffpct.phv <= 200 | swediffpct.phv >= -200),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,swediffpct.phv > 200 | swediffpct.phv < -200),aes(x=long,y=lat,shape=extreme),size=3,fill='black',alpha=0.75)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
##     scale_shape_manual(values=c(25,24))+
##     guides(shape=guide_legend('Extremes'),
##            alpha=F,
##            size=guide_legend('SWE Differences (%)'),
##            colour=guide_legend('SWE Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
##     ggtitle('April 1 Differences with Observed SNOTEL SWE\nRegression with only Physiographics')+facet_wrap(~yr)
## show(gphv)

## ## ggsave(plot=gphv,file='plots/swediffpct.phv.map.pdf',height=12,width=12)

## dev.set(3)
## ## Map regression with rcn % swe differences with observed facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr & yr==2002)
## swediffmapplot$swediffpct.phvrcn=swediffmapplot$swediffpct.phvrcn*100
## swediffmapplot$cuts=cut(swediffmapplot$swediffpct.phvrcn,breaks=seq(-200,200,length.out=17))
## #
## swediffmapplot$extreme=NA
## swediffmapplot$extreme[swediffmapplot$swediffpct.phvrcn > 200] ='> 200%'
## swediffmapplot$extreme[swediffmapplot$swediffpct.phvrcn < -200] = '< -200%'
## #
## extremesz=2
## szmin=2
## szstep=1
## szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts))-2)/2)
## #
## gphvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,swediffpct.phvrcn <= 200 | swediffpct.phvrcn >= -200),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,swediffpct.phvrcn > 200 | swediffpct.phvrcn < -200),aes(x=long,y=lat,shape=extreme),size=3,fill='black',alpha=0.75)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
##     scale_shape_manual(values=c(25,24))+
##     guides(shape=guide_legend('Extremes'),
##            alpha=F,
##            size=guide_legend('SWE Differences (%)'),
##            colour=guide_legend('SWE Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
##     ggtitle('April 1 Differences with Observed SNOTEL SWE\nRegression with Physiographics and RCN')+facet_wrap(~yr)
## show(gphvrcn)

## ## ggsave(plot=gphvrcn,file='plots/swediffpct.phvrcn.2001.map.pdf',height=12,width=10)




## dev.set(4)
## ## Map reconstruction % swe differences with observed. facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr')
## swediffmapplot$swediffpct.recon=swediffmapplot$swediffpct.recon*100
## swediffmapplot$cuts=cut(swediffmapplot$swediffpct.recon,breaks=seq(-200,200,length.out=17))
## #
## swediffmapplot$extreme=NA
## swediffmapplot$extreme[swediffmapplot$swediffpct.recon > 200] ='> 200%'
## swediffmapplot$extreme[swediffmapplot$swediffpct.recon < -200] = '< -200%'
## #
## extremesz=2
## szmin=2
## szstep=1
## szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
## #
## grecon=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,swediffpct.recon <= 200 | swediffpct.recon >= -200),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,swediffpct.recon > 200 | swediffpct.recon < -200),aes(x=long,y=lat,shape=extreme),size=3,fill='black',alpha=0.75)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
##     scale_shape_manual(values=c(25,24),drop=F)+
##     guides(shape=guide_legend('Extremes'),
##            alpha=F,
##            size=guide_legend('SWE Differences (%)'),
##            colour=guide_legend('SWE Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),
##           legend.background=element_rect(fill='grey80'))+
##     ggtitle('April 1 Differences with Observed SNOTEL SWE\nReconstruction')+
##     facet_wrap(~yr)
## show(grecon)

## ## ggsave(plot=grecon,file='plots/swediffpct.recon.map.pdf',height=12,width=12)


## facet1_names=list(
##     'swephvpct.phv'='Regression w/o RCN',
##     'swephvpct.phvrcn'='Regression w/ RCN',
##     'swephvpct.recon'='Reconstruction')
## plot_labeller <- function(variable,value){
##     if(variable=='yr'){
##         return(value)
##     } else {
##         return(facet1_names[value])
##     }
##  }

## #### Combined mapping - rows = modeltype, columns 4 years.
## numf=11
## outlier=250
## swediffmapplot=melt(swediffmap,id=c('mth','yr','reconyear','lat','long','swediff.phv','swediff.phvrcn','swediff.recon'))
## #
## mnth='Apr'
## mthcol=which(names(bestmodel)==mnth)
## swediffmapplot=subset(swediffmapplot,(yr==2002 & reconyear==bestmodel[2002-1999,mthcol]) | (yr==2000 & reconyear==bestmodel[2000-1999,mthcol]) | (yr==2011 & reconyear==bestmodel[2011-1999,mthcol]))
## swediffmapplot=subset(swediffmapplot,as.character(mth)==mnth)
## swediffmapplot$value=swediffmapplot$value*100
## swediffmapplot$value[swediffmapplot$value > outlier]=outlier-1
## swediffmapplot$value[swediffmapplot$value < -outlier]=-outlier+1
## swediffmapplot$cuts=cut(swediffmapplot$value,breaks=seq(-outlier,outlier,length.out=numf),right=F,include.lowest=T)
## #
## legendlabels=levels(swediffmapplot$cuts)
## outlabel=outlier*2/(numf-1)*((numf-1)/2-1)
## legendlabels[1]=paste('< ',-outlabel,sep='')
## legendlabels[length(levels(swediffmapplot$cuts))]=paste('>= ',outlabel,sep='')
## #
## szmin=2
## szstep=.5
## szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2)
## #
## gphvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,value < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,value >= 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F,labels=legendlabels)+
##     scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,labels=legendlabels)+
##     labs(x='Longitude',y='Latitude')+
##     guides(alpha=F,
##            size=guide_legend('SWE Differences (%)'),
##            colour=guide_legend('SWE Differences (%)'),
##            shape=guide_legend('SWE Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),
##           legend.key.size=unit(1.5,'lines'),
##           legend.text=element_text(size=14),
##           legend.title=element_text(size=16),
##           ## legend.justification='right',
##           legend.background=element_rect(fill='grey80'),
##           legend.position='right',
##           axis.text=element_text(size=14),
##           plot.title=element_text(size=18),
##           axis.title=element_text(size=16),
##           strip.text=element_text(size=14,face='bold'))+
##     ggtitle('April 1 % Differences with Observed SNOTEL SWE\n')+
##     facet_grid(variable~yr,labeller=plot_labeller)
## show(gphvrcn)
## ggsave(plot=gphvrcn,'plots/allmodels.pctdiff.2002.2000.2011.png',dpi=300,width=11,height=10)


## #********************************
## ## Map differences between %differences of the regression models facet by year
## mnth='Apr'
## mthcol=which(names(bestmodel)==mnth)
## swediffmapplot=subset(swediffmap,(yr==2002 & reconyear==bestmodel[2002-1999,mthcol]) | (yr==2000 & reconyear==bestmodel[2000-1999,mthcol]) | (yr==2011 & reconyear==bestmodel[2011-1999,mthcol]))
## swediffmapplot=subset(swediffmapplot,as.character(mth)==mnth)
##      #
## numf=11
## #
## swediffmapplot$mdldiffpct.phv.phvrcn=(abs(swediffmapplot$swediffpct.phv)-abs(swediffmapplot$swediffpct.phvrcn))*100
## swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn > outlier]=outlier-1
## swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn < -outlier]=-outlier+1
## swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.phv.phvrcn,breaks=seq(-outlier,outlier,length.out=numf),right=F, include.lowest=T)
## #
## legendlabels=levels(swediffmapplot$cuts)
## outlabel=outlier*2/(numf-1)*(numf-1)/2-outlier*2/(numf-1)
## legendlabels[1]=paste('< ',-outlabel,sep='')
## legendlabels[length(levels(swediffmapplot$cuts))]=paste('>= ',outlabel,sep='')
## #
## extremesz=4
## szmin=3
## szstep=1.5
## szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2) 
## #
## gdiff.phv.phvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn >= 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F, labels=legendlabels)+
##     scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,
##                       labels=legendlabels)+
##     labs(x='Longitude', y='Latitude')+
##     guides(alpha=F,
##            size=guide_legend('Error Differences (%)'),
##            colour=guide_legend('Error Differences (%)'),
##            shape=guide_legend('Error Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),
##           legend.key.size=unit(1.5,'lines'),
##           legend.background=element_rect(fill='grey80'),
##           legend.text=element_text(size=14),
##           legend.title=element_text(size=16),
##           axis.text=element_text(size=14),
##           plot.title=element_text(size=18),
##           axis.title=element_text(size=16),
##           strip.text=element_text(size=14,face='bold'))+
##     facet_wrap(~yr)
##     ## ggtitle(paste('April 1 Differences between % Error of Regression w/o RCN and Regression w/ RCN\nRelative to Observed SNOTEL SWE\nBlue triangles indicate a smaller error with Regression w/ RCN',sep=''))+
##    #
## show(gdiff.phv.phvrcn)
##  ## ggsave(plot=gdiff.phv.phvrcn,file='plots/mdldiff.phv.phvrcn.2000.2002.2011.png',dpi=300,width=14, height=6.5)
## ## }




## ## Map differences between %differences of phvrcn and recon models. facet by year
## outliers=100
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## swediffmapplot$mdldiffpct.recon.phvrcn=(abs(swediffmapplot$swediffpct.recon)-abs(swediffmapplot$swediffpct.phvrcn))*100
## swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.recon.phvrcn,breaks=seq(-outliers,outliers,length.out=21))
## #
## swediffmapplot$extreme=NA
## swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phvrcn > outliers] = paste('> ',outliers,'%',sep='')
## swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phvrcn < -outliers] = paste('< ',-outliers,'%',sep='')
## #
## extremesz=1
## szmin=1
## szstep=.5
## szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
## #
## gdiff.recon.phvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phvrcn <= outliers | mdldiffpct.recon.phvrcn >= -outliers),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phvrcn > outliers | mdldiffpct.recon.phvrcn < -outliers),aes(x=long,y=lat,shape=extreme),size=extremesz,fill='black',alpha=0.75)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
##     scale_shape_manual(values=c(25,24),drop=F)+
##     guides(shape=guide_legend('Extremes'),
##            alpha=F,
##            size=guide_legend('Error Differences (%)'),
##            colour=guide_legend('Error Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
##     ggtitle('April 1 Differences between Reconstruction and Regression with PHV/RCN\nRelative to Observed SNOTEL SWE\nBlue points indicate a smaller error with Regression w/ RCN')+
##     facet_wrap(~yr)
## show(gdiff.recon.phvrcn)

## ## ggsave(plot=gdiff.recon.phvrcn,file='plots/modeldiffpct.recon.phvrcn.map.pdf',height=12,width=12)

## dev.new()
## ## Map differences between %differences of phv and recon models. facet by year
## outliers=100
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## swediffmapplot$mdldiffpct.recon.phv=(abs(swediffmapplot$swediffpct.recon)-abs(swediffmapplot$swediffpct.phv))*100
## swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.recon.phv,breaks=seq(-outliers,outliers,length.out=21))
## #
## swediffmapplot$extreme=NA
## swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phv > outliers] = paste('> ',outliers,'%',sep='')
## swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phv < -outliers] = paste('< ',-outliers,'%',sep='')
## #
## extremesz=1
## szmin=1
## szstep=.5
## szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
## #
## gdiff.recon.phv=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phv <= outliers | mdldiffpct.recon.phv >= -outliers),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phv > outliers | mdldiffpct.recon.phv < -outliers),aes(x=long,y=lat,shape=extreme),size=extremesz,fill='black',alpha=0.75)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
##     scale_shape_manual(values=c(25,24),drop=F)+
##     guides(shape=guide_legend('Extremes'),
##            alpha=F,
##            size=guide_legend('Error Differences (%)'),
##            colour=guide_legend('Error Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
##     ggtitle('April 1 Differences between Reconstruction and Regression with only PHV\nRelative to Observed SNOTEL SWE\nBlue points indicate a smaller error with Regression with PHV')+
##     facet_wrap(~yr)
## show(gdiff.recon.phv)

## ## ggsave(plot=gdiff.recon.phvrcn,file='plots/modeldiffpct.recon.phvrcn.map.pdf',height=12,width=12)







## str(swediffmap)

## swehist=swediffmap
## swehist[,c('swediffpct.phv','swediffpct.phvrcn','swediffpct.recon')]=swehist[,c('swediffpct.phv','swediffpct.phvrcn','swediffpct.recon')]*100
## #
## swehist=subset(swehist,yr==reconyear & mth=='Apr')
## avgs=ddply(swehist[,c('mth','yr','swediffpct.phv','swediffpct.phvrcn','swediffpct.recon')],.(yr,mth),summarize,avg.phv=mean(swediffpct.phv,na.rm=T),avg.phvrcn=mean(swediffpct.phvrcn,na.rm=T),avg.recon=mean(swediffpct.recon,na.rm=T))
## #
## avgs=subset(avgs,mth='Apr')
## #
## avg.phv=NULL
## avg.phvrcn=NULL
## avg.recon=NULL
## for(i in 1:12){
##     avg.phv=rbind(avg.phv,rep(avgs$avg.phv[i],237))
##     avg.phvrcn=rbind(avg.phvrcn,rep(avgs$avg.phvrcn[i],237))
##     avg.recon=rbind(avg.recon,rep(avgs$avg.recon[i],237))
## }

## swehist=cbind(swehist,avg=c(as.vector(t(avg.phv)),as.vector(t(avg.phvrcn)),as.vector(t(avg.recon))))
## swehistplot=melt(swehist,id=c('mth','yr','reconyear','lat','long','swediff.phv','swediff.phvrcn','swediff.recon','avg'))
## swediffhistplot1=subset(swehistplot, yr<2004)

## ggplot(data=swediffhistplot1)+geom_histogram(aes(x=value),binwidth=1000)
## +
##     facet_grid(variable~yr)



## ggplot()+geom_histogram(data=swediffmap,aes(x=swediffpct.phvrcn),binwidth=0.05)+facet_wrap(~yr)

## ggplot()+geom_histogram(data=swediffmap,aes(x=swediffpct.phv),binwidth=0.05)+facet_wrap(~yr)

## ggplot()+geom_histogram(data=swediffmap,aes(x=swediffpct.recon),binwidth=0.05)+facet_wrap(~yr)






## #Bias
## bias.phvrcn=mean(swediffmap$swediff.phvrcn,na.rm=T)
## bias.phv=mean(swediffmap$swediff.phv,na.rm=T)
## bias.recon=mean(swediffmap$swediff.recon,na.rm=T)

## #MAE
## mae.phvrcn=sum(abs(swediffmap$swediff.phvrcn),na.rm=T)/nrow(swediffmap)#0.11
## mae.phv=sum(abs(swediffmap$swediff.phv),na.rm=T)/nrow(swediffmap)#0.13
## mae.recon=sum(abs(swediffmapplot$swediff.recon),na.rm=T)/nrow(swediffmapplot)#0.20

## # phv - phvrcn
## mdldiff.phv.phvrcn=abs(swediffmap$swediff.phv)-abs(swediffmap$swediff.phvrcn)
## mean(mdldiff.phv.phvrcn,na.rm=T)
## ggplot()+geom_histogram(aes(x=mdldiff.phv.phvrcn),binwidth=0.05)
## #
## # recon - phvrcn
## mdldiff.recon.phvrcn=abs(swediffmap$swediff.recon)-abs(swediffmap$swediff.phvrcn)
## mean(mdldiff.recon.phvrcn,na.rm=T)
## ggplot()+geom_histogram(aes(x=mdldiff.recon.phvrcn),binwidth=0.05)
## #
## # recon - phv
## mdldiff.recon.phv=abs(swediffmap$swediff.recon)-abs(swediffmap$swediff.phv)
## mean(mdldiff.recon.phv,na.rm=T)
## ggplot()+geom_histogram(aes(x=mdldiff.recon.phv),binwidth=0.05)
