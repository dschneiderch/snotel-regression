snotelQAQC=function(snotelrec,precipmeth){
     
     snotelrec=ddply(snotelrec,.inform=F,.(Station_ID),function(dF){
     
          #save raw snotelrec for plotting later
          dFraw=dF
          dFraw$version='raw'
          
          ##---------- temperature
          #only 1 sensor so if any measurement is bad, flag all of them.
          #QA temperature
          tempind=with(dF,
                       which(tmin < -40 | tmax < -40 | tavg < -20 | tmin > 40 | tmax > 40 | tavg > 20))
          tempind2=with(dF,
                        which(is.na(tmin) | is.na(tmax) | is.na(tavg)))
          dF$tmin[c(tempind,tempind2)]=NA
          dF$tmax[c(tempind,tempind2)]=NA
          dF$tavg[c(tempind,tempind2)]=NA
          #QC temperature
          dF$tmin=fillMissing(dF$tmin,max.fill=14) #I probably don't need a complete series
          dF$tmax=fillMissing(dF$tmax,max.fill=14)
          dF$tavg=fillMissing(dF$tavg,max.fill=14)
          
          ## ----------- swe and precip
          #set summer SWE NAs to 0
          summerNA=which(dF$mth>=7 & dF$mth<11 & is.na(dF[,'swe']))
          dF[summerNA,'swe']=0
          
          #check for and convert neg swe to 0
          negs=which(dF$swe<0)
          if(length(negs)>0)  dF[negs,'swe']=0
          
          #make continuous swe timeseries 
          print('----------------qc1---------------')
          dFqc1=ddply(dF,.(wy),fixSWE,.inform=T)#does some minor precip QA, then uses precip to fill swe. should never see "end of ts not filled if summerNA sets Sept NAs to 0
          
          #adjust precip to fix some of the undercatch compared to the pillow.
          
          dFqc1=precip_adjust(dFqc1,precipmeth)#script transcribed from J.Barsugli's matlab script. adjusts precip for undercatch based on pillow. several methods available
          dFqc1$version='qc1'
          #write new record to file
          dirname=paste0('data/snotelQAQC')
          fn=file.path(dirname,'qc',paste0(unique(dFqc1$Station_ID),'-qc1_precipmeth',precipmeth,'.txt'))
          write.table(dFqc1[,c('snoteldate','swe','apcp','tmax','tmin','tavg','precip')],fn,sep='\t',row.names=F,quote=F)
          
          ##CAVEAT: precip and aprecip are used to fill SWE. precip_adjust then increases apcp and precip based on max(dswe,precip). since there are no NA's in swe, rerunning fixSWE on dFqc1 wouldn't do anything. rerun fixSWE with raw swe but updated apcp and precip? -- there are so few days where SWE =NA. plus with precipmeth=1 swe needs to increase for precip to be adjusted so it might not change swe much to rerun fixSWE...? maybe more with precipmeth=2
          dFqc2=dFqc1
#                print('----------------------qc2----------------')
#                dFqc2=dF
#                dFqc2$precip=dFqc1$precip
#                dFqc2$apcp=dFqc1$apcp
#                dFqc2=ddply(dFqc2,.(wy),fixSWE,.inform=T)
#                dFqc2$version='qc2'
#                #write new record to file
#                dirname=paste0('data/snotelQAQC')
#                fn=file.path(dirname,'qc',paste0(unique(dFqc2$Station_ID),'-qc2_precipmeth',precipmeth,'.txt'))
#                write.table(dFqc1[,c('snoteldate','swe','apcp','tmax','tmin','tavg','precip')],fn,sep='\t',row.names=F,quote=F)
          
          #combine raw and qc snotelrec and plot together.
          dFall=rbind(dFraw,dFqc1,dFqc2)
          dFplot=melt(dFall[,c('date','swe','apcp','version')],.(date,version))
          g=ggplot(dFplot)+
               geom_line(aes(x=date,y=value,colour=version,linetype=variable),alpha=0.5)+
               labs(title=paste(unique(dF$Station_ID),'\nprecip method',precipmeth))+
               theme_minimal()
          ggsave(plot=g,filename=paste0(dirname,'/plots/snotelrec_',unique(dF$Station_ID),'_precipmeth',precipmeth,'.pdf'),height=3,width=6)
          ## ggplot is giving a warning that 23 rows contain missing values but snotelrec doesn't have any missing values so something in dFplot after melt?
          return(dFqc2[,-ncol(dFqc2)])
     })
     return(snotelrec)
}