generate_swedata=function(snotelrecon){
          
# scale the x variables -----------------------------------------
phv_scatt=scale(phvrec[,-ncol(phvrec)])
phvrec.sc=data.frame(phv_scatt,Station_ID=phvrec$Station_ID)
phvrec.sc=phvrec.sc[order(phvrec$Station_ID),]

# create temporally continuous snotelrecon record ----------------
temp1=snotelrec[,c('swe','date','Station_ID')]
temp2=snotelrecon[,c('recon','date','Station_ID')]
swedata=merge(temp1,temp2,all=T)
swedata$yrdoy=as.numeric(strftime(swedata$date,'%Y%j'))
swedata$dy=as.numeric(strftime(swedata$date,'%d'))
swedata$mth=as.numeric(strftime(swedata$date,'%m'))
swedata$yr=as.numeric(strftime(swedata$date,'%Y'))
#
originaldata=merge(swedata,phvrec,by='Station_ID',all=T)
colnames(originaldata) <- gsub(pattern='swe',replacement='snotel',x=colnames(originaldata))
#
swedata=merge(swedata,phvrec.sc,by='Station_ID',all=T)
colnames(swedata) <- gsub(pattern='swe',replacement='snotel',x=colnames(swedata))
#
#difference between swedata and original data is phvrec was scaled in swedata
swedata$mnth=factor(strftime(swedata$date,'%b'),levels=unique(strftime(swedata$date[order(swedata$date)],'%b')))
originaldata$mnth=swedata$mnth

return(swedata)
}