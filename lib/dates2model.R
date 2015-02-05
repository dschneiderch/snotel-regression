dates2model=function(opt){
     
     # Option A will only model dates for dates selected from modscag images.
     # Option B will model selected dates from modscag for months prior to March and then 1st and 15th of March, April, May

     if(opt=='A'){
          ## option A
          dateselect=read.table('data/selectdates/alldates.txt')#dates selected as clear from modscag.
          dateselect=as.POSIXct(strptime(dateselect$V1,'%d%b%Y',tz='MST'))
          
     } else {
          
          ## option B
          dateselect=read.table('data/selectdates/alldates.txt')#dates selected as clear from modscag.
          dts=data.frame(date=as.POSIXct(strptime(dateselect$V1,'%d%b%Y',tz='MST')))
          dts=mutate(dts,
                     yr=as.numeric(strftime(date,'%Y')),
                     mth=as.numeric(strftime(date,'%m')),
                     dy=as.numeric(strftime(date,'%d')))
          dts=filter(dts,mth>=3 & mth<6)
          make_gendts=function(){
               yr=seq(min(dts$yr),max(dts$yr))
               mth=c(3,4,5)
               dy=c(1,15)
               dF=expand.grid(yr,mth,dy)
               colnames(dF)=c('yr','mth','dy')
               dF$date=strptime(with(dF,paste(yr,mth,dy)),'%Y %m %d',tz='MST')
               return(dF)
          }
          gendts=make_gendts()
          dts=rbind(dts,gendts)
          dateselect=dts$date
     }
     return(dateselect)
}