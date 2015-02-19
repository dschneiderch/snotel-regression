fixNA = function(dF,colnm){
     dF=dF[order(dF$Station_ID,dF$date),]
     dF$mth=as.numeric(strftime(dF$date,'%m'))
     summerNA=which(dF$mth>=7 & dF$mth<11 & is.na(dF[,colnm]))
     dF[summerNA,colnm]=0
     #print(dF[summerNA,c('date','Station_ID','snoteldate','swe')])
     newdF=dF
     ind=which((dF$mth<7 | dF$mth>=11 ) & is.na(dF[,colnm]))
     print(dF[ind,c('date','Station_ID','snoteldate','swe')])
     if(length(ind)>0){
          winsize=15#needs to be >=10
          dFfix=ldply(as.list(ind),function(i){     
               tmp=dF[(i-winsize):(i+winsize),colnm]
               tmpna=which(is.na(tmp))
               #           print(tmpna)
               y=tmp
               #           print(y)
               if(length(tmpna)>0){
                    mind=diff(which(!is.na(tmp)))
                    i1=which(diff(mind)>(length(tmpna)-1))+1# +1 to account for diff being 1 element shorter
                    i2=i1+length(tmpna)+1# add number of na values to index
                    iv1=mean(tmp[(i1-1):i1])
                    iv2=mean(tmp[(i2+1):i2])
                    m=(iv1-iv2)/(i1-i2)
                    x=tmpna-min(tmpna)+1
                    y[is.na(y)]=m*x+iv1
               }
               data.frame(colnew=y,ind=seq((i-winsize),(i+winsize)))
               #dF[(i-winsize):(i+winsize),colnm]=y
               #           return(dF)
          })
          newdF[dFfix$ind,colnm]=dFfix[,'colnew']
     }
     return(newdF)
     
}