fixNA = function(dF,colnm){
     dF=dF[order(dF$Station_ID,dF$date),]
     ind=which(is.na(dF[,colnm]))
     winsize=10#needs to be >=10
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
          data.frame(snotel=y,ind=seq((i-winsize),(i+winsize)))
          #dF[(i-winsize):(i+winsize),colnm]=y
#           return(dF)
     })
     newdF=dF
     newdF[dFfix$ind,'snotel']=dFfix$snotel
     return(newdF)
}