fixSWE = function(dF,colnm){
print(unique(dF$Station_ID))
     colnm='swe'#keep for prosperity
     # user adjustments
     if(unique(dF$Station_ID)=='05P04S' & unique(dF$wy)==2012){
          dF$precip[183:184]=0
          dF$precip[185]=0.01016
          dF$apcp[183:185]=0.24892
     }
     
     #fix missing precip end of every water year
     dF[dF$mth==9 & dF$dy==30,'precip']=0
     
     #check for negative/positive pairs in precip 
     negind=which(dF$precip<0)
     for( d in negind){
          #           print(d)
          negp=dF$precip[d]
          dF$precip[d]=0
          #dF$apcp[(d-1):nrow(dF)]=dF$apcp[(d-1):nrow(dF)]
          for( dcount in seq(-7,7,1)){
               di=d+dcount
               #                print(paste('---',di))
               if(identical(abs(negp),dF[di,'precip'])){
                    dF$precip[di]=0
                    dF$apcp=dF$apcp[di-1]
                    break
               }
               if(dcount==7){
                    dF$apcp[(d+1):nrow(dF)]=dF$apcp[(d+1):nrow(dF)]+abs(negp) #if doesn't find a positive counterpart, then add the negp back into ap     
               }
          }     
     }
     
     #check for and fix swe NA
     z=dF[,colnm]
     print(unique(dF$wy))
     #check beginning of record for NA. find first non-NA value
     if(is.na(z[1])) {
          for(k in 2:nrow(dF)){
               if(!is.na(z[k])){
                    zref1=z[k]
                    print('warning: missing data at beginning of ts')
                    break
               }
          }
     }
     #step through vector of colnm. for each NA in SWE, fill with swe+precip recorded if in accumulation season
     if(!exists('k')) k=2
     #      print(k)
     for(i in k:nrow(dF)){
          if(is.na(z[i])){#stop where z is na
               zref1=z[i-1]#set 1st endpoint
               print(i)
               #print(zref1)
               for(n in 1:nrow(dF)){ #find first non-na value.
                    #                     print(n)
                    #print(i+n)
                    if(i+n > length(z)){#break out of loop if the end of record
                         zref2=NA
                         print('warning: missing data at end of time series not filled')
                         break
                    }
                    if(!is.na(z[i+n])){#stop where z is no longer NA
                         zref2=z[i+n]#set second endpoint
                         nref=n
                         dswe=dF$precip[(i):(i+n-1)]#get the precip for days where swe was NA
                         #check for accumulation season. in melt season NA can be linearly interpolated with precip added in. in accum season, swe is adjusted by precip. if precip is missing then the change in apcp is used.
                         if(i<=14) {
                              accumseason=TRUE
                         } else if( (z[i-14]-z[i-1]) < 0 ){
                              accumseason=TRUE
                         } else { accumseason=FALSE }
                         if( accumseason ){#accumulation season
                              if( length(which(is.na(dswe)))==length(dswe) ){#sometimes the whole snotel station goes down, everything is NA so the change in apcp is used.
                                   dswe=dF$apcp[i+n]-dF$apcp[i-1]
                              }
                              z[i:(i+n-1)]=zref1+dswe 
                         } else {#ablation season
                              for (n in 1:nref){
                                   if( length(which(is.na(dswe)))==length(dswe) ){
                                        z[i+n-1]=zref1+(n/(nref+1))*(zref2-zref1)
                                   } else {
                                        z[i+n-1]=zref1+(n/(nref+1))*(zref2-zref1)+dswe[n]
                                   }
                              }
                         }
                         break
                    }
               }
          }     
     }
     dF[,colnm]=z
     #
     print(paste0('precip < 0: ',length(which(dF$precip<0))))
     print(paste0('swe NA: ',length(which(is.na(dF$swe)))))
     return(dF)
     
}