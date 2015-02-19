precip_adjust=function(dF,meth){
     
     #      % Assumes that the data have been cleaned up so that spurious "negative-positive" pairs of precipitaiton data,  and other obvious artifacts are cleaned up
     #      %
     #      
     #      %  USAGE
     #      % y=dlmread('jackwhacker_gulch.txt','',1,0);
     #      % swe = y(:,2);ap=y(:,3);pp=y(:,7);
     #      %  phat=precip_adjust(swe,ap,pp,1);
     #      %
     #      %
     #      %  swe is a timeseries of SWE, with no missing values
     #      %  ap  is a timeseries of accumulate precipitaiton for the same time interval
     #      %  pp is the timeseries of precipitation 
     #      %
     #      %  meth   method
     #      %    meth = 1     use  maximum of either precipitation or dSWE for each day only for days on which precip and swe both increase.
     #      %    meth = 2     use  maximum of either precipitation or dSWE for each day for all days
     #      %    meth = 3     use  maximum of either precipitation or dSWE for each day averaged over the n-dayperiod.  Then add back in daily SWE/p adjustment
     #      %    meth = 4     use Algorithm  use dSWE if SWE >2 and is increasing
     #      
     
     #fix missing precip end of every water year
     dF[dF$mth==9 & dF$dy==30,'precip']=0
     
     #get variables
     swe=dF$swe
     ap=dF$apcp
     pp=dF$precip
     #numer of whole years
     ny=length(unique(dF$wy))
     
     dswe=c(0,diff(swe)) #calculate delta-SWE and add back in one day at the end to make it the same as the precip record. 
     dp=c(0,diff(ap)) #Same for accumulated precip
     
     # Quality control check that accumulated precipitation is always increasing. (except at end of water year when it is reset to zero. 
     if(length(which(dp<0)) > ny){
          print('Negative values in delta-accumulated-precip -- more than simply at the end of each year.  Please clean up data') 
          stop()
     }          
     phat=rep(0,nrow(dF))
     swe_thresh=0
     s0 = swe>=swe_thresh # %  Find days when SWE >= swe_thresh inches
     
     #      %  calculate everything in terms of precipitation (as input)
     
     if (meth == 1){
          f0=which(!s0)
          phat[f0]=pp[f0]#  %  SWE < swe_thresh  use change in accumulated precipitation.
          f1=which(s0 & dswe<=0 & pp>0)
          phat[f1]=pp[f1]#          % SWE decreasing and Precip >0  =>  use Precip
          f2=which(s0 & dswe>0 & pp>0)
          phat[f2]=apply(matrix(c(dswe[f2],pp[f2]),ncol=2),1,max)#   %  SWE > swe_thresh and both SWE and Precip increase, then use the maximum of the two.
     }
     if (meth == 2){     
          f0=which(!s0)
          phat[f0]=pp[f0]    #  SWE < swe_thresh  use change in accumulated precipitation (same as precip)
          f1=which(s0 & dswe<=0 & pp>0)
          phat[f1]=pp[f1]#          % SWE decreasing and Precip >0  =>  use Precip
          f2=which(s0)
          phat[f2]=apply(matrix(c(dswe[f2],pp[f2]),ncol=2),1,max)#   %  SWE > swe_thresh  then use the maximum of the two.
     }
     if(meth == 3){
          nseg=7#   % chunks over which SWE vs P is analyzed and maximum taken
          phat=pp#   % set adjusted precip to observed precip to start
          
          for (k in seq(1,(length(phat)-nseg),10)){#       % loop over timeseries one nseg-chunk of days at at time
               kseg = k+seq(1,nseg)-1
               if ( sum(dswe[kseg]) > sum(pp[kseg]) ){#   % if the total change in SWE is greater than the total precipitation during that period, use the SWE record
                    phat[kseg]=dswe[kseg]
               }
          }
          # daily corrections 
          
          f1=which(s0 & dswe<=0 & pp>0)
          phat[f1]=pp[f1] #  SWE decreasing and Precip >0  =>  use Precip       
          f2=which(s0 & dswe>0 & phat>0)
          phat[f2]=apply(matrix(c(dswe[f2],phat[f2]),ncol=2),1,max) # SWE > swe_thresh and both SWE and Precip increase, then use the maximum of the two.
          f3=which(phat<0)
          phat[f3]=0# phat could still by <0 if it was taken from a dSWE "chunk" and had decreasing swe without precipitation >0.  Correct for this. 
     }
     if( meth == 4){
          f0=which(!s0  & pp>0)
          phat[f0]=pp[f0]#   SWE < swe_thresh  and precip > 0  use change in accumulated precipitation.
          f1=which(s0 & dswe<=0 & pp>0)
          phat[f1]=pp[f1]# SWE decreasing and Precip >0  =>  use Precip
          f2=which(s0 & dswe>0 & pp>0)
          phat[f2]=dswe[f2]#  SWE > swe_thresh and both SWE and Precip increase, then use SWE.
     }
     
     #           %  Calculate everything in terms of accumulated precip
     #           %f0=find (~s0  & dp>0);  phat(f0)=dp(f0);    %  SWE < 1 inch +> use change in accumulated precipitation.
     #           %f1=find (s0 & dswe<=0 & dp>0);phat(f1)=dp(f1);          % SWE decreasing and Precip >0  =>  use Precip
     #           %%f2=find (s0 & dswe>0 & dp>0);phat(f2)=max(dswe(f2),dp(f2));   %  SWE > 1 and both SWE and Precip increase, then use the maximum of the two.
     #           %f2=find (s0 & dswe>0 & dp>0);phat(f2)=dswe(f2);   %  SWE > 1 and both SWE and Precip increase, then use the maximum of the two.
     #           
     dF$precip=phat
     dF=ddply(dF,.(wy),function(dF) {
          dF$apcp=cumsum(dF$precip) 
          return(dF)
     })
          
     return(dF)
}