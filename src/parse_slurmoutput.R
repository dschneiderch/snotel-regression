setwd('~/Documents/snotel-regression_project')
library(plyr)

recon.version='v3.1'
covrange='idp1'
config='v2.4'

f2r=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange)
fns=list.files(f2r,pattern='.Sout$',full.names=T)
fn=fns[1]
for(fn in fns){
	print(fn)
	temp=readLines(fn)	
	# print(head(temp))
	for(i in 1:length(temp)){
		if(grepl('mdlyr:',temp[i])){
			#print(temp[i])
			starti=i
			oi=starti+3
			if(substr(temp[oi],6,7)=='r2'){ 
				cost=substr(temp[oi], 6,7)
			} else { 
				cost=substr(temp[oi],6,9)
			}
			sci=starti+6
			if(substr(temp[sci],6,7)=='sc'){ 
				scalesnotel=substr(temp[sci], 6,10)
			} else { 
				scalesnotel=substr(temp[sci],6,12)
			}
			fi=starti+7
			if(substr(temp[fi],6,7)=='fs'){
				fscaMatch=substring(temp[fi],6,9)
			} else {
				fscaMatch=substring(temp[fi],6,11)
			}
			# next
		}
		if(grepl('size',temp[i])){
			keep=temp[starti:(i-1)]
			writeLines(keep,paste0(fn,'-',scalesnotel,'-',fscaMatch,'-',cost,'.txt'))
		}	
	}
}