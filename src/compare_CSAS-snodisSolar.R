library('ProjectTemplate')
setwd('~/Documents/snotel-regression_project')
load.project()

sbsp=read.fwf('data/metvalidation/SenatorBeck/SBSP-forcing_seriallycomplete.txt',
		widths=c(4,2,2,2,2,2,6,14,10,10,10,10,10,10,10,10,10,10,10,10,4,4,4,4,4,4),
		sep='')
sasp=read.fwf('data/metvalidation/SenatorBeck/SASP-forcing_seriallycomplete.txt',
		widths=c(4,2,2,2,2,2,6,14,10,10,10,10,10,10,10,10,10,10,10,10,4,4,4,4,4,4),
		sep='')
# colnames(sbsp) <- colnames(sasp) <- c('yr','mth','dy','hour','minute','second','precip','SWin','LWin','airtemp','wind','relhum','pressure','spechum','dewtemp','precip-wmo','at-kent','at-anderson','at-nakamura','at-huwald')#,'qc-precip','qc-SWin','qc-LWin','qc-at','qc-wind','qc-relhum')
valdata$date=as.POSIXct(strptime(paste0(valdata$yr,sprintf('%02d',valdata$mth),sprintf('%02d',valdata$dy),sprintf('%02d',valdata$hour),sprintf('%02d',valdata$minute),sprintf('%02d',valdata$second)),'%Y%m%d%H%M%S'),tz='MST'))
valdata$date=as.POSIXct(strptime(valdata$hour,'%Y%j%H'))

