read_surveydates=function(fn){

tmp=scan(fn,comment.char='#',what=character(),sep='\t')
dte=tmp[seq(1,length(tmp),2)]
site=tmp[seq(2,length(tmp),2)]
data.frame(date=as.Date(dte),site)

}