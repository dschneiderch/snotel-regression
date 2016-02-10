lm_eqn = function(m) {
	yobs=m$model[[1]]
	yhat=fitted(m)
	l <- list(a = format(coef(m)[1], digits = 2),
						b = format(abs(coef(m)[2]), digits = 2),
						r2 = format(cor(yobs,yhat)^2, digits = 3),
						pval = format(anova(m,test='F')$P[2],digits=2));
	
	if (coef(m)[2] >= 0)  {
		eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pval,l)
	} else {
		eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(p)~"="~pval,l)    
	}
		as.character(as.expression(eq));                 
}
