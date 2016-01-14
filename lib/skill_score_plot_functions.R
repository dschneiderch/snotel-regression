monthlybarplot=function(monthskillplot,yvar){
ggmnth=ggplot(monthskillplot,aes_string(y=yvar))+
    geom_bar(aes(x=variable,fill=modeltype),stat='identity',width=.8,position=position_dodge(width=.8))+
    # facet_grid(.~modeltype)+
    labs(x='Month',y=costlong)+
    theme_minimal(base_size=8)+
    theme(
            strip.text=element_text(size=10,face='bold'),
            axis.line=element_line(),
            legend.position=c(.1,.85),
            legend.justification='left',
            legend.key.size=unit(.5,'lines'))
}


myboxplot=function(skill.m,yvar){
skill.mplot=subset(skill.m,as.character(set)=='data')
ggplot(skill.mplot,aes_string(y=yvar))+
    geom_boxplot(aes(x=as.factor(yr),colour=modeltype,linetype=residblend),outlier.size=1)+
    labs(x='Year',y=costlong)+
    scale_linetype_manual(name='Residual Blending',values=c(2,1))+
    scale_colour_brewer(palette='Set1',name='Model',labels=c('PHV Regression' ,'PHV+RCN Regression'))+
    theme_minimal()+
    theme(axis.line=element_line(colour='grey10'),
            legend.title=element_text(size=8),
           legend.text=element_text(size=6),
           legend.key=element_rect(colour='white'),
           legend.key.size=unit(0.5,'lines'),
           legend.position = 'bottom', 
           legend.box = "horizontal")
}


## create a facet version so we can steal the legend.
multifacetplot=function(skill.m,yvar,shapeguide){
  ggplot(skill.m,aes_string(y=yvar))+      
     geom_point(aes(x=as.numeric(doy),colour=modeltype,shape=residblend),size=1)+
     facet_grid(set~yr,scale='free_y',space='free')+
     labs(y="")+
     guides(colour=guide_legend(ncol=3, byrow=F,override.aes=list(alpha=1,size=2)))+#,
               # shape=guide_legend(override.aes=list(size=2)))+
     #scale_alpha_manual(name='Residual Blending',values=c(1,.5))+
     scale_shape_manual(guide=as.logical(shapeguide),name='Residual Blending',values=c(4,1))+
     scale_colour_manual(name='Model',values=c('blue','red','grey40'),labels=c('PHV-baseline','PHV-RCN','Skill Difference'))+
     scale_x_continuous('day of year',limits=c(1,160),labels=c(1,seq(30,160,30)),breaks=c(1,seq(30,160,30)))+
     scale_y_continuous(breaks=seq(0,.7,.1),labels=seq(0,.7,.1))+
     theme_bw()+
     theme(axis.line=element_line(colour='grey10'),
           strip.background=element_rect(fill=NA,colour=NA),
           strip.text=element_blank(),#element_text(face='bold',size='12'),
           axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=.5),
           axis.text.y=element_text(size=6,angle=0),
           legend.title=element_text(size=6),
           legend.text=element_text(size=6),
           legend.key=element_rect(colour='white'),
           legend.key.size=unit(.2,'lines'),
           legend.position = 'bottom', 
           legend.box = "horizontal")
}


## extract legend from plot
glegend<-function(a.ggplot){
  tmp <- ggplot_gtable(ggplot_build(a.ggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


## 
dataplot=function(skill.mplot,yvar){
ggplot(skill.mplot,aes_string(y=yvar))+      
      geom_point(aes(x=as.numeric(doy),colour=modeltype,shape=residblend),size=1,alpha=0.6)+
      facet_grid(~yr,scale='free',space='free')+
      labs(y=costlong,x='')+
#      scale_alpha_manual(guide=F,name='',values=c(1,0.5))+
      scale_shape_manual(guide=F,name='',values=c(4,1))+
      scale_colour_manual(guide=F,name='',values=c('blue','red'),labels=c('PHV','PHV+RCN'))+
      scale_x_continuous(limits=c(1,160),labels=c(1,seq(30,180,30)),breaks=c(1,seq(30,180,30)))+
     # coord_cartesian(ylim=c(0,0.7))+
     #guides(colour=guide_legend(ncol=3, byrow=F,override.aes=list(alpha=1,size=2)))+
      theme_bw()+
      theme(axis.line=element_line(colour='grey10'),
           strip.background=element_rect(fill=NA,colour=NA),
           strip.text=element_text(face='bold',size='8'),
           axis.line.x=element_blank(),#element_text(size=6,angle=90,hjust=1,vjust=.5),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank(),
           axis.text.y=element_text(size=6,angle=0),
            axis.title=element_text(size=6),
           legend.key=element_rect(colour='white'),
           legend.position = c(.1,-1), 
           legend.box = "horizontal",
           plot.margin=unit(c(1,1,-.4,1),"cm"))
}

## 2nd plot in combined figure
diffplot=function(skill.mplot,yvar){
ggplot(skill.mplot,aes_string(y=yvar))+      
     geom_point(aes(x=as.numeric(doy),colour=variable),alpha=1,size=1)+
     facet_grid(~yr,scale='free',space='free')+
     labs(y=diffcostlong)+
     scale_colour_manual(guide=F,name='',values='grey40',labels=c('Skill Difference'))+
     scale_x_continuous('day of year',limits=c(1,160),labels=c(1,seq(30,180,30)),breaks=c(1,seq(30,180,30)))+ 
     # guides(colour=guide_legend(ncol=3, byrow=F,override.aes=list(alpha=1,size=2)))+
     theme_bw()+
     theme(axis.line=element_line(colour='grey10'),
           strip.background=element_rect(fill=NA,colour=NA),
           strip.text=element_blank(),#element_text(face='bold',size='12'),
            axis.title=element_text(size=6),
           axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=.5),
           axis.text.y=element_text(size=6,angle=0),
           legend.key=element_rect(),#element_rect(colour='white'),
           legend.position = 'bottom', 
           legend.box = "horizontal",
           plot.margin=unit(c(-.4,1,-.2,1),"cm"))
}

