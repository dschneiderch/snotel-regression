
setup_plotdata=function(cost,recon.version,covrange,snotelscale,dateflag,fscaMatch,style){
     #skill=read.table(paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/fassnacht/',style,'_recondate_selection_',cost,'.txt'),sep='\t',header=T)
     fnpath=paste0('diagnostics/rswe_',recon.version,'/covrange',covrange,'/snotel',snotelscale,'/',dateflag)
     fnpattern=paste0(style,'_recondate_selection_',fscaMatch,'_',cost)
     fns=list.files(path=fnpath,pattern=fnpattern,full.names=T)
     skill=ldply(fns,read.table,sep='\t',header=T,stringsAsFactors=F)
     skill$date=as.POSIXct(skill$date,tz='MST')
     skill$phvrcn_recondate=as.POSIXct(skill$phvrcn_recondate,tz='MST')
     skill$recon_costdate=as.POSIXct(skill$recon_costdate,tz='MST')
     # str(skill)
     ##
     skill$fulldiff=skill$skill_phvrcnfull-skill$skill_phvfull
     skill$diff=skill$skill_phvrcn-skill$skill_phv
     skill.m=melt(skill[,c('date','yr','skill_phv','skill_phvrcn','skill_phvfull','skill_phvrcnfull','fulldiff','diff','yrdoy')],.(date,yr,yrdoy))
     skill.m$set=ifelse(grepl('skill',as.character(skill.m$variable)),'data','diff')
     skill.m$set=as.factor(skill.m$set)
     skill.m$modeltype=ifelse(grepl('phvrcn',as.character(skill.m$variable)),'phvrcn','phv')
     skill.m$modeltype=ifelse(grepl('diff',as.character(skill.m$variable)),'diff',skill.m$modeltype)
     skill.m$modeltype=as.factor(skill.m$modeltype)
     skill.m$residblend=ifelse(grepl('full',as.character(skill.m$variable)),'blended','unblended')
     skill.m$residblend=factor(skill.m$residblend,levels=c('unblended','blended'))
     skill.m$mnthdy=strftime(skill.m$date,'%b%d')
     skill.m$mnth=factor(strftime(skill.m$date,'%b'),levels=c('Jan','Feb','Mar','Apr','May','Jun'))
     skill.m$doy=as.numeric(strftime(skill.m$date,'%j'))
     
     ## -- calculate daily average snotel for relative calcs below.
     if(snotelscale=='noscale'){
          dailyavg=ddply(snotelrec[snotelrec$date %in% skill$date,],.(date),function(x){
               summarise(x,
                         avgsnotel=mean(swe)
               )})
     } else if(snotelscale=='scale'){
          dailyavg=dailyavg_scaled
     }
     # write.table(x=dailyavg,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_dailysnotelswe.txt'),quote=F,row.names=F)
     
     ## --- set up labels based on skill metric and compute relative rmse if applicable
     if(cost=='rmse'){
          
          costshort='rmse'
          costlong=expression('RMSE')
          diffcostlong=expression(Delta * RMSE)
          
          skill.m=cbind(skill.m,avgsnotel=dailyavg$avgsnotel)
          skill.m$Rrmse=skill.m$value/skill.m$avgsnotel
          
          summrel=ddply(skill.m,.(variable,yr),summarise,
                        avg=mean(Rrmse,na.rm=T),
                        med=median(Rrmse,na.rm=T),
                        min=min(Rrmse,na.rm=T),
                        max=max(Rrmse,na.rm=T))
          # write.table(x=summrel,file=paste0('diagnostics/rswe_',recon.version,'/',style,'_summary_',costshort,'.txt'),quote=F,row.names=F)
          
          # ddply(skill.m,.(variable,yr),function(x){
          #     summarise(x,
          #     min=doy[which.min(Rrmse)],
          #     max=doy[which.max(Rrmse)]
          #     )
          # })
     } else if(cost=='r2'){
          costshort='r2'
          costlong=expression(r^2)
          diffcostlong=expression(Delta * r^2)
          # bquote(r^2)
     }
     return(list(costshort,costlong,diffcostlong,skill.m))
}

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
          theme_bw(base_size=14)+
          theme(axis.line=element_line(colour='grey10'),
                strip.background=element_rect(fill=NA,colour=NA),
                strip.text=element_blank(),#element_text(face='bold',size='12'),
                axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=.5),
                axis.text.y=element_text(size=6,angle=0),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8),
                legend.text.align=0,
                legend.key=element_rect(colour='white'),
                legend.key.width=unit(.7,'lines'),
                legend.position = 'bottom', 
                legend.box = "horizontal")
}


## extract legend from plot
# glegend<-function(a.ggplot){
#   tmp <- ggplot_gtable(ggplot_build(a.ggplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)
# }

##
unify_grobwidth=function(ggobj,maxWidth){
     gggrob=ggplotGrob(ggobj)
     gggrob$widths=maxWidth
     return(gggrob)
}

## 
dataplot=function(skill.mplot,yvar){
     ggplot(skill.mplot,aes_string(y=yvar))+      
          geom_point(aes(x=as.numeric(doy),colour=modeltype,shape=residblend),size=1,alpha=0.6)+
          # facet_wrap(~yr)+
          facet_grid(~yr,scale='free',space='free')+
          labs(y=costlong,x='')+
          #      scale_alpha_manual(guide=F,name='',values=c(1,0.5))+
          scale_shape_manual(guide=F,name='',values=c(4,1))+
          scale_colour_manual(guide=F,name='',values=c('blue','red'),labels=c('PHV-baseline','PHV-RCN'))+
          scale_x_continuous(limits=c(1,160),labels=c(1,seq(30,180,30)),breaks=c(1,seq(30,180,30)))+
          # coord_cartesian(ylim=c(0,0.7))+
          # guides(colour=guide_legend(ncol=3, byrow=F,override.aes=list(alpha=1,size=2)))+
          theme_bw(base_size = 12)+
          theme(axis.line=element_line(colour='grey10'),
                strip.background=element_rect(fill=NA,colour=NA),
                strip.text=element_text(face='bold',size=13),
                axis.line.x=element_blank(),#element_text(size=6,angle=90,hjust=1,vjust=.5),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(size=7,angle=0),
                axis.text.x=element_text(size=7),
                axis.title=element_text(),
                # legend.key=element_rect(colour='white'),
                # legend.position = c(.1,-1), 
                # legend.box = "horizontal",
                plot.margin=unit(c(0,0.5,-0.25,0),"cm")
          )
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
          theme_bw(base_size=12)+
          theme(axis.line=element_line(colour='grey10'),
                strip.background=element_rect(fill=NA,colour=NA),
                strip.text=element_blank(),#element_text(face='bold',size='12'),
                axis.title=element_text(),
                axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=7),
                axis.text.y=element_text(size=7,angle=0),
                # legend.key=element_rect(),#element_rect(colour='white'),
                # legend.position = 'bottom', 
                # legend.box = "horizontal",
                plot.margin=unit(c(-0.01,0.5,0,0),"cm")
          )
}

