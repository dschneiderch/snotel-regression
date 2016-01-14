setwd('~/Documents/snotel-regression_project')
library(raster)
library(dplyr)
library(tidyr)


r=raster('data/forestdensity/nlcd_2011_landcover_500m.tif')

dF=as.data.frame(r,xy=T) %>% tbl_df
colnames(dF)=c('x','y','nlcd')

dF=dF %>% mutate(
	landcover=ifelse(nlcd==11,'water',
		ifelse(nlcd==12,'ice',
		ifelse(nlcd>=21 & nlcd<= 24,'developed',
		ifelse(nlcd==31,'barren',
		ifelse(nlcd==41,'deciduous',
		ifelse(nlcd==42,'evergreen',
			ifelse(nlcd==43,'mixed forest',
				ifelse(nlcd==52 | nlcd==71 | nlcd==72 | nlcd==72,'herbaceous',
					ifelse(nlcd==81 |nlcd==82 ,'ag',
						ifelse(nlcd>=90,'wetlands',NA)))))))))))

dF %>%
	group_by(landcover) %>%
	summarise(
		per=n()/nrow(dF)
		) 