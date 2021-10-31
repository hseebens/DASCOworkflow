
graphics.off()
rm(list=ls())

library(sf)
library(shape)
library(data.table)
library(RColorBrewer)

setwd("/home/hanno/Bioinvasion/DASCO/DASCOworkflow/")

source("../../SpatialSpread/owncolleg.r") # adjusted color legend
source(file.path("R","standardise_location_names.R")) # standardise location names (for matching with shapefile)

## set max limits for species numbers
max_terr <- 2500
max_mar  <- 140


### load data ###############################################################

all_regspec_fr <- fread("Data/Output/DASCO_AlienRegions_SInAS_2.4.1.csv")

## Polygon file of marine and terrestrial regions
regions <- st_read(dsn="Data/Input/Shapefiles",layer="RegionsTerrMarine_160621",stringsAsFactors = F)
colnames(regions)[colnames(regions)=="Location"] <- "location"

## standardise location names 
file_name_extension <- "SInAS_2"
newLocNames <- standardise_location_names(regions$location,file_name_extension,data_set="Shapefile")
if (nrow(newLocNames)!=nrow(regions)){
  stop("\n Standardisation of location names went wrong. Check standardise_location_names.R in coords_to_regions.R \n")
} 
regions$location <- newLocNames$location

## change projection to Robinson
crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m" # Robinson projection
regions <- st_wrap_dateline(regions, options = c("WRAPDATELINE=YES"))
regions <- st_transform(regions,crs)


## world map ###################################

nspec_reg <- aggregate(scientificName ~ location,data=all_regspec_fr,FUN=length)
colnames(nspec_reg)[dim(nspec_reg)[2]] <- "nSpec"
# colnames(nspec_reg)[1] <- "Region"

nspec_reg$nSpec[nspec_reg$nSpec>max_terr] <- max_terr

spatial_nspec <- merge(regions,nspec_reg,by=c("location"),all.x=T)


### colour coding ######

ind_terr <- !grepl("_MEOW",spatial_nspec$location) 
spatial_nspec$col_norm[ind_terr] <- round(spatial_nspec$nSpec[ind_terr]/max(spatial_nspec$nSpec[ind_terr],na.rm=T)*99)+1
cols <- colorRampPalette((brewer.pal(n=9,name="YlOrRd")[1:9]))(max(spatial_nspec$col_norm[ind_terr],na.rm=T))

# cols <- rev(colorRampPalette(c("brown4","brown3","orange","yellow","darkolivegreen3","lightblue"))(max(data_regs$col_norm,na.rm=T)))
# cols <- colorRampPalette(c("yellow","yellow","yellow","orange","darkred"))(max(n_sp$Freq_norm)) # 
# cols <- (colorRampPalette(cols)(max(n_sp$Freq_norm)))

spatial_nspec$col <- NA
spatial_nspec[ind_terr,]$col <- cols[spatial_nspec[ind_terr,]$col_norm]
spatial_nspec[ind_terr,]$col[is.na(spatial_nspec[ind_terr,]$col_norm)] <- grey(0.95)


### marine colour coding #######################
ind_marine <- grepl("_MEOW",spatial_nspec$location) 

spatial_nspec$nSpec[ind_marine][spatial_nspec$nSpec[ind_marine]>150] <- 150

spatial_nspec$col_norm[ind_marine] <- round(spatial_nspec$nSpec[ind_marine]/max(spatial_nspec$nSpec[ind_marine],na.rm=T)*99)+1
cols_mar <- colorRampPalette((brewer.pal(n=9,name="YlGnBu")[4:9]))(max(spatial_nspec$col_norm[ind_marine],na.rm=T))

spatial_nspec$col[ind_marine] <- cols_mar[spatial_nspec$col_norm[ind_marine]]



x11(width=9,height=3.3)
# png("../Figures/Worldmap_NumberTaxa_TerrMar.png",unit="in",width=9,height=3.3,res=300)
# png("../Figures/Worldmap_NumberTaxa_GRIISInvasive.png",unit="in",width=8,height=3.3,res=300)
layout(matrix(1:3,nc=3),widths=c(0.7,0.15,0.15))
op <- par(mar=c(0,0,0,0),las=1,cex=0.9,tck=-0.02,mgp=c(2,0.3,0))

plot(st_geometry(spatial_nspec),col=spatial_nspec$col,lwd=0.5) # for some reason, this needs to be plotted first

plot(1:10,axes=F,type="n")
owncolleg(posx=c(0.1,0.2),posy=c(0.1,0.85),col=cols,cex=1.2,digit=0,zlim=c(1,max_terr))#,">2000"
text(x=0.3,y=1.5,labels="Terrestrial",cex=1.3,xpd=T)
mtext(">",side=3,line=-2.8,adj=0.27)

## marine ######

plot(1:10,axes=F,type="n")
owncolleg(posx=c(0.1,0.2),posy=c(0.1,0.85),col=cols_mar,cex=1.2,digit=0,zlim=c(1,max_mar))#,">2000"
text(x=0.3,y=1.5,labels="Marine",cex=1.3,xpd=T)
mtext(">",side=3,line=-2.8,adj=0.27)
# rect(0.3,-2,1,2,col="white",border="white")
# text(x=0.4,y=c(-0.65,0.3,0.95,1.25,1.53),labels=c(1,10,50,100,200))
text(x=0.7,y=0.6,labels="Number of alien taxa",srt=270,cex=1.2)


par(op)
# dev.off()

