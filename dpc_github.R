#################
# HPSA and DPCs
# Citation: Goldstein ND, Yerkes P. Are direct primacy care practices located in health professional shortage areas? Manuscript in preparation.
# 10/9/23 -- Neal Goldstein
#################


### FUNCTIONS ###

library(rgdal) #shapefile read
library(sp) #spatial data support
library(maptools) #manipulate maps


### READ DATA ###

#HRSA Health Professional Shortage Areas primary care: https://data.hrsa.gov/data/download (Oct 9 2023)
hpsa = readOGR("HPSA_PLYPC_SHP")

#Direct Primary Care Practices: https://mapper.dpcfrontier.com (Oct 9 2023)
#lat/lon embedded in web page source code
dpc = read.csv("DPC.csv", stringsAsFactors=F, as.is=T)

#state cartographic boundaries: https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
us_carto = readOGR("cb_2018_us_state_20m/", "cb_2018_us_state_20m")


### ANALYSIS ###

#set as spatial points type
dpc_sp = dpc
coordinates(dpc_sp) = ~ longitude + latitude

#assign common CRS
proj4string(dpc_sp) = proj4string(hpsa)

#calculate overlap in DPCs with HPSAs
dpc_in_hpsa = cbind(dpc, over(dpc_sp, hpsa))
dpc_in_hpsa$in_hpsa = ifelse(is.na(dpc_in_hpsa$objectid_1), F, T)

dpc_in_hpsa_low = cbind(dpc, over(dpc_sp, hpsa[hpsa$HpsScore<=13, ]))
dpc_in_hpsa_low$in_hpsa = ifelse(is.na(dpc_in_hpsa_low$objectid_1), F, T)

dpc_in_hpsa_med = cbind(dpc, over(dpc_sp, hpsa[hpsa$HpsScore>=14 & hpsa$HpsScore<=17, ]))
dpc_in_hpsa_med$in_hpsa = ifelse(is.na(dpc_in_hpsa_med$objectid_1), F, T)

dpc_in_hpsa_high = cbind(dpc, over(dpc_sp, hpsa[hpsa$HpsScore>=18, ]))
dpc_in_hpsa_high$in_hpsa = ifelse(is.na(dpc_in_hpsa_high$objectid_1), F, T)

dpc_in_hpsa_rural = cbind(dpc, over(dpc_sp, hpsa[hpsa$RurStatDes=="Rural", ]))
dpc_in_hpsa_rural$in_hpsa = ifelse(is.na(dpc_in_hpsa_rural$objectid_1), F, T)

dpc_in_hpsa_partrural = cbind(dpc, over(dpc_sp, hpsa[hpsa$RurStatDes=="Partially Rural", ]))
dpc_in_hpsa_partrural$in_hpsa = ifelse(is.na(dpc_in_hpsa_partrural$objectid_1), F, T)

dpc_in_hpsa_nonrural = cbind(dpc, over(dpc_sp, hpsa[hpsa$RurStatDes=="Non-Rural", ]))
dpc_in_hpsa_nonrural$in_hpsa = ifelse(is.na(dpc_in_hpsa_nonrural$objectid_1), F, T)

#calculate overlaps in HPSAs with DPCs 
hpsa_with_dpc = cbind(hpsa@data, "dpc_id"=over(hpsa, dpc_sp))
hpsa_with_dpc$with_dpc = ifelse(is.na(hpsa_with_dpc$dpc_id), F, T)

#metrics, overall and stratified (note that there are duplicate HPSA geographies, hence % may be greater than 100, e.g. 1126213970 and 1125220326)
sum(dpc_in_hpsa$in_hpsa)/nrow(dpc_in_hpsa)*100
sum(hpsa_with_dpc$with_dpc)/nrow(hpsa_with_dpc)*100

sum(hpsa_with_dpc$HpsScore<=13)
sum(hpsa$HpsScore>=14 & hpsa$HpsScore<=17)
sum(hpsa_with_dpc$HpsScore>=18)

sum(dpc_in_hpsa_low$in_hpsa)/sum(dpc_in_hpsa$in_hpsa)*100
sum(dpc_in_hpsa_med$in_hpsa)/sum(dpc_in_hpsa$in_hpsa)*100
sum(dpc_in_hpsa_high$in_hpsa)/sum(dpc_in_hpsa$in_hpsa)*100

sum(hpsa_with_dpc$RurStatDes=="Rural")
sum(hpsa_with_dpc$RurStatDes=="Partially Rural")
sum(hpsa_with_dpc$RurStatDes=="Non-Rural")

sum(dpc_in_hpsa_rural$in_hpsa)/sum(dpc_in_hpsa$in_hpsa)*100
sum(dpc_in_hpsa_partrural$in_hpsa)/sum(dpc_in_hpsa$in_hpsa)*100
sum(dpc_in_hpsa_nonrural$in_hpsa)/sum(dpc_in_hpsa$in_hpsa)*100

#dpc in hpsa stratified
#sum(dpc_in_hpsa$in_hpsa[dpc_in_hpsa$HpsScore<=13])/sum(dpc_in_hpsa$HpsScore<=13)*100


### CHOROPLETH US MAP ###

#retain only 50 U.S. states
hpsa_carto = hpsa[hpsa$CStNM %in% state.name, ]
us_carto = us_carto[us_carto$NAME %in% state.name, ]

#set projected coordinate system for U.S.
hpsa_carto_proj = spTransform(hpsa_carto,CRS("+init=epsg:2163"))
us_carto_proj = spTransform(us_carto,CRS("+init=epsg:2163"))

#need to transform AK, HI to fit under U.S. for single map; see https://stackoverflow.com/questions/13757771/relocating-alaska-and-hawaii-on-thematic-map-of-the-usa-with-ggplot2
fixup <- function(usa,alaskaFix,hawaiiFix){
  
  alaska=usa[usa$NAME=="Alaska",]
  alaska = fix1(alaska,alaskaFix)
  proj4string(alaska) <- proj4string(usa)
  
  hawaii = usa[usa$NAME=="Hawaii",]
  hawaii = fix1(hawaii,hawaiiFix)
  proj4string(hawaii) <- proj4string(usa)
  
  usa = usa[! usa$NAME %in% c("Alaska","Hawaii"),]
  usa = rbind(usa,alaska,hawaii)
  
  return(usa)
  
}

fix1 <- function(object,params){
  r=params[1];scale=params[2];shift=params[3:4]
  object = elide(object,rotate=r)
  size = max(apply(bbox(object),1,diff))/scale
  object = elide(object,scale=size)
  object = elide(object,shift=shift)
  object
}

hpsa_carto_proj$NAME = hpsa_carto_proj$CStNM
hpsa_map = fixup(hpsa_carto_proj,c(-35,2,-2500000,-2500000),c(-35,1,5500000,-1600000))
us_map = fixup(us_carto_proj,c(-35,2,-2500000,-2500000),c(-35,1,5500000,-1600000))

#flag DPC practices not in contiguous US
dpc$NAME = ifelse(dpc$latitude>=50, "Alaska", ifelse(dpc$longitude<=-130, "Hawaii", "Contiguous"))

#trim Alaska and Hawaii from map, will need to manually enter since transformation does not work
dpc = subset(dpc, dpc$NAME!="Alaska" & dpc$NAME!="Hawaii")

#place points on a projected map: https://gis.stackexchange.com/questions/194474/put-points-on-a-projected-world-map
dpc_sp = dpc
coordinates(dpc_sp) = ~ longitude + latitude
proj4string(dpc_sp) <- CRS(paste("+init=epsg:4326"))
spatial_pt_proj = spTransform(dpc_sp, CRS=CRS("+init=epsg:2163"))

#fix AK, HI points (doesn't work for AK and HI right now)
#dpc_map = fixup(spatial_pt_proj,c(-35,2,-2500000,-2500000),c(-35,1,5500000,-1600000))

rm(fix1,fixup,us_carto,us_carto_proj,hpsa_carto,hpsa_carto_proj,dpc_sp)

#choropleth shading for IMR
map_col = c("#26577C","#B4B4B3","#EBE4D1") #https://www.colorhunt.co/palettes/popular
map_col_index = ifelse(hpsa_map$HpsScore<=13, 3, ifelse(hpsa_map$HpsScore>=18, 1, 2))

#draw choropleth map
par(mar=rep(0.1,4))
plot(hpsa_map, col=map_col[map_col_index], border=NA)
plot(us_map, add=T)
points(spatial_pt_proj, pch=16, col="#E55604", cex=0.5)
legend("bottomright", legend=c("None", "Low", "Medium", "High"), title="HPSA Priority Score", fill=c("#FFFFFF",rev(map_col)), cex=0.8, horiz=T)

#manually add AK, HI
ak_pts = matrix(c(-1350000,-2150000,-1300000,-1950000), ncol=2, byrow=T)
hi_pts = matrix(c(-400000,-2200000,-100000,-2400000,-140000,-2450000,-60000,-2450000), ncol=2, byrow=T)
points(ak_pts, pch=16, col="#E55604", cex=0.5)
points(hi_pts, pch=16, col="#E55604", cex=0.5)

