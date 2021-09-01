# File for processing eBird ERD checklists downloaded on July 19, 2021.
# Created by Erin Conlisk from April-August 2021

# This code downloads all the eBird relative abundance rasters within the East Bay 
# Stewardship Network (EBSN) for each week for a specific species and averages them across 
# the breeding season. Relative abundance rasters were created by Cornell to use all of 
# the North American eBird data to model locations where abundances are expected to be 
# highest or lowest for each week across the year (they do not include trends across years).  

##### These rasters are not necessary for the EBSN Eco-Health analysis of trends. ##### 

# Instead, these rasters were used as a background check for biases in the EBSN eBird data.  
# For example, observations made in 2010 could theoretically occur only in locations where 
# Cornell predicted the highest relative abundance whereas observations made in 2020 could 
# occur only in locations with the lowest relative abundance.  If this were the case then 
# the 2010 data would be occurring in locations with high habitat suitability whereas the 
# 2020 data would be in locations with low habitat suitability.  This would introduce 
# significant bias.  

# Unfortunately, the eBird relative abundance rasters were not significantly resolved enough 
# to expose biases.  While we saw little bias in the eBird relative abundance values across 
# locations and years, we did see bias when we looked at the land cover types underlying the 
# eBird observations.  Because the land cover background check exposed more bias, we focused 
# on land cover as the primary means to mitigate biased sampling in the EBSN data.


.libPaths("C:/R_packages")
library(data.table)
library(ggplot2)
library(sp)
library(raster)
library(rgdal)
library(GISTools)
library(ggmap)
library(plyr)
library(sf)
library(lme4)
library(glmmTMB)
library(spatstat)
library(lattice)
library(maptools)
library(DHARMa)
library(pscl)

# The package "auk" is specifically for downloading and using eBird data.
# I do not use all of the functionality of the package. Thus, the description
# of the package might be very useful:
# https://cornelllabofornithology.github.io/auk/
library(auk)

# The next package downloads spatial layers of eBird relative population trends
# More info can be found here: https://cornelllabofornithology.github.io/ebirdst/
library(ebirdst)

# As the eBird dataset evolves, both of these packages could be subject to changes,
# which may require changes to the code.

# This code depends on two raster files:

# 1) focal_rast - for all species with max territory size less than 4 ha, we use
#                 a 100 m resolution
# 2) foc_rast200m - for all species with max territory size more than 4 ha, we
#                   use a 200 m resolution

# In the code, I delineate sections with "########" and section headings.

##########################################################
# Multipanel figure function:
# The following is not critical for downloading and using 
# the eBird data. It is just for drawing figures.
##########################################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##########################################################

################
# DEFINE SPECIES
# This is probably not the most elegant way to do this, but it works.
# I need to assign the species and all of the "synonyms" for the species
# to accurately assign presence and absence.  So I list all the species 
# and their potential synonyms.  

# Only the species you are looking at should be activated (as opposed to greened out with "#")
#sppx<-"acowoo" #"Melanerpes formicivorus" # 200m resolution
#sppx<-"astfly" #"Myiarchus cinerascens" # 200m resolution
#sppx<-"bkcspa" #"Spizella atrogularis"
#sppx<-"belkin1" #"Megaceryle alcyon"
#sppx<-"buggna" #"Polioptila caerulea"
#sppx<-"bkhgro" #"Pheucticus melanocephalus"
#sppx<-"cowscj1" #"Aphelocoma californica" # 200m resolution
#sppx<-"calthr" #"Toxostoma redivivum" # 200m resolution
#sppx<-"dowwoo" #"Dryobates pubescens" # 200m resolution
#sppx<-"graspa" #"Ammodramus savannarum"
#sppx<-"horlar" #"Eremophila alpestris" # 200m resolution
#sppx<-"larspa" #"Chondestes grammacus"
#sppx<-"logshr" #"Lanius ludovicianus" # 200m resolution
#sppx<-"norhar2" #"Circus hudsonius"
#sppx<-"nutwoo" #"Dryobates nuttallii"
#sppx<-"oaktit" #"Baeolophus inornatus" # 200m resolution
sppx<-"rucspa" #"Aimophila ruficeps"
#sppx<-"savspa" #"Passerculus sandwichensis"
#sppx<-"sospti" #"Melospiza melodia"
#sppx<-"spotow" #"Pipilo maculatus"
#sppx<-"treswa" #"Tachycineta bicolor"
#sppx<-"waviti" #"Vireo gilvus"
#sppx<-"whbnut" #"Sitta carolinensis" # 200m resolution
#sppx<-"wesblu" #"Sialia mexicana"
#sppx<-"wesmea" #"Sturnella neglecta"  # 200m resolution
#sppx<-"wlswar" #"Cardellina pusilla"
#sppx<-"wesmea" #"Chamaea fasciata"
#sppx<-"whtkit" #"Elanus leucurus"
#sppx<-"yelwar" #"Setophaga petechia"


#####################################
# Load eBird abundance rasters
#####################################

# Upload the raster file for the east bay extent:
if (sppx=="acowoo") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="astfly") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="cowscj1") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="calthr") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="dowwoo") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="horlar") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="logshr") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="oaktit") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="whbnut") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else if (sppx=="wesmea") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
} else {
  eb_extent<-raster("G:/EcoHealthEastBay/focal_rast") # Chose a 0.001 decimal degree resolution (approximates Tom's recommended hectare)
  #eb_extent<-raster("G:/EcoHealthEastBay/rast_001dd") # Chose a 0.001 decimal degree resolution (approximates Tom's recommended hectare)
}

# Set the file path
ebirdpath<-"G:/EcoHealthEastBay/"

# First time, create a folder to house files:
dir.create(paste(ebirdpath, "eBird_MAPS/", sppx, sep=""))

# Start downloading the abundance rasters
sp_path <- ebirdst_download(species = sppx)
# see other species with xx<-ebirdst_runs  # https://www.rdocumentation.org/packages/ebirdst/versions/0.2.2/topics/ebirdst_runs
# load the abundance data
# this automatically labels layers with their dates
abd <- load_raster("abundance", path = sp_path)

# Set the spatial extent of the EBSN
ebay_extent<-extent(-10793602, -10659833, 4109000, 4237029)

# For each week starting in the first week of April and ending in the third week of July, 
# (i.e. the breeding season) save the eBird relative abundance rasters.  Label them with 
# the appropriate date.
april_week1<-abd[[14]]
april_wk1<-crop(april_week1, ebay_extent)
writeRaster(april_wk1, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "04052018.tif", sep=""))

april_week2<-abd[[15]]
april_wk2<-crop(april_week2, ebay_extent)
writeRaster(april_wk2, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "04122018.tif", sep=""))

april_week3<-abd[[16]]
april_wk3<-crop(april_week3, ebay_extent)
writeRaster(april_wk3, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "04192018.tif", sep=""))

april_week4<-abd[[17]]
april_wk4<-crop(april_week4, ebay_extent)
writeRaster(april_wk4, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "04262018.tif", sep=""))

may_week1<-abd[[18]]
may_wk1<-crop(may_week1, ebay_extent)
writeRaster(may_wk1, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "05032018.tif", sep=""))

may_week2<-abd[[19]]
may_wk2<-crop(may_week2, ebay_extent)
writeRaster(may_wk2, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "05102018.tif", sep=""))

may_week3<-abd[[20]]
may_wk3<-crop(may_week3, ebay_extent)
writeRaster(may_wk3, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "05172018.tif", sep=""))

may_week4<-abd[[21]]
may_wk4<-crop(may_week4, ebay_extent)
writeRaster(may_wk4, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "05242018.tif", sep=""))

may_week5<-abd[[22]]
may_wk5<-crop(may_week5, ebay_extent)
writeRaster(may_wk5, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "05312018.tif", sep=""))

jun_week1<-abd[[23]]
jun_wk1<-crop(jun_week1, ebay_extent)
writeRaster(jun_wk1, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "06072018.tif", sep=""))

jun_week2<-abd[[24]]
jun_wk2<-crop(jun_week2, ebay_extent)
writeRaster(jun_wk2, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "06142018.tif", sep=""))

jun_week3<-abd[[25]]
jun_wk3<-crop(jun_week3, ebay_extent)
writeRaster(jun_wk3, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "06212018.tif", sep=""))

jun_week4<-abd[[26]]
jun_wk4<-crop(jun_week4, ebay_extent)
writeRaster(jun_wk4, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "06282018.tif", sep=""))

jul_week1<-abd[[27]]
jul_wk1<-crop(jul_week1, ebay_extent)
writeRaster(jul_wk1, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "07062018.tif", sep=""))

jul_week2<-abd[[28]]
jul_wk2<-crop(jul_week2, ebay_extent)
writeRaster(jul_wk2, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "07132018.tif", sep=""))

# Average the loaded rasters across April to mid-July:
ave_spp<-mean(stack(april_wk1, april_wk2, april_wk3, april_wk4, may_wk1, may_wk2, may_wk3, may_wk4, 
                     jun_wk1, jun_wk2, jun_wk3, jun_wk4, jul_wk1, jul_wk2))

# Plot the output:
plot(ave_spp)

# Write the breeding season average relative abundance raster.
writeRaster(ave_spp, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "_ave0401.07152018.tif", sep=""))

# Take the standard deviation in relative abundance and save that raster also. 
stdev_spp<-calc(stack(april_wk1, april_wk2, april_wk3, april_wk4, may_wk1, may_wk2, may_wk3, may_wk4, 
                       jun_wk1, jun_wk2, jun_wk3, jun_wk4, jul_wk1, jul_wk2), fun=sd)
writeRaster(stdev_spp, paste(ebirdpath, "eBird_MAPS/", sppx, "/", sppx, "_stdev0401.07152018.tif", sep=""))

# To use the files created above, we need to change the resolution (i.e. resample) and 
# projection to match the other files we are working with.  Theoretically this could be 
# done in R, however, I never was able to do that.  Thus, I write instructions below for 
# doing it in ArcGIS. If you cannot make the raster file because ArcGIS is not available, 
# you can still continue. You just won't have the summary statistics for the eBird 
# relative abundance value at each site.  In the file, "eBird_grid_analyze_and_figures.R" 
# you will set ArcGIS <- 0 and then you can proceed.

# Instructions for creating the resampled and reprojected file in ESRI: 
# 1. Load the ave file you just wrote.
# 2. Project raster using this file.
# 3. Change raster extent to match "rast_001dd"
# 4. Snap raster to match "rast_001dd"
# 5. Make resolution match "rast_001dd"
# 6. File path and name should be: ("G:/EcoHealthEastBay/eBird_MAPS/xxxxxx/avexxxxxxrprj")


# What is below attempts to do the task above in R. However, it produces a slight mismatch 
# in resolution with the eb_extent file. Further, there are slight differences in how ESRI 
# and R do "resampling". I have not vetted the code below. 

#resampleresolution<-raster("G:/EcoHealthEastBay/yearly_webl/EBrsmplprjxxxx") 
#avexxxxrsmpl<-resample(ave_xxxx, resampleresolution, method='bilinear') 
# Note that bilinear resampling interpolates across space; to preserve EXACT eBird values use "ngb"
#crs(avexxxxrsmpl)<-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs "
#avexxxxrprj<-projectRaster(avexxxxxxrsmpl, crs=crs(xxxx2010_subset))
#writeRaster(avexxxxrprj, "G:/EcoHealthEastBay/eBird_MAPS/xxxxxx/avexxxxxxrprj_R.tif")
