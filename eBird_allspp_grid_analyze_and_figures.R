# File for processing eBird ERD checklists downloaded on July 19, 2021.
# Created by Erin Conlisk from April-August 2021

# This code takes the output of "eBird_text_to_shapefile_allapp.R"
# and averages the observations within a grid cell on a specific day.
# Then the resulting daily grids are averaged across the breeding season.
# The result is one grid per year of average counts (across the grid cells 
# and across days) that is used in trend analyses. 

# After filtering the data, we analyzed two trends.  The trends were meant to bracket two possible approaches to analysis: 

# (i)  A conservative response variable (presence-absence) with an inclusive dataset that included all observations from 
#      2010-2020.  By including all the data, the number of observations grew annually over three-fold from 2010 to 2020.  
#      By limiting the analysis to presence-absence data, we removed the possibility that a few observations with high 
#      abundance would have a big influence on the results, where we might expect more high abundance observations with 
#      increased sampling.  While choosing the right distributional form for the analysis can mitigate problems arising 
#      from a few very high bird counts, real data often deviate from statistical distributions.

# (ii) A highly resolved response variable (observed abundance or counts) with a conservative dataset (only locations that 
#      had repeat observations).  By "repeat observations", we mean any location (grid cell) where an observation was made 
#      in at least 8 out of 11 years from 2010-2020, resulting in ~57 observation locations (where the number of observation 
#      locations ranges from 55-58 across species). (Given that each species had slightly different "incidental" observations, 
#      some species had more or fewer than 57 observations.)  We also excluded observations where data were missing from all 
#      years within 2010-2012 because we wanted to have at least one observation within the first three years of the 2010-2020 
#      time series.

# Additional auxilliary analyses are described but not included in the analysis. 

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

# This code depends on the following files:

# 1) "focal_buffer_dissolve" - creates a shapefile for use in point density
#                              calculations. NOT ESSENTIAL.

eb_poly<-readOGR("G:/EcoHealthEastBay", layer="focal_buffer_dissolve")

# 2) "focal_rast" OR "foc_rast200m" - the former is for species with max territory size
#                   less than 4 ha (which uses a grid cell resolution of 100m) and 
#                   the latter is for all species with max territory size greater than 
#                   4 ha (which uses a grid cell resolution of 200 m). 

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
#sppx<-"Melanerpes formicivorus" # 200m resolution
#sppx<-"Myiarchus cinerascens" # 200m resolution
#sppx<-"Spizella atrogularis"
#sppx<-"Megaceryle alcyon"
#sppx<-"Polioptila caerulea"
#sppx<-"Pheucticus melanocephalus"
#sppx<-"Aphelocoma californica" # 200m resolution
#sppx<-"Toxostoma redivivum" # 200m resolution
#sppx<-"Dryobates pubescens" # 200m resolution
#sppx<-"Ammodramus savannarum"
#sppx<-"Eremophila alpestris" # 200m resolution
#sppx<-"Chondestes grammacus"
#sppx<-"Lanius ludovicianus" # 200m resolution
#sppx<-"Circus hudsonius"
#sppx<-"Dryobates nuttallii"
#sppx<-"Baeolophus inornatus" # 200m resolution
sppx<-"Aimophila ruficeps"
#sppx<-"Passerculus sandwichensis"
#sppx<-"Melospiza melodia"
#sppx<-"Pipilo maculatus"
#sppx<-"Tachycineta bicolor"
#sppx<-"Vireo gilvus"
#sppx<-"Sitta carolinensis" # 200m resolution
#sppx<-"Sialia mexicana"
#sppx<-"Sturnella neglecta"  # 200m resolution
#sppx<-"Cardellina pusilla"
#sppx<-"Chamaea fasciata"
#sppx<-"Elanus leucurus"
#sppx<-"Setophaga petechia"

# There is no compatible (i.e. resampled and reprojected) eBird relative abundance 
# file, so set "ArcGIS <- 0".  Set ArcGIS <- 1, if there is a file.
ArcGIS <- 0

if (ArcGIS==1) { # Naming the species with the eBird relative abundance name is only 
  # necessary if there is a compatible relative abundance raster. 

# Only the species you are looking at should be activated (as opposed to greened out with "#")
#sppx2<-"acowoo" #"Melanerpes formicivorus" # 200m resolution
#sppx2<-"astfly" #"Myiarchus cinerascens" # 200m resolution
#sppx2<-"bkcspa" #"Spizella atrogularis"
#sppx2<-"belkin1" #"Megaceryle alcyon"
#sppx2<-"buggna" #"Polioptila caerulea"
#sppx2<-"bkhgro" #"Pheucticus melanocephalus"
#sppx2<-"cowscj1" #"Aphelocoma californica" # 200m resolution
#sppx2<-"calthr" #"Toxostoma redivivum" # 200m resolution
#sppx2<-"dowwoo" #"Dryobates pubescens" # 200m resolution
#sppx2<-"graspa" #"Ammodramus savannarum"
#sppx2<-"horlar" #"Eremophila alpestris" # 200m resolution
#sppx2<-"larspa" #"Chondestes grammacus"
#sppx2<-"logshr" #"Lanius ludovicianus" # 200m resolution
#sppx2<-"norhar2" #"Circus hudsonius"
#sppx2<-"nutwoo" #"Dryobates nuttallii"
#sppx2<-"oaktit" #"Baeolophus inornatus" # 200m resolution
sppx2<-"rucspa" #"Aimophila ruficeps"
#sppx2<-"savspa" #"Passerculus sandwichensis"
#sppx2<-"sospti" #"Melospiza melodia"
#sppx2<-"spotow" #"Pipilo maculatus"
#sppx2<-"treswa" #"Tachycineta bicolor"
#sppx2<-"waviti" #"Vireo gilvus"
#sppx2<-"whbnut" #"Sitta carolinensis" # 200m resolution
#sppx2<-"wesblu" #"Sialia mexicana"
#sppx2<-"wesmea" #"Sturnella neglecta"  # 200m resolution
#sppx2<-"wlswar" #"Cardellina pusilla"
#sppx2<-"wesmea" #"Chamaea fasciata"
#sppx2<-"whtkit" #"Elanus leucurus"
#sppx2<-"yelwar" #"Setophaga petechia"

}

# Load the raster file that will be needed to grid the eBird
# observations.

if (sppx=="Melanerpes formicivorus") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Myiarchus cinerascens") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Aphelocoma californica") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Toxostoma redivivum") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Dryobates pubescens") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m") 
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Eremophila alpestris") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Lanius ludovicianus") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m") 
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Baeolophus inornatus") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Sitta carolinensis") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else if (sppx=="Sturnella neglecta") {
  eb_extent<-raster("G:/EcoHealthEastBay/foc_rast200m")  
  veg2<-raster("G:/EcoHealthEastBay/veg200m")
} else {
  eb_extent<-raster("G:/EcoHealthEastBay/focal_rast") 
  veg2<-raster("G:/EcoHealthEastBay/veg_rast")
}


##############################################################################
# Define whether to restrict the data to certain times of day, distances
# travelled, maximum number of counts, etc.
# I do not employ most of these constraints because I retained effort 
# variables across grids and used these covariates in the statistical 
# models. 
##############################################################################

# From: https://www.sciencedirect.com/science/article/pii/S0006320720307114

# Finally, following eBird "best practices" (Johnston et al., 2019), we removed checklists 
# that took place over >5 h, that traveled distances of >5 km, and that involved >10 observers. 
# Additionally, we removed any checklists lacking data for the distance, duration, and number 
# of observers covariates. Our final dataset contained >210 million records.

use_time_contraints<-1
minhr<-5 # Discard points before 9AM
maxhr<-22 # Discard points after 10PM

use_distance_contraints<-0
maxdist<-10 # Discard points moving > 10km

cap_counts<-0
# For travelling counts, discard points with counts > (distance)/(territory size = 100 km)
max_stat<-10 # For all other counts, discard points with counts > max_stat

use_obstime_contraints<-0
maxtime<-180 # No surveys that lasted greater than 3 hours

# Are we analyzing ALL the data as presence-absence:
presence_only<-0

# Set the path:
ebirdpath<-"G:/EcoHealthEastBay/"

# Get the species file name abbreviation:
spaceis<-regexpr(" ", sppx)[1]
sppname<-paste(substr(sppx, 1,2), substr(sppx, (spaceis+1),(spaceis+2)), sep="") 


#######################
# 2010
#######################

# Load the file:
sp2010_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2010", sep="")) 
# Make sure counts are numeric: 
sp2010_subset@data$counts<-as.numeric(sp2010_subset@data$counts)
hist(sp2010_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2010_subset@data$DATE<-paste(sp2010_subset@data$yr, sp2010_subset@data$month, sp2010_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2010_subset@data$EFFORT_K[which(sp2010_subset@data$PROTOCO=="Incidental")]<-0
sp2010_subset@data$EFFORT_K[which(sp2010_subset@data$PROTOCO=="Stationary")]<-0
sp2010_subset@data$EFFORT_K[which(sp2010_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2010_subset<-sp2010_subset[which(sp2010_subset@data$hour>=minhr),]
  sp2010_subset<-sp2010_subset[which(sp2010_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2010_subset<-sp2010_subset[which(sp2010_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
sp2010_subset@data$max_territories<-sp2010_subset@data$EFFORT_K/0.1
traveli<-intersect(which(sp2010_subset@data$PROTCL_=="P22"), which(sp2010_subset@data$counts>sp2010_subset@data$max_territories))
if (length(traveli)>0){
  sp2010_subset<-sp2010_subset[-traveli,]
}
travelii<-intersect(which(sp2010_subset@data$PROTCL_!="P22"), which(sp2010_subset@data$counts>max_stat))
if (length(travelii)>0){
  sp2010_subset<-sp2010_subset[-travelii,]
}
}
if (use_obstime_contraints==1) {
  sp2010_subset@data<-subset(sp2010_subset@data, sp2010_subset@data$DUR_MIN<=maxtime)
  sp2010_subset<-sp2010_subset[which(sp2010_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2010_subset<-sp2010_subset[which(sp2010_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2010_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2010_subset[which(sp2010_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
  print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season.
# The average of each grid cell will be used to determine the trend
ave2010<-mean(day_stack, na.rm=TRUE)
ave_time_2010<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2010<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2010<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2010) # 65 total points
#length(which(sp2010_subset@data$counts>0))/nrow(sp2010_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
# Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
one2010<-ave2010
one2010[one2010>=0]<-1
ebird_pred<-crop(ebird_pred, extent(one2010))
maskvalues2010<-ebird_pred*one2010
#hist(na.omit(getValues(maskvalues2010)), breaks=20, xlab="eBird Relative Abundance", main="Year 2010")
} else {
  one2010<-ave2010 
  one2010[one2010>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2010))
maskveg<-veg*one2010
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2010<-vv
#ggplot(data=vv2010, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2010")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2010_subset)
points2010 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2010 <- density(points2010) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2010)
points(points2010)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2010<-stack(ave2010, ave_time_2010, ave_dist_2010, ave_obs_2010, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2010<-rasterToPoints(data2010, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2010@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

# The following allows you to save this "processed" shapefile for each species-year
#sppt2010_v2<-sppt2010
#names(sppt2010@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2010_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2010_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2010_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2010_processed.kml", sep=""), layer=paste(sppname, "_2010_processed", sep=""), driver='KML')

#######################
# 2011
#######################

# Load the file:
sp2011_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2011", sep="")) 
# Make sure counts are numeric: 
sp2011_subset@data$counts<-as.numeric(sp2011_subset@data$counts)
hist(sp2011_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2011_subset@data$DATE<-paste(sp2011_subset@data$yr, sp2011_subset@data$month, sp2011_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2011_subset@data$EFFORT_K[which(sp2011_subset@data$PROTOCO=="Incidental")]<-0
sp2011_subset@data$EFFORT_K[which(sp2011_subset@data$PROTOCO=="Stationary")]<-0
sp2011_subset@data$EFFORT_K[which(sp2011_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2011_subset<-sp2011_subset[which(sp2011_subset@data$hour>=minhr),]
  sp2011_subset<-sp2011_subset[which(sp2011_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2011_subset<-sp2011_subset[which(sp2011_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2011_subset@data$max_territories<-sp2011_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2011_subset@data$PROTCL_=="P22"), which(sp2011_subset@data$counts>sp2011_subset@data$max_territories))
  if (length(traveli)>0){
    sp2011_subset<-sp2011_subset[-traveli,]
  }
  travelii<-intersect(which(sp2011_subset@data$PROTCL_!="P22"), which(sp2011_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2011_subset<-sp2011_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2011_subset@data<-subset(sp2011_subset@data, sp2011_subset@data$DUR_MIN<=maxtime)
  sp2011_subset<-sp2011_subset[which(sp2011_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2011_subset<-sp2011_subset[which(sp2011_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2011_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2011_subset[which(sp2011_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2011<-mean(day_stack, na.rm=TRUE)
ave_time_2011<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2011<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2011<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2011) # 65 total points
#length(which(sp2011_subset@data$counts>0))/nrow(sp2011_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2011<-ave2011
  one2011[one2011>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2011))
  maskvalues2011<-ebird_pred*one2011
  #hist(na.omit(getValues(maskvalues2011)), breaks=20, xlab="eBird Relative Abundance", main="Year 2011")
} else {
  one2011<-ave2011 
  one2011[one2011>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2011))
maskveg<-veg*one2011
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2011<-vv
#ggplot(data=vv2011, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2011")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2011_subset)
points2011 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2011 <- density(points2011) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2011)
points(points2011)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2011<-stack(ave2011, ave_time_2011, ave_dist_2011, ave_obs_2011, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2011<-rasterToPoints(data2011, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2011@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

# The following allows you to save this "processed" shapefile for each species-year
#sppt2011_v2<-sppt2011
#names(sppt2011@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2011_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2011_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2011_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2011_processed.kml", sep=""), layer=paste(sppname, "_2011_processed", sep=""), driver='KML')

#######################
# 2012
#######################

# Load the file:
sp2012_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2012", sep="")) 
# Make sure counts are numeric: 
sp2012_subset@data$counts<-as.numeric(sp2012_subset@data$counts)
hist(sp2012_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2012_subset@data$DATE<-paste(sp2012_subset@data$yr, sp2012_subset@data$month, sp2012_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2012_subset@data$EFFORT_K[which(sp2012_subset@data$PROTOCO=="Incidental")]<-0
sp2012_subset@data$EFFORT_K[which(sp2012_subset@data$PROTOCO=="Stationary")]<-0
sp2012_subset@data$EFFORT_K[which(sp2012_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2012_subset<-sp2012_subset[which(sp2012_subset@data$hour>=minhr),]
  sp2012_subset<-sp2012_subset[which(sp2012_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2012_subset<-sp2012_subset[which(sp2012_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2012_subset@data$max_territories<-sp2012_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2012_subset@data$PROTCL_=="P22"), which(sp2012_subset@data$counts>sp2012_subset@data$max_territories))
  if (length(traveli)>0){
    sp2012_subset<-sp2012_subset[-traveli,]
  }
  travelii<-intersect(which(sp2012_subset@data$PROTCL_!="P22"), which(sp2012_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2012_subset<-sp2012_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2012_subset@data<-subset(sp2012_subset@data, sp2012_subset@data$DUR_MIN<=maxtime)
  sp2012_subset<-sp2012_subset[which(sp2012_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2012_subset<-sp2012_subset[which(sp2012_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2012_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2012_subset[which(sp2012_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2012<-mean(day_stack, na.rm=TRUE)
ave_time_2012<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2012<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2012<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2012) # 65 total points
#length(which(sp2012_subset@data$counts>0))/nrow(sp2012_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2012<-ave2012
  one2012[one2012>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2012))
  maskvalues2012<-ebird_pred*one2012
  #hist(na.omit(getValues(maskvalues2012)), breaks=20, xlab="eBird Relative Abundance", main="Year 2012")
} else {
  one2012<-ave2012 
  one2012[one2012>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2012))
maskveg<-veg*one2012
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2012<-vv
#ggplot(data=vv2012, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2012")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2012_subset)
points2012 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2012 <- density(points2012) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2012)
points(points2012)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2012<-stack(ave2012, ave_time_2012, ave_dist_2012, ave_obs_2012, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2012<-rasterToPoints(data2012, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2012@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

# The following allows you to save this "processed" shapefile for each species-year
#sppt2012_v2<-sppt2012
#names(sppt2012@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2012_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2012_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2012_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2012_processed.kml", sep=""), layer=paste(sppname, "_2012_processed", sep=""), driver='KML')

#######################
# 2013
#######################

# Load the file:
sp2013_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2013", sep="")) 
# Make sure counts are numeric: 
sp2013_subset@data$counts<-as.numeric(sp2013_subset@data$counts)
hist(sp2013_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2013_subset@data$DATE<-paste(sp2013_subset@data$yr, sp2013_subset@data$month, sp2013_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2013_subset@data$EFFORT_K[which(sp2013_subset@data$PROTOCO=="Incidental")]<-0
sp2013_subset@data$EFFORT_K[which(sp2013_subset@data$PROTOCO=="Stationary")]<-0
sp2013_subset@data$EFFORT_K[which(sp2013_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2013_subset<-sp2013_subset[which(sp2013_subset@data$hour>=minhr),]
  sp2013_subset<-sp2013_subset[which(sp2013_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2013_subset<-sp2013_subset[which(sp2013_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2013_subset@data$max_territories<-sp2013_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2013_subset@data$PROTCL_=="P22"), which(sp2013_subset@data$counts>sp2013_subset@data$max_territories))
  if (length(traveli)>0){
    sp2013_subset<-sp2013_subset[-traveli,]
  }
  travelii<-intersect(which(sp2013_subset@data$PROTCL_!="P22"), which(sp2013_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2013_subset<-sp2013_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2013_subset@data<-subset(sp2013_subset@data, sp2013_subset@data$DUR_MIN<=maxtime)
  sp2013_subset<-sp2013_subset[which(sp2013_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2013_subset<-sp2013_subset[which(sp2013_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2013_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2013_subset[which(sp2013_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2013<-mean(day_stack, na.rm=TRUE)
ave_time_2013<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2013<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2013<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2013) # 65 total points
#length(which(sp2013_subset@data$counts>0))/nrow(sp2013_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2013<-ave2013
  one2013[one2013>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2013))
  maskvalues2013<-ebird_pred*one2013
  #hist(na.omit(getValues(maskvalues2013)), breaks=20, xlab="eBird Relative Abundance", main="Year 2013")
} else {
  one2013<-ave2013 
  one2013[one2013>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2013))
maskveg<-veg*one2013
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2013<-vv
#ggplot(data=vv2013, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2013")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2013_subset)
points2013 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2013 <- density(points2013) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2013)
points(points2013)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2013<-stack(ave2013, ave_time_2013, ave_dist_2013, ave_obs_2013, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2013<-rasterToPoints(data2013, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2013@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

# The following allows you to save this "processed" shapefile for each species-year
#sppt2013_v2<-sppt2013
#names(sppt2013@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2013_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2013_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2013_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2013_processed.kml", sep=""), layer=paste(sppname, "_2013_processed", sep=""), driver='KML')

#######################
# 2014
#######################

# Load the file:
sp2014_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2014", sep="")) 
# Make sure counts are numeric: 
sp2014_subset@data$counts<-as.numeric(sp2014_subset@data$counts)
hist(sp2014_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2014_subset@data$DATE<-paste(sp2014_subset@data$yr, sp2014_subset@data$month, sp2014_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2014_subset@data$EFFORT_K[which(sp2014_subset@data$PROTOCO=="Incidental")]<-0
sp2014_subset@data$EFFORT_K[which(sp2014_subset@data$PROTOCO=="Stationary")]<-0
sp2014_subset@data$EFFORT_K[which(sp2014_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2014_subset<-sp2014_subset[which(sp2014_subset@data$hour>=minhr),]
  sp2014_subset<-sp2014_subset[which(sp2014_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2014_subset<-sp2014_subset[which(sp2014_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2014_subset@data$max_territories<-sp2014_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2014_subset@data$PROTCL_=="P22"), which(sp2014_subset@data$counts>sp2014_subset@data$max_territories))
  if (length(traveli)>0){
    sp2014_subset<-sp2014_subset[-traveli,]
  }
  travelii<-intersect(which(sp2014_subset@data$PROTCL_!="P22"), which(sp2014_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2014_subset<-sp2014_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2014_subset@data<-subset(sp2014_subset@data, sp2014_subset@data$DUR_MIN<=maxtime)
  sp2014_subset<-sp2014_subset[which(sp2014_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2014_subset<-sp2014_subset[which(sp2014_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2014_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2014_subset[which(sp2014_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2014<-mean(day_stack, na.rm=TRUE)
ave_time_2014<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2014<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2014<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2014) # 65 total points
#length(which(sp2014_subset@data$counts>0))/nrow(sp2014_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2014<-ave2014
  one2014[one2014>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2014))
  maskvalues2014<-ebird_pred*one2014
  #hist(na.omit(getValues(maskvalues2014)), breaks=20, xlab="eBird Relative Abundance", main="Year 2014")
} else {
  one2014<-ave2014 
  one2014[one2014>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2014))
maskveg<-veg*one2014
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2014<-vv
#ggplot(data=vv2014, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2014")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2014_subset)
points2014 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2014 <- density(points2014) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2014)
points(points2014)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2014<-stack(ave2014, ave_time_2014, ave_dist_2014, ave_obs_2014, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2014<-rasterToPoints(data2014, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2014@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

# The following allows you to save this "processed" shapefile for each species-year
#sppt2014_v2<-sppt2014
#names(sppt2014@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2014_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2014_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2014_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2014_processed.kml", sep=""), layer=paste(sppname, "_2014_processed", sep=""), driver='KML')

#######################
# 2015
#######################

# Load the file:
sp2015_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2015", sep="")) 
# Make sure counts are numeric: 
sp2015_subset@data$counts<-as.numeric(sp2015_subset@data$counts)
hist(sp2015_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2015_subset@data$DATE<-paste(sp2015_subset@data$yr, sp2015_subset@data$month, sp2015_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2015_subset@data$EFFORT_K[which(sp2015_subset@data$PROTOCO=="Incidental")]<-0
sp2015_subset@data$EFFORT_K[which(sp2015_subset@data$PROTOCO=="Stationary")]<-0
sp2015_subset@data$EFFORT_K[which(sp2015_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2015_subset<-sp2015_subset[which(sp2015_subset@data$hour>=minhr),]
  sp2015_subset<-sp2015_subset[which(sp2015_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2015_subset<-sp2015_subset[which(sp2015_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2015_subset@data$max_territories<-sp2015_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2015_subset@data$PROTCL_=="P22"), which(sp2015_subset@data$counts>sp2015_subset@data$max_territories))
  if (length(traveli)>0){
    sp2015_subset<-sp2015_subset[-traveli,]
  }
  travelii<-intersect(which(sp2015_subset@data$PROTCL_!="P22"), which(sp2015_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2015_subset<-sp2015_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2015_subset@data<-subset(sp2015_subset@data, sp2015_subset@data$DUR_MIN<=maxtime)
  sp2015_subset<-sp2015_subset[which(sp2015_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2015_subset<-sp2015_subset[which(sp2015_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2015_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2015_subset[which(sp2015_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2015<-mean(day_stack, na.rm=TRUE)
ave_time_2015<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2015<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2015<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2015) # 65 total points
#length(which(sp2015_subset@data$counts>0))/nrow(sp2015_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2015<-ave2015
  one2015[one2015>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2015))
  maskvalues2015<-ebird_pred*one2015
  #hist(na.omit(getValues(maskvalues2015)), breaks=20, xlab="eBird Relative Abundance", main="Year 2015")
} else {
  one2015<-ave2015 
  one2015[one2015>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2015))
maskveg<-veg*one2015
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2015<-vv
#ggplot(data=vv2015, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2015")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2015_subset)
points2015 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2015 <- density(points2015) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2015)
points(points2015)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2015<-stack(ave2015, ave_time_2015, ave_dist_2015, ave_obs_2015, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2015<-rasterToPoints(data2015, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2015@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
#sppt2015_v2<-sppt2015
#names(sppt2015@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
#writeOGR(sppt2015_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2015_processed", sep=""), driver='ESRI Shapefile')

# The following allows you to save this "processed" shapefile for each species-year
#sppt2015_v2<-sppt2015
#names(sppt2015@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2015_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2015_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2015_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2015_processed.kml", sep=""), layer=paste(sppname, "_2015_processed", sep=""), driver='KML')

#######################
# 2016
#######################

# Load the file:
sp2016_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2016", sep="")) 
# Make sure counts are numeric: 
sp2016_subset@data$counts<-as.numeric(sp2016_subset@data$counts)
hist(sp2016_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2016_subset@data$DATE<-paste(sp2016_subset@data$yr, sp2016_subset@data$month, sp2016_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2016_subset@data$EFFORT_K[which(sp2016_subset@data$PROTOCO=="Incidental")]<-0
sp2016_subset@data$EFFORT_K[which(sp2016_subset@data$PROTOCO=="Stationary")]<-0
sp2016_subset@data$EFFORT_K[which(sp2016_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2016_subset<-sp2016_subset[which(sp2016_subset@data$hour>=minhr),]
  sp2016_subset<-sp2016_subset[which(sp2016_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2016_subset<-sp2016_subset[which(sp2016_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2016_subset@data$max_territories<-sp2016_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2016_subset@data$PROTCL_=="P22"), which(sp2016_subset@data$counts>sp2016_subset@data$max_territories))
  if (length(traveli)>0){
    sp2016_subset<-sp2016_subset[-traveli,]
  }
  travelii<-intersect(which(sp2016_subset@data$PROTCL_!="P22"), which(sp2016_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2016_subset<-sp2016_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2016_subset@data<-subset(sp2016_subset@data, sp2016_subset@data$DUR_MIN<=maxtime)
  sp2016_subset<-sp2016_subset[which(sp2016_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2016_subset<-sp2016_subset[which(sp2016_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2016_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2016_subset[which(sp2016_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2016<-mean(day_stack, na.rm=TRUE)
ave_time_2016<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2016<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2016<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2016) # 65 total points
#length(which(sp2016_subset@data$counts>0))/nrow(sp2016_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2016<-ave2016
  one2016[one2016>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2016))
  maskvalues2016<-ebird_pred*one2016
  #hist(na.omit(getValues(maskvalues2016)), breaks=20, xlab="eBird Relative Abundance", main="Year 2016")
} else {
  one2016<-ave2016 
  one2016[one2016>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2016))
maskveg<-veg*one2016
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2016<-vv
#ggplot(data=vv2016, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2016")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2016_subset)
points2016 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2016 <- density(points2016) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2016)
points(points2016)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2016<-stack(ave2016, ave_time_2016, ave_dist_2016, ave_obs_2016, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2016<-rasterToPoints(data2016, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2016@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
#sppt2016_v2<-sppt2016
#names(sppt2016@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
#writeOGR(sppt2016_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2016_processed", sep=""), driver='ESRI Shapefile')

# The following allows you to save this "processed" shapefile for each species-year
#sppt2016_v2<-sppt2016
#names(sppt2016@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2016_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2016_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2016_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2016_processed.kml", sep=""), layer=paste(sppname, "_2016_processed", sep=""), driver='KML')

#######################
# 2017
#######################

# Load the file:
sp2017_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2017", sep="")) 
# Make sure counts are numeric: 
sp2017_subset@data$counts<-as.numeric(sp2017_subset@data$counts)
hist(sp2017_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2017_subset@data$DATE<-paste(sp2017_subset@data$yr, sp2017_subset@data$month, sp2017_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2017_subset@data$EFFORT_K[which(sp2017_subset@data$PROTOCO=="Incidental")]<-0
sp2017_subset@data$EFFORT_K[which(sp2017_subset@data$PROTOCO=="Stationary")]<-0
sp2017_subset@data$EFFORT_K[which(sp2017_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2017_subset<-sp2017_subset[which(sp2017_subset@data$hour>=minhr),]
  sp2017_subset<-sp2017_subset[which(sp2017_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2017_subset<-sp2017_subset[which(sp2017_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2017_subset@data$max_territories<-sp2017_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2017_subset@data$PROTCL_=="P22"), which(sp2017_subset@data$counts>sp2017_subset@data$max_territories))
  if (length(traveli)>0){
    sp2017_subset<-sp2017_subset[-traveli,]
  }
  travelii<-intersect(which(sp2017_subset@data$PROTCL_!="P22"), which(sp2017_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2017_subset<-sp2017_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2017_subset@data<-subset(sp2017_subset@data, sp2017_subset@data$DUR_MIN<=maxtime)
  sp2017_subset<-sp2017_subset[which(sp2017_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2017_subset<-sp2017_subset[which(sp2017_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2017_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2017_subset[which(sp2017_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2017<-mean(day_stack, na.rm=TRUE)
ave_time_2017<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2017<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2017<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2017) # 65 total points
#length(which(sp2017_subset@data$counts>0))/nrow(sp2017_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2017<-ave2017
  one2017[one2017>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2017))
  maskvalues2017<-ebird_pred*one2017
  #hist(na.omit(getValues(maskvalues2017)), breaks=20, xlab="eBird Relative Abundance", main="Year 2017")
} else {
  one2017<-ave2017 
  one2017[one2017>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2017))
maskveg<-veg*one2017
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2017<-vv
#ggplot(data=vv2017, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2017")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2017_subset)
points2017 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2017 <- density(points2017) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2017)
points(points2017)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2017<-stack(ave2017, ave_time_2017, ave_dist_2017, ave_obs_2017, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2017<-rasterToPoints(data2017, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2017@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
#sppt2017_v2<-sppt2017
#names(sppt2017@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
#writeOGR(sppt2017_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2017_processed", sep=""), driver='ESRI Shapefile')

# The following allows you to save this "processed" shapefile for each species-year
#sppt2017_v2<-sppt2017
#names(sppt2017@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2017_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2017_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2017_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2017_processed.kml", sep=""), layer=paste(sppname, "_2017_processed", sep=""), driver='KML')

#######################
# 2018
#######################

# Load the file:
sp2018_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2018", sep="")) 
# Make sure counts are numeric: 
sp2018_subset@data$counts<-as.numeric(sp2018_subset@data$counts)
hist(sp2018_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2018_subset@data$DATE<-paste(sp2018_subset@data$yr, sp2018_subset@data$month, sp2018_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2018_subset@data$EFFORT_K[which(sp2018_subset@data$PROTOCO=="Incidental")]<-0
sp2018_subset@data$EFFORT_K[which(sp2018_subset@data$PROTOCO=="Stationary")]<-0
sp2018_subset@data$EFFORT_K[which(sp2018_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2018_subset<-sp2018_subset[which(sp2018_subset@data$hour>=minhr),]
  sp2018_subset<-sp2018_subset[which(sp2018_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2018_subset<-sp2018_subset[which(sp2018_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2018_subset@data$max_territories<-sp2018_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2018_subset@data$PROTCL_=="P22"), which(sp2018_subset@data$counts>sp2018_subset@data$max_territories))
  if (length(traveli)>0){
    sp2018_subset<-sp2018_subset[-traveli,]
  }
  travelii<-intersect(which(sp2018_subset@data$PROTCL_!="P22"), which(sp2018_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2018_subset<-sp2018_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2018_subset@data<-subset(sp2018_subset@data, sp2018_subset@data$DUR_MIN<=maxtime)
  sp2018_subset<-sp2018_subset[which(sp2018_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2018_subset<-sp2018_subset[which(sp2018_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2018_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2018_subset[which(sp2018_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2018<-mean(day_stack, na.rm=TRUE)
ave_time_2018<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2018<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2018<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2018) # 65 total points
#length(which(sp2018_subset@data$counts>0))/nrow(sp2018_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2018<-ave2018
  one2018[one2018>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2018))
  maskvalues2018<-ebird_pred*one2018
  #hist(na.omit(getValues(maskvalues2018)), breaks=20, xlab="eBird Relative Abundance", main="Year 2018")
} else {
  one2018<-ave2018 
  one2018[one2018>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2018))
maskveg<-veg*one2018
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2018<-vv
#ggplot(data=vv2018, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2018")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2018_subset)
points2018 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2018 <- density(points2018) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2018)
points(points2018)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2018<-stack(ave2018, ave_time_2018, ave_dist_2018, ave_obs_2018, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2018<-rasterToPoints(data2018, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2018@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

# The following allows you to save this "processed" shapefile for each species-year
#sppt2018_v2<-sppt2018
#names(sppt2018@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2018_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2018_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2018_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2018_processed.kml", sep=""), layer=paste(sppname, "_2018_processed", sep=""), driver='KML')

#######################
# 2019
#######################

# Load the file:
sp2019_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2019", sep="")) 
# Make sure counts are numeric: 
sp2019_subset@data$counts<-as.numeric(sp2019_subset@data$counts)
hist(sp2019_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2019_subset@data$DATE<-paste(sp2019_subset@data$yr, sp2019_subset@data$month, sp2019_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2019_subset@data$EFFORT_K[which(sp2019_subset@data$PROTOCO=="Incidental")]<-0
sp2019_subset@data$EFFORT_K[which(sp2019_subset@data$PROTOCO=="Stationary")]<-0
sp2019_subset@data$EFFORT_K[which(sp2019_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2019_subset<-sp2019_subset[which(sp2019_subset@data$hour>=minhr),]
  sp2019_subset<-sp2019_subset[which(sp2019_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2019_subset<-sp2019_subset[which(sp2019_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2019_subset@data$max_territories<-sp2019_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2019_subset@data$PROTCL_=="P22"), which(sp2019_subset@data$counts>sp2019_subset@data$max_territories))
  if (length(traveli)>0){
    sp2019_subset<-sp2019_subset[-traveli,]
  }
  travelii<-intersect(which(sp2019_subset@data$PROTCL_!="P22"), which(sp2019_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2019_subset<-sp2019_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2019_subset@data<-subset(sp2019_subset@data, sp2019_subset@data$DUR_MIN<=maxtime)
  sp2019_subset<-sp2019_subset[which(sp2019_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2019_subset<-sp2019_subset[which(sp2019_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2019_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2019_subset[which(sp2019_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2019<-mean(day_stack, na.rm=TRUE)
ave_time_2019<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2019<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2019<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2019) # 65 total points
#length(which(sp2019_subset@data$counts>0))/nrow(sp2019_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2019<-ave2019
  one2019[one2019>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2019))
  maskvalues2019<-ebird_pred*one2019
  #hist(na.omit(getValues(maskvalues2019)), breaks=20, xlab="eBird Relative Abundance", main="Year 2019")
} else {
  one2019<-ave2019 
  one2019[one2019>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2019))
maskveg<-veg*one2019
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2019<-vv
#ggplot(data=vv2019, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2019")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2019_subset)
points2019 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2019 <- density(points2019) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2019)
points(points2019)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2019<-stack(ave2019, ave_time_2019, ave_dist_2019, ave_obs_2019, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2019<-rasterToPoints(data2019, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2019@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

# The following allows you to save this "processed" shapefile for each species-year
#sppt2019_v2<-sppt2019
#names(sppt2019@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
# Save as a shapefile: 
#writeOGR(sppt2019_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2019_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file: 
#writeOGR(sppt2019_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2019_processed.kml", sep=""), layer=paste(sppname, "_2019_processed", sep=""), driver='KML')

#######################
# 2020
#######################

# Load the file:
sp2020_subset<-readOGR(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "2020", sep="")) 
# Make sure counts are numeric: 
sp2020_subset@data$counts<-as.numeric(sp2020_subset@data$counts)
hist(sp2020_subset@data$counts)
# Create a "DATE" column for sorting the data by day in the breeding season:
sp2020_subset@data$DATE<-paste(sp2020_subset@data$yr, sp2020_subset@data$month, sp2020_subset@data$day, sep=".")

# Assign an effort distance of 0 to the observations where observers were stationary.
# The assignment of stationary counts is articulated in the "PROTOCOL" column:
sp2020_subset@data$EFFORT_K[which(sp2020_subset@data$PROTOCO=="Incidental")]<-0
sp2020_subset@data$EFFORT_K[which(sp2020_subset@data$PROTOCO=="Stationary")]<-0
sp2020_subset@data$EFFORT_K[which(sp2020_subset@data$PROTOCO=="Area")]<-0

# We do not remove "historical" observations here, but they are typically removed
# in the analysis because they lack informaiton on the effort covariates.

# Impose constraints that were defined above:
if (use_time_contraints==1) {
  sp2020_subset<-sp2020_subset[which(sp2020_subset@data$hour>=minhr),]
  sp2020_subset<-sp2020_subset[which(sp2020_subset@data$hour<maxhr),]
}
if (use_distance_contraints==1) {
  sp2020_subset<-sp2020_subset[which(sp2020_subset@data$EFFORT_K<=maxdist),]
}
if (cap_counts==1) {
  sp2020_subset@data$max_territories<-sp2020_subset@data$EFFORT_K/0.1
  traveli<-intersect(which(sp2020_subset@data$PROTCL_=="P22"), which(sp2020_subset@data$counts>sp2020_subset@data$max_territories))
  if (length(traveli)>0){
    sp2020_subset<-sp2020_subset[-traveli,]
  }
  travelii<-intersect(which(sp2020_subset@data$PROTCL_!="P22"), which(sp2020_subset@data$counts>max_stat))
  if (length(travelii)>0){
    sp2020_subset<-sp2020_subset[-travelii,]
  }
}
if (use_obstime_contraints==1) {
  sp2020_subset@data<-subset(sp2020_subset@data, sp2020_subset@data$DUR_MIN<=maxtime)
  sp2020_subset<-sp2020_subset[which(sp2020_subset@data$DUR_MIN<=maxtime),]
}
if (presence_only==1) {
  sp2020_subset<-sp2020_subset[which(sp2020_subset@data$counts>0),]
}


# For each day within the breeding season of a year, we average all points within the same grid cell.
# Make a list of the existing days:
days_to_ave<-unique(sp2020_subset@data$DATE)
# Keep the following line which defines how daily lists are made...
make_list<-1

# Start running through the days:
for (ii in days_to_ave) {
  # subset the points on a given day
  thesepoints<-sp2020_subset[which(sp2020_subset@data$DATE==ii),]
  
  # analyze the counts on that day
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'counts', fun=mean)", sep="")))
  eval(parse(text=paste("day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count<-rasterize(thesepoints, eb_extent, 'counts', fun='count')", sep="")))
  eval(parse(text=paste("cc<-cellStats(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "_count, stat=max)", sep="")))
  
  # analyze the effort on that day: "DUR_MIN", "EFFORT_K", "NUM_OBS"
  eval(parse(text=paste("day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'DUR_MIN', fun=mean)", sep="")))
  eval(parse(text=paste("day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'EFFORT_K', fun=mean)", sep="")))
  eval(parse(text=paste("day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), "<-rasterize(thesepoints, eb_extent, 'NUM_OBS', fun=mean)", sep="")))
  
  if (cc>1) {
    print(paste("More than two points in a cell in day ", ii, sep=""))
  }
  if (make_list==1) {
    eval(parse(text=paste("day_stack<-stack(day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    make_list<-2
  } else {
    eval(parse(text=paste("day_stack<-stack(day_stack, day_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_time_stack<-stack(day_time_stack, day_time_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_dist_stack<-stack(day_dist_stack, day_dist_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
    eval(parse(text=paste("day_obs_stack<-stack(day_obs_stack, day_obs_", substr(ii, 1, 4), substr(ii, 6,7), substr(ii, 9, 10), ")", sep="")))
  }
}
# Next we average each location across the breeding season - this could be revisited!
# The average of each grid cell will be used to determine the trend
ave2020<-mean(day_stack, na.rm=TRUE)
ave_time_2020<-mean(day_time_stack, na.rm=TRUE)
ave_dist_2020<-mean(day_dist_stack, na.rm=TRUE)
ave_obs_2020<-mean(day_obs_stack, na.rm=TRUE)
#freq(ave2020) # 65 total points
#length(which(sp2020_subset@data$counts>0))/nrow(sp2020_subset@data)

# Code for using the eBird abundance files. Only works if ArcGIS==1.
if (ArcGIS==1) {
  # Range map downloaded from: https://ebird.org/science/status-and-trends/download-data?forceLogin=true
  ebird_pred<-raster(paste(ebirdpath, "eBird_MAPS/", sppx2,"/ave", sppx2, "rprj", sep=""))
  one2020<-ave2020
  one2020[one2020>=0]<-1
  ebird_pred<-crop(ebird_pred, extent(one2020))
  maskvalues2020<-ebird_pred*one2020
  #hist(na.omit(getValues(maskvalues2020)), breaks=20, xlab="eBird Relative Abundance", main="Year 2020")
} else {
  one2020<-ave2020 
  one2020[one2020>=0]<-1
}

# Now upload the East Bay vegetation maps and sample the cover types where bird observations occurred.
veg<-extend(veg2, extent(one2020))
maskveg<-veg*one2020
vv<-as.data.frame(freq(maskveg))
vv$veg_type<-rep(NA, nrow(vv))
vv$count[which(vv$value==10)]<-vv$count[which(vv$value==10)]+ifelse(length(which(vv$value==2))>0, vv$count[which(vv$value==2)], 0)
vv$veg_type[which(vv$value==3)]<-"riparian"
vv$veg_type[which(vv$value==4)]<-"oak"
vv$veg_type[which(vv$value==5)]<-"grass"
vv$veg_type[which(vv$value==6)]<-"develop"
vv$veg_type[which(vv$value==7)]<-"forest"
vv$veg_type[which(vv$value==8)]<-"chaparral"
vv$veg_type[which(vv$value==9)]<-"scrub"
vv$veg_type[which(vv$value==10)]<-"water"
vv$veg_type[which(vv$value==11)]<-"ag"
vv$veg_type[which(vv$value==13)]<-"barren"
vv$veg_type[which(vv$value==14)]<-"wetland"
vv$veg_type[which(vv$value==15)]<-"redwood"
vv$veg_type[which(vv$value==16)]<-"non-native"
vv$veg_type[which(vv$value==17)]<-"conifer"
vv<-vv[-which(is.na(vv$veg_type)==TRUE),]
vv2020<-vv
#ggplot(data=vv2020, aes(x=veg_type, y=count)) +
#  geom_bar(stat="identity") + ggtitle("Year 2020")

# Draw a point density figure (point density of original data; 
# NOT the resulting points averaged acorss grid cells and the 
# breeding season:
extwind<-owin(xrange=c(-122.4,-121.4), yrange=c(37.4,38.1))
pts <- coordinates(sp2020_subset)
points2020 <- ppp(pts[,1], pts[,2], window=extwind)
Point_Density_2020 <- density(points2020) # Using the default bandwidth
shapefile_df <- fortify(eb_poly)
plot(Point_Density_2020)
points(points2020)
lines(shapefile_df)

# Stack all the rasters for:
# (1) average counts, 
# (2) average sampling duration, 
# (3) average sampling distance,
# (4) average number of observers, and
# (5) vegetation type
data2020<-stack(ave2020, ave_time_2020, ave_dist_2020, ave_obs_2020, maskveg)
# Create points for these raster values. Now there is a shapefile with spatial
# information for points and a data table with the important attributes:
sppt2020<-rasterToPoints(data2020, fun=NULL, spatial=TRUE)
# Clean up the attribute table names:
names(sppt2020@data)<-c("counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")

#sppt2020_v2<-sppt2020
#names(sppt2020@data)<-c("counts", "mins", "dist", "num_obs", "vegid")
#writeOGR(sppt2020_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""), layer=paste(sppname, "_2020_processed", sep=""), driver='ESRI Shapefile')
# Here is an option to save as a kml file
#writeOGR(sppt2020_v2, paste(ebirdpath, sppname, "_eBird_data_by_yr/", sppname, "_2020_processed.kml", sep=""), layer=paste(sppname, "_2020_processed", sep=""), driver='KML')

################################


#################################################
# Histograms of veg and eBird relative abundance
#################################################

# Plots vegetation in the larger landscape for comparison to 
# the vegetation cover present at locations of bird observations.
vegland<-as.data.frame(freq(veg))
vegland$veg_type<-rep(NA, nrow(vegland))
vegland$veg_type[which(vegland$value==3)]<-"riparian"
vegland$veg_type[which(vegland$value==4)]<-"oak"
vegland$veg_type[which(vegland$value==5)]<-"grass"
vegland$veg_type[which(vegland$value==6)]<-"develop"
vegland$veg_type[which(vegland$value==7)]<-"forest"
vegland$veg_type[which(vegland$value==8)]<-"chaparral"
vegland$veg_type[which(vegland$value==9)]<-"scrub"
vegland$veg_type[which(vegland$value==10)]<-"water"
vegland$veg_type[which(vegland$value==11)]<-"ag"
vegland$veg_type[which(vegland$value==13)]<-"barren"
vegland$veg_type[which(vegland$value==14)]<-"wetland"
vegland$veg_type[which(vegland$value==15)]<-"redwood"
vegland$veg_type[which(vegland$value==16)]<-"non-native"
vegland$veg_type[which(vegland$value==17)]<-"conifer"

ggplot(data=vegland[c(2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15),], aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Landscape Vegetation")


# If there are eBird relative abundance maps, this code plots the
# predicted relative abundance at locations where birds were 
# observed.
if (ArcGIS==1) {
par(mfrow=c(3, 4))
hist(na.omit(getValues(maskvalues2010)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2010))), 2), sep=""), main="Year 2010")
hist(na.omit(getValues(maskvalues2011)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2011))), 2), sep=""), main="Year 2011")
hist(na.omit(getValues(maskvalues2012)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2012))), 2), sep=""), main="Year 2012")
hist(na.omit(getValues(maskvalues2013)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2013))), 2), sep=""), main="Year 2013")
hist(na.omit(getValues(maskvalues2014)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2014))), 2), sep=""), main="Year 2014")
hist(na.omit(getValues(maskvalues2015)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2015))), 2), sep=""), main="Year 2015")
hist(na.omit(getValues(maskvalues2016)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2016))), 2), sep=""), main="Year 2016")
hist(na.omit(getValues(maskvalues2017)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2017))), 2), sep=""), main="Year 2017")
hist(na.omit(getValues(maskvalues2018)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2018))), 2), sep=""), main="Year 2018")
hist(na.omit(getValues(maskvalues2019)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2019))), 2), sep=""), main="Year 2019")
hist(na.omit(getValues(maskvalues2020)), breaks=20, xlim=c(0, 5), xlab=paste("eBird Rel Abund - mean = ", round(mean(na.omit(getValues(maskvalues2020))), 2), sep=""), main="Year 2020")
}


# Plots vegetation in each year at locations of bird observations.
par(mfrow=c(3, 4))
y10<-ggplot(data=vv2010, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2010")
y11<-ggplot(data=vv2011, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2011")
y12<-ggplot(data=vv2012, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2012")
y13<-ggplot(data=vv2013, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2013")
y14<-ggplot(data=vv2014, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2014")
y15<-ggplot(data=vv2015, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2015")
y16<-ggplot(data=vv2016, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2016")
y17<-ggplot(data=vv2017, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2017")
y18<-ggplot(data=vv2018, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2018")
y19<-ggplot(data=vv2019, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2019")
y20<-ggplot(data=vv2020, aes(x=veg_type, y=count)) +
  geom_bar(stat="identity") + ggtitle("Year 2020")
multiplot(y10, y14, y18, y11, y15, y19, y12, y16, y20, y13, y17, cols=4)

#############################################################
# Create table with all of the annual data
#############################################################

sppt2010@data$year<-rep(2010, nrow(sppt2010@data))
sppt2011@data$year<-rep(2011, nrow(sppt2011@data))
sppt2012@data$year<-rep(2012, nrow(sppt2012@data))
sppt2013@data$year<-rep(2013, nrow(sppt2013@data))
sppt2014@data$year<-rep(2014, nrow(sppt2014@data))
sppt2015@data$year<-rep(2015, nrow(sppt2015@data))
sppt2016@data$year<-rep(2016, nrow(sppt2016@data))
sppt2017@data$year<-rep(2017, nrow(sppt2017@data))
sppt2018@data$year<-rep(2018, nrow(sppt2018@data))
sppt2019@data$year<-rep(2019, nrow(sppt2019@data))
sppt2020@data$year<-rep(2020, nrow(sppt2020@data))
all_counts<-as.data.frame(sppt2010@data)
all_counts<-rbind(all_counts, sppt2011@data)
all_counts<-rbind(all_counts, sppt2012@data)
all_counts<-rbind(all_counts, sppt2013@data)
all_counts<-rbind(all_counts, sppt2014@data)
all_counts<-rbind(all_counts, sppt2015@data)
all_counts<-rbind(all_counts, sppt2016@data)
all_counts<-rbind(all_counts, sppt2017@data)
all_counts<-rbind(all_counts, sppt2018@data)
all_counts<-rbind(all_counts, sppt2019@data)
all_counts<-rbind(all_counts, sppt2020@data)
all_counts$yr<-all_counts$year-2010
# Create PA data to use in analysis below:
all_counts$pres<-ifelse(all_counts$count>0, 1, 0)

# Change the numeric vegetation codes to text descriptions:
all_counts$veg_type<-rep(NA, nrow(all_counts))
all_counts$veg_type[which(all_counts$veg_code==2)]<-"water"
all_counts$veg_type[which(all_counts$veg_code==3)]<-"riparian"
all_counts$veg_type[which(all_counts$veg_code==4)]<-"oak"
all_counts$veg_type[which(all_counts$veg_code==5)]<-"grass"
all_counts$veg_type[which(all_counts$veg_code==6)]<-"develop"
all_counts$veg_type[which(all_counts$veg_code==7)]<-"forest"
all_counts$veg_type[which(all_counts$veg_code==8)]<-"chaparral"
all_counts$veg_type[which(all_counts$veg_code==9)]<-"scrub"
all_counts$veg_type[which(all_counts$veg_code==10)]<-"water"
all_counts$veg_type[which(all_counts$veg_code==11)]<-"ag"
all_counts$veg_type[which(all_counts$veg_code==13)]<-"barren"
all_counts$veg_type[which(all_counts$veg_code==14)]<-"wetland"
all_counts$veg_type[which(all_counts$veg_code==15)]<-"redwood"
all_counts$veg_type[which(all_counts$veg_code==16)]<-"non-native"
all_counts$veg_type[which(all_counts$veg_code==17)]<-"conifer"

# Plot all the observations 
# THE RESULTING PLOT OF ALL THE ABUNDANCE DATA IS NOT USED IN ANALYSIS
# This plot is just for review.
if (presence_only==1) {
ggplot(all_counts, aes(x=year, y=counts)) + 
  geom_point(alpha=0.1, size=5)+
  scale_x_continuous(breaks=seq(2010, 2020, 2),limits=c(2010, 2020)) +
  theme_classic() + xlab("Year") + ylab("Bird Counts") + ggtitle(paste(sppx, " Trend - Presence Only", sep=""))
} else  {
  ggplot(all_counts, aes(x=year, y=counts)) + 
    geom_point(alpha=0.1, size=5)+
    scale_x_continuous(breaks=seq(2010, 2020, 2),limits=c(2010, 2020)) +
    theme_classic() + xlab("Year") + ylab("Bird Counts") + ggtitle(paste(sppx, " Trend - Presences and Absences", sep=""))
}

# Save the full dataset:
all_counts_with_dev<-all_counts
# Take out any points in developed areas
all_counts<-subset(all_counts, all_counts$veg_type!="develop")

# List of number of presences by year:
#tapply(all_counts$pres, all_counts$year, sum)
# List of average occupnacy by year:
#tapply(all_counts$pres, all_counts$year, sum)/tapply(all_counts$pres, all_counts$year, length)

##############################################################
# Trend analyses

# We consider two approaches that bracket the possiblities:
# (1) In the first approach, we consider ALL the data, but
#     reduce the potential biases in this data by only
#     considering presence-absence.
# (2) In the second approach we consider a restricted 
#     dataset, but retain all the information provided by
#     the counts (or abundances).

# We also performed some auxilliary tests to check for 
# consistent answers across approaches.  The auxilliary 
# tests were not contained in the final report and are 
# flagged as auxilliary approaches below.
##############################################################

#######################################
# APPROACH 1
# Full Dataset, Presence-Absence
#######################################

b4 <- glm(pres ~ yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs, family = "binomial", data = all_counts)
summary(b4)

yearlyPA<-as.data.frame(tapply(all_counts$pres, all_counts$year, mean))
names(yearlyPA)<-"PA"
yearlyPA$yr<-row.names(yearlyPA)
yearlyPA$Number_Observations<-tapply(all_counts$pres, all_counts$year, length)
yearlyPA$SE<-sqrt((yearlyPA$PA*(1-yearlyPA$PA))/yearlyPA$Number_Observations)

ggplot(yearlyPA, aes(x=yr, y=PA, group=1)) + 
  geom_line(colour = "black") + 
  geom_point(aes(size=Number_Observations), show.legend=FALSE)+
  geom_errorbar(aes(ymin=PA-SE, ymax=PA+SE), width=.2,
                position=position_dodge(0.05)) + 
  #scale_x_continuous(breaks=seq(2010, 2020, 2),limits=c(2010, 2020)) +
  #coord_cartesian(ylim=c(0,0.25)) + 
  scale_x_discrete(breaks=seq(2010, 2020, 2)) + 
  theme_classic() + xlab("Year") + ylab(paste("Fraction Sites ", sppx, " Present", sep=""))


#######################################################
# Full Dataset, Counts

# THIS ANALYSIS WAS NOT USED IN THE FINAL REPORT,
# ONLY USED TO 'CHECK' FOR DISCREPANCIES WITH
# PRESENCE-ABSENCE RESULTS...
# Thus, I have removed the outputs of this code 
# with "#". If you want to see summaries, remove 
# the "#" before the command summary().
#######################################################

# Round all data to nearest integer for count analysis
all_counts$roundcount<-round(all_counts$count)
all_counts$notcount<-all_counts$roundcount-all_counts$count #261/1759 not counts...14.8%

# Poisson model:
m3<-glm(roundcount~yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs, family = poisson, data=all_counts)
#summary(m3)

# Negative binomial model:
n3<-glm.nb(roundcount~yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs, data=all_counts)
#summary(n3)
#pchisq(2 * (logLik(n3) - logLik(m3)), df = 1, lower.tail = FALSE)

# Zero-inflated, Poisson model (typically less likely to show signficant results):
m4 <- zeroinfl(roundcount~yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs | 1, data = all_counts)
#summary(m4)
#mnull <- update(m4, . ~ 1)
#pchisq(2 * (logLik(m4) - logLik(mnull)), df = 1, lower.tail = FALSE) # Tests whether the model is a signicant improvement over null model.

# Zero-inflated, Negative binomialmodel (typically less likely to show signficant results):
n4 <- zeroinfl(roundcount~yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs | 1, dist = "negbin", data = all_counts)
#summary(n4)
#nnull <- update(n4, . ~ 1)
#pchisq(2 * (logLik(n4) - logLik(nnull)), df = 1, lower.tail = FALSE) # Tests whether the model is a signicant improvement over nothing.

#vuong(m3, m4) # Vuong test to show whether the zero-inflated poisson model is better.
#vuong(n3, n4) # Vuong test to show whether the zero-inflated negative binomial model is better.

# Code to take out individual years, or groups of years
# The example below takes out year 2010
after2010<-subset(all_counts, all_counts$year>2010)
an4 <- zeroinfl(roundcount~yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs | 1, dist = "negbin", data = after2010)
#summary(an4)


#######################################################
# Full Dataset, Subsetted to 40 points in each year
# Both Presence-Absence and Counts were analyzed.

# THIS ANALYSIS WAS NOT USED IN THE FINAL REPORT,
# ONLY USED TO 'CHECK' FOR DISCREPANCIES WITH
# PRESENCE-ABSENCE RESULTS...
#######################################################

# Subsample down to 40 points in each year.  Do this a bunch
sub_pts<-40 # Number of points in each year
num_runs<-1000 # Number of times to subsample
run_num<-1 # Running index for summary tables

# Create summary tables to store results.
summary_subsets<-as.data.table(matrix(NA, num_runs, 3))
names(summary_subsets)<-c("run_num", "coefficient", "pvalue")
summary_bin<-as.data.table(matrix(NA, num_runs, 3))
names(summary_bin)<-c("run_num", "coefficient", "pvalue")

# Code to take out individual years, or groups of years
# The example below takes out year 2010
#all_counts_save<-all_counts # First save the original file so it can be recovered.
#all_counts<-subset(all_counts, all_counts$year>2010) 
# If you use the code above, set: take_out<-1
take_out<-0 # This adjusts the for-loop below.  If a year other than
# 2010 is chosen than the code below will have to be modified to 
# remove the appropriate years...

for (iii in 1:num_runs) {
  if (take_out==0) {
  pt2010<-subset(all_counts, all_counts$year==2010)
  pt2010$rand<-runif(nrow(pt2010))
  pt2010sub<-pt2010[order(pt2010$rand),]
  pt2010sub<-pt2010sub[1:sub_pts,]
  }
  pt2011<-subset(all_counts, all_counts$year==2011)
  pt2011$rand<-runif(nrow(pt2011))
  pt2011sub<-pt2011[order(pt2011$rand),]
  pt2011sub<-pt2011sub[1:sub_pts,]
  pt2012<-subset(all_counts, all_counts$year==2012)
  pt2012$rand<-runif(nrow(pt2012))
  pt2012sub<-pt2012[order(pt2012$rand),]
  pt2012sub<-pt2012sub[1:sub_pts,]
  pt2013<-subset(all_counts, all_counts$year==2013)
  pt2013$rand<-runif(nrow(pt2013))
  pt2013sub<-pt2013[order(pt2013$rand),]
  pt2013sub<-pt2013sub[1:sub_pts,]
  pt2014<-subset(all_counts, all_counts$year==2014)
  pt2014$rand<-runif(nrow(pt2014))
  pt2014sub<-pt2014[order(pt2014$rand),]
  pt2014sub<-pt2014sub[1:sub_pts,]
  pt2015<-subset(all_counts, all_counts$year==2015)
  pt2015$rand<-runif(nrow(pt2015))
  pt2015sub<-pt2015[order(pt2015$rand),]
  pt2015sub<-pt2015sub[1:sub_pts,]
  pt2016<-subset(all_counts, all_counts$year==2016)
  pt2016$rand<-runif(nrow(pt2016))
  pt2016sub<-pt2016[order(pt2016$rand),]
  pt2016sub<-pt2016sub[1:sub_pts,]
  pt2017<-subset(all_counts, all_counts$year==2017)
  pt2017$rand<-runif(nrow(pt2017))
  pt2017sub<-pt2017[order(pt2017$rand),]
  pt2017sub<-pt2017sub[1:sub_pts,]
  pt2018<-subset(all_counts, all_counts$year==2018)
  pt2018$rand<-runif(nrow(pt2018))
  pt2018sub<-pt2018[order(pt2018$rand),]
  pt2018sub<-pt2018sub[1:sub_pts,]
  pt2019<-subset(all_counts, all_counts$year==2019)
  pt2019$rand<-runif(nrow(pt2019))
  pt2019sub<-pt2019[order(pt2019$rand),]
  pt2019sub<-pt2019sub[1:sub_pts,]
  pt2020<-subset(all_counts, all_counts$year==2020)
  pt2020$rand<-runif(nrow(pt2020))
  pt2020sub<-pt2020[order(pt2020$rand),]
  pt2020sub<-pt2020sub[1:sub_pts,]
  if (take_out==1) {
    pt1011<-pt2010sub
  } else {
  pt1011<-rbind(pt2010sub, pt2011sub)
  }
  pt1213<-rbind(pt2012sub, pt2013sub)
  pt1415<-rbind(pt2014sub, pt2015sub)
  pt1617<-rbind(pt2016sub, pt2017sub)
  pt1819<-rbind(pt2018sub, pt2019sub)
  pt1013<-rbind(pt1011, pt1213)
  pt1417<-rbind(pt1415, pt1617)
  pt1820<-rbind(pt1819, pt2020sub)
  pt1017<-rbind(pt1013, pt1417)
  pt1020<-rbind(pt1017, pt1820)
  #modsub<-lm(count~yr, data=pt1020)
  modsub <- zeroinfl(roundcount~yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs | 1, dist = "negbin", data = pt1020)
  sumsub<-summary(modsub)
  summary_subsets$run_num[run_num]<-run_num 
  summary_subsets$coefficient[run_num]<-sumsub$coefficients$count[2,1]
  summary_subsets$pvalue[run_num]<-sumsub$coefficients$count[2,4]
  
  b4 <- glm(pres ~ yr+minutes_effort+distance_effort+I(distance_effort^2)+num_obs, family = "binomial", data = pt1020)
  binsub <- summary(b4)
  summary_bin$run_num[run_num]<-run_num 
  summary_bin$coefficient[run_num]<-binsub$coefficients[2,1]
  summary_bin$pvalue[run_num]<-binsub$coefficients[2,4]
  
  run_num<-run_num+1
}

# Take out runs when the model did not converge
summary_subsets<-subset(summary_subsets, is.na(summary_subsets$pvalue)==FALSE)
summary_bin<-subset(summary_bin, is.na(summary_bin$pvalue)==FALSE)

# Plots for the Count-based Analysis
ggplot(summary_subsets, aes(x=pvalue)) + 
  geom_histogram(binwidth=0.02, color="black", fill="gray") + theme_classic() + 
  xlab("Year P-values") + coord_cartesian(xlim=c(0,1)) +
  ggtitle("P-values of zi-NegBin models with 40 points/year")

ggplot(summary_subsets, aes(x=coefficient)) + 
  geom_histogram(color="black", fill="gray") + theme_classic() + 
  xlab("Coefficients") + 
  ggtitle("Year Coefficients of zi-NegBin models with 40 points/year")

# Plots for the Presence-Absence Analysis
ggplot(summary_bin, aes(x=pvalue)) + 
  geom_histogram(binwidth=0.02, color="black", fill="gray") + theme_classic() + 
  xlab("Year P-values") + coord_cartesian(xlim=c(0,1)) +
  ggtitle("P-values of binomial models with 40 points/year")

ggplot(summary_bin, aes(x=coefficient)) + 
  geom_histogram(color="black", fill="gray") + theme_classic() + 
  xlab("Year Coefficients") + 
  ggtitle("Coefficients of binomial models with 40 points/year")


#########################################################################################################
# APPROACH 2
# Now we reduce the amount of data and only consider locations with repeat measures
# However we retain the precise abundance estimate contained in the count-based data
# (instead of converting to presence-absence daat).

# First we organize the data, then we look at statistical models.
#########################################################################################################

######################################
# Data organization
######################################

# See how many grid cells have repeat measures. Are there enough to do an analysis with a random effect?
# For each year we create a raster with a 1 if there was an observation and a zero if there was no obervation
isBecomes <- cbind(c(NA, 1), c(0, 1)) 
ones2010<-ave2010
ones2010[ones2010>=0]<-1
ones2010<-reclassify(ones2010, rcl=isBecomes)
ones2011<-ave2011
ones2011[ones2011>=0]<-1
ones2011<-reclassify(ones2011, rcl=isBecomes)
ones2012<-ave2012
ones2012[ones2012>=0]<-1
ones2012<-reclassify(ones2012, rcl=isBecomes)
ones2013<-ave2013
ones2013[ones2013>=0]<-1
ones2013<-reclassify(ones2013, rcl=isBecomes)
ones2014<-ave2014
ones2014[ones2014>=0]<-1
ones2014<-reclassify(ones2014, rcl=isBecomes)
ones2015<-ave2015
ones2015[ones2015>=0]<-1
ones2015<-reclassify(ones2015, rcl=isBecomes)
ones2016<-ave2016
ones2016[ones2016>=0]<-1
ones2016<-reclassify(ones2016, rcl=isBecomes)
ones2017<-ave2017
ones2017[ones2017>=0]<-1
ones2017<-reclassify(ones2017, rcl=isBecomes)
ones2018<-ave2018
ones2018[ones2018>=0]<-1
ones2018<-reclassify(ones2018, rcl=isBecomes)
ones2019<-ave2019
ones2019[ones2019>=0]<-1
ones2019<-reclassify(ones2019, rcl=isBecomes)
ones2020<-ave2020
ones2020[ones2020>=0]<-1
ones2020<-reclassify(ones2020, rcl=isBecomes)

# Sum the observations across years
num_repeats2<-ones2010+ones2011+ones2012+ones2013+ones2014+ones2015+ones2016+ones2017+ones2018+ones2019+ones2020
#freq(num_repeats2)

# Take out locations where there are missing data for 2010-2012.
first3_repeats<-ones2010+ones2011+ones2012
first3_repeats[first3_repeats>0]<-1
num_repeats<-first3_repeats*num_repeats2

# See how many grid cells have >=8 points in them:
freq(num_repeats)

# Create point file to extract bird data from stack:
point8<-rasterToPoints(num_repeats, fun=function(x){x>7})
# ggplot(aes(x,y, color=layer), data=as.data.frame(point8))+geom_point() # plot points
point8<-as.data.frame(point8)
coords8<-cbind(point8$x, point8$y)
sp8 <- SpatialPointsDataFrame(coords = coords8, data=point8, proj4string = CRS('+proj=longlat +ellps=GRS80 +no_defs'))
names(sp8@data)<-c("x", "y", "num_repeats")
#plot(eb_extent)
#plot(sp8)


# Now go through each year and extract the relevant data.
# Note that attributes are being extracted at the same 
# sp8 points in each year. These points are given an ID
# at the time of extraction. This unique location ID will
# be used as a random effect.
repeat2010<-extract(data2010, sp8, df=T)
names(repeat2010)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2010<-repeat2010$counts
sp8@data$min2010<-repeat2010$minutes_effort
sp8@data$dist2010<-repeat2010$distance_effort
sp8@data$obs2010<-repeat2010$num_obs

repeat2011<-extract(data2011, sp8, df=T)
names(repeat2011)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2011<-repeat2011$counts
sp8@data$min2011<-repeat2011$minutes_effort
sp8@data$dist2011<-repeat2011$distance_effort
sp8@data$obs2011<-repeat2011$num_obs

repeat2012<-extract(data2012, sp8, df=T)
names(repeat2012)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2012<-repeat2012$counts
sp8@data$min2012<-repeat2012$minutes_effort
sp8@data$dist2012<-repeat2012$distance_effort
sp8@data$obs2012<-repeat2012$num_obs

repeat2013<-extract(data2013, sp8, df=T)
names(repeat2013)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2013<-repeat2013$counts
sp8@data$min2013<-repeat2013$minutes_effort
sp8@data$dist2013<-repeat2013$distance_effort
sp8@data$obs2013<-repeat2013$num_obs

repeat2014<-extract(data2014, sp8, df=T)
names(repeat2014)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2014<-repeat2014$counts
sp8@data$min2014<-repeat2014$minutes_effort
sp8@data$dist2014<-repeat2014$distance_effort
sp8@data$obs2014<-repeat2014$num_obs

repeat2015<-extract(data2015, sp8, df=T)
names(repeat2015)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2015<-repeat2015$counts
sp8@data$min2015<-repeat2015$minutes_effort
sp8@data$dist2015<-repeat2015$distance_effort
sp8@data$obs2015<-repeat2015$num_obs

repeat2016<-extract(data2016, sp8, df=T)
names(repeat2016)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2016<-repeat2016$counts
sp8@data$min2016<-repeat2016$minutes_effort
sp8@data$dist2016<-repeat2016$distance_effort
sp8@data$obs2016<-repeat2016$num_obs

repeat2017<-extract(data2017, sp8, df=T)
names(repeat2017)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2017<-repeat2017$counts
sp8@data$min2017<-repeat2017$minutes_effort
sp8@data$dist2017<-repeat2017$distance_effort
sp8@data$obs2017<-repeat2017$num_obs

repeat2018<-extract(data2018, sp8, df=T)
names(repeat2018)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2018<-repeat2018$counts
sp8@data$min2018<-repeat2018$minutes_effort
sp8@data$dist2018<-repeat2018$distance_effort
sp8@data$obs2018<-repeat2018$num_obs

repeat2019<-extract(data2019, sp8, df=T)
names(repeat2019)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2019<-repeat2019$counts
sp8@data$min2019<-repeat2019$minutes_effort
sp8@data$dist2019<-repeat2019$distance_effort
sp8@data$obs2019<-repeat2019$num_obs

repeat2020<-extract(data2020, sp8, df=T)
names(repeat2020)<-c("ID", "counts", "minutes_effort", "distance_effort", "num_obs", "veg_code")
sp8@data$y2020<-repeat2020$counts
sp8@data$min2020<-repeat2020$minutes_effort
sp8@data$dist2020<-repeat2020$distance_effort
sp8@data$obs2020<-repeat2020$num_obs
sp8@data$veg_code2018<-repeat2018$veg_code


# Add vegetation text to match vegetation code:
sp8@data$veg_type<-rep(NA, nrow(sp8@data))
sp8@data$veg_type[which(sp8@data$veg_code2018==3)]<-"riparian"
sp8@data$veg_type[which(sp8@data$veg_code2018==4)]<-"oak"
sp8@data$veg_type[which(sp8@data$veg_code2018==5)]<-"grass"
sp8@data$veg_type[which(sp8@data$veg_code2018==6)]<-"develop"
sp8@data$veg_type[which(sp8@data$veg_code2018==7)]<-"forest"
sp8@data$veg_type[which(sp8@data$veg_code2018==8)]<-"chaparral"
sp8@data$veg_type[which(sp8@data$veg_code2018==9)]<-"scrub"
sp8@data$veg_type[which(sp8@data$veg_code2018==10)]<-"water"
sp8@data$veg_type[which(sp8@data$veg_code2018==11)]<-"ag"
sp8@data$veg_type[which(sp8@data$veg_code2018==13)]<-"barren"
sp8@data$veg_type[which(sp8@data$veg_code2018==14)]<-"wetland"
sp8@data$veg_type[which(sp8@data$veg_code2018==15)]<-"redwood"
sp8@data$veg_type[which(sp8@data$veg_code2018==16)]<-"non-native"
sp8@data$veg_type[which(sp8@data$veg_code2018==17)]<-"conifer"

# Extract the relative abundance data if available:
if (ArcGIS==1) {
eba<-extract(ebird_pred, sp8, df=T)
sp8@data$eBird_abnd_pred<-eba$averxxxxxxrprj
}

# Now go through each year and instead of having the rows be the points
# and the columns be the yearly counts, make each count it's own row
# with a count, year, ID, and effort variables associated with it.
repeat2010$year<-rep(2010, nrow(repeat2010))
countbyyr<-repeat2010
repeat2011$year<-rep(2011, nrow(repeat2011))
countbyyr<-rbind(countbyyr, repeat2011)
repeat2012$year<-rep(2012, nrow(repeat2012))
countbyyr<-rbind(countbyyr, repeat2012)
repeat2013$year<-rep(2013, nrow(repeat2013))
countbyyr<-rbind(countbyyr, repeat2013)
repeat2014$year<-rep(2014, nrow(repeat2014))
countbyyr<-rbind(countbyyr, repeat2014)
repeat2015$year<-rep(2015, nrow(repeat2015))
countbyyr<-rbind(countbyyr, repeat2015)
repeat2016$year<-rep(2016, nrow(repeat2016))
countbyyr<-rbind(countbyyr, repeat2016)
repeat2017$year<-rep(2017, nrow(repeat2017))
countbyyr<-rbind(countbyyr, repeat2017)
repeat2018$year<-rep(2018, nrow(repeat2018))
countbyyr<-rbind(countbyyr, repeat2018)
repeat2019$year<-rep(2019, nrow(repeat2019))
countbyyr<-rbind(countbyyr, repeat2019)
repeat2020$year<-rep(2020, nrow(repeat2020))
countbyyr<-rbind(countbyyr, repeat2020)
countbyyr$year<-as.numeric(countbyyr$year) # make sure that year is numeric

# Add vegetation text to match vegetation code:
countbyyr$veg_type<-rep(NA, nrow(countbyyr))
countbyyr$veg_type[which(countbyyr$veg_code==3)]<-"riparian"
countbyyr$veg_type[which(countbyyr$veg_code==4)]<-"oak"
countbyyr$veg_type[which(countbyyr$veg_code==5)]<-"grass"
countbyyr$veg_type[which(countbyyr$veg_code==6)]<-"develop"
countbyyr$veg_type[which(countbyyr$veg_code==7)]<-"forest"
countbyyr$veg_type[which(countbyyr$veg_code==8)]<-"chaparral"
countbyyr$veg_type[which(countbyyr$veg_code==9)]<-"scrub"
countbyyr$veg_type[which(countbyyr$veg_code==10)]<-"water"
countbyyr$veg_type[which(countbyyr$veg_code==11)]<-"ag"
countbyyr$veg_type[which(countbyyr$veg_code==13)]<-"barren"
countbyyr$veg_type[which(countbyyr$veg_code==14)]<-"wetland"
countbyyr$veg_type[which(countbyyr$veg_code==15)]<-"redwood"
countbyyr$veg_type[which(countbyyr$veg_code==16)]<-"non-native"
countbyyr$veg_type[which(countbyyr$veg_code==17)]<-"conifer"

countbyyr$yr<-countbyyr$year-2010 # Change the year data to be 0-10 instead of 2010-2020.
countbyyr$uniqueID<-countbyyr$ID
countbyyr<-countbyyr[which(is.na(countbyyr$counts)==FALSE),] # Take out NA count data (remember that up to 3/11 years of bird observations can be missing from a site)

######################################
# Plotting
######################################

# Save a copy of the main data table:
countbyyr_saved<-countbyyr

# Plot Every Count
ggplot(countbyyr, aes(x=year, y=counts)) + 
  geom_point(alpha=0.1, size=5)+
  scale_x_continuous(breaks=seq(2010, 2020, 3),limits=c(2010, 2020)) +
  geom_line(data=xx, color="red", size=2)+
  coord_cartesian(ylim=c(0,1.5)) +
  theme_classic() + xlab("Year") + ylab(paste(sppx, " Counts", sep="")) #

# Now Plot Average of Averaged Counts and Rounded Counts:

# Create a summary table of average counts for plotting:
countbyyr2<-as.data.frame(tapply(countbyyr$counts, countbyyr$year, mean, na.rm=TRUE))
names(countbyyr2)<-"mean"
countbyyr2$stdev<-tapply(countbyyr$counts, countbyyr$year, sd, na.rm=TRUE)/sqrt(tapply(countbyyr$counts, countbyyr$year, length))
countbyyr2$year<-row.names(countbyyr2)
countbyyr2$yr<-as.numeric(countbyyr2$year)-2010
countbyyr2$year2<-as.numeric(countbyyr2$year)

pave<-ggplot(countbyyr2, aes(x=year2, y=mean, group=1)) + 
  geom_point(size=3)+
  geom_line() +
  geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), width=.2,
                position=position_dodge(0.05))+
  scale_x_continuous(breaks=seq(2010, 2020, 2)) +
  theme_classic() + xlab("Year") + ylab(paste(sppx, " Abundance", sep=""))
pave

# Now round the data to integer:
countbyyr$countround<-countbyyr$counts

# All values > 0.5 will be rounded to nearest integer.  
# For values <0.5 we can set a lower bound lower than 0.5. This would make values between
# the cut-off and 0.5 equal to 1. 
cutoff<-0.5
countbyyr$countround[intersect(which(countbyyr$countround>cutoff), which(countbyyr$countround<0.50))]<-1
countbyyr$countround<-round(countbyyr$countround, 0)
countbyyr$pres<-ifelse(countbyyr$countround>0, 1, 0)

# Create a summary table of average ROUNDED counts for plotting:
countbyyr2$roundmean<-tapply(countbyyr$countround, countbyyr$year, mean, na.rm=TRUE)
countbyyr2$roundsd<-tapply(countbyyr$countround, countbyyr$year, sd, na.rm=TRUE)/sqrt(nrow(sp8@data))
countbyyr2$year3<-0.1+countbyyr2$year2

pave + geom_line(aes(x=year3, y=roundmean), data=countbyyr2, color="red") + 
  geom_point(aes(x=year3, y=roundmean), data=countbyyr2, size=3, color="red") +
  geom_errorbar(aes(x=year3, ymin=roundmean-roundsd, ymax=roundmean+roundsd), data=countbyyr2, width=.2,
                position=position_dodge(0.05), color="red") +
  scale_x_continuous(breaks=seq(2010, 2020, 2)) +
  theme_classic() + xlab("Year") + ylab(paste(sppx, " Abundance", sep="")) + 
  ggtitle("Black = averaged counts, Red = Averaged ROUNDED counts \n       Ideally, the pattern across years would match") 
  
# If the pattern across years matches for the red and black line, move on.  If not, fuss with the cut-off value.

######################################
# Statsitical Models
######################################

# Mannually, I compare:
# 1) A Poisson Model to a negative binomial model, choosing the one with lower AIC.
# 2) No zero-inflation, intercept-only zero inflation, and intercept + effort minutes 
#    zero inflation. Again, I decide based on AIC.

# Try a Poisson model first:
remod2<-glmmTMB(countround ~ yr + minutes_effort + distance_effort + I(distance_effort^2) + num_obs +
                  (1|uniqueID), family='poisson', ziformula = ~1, #1 + minutes_effort, 
                na.action = na.omit, data=countbyyr) #[countbyyr$year>2013,])
summary(remod2)
# Check residuals:
counts <- simulateResiduals(remod2) 
plot(counts)  #1 - pchisq(2809.7, 941)


# Try a Negative Binomial model next:
remodnbd<-glmmTMB(countround ~ yr + minutes_effort + distance_effort + I(distance_effort^2) + num_obs + 
                  (1|uniqueID), family='nbinom2', ziformula = ~1, na.action = na.omit, data=countbyyr, 
                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))) #[countbyyr$year>2013,])
summary(remodnbd)
# Check residuals:
counts <- simulateResiduals(remodnbd) 
plot(counts)  #1 - pchisq(2809.7, 941)


# Look at landscape characteristics of repeat locations:
ggplot(sp8@data, aes(x=veg_type)) + 
  geom_histogram(stat="count", color="black", fill="gray") + theme_classic() + 
  xlab("Land cover") + 
  ggtitle("Land cover for locations with repeat points")

