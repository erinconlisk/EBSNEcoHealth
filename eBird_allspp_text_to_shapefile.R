# File for processing eBird ERD checklists downloaded on July 19, 2021.
# Created by Erin Conlisk from April-August 2021

# This code takes an eBird text file and creates GIS shapefiles for the 
# observations that lie within the East Bay stewardship Network (EBSN) 
# agency boundaries for a specific species in a specific year.

# The folder that is created for the files will be made in this script 
# and named as a four-letter code made up of the first two letters of 
# the genus combined with the first two letters of the species.  
# Within the folder will be files for each year, labeled with the same
# four-letter code followed by the year.

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
# I barely use this package. However, it might be useful with larger eBird
# datasets. Thus, I list the link that describes the package:
# https://cornelllabofornithology.github.io/auk/
library(auk)

# The next package downloads spatial layers of eBird relative population trends
# More info can be found here: https://cornelllabofornithology.github.io/ebirdst/
library(ebirdst)

# As the eBird dataset evolves, both of these packages could be subject to changes,
# which may require changes to the code.

# This code depends on one file: "AgencyBoundary_dissolve_rprj" which defines the
# EBSN agency boundaries.  You must make sure that the projection for the eBird 
# files and the AgencyBoundary_dissolve_rprj are the same.

# This code requires that you choose a species to focus on (the list of species is
# at ~line 150).

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

ebirdabundrast<-0

######################################
# Extract most relevant eBird data
######################################

# Set the path:
ebirdpath<-"G:/EcoHealthEastBay/"

# Load new ebird data.  The title of this file is reflects that I extracted a small portion of the 
# eBird dataset.  After "ebd_US-" is "CA", the state I used to extract data and then the dates
# are provided: 201001 for Jan 2010 to 202012 or Dec 2020:
dat <- fread(paste(ebirdpath, "ebd_US-CA_201001_202012_relJan-2021/ebd_US-CA_201001_202012_relJan-2021.txt", sep=""),header=T) #select=colsToKeep) 
# To download other eBird data you will need to "register", sign in, and request access.  The process starts here:
# https://ebird.org/data/download

# Trim to focal region within CA
CA_lat<-subset(dat, dat$LATITUDE>36.9 & dat$LATITUDE<38.1)
EB<-subset(CA_lat, CA_lat$LONGITUDE<-121.2 & CA_lat$LONGITUDE>-122.4)

# Re-format date variables so that you can parse data by day, month, and year:
EB$yr<-substr(EB$'OBSERVATION DATE', 1, 4)
EB$month<-substr(EB$'OBSERVATION DATE', 6, 7)
EB$day<-substr(EB$'OBSERVATION DATE', 9, 10)

# Re-format time variable as a number:
EB$hour<-as.numeric(substr(EB$'TIME OBSERVATIONS STARTED', 1, 2))

# Some species can have synonyms or imprecise identification (e.g. identified
# to the genus instead of the species). For each species, we checked for this.
# The first step was to find all the unique species names and how many observations
# there were for each species:
#xx<-as.data.frame(tapply(EB$'SCIENTIFIC NAME', EB$'SCIENTIFIC NAME', length))
#names(xx)[1]<-"n_obs"
#xx$species<-rownames(xx)
#  write.csv(xx, "G:/EcoHealthEastBay/species_in_eBird.csv", sep="")

# eBird creates unique locationIDs.  If we combine the location with the date, we 
# can remove data points occurring at the same spot on the same day (e.g. because
# a school group was looking at birds and all students reported their checklists.)
# Create a column in the dataset for unique location-date observations 
EB$date_locality<-paste(EB$'OBSERVATION DATE', EB$'LOCALITY ID', sep="-")
# In the script "eBird_allspp_grid_analyze_and_figures.R", we will trim this further 
# so that no observations can happen at sites within 100 meters of one another.

# Save the file for later use because loading large datasets into R can be time consuming.
#save.image(paste(ebirdpath, "eBird_loaded.RData", sep=""))

# Then you can load the data from here:
#load(paste(ebirdpath, "eBird_loaded.RData", sep=""))

################
# DEFINE SPECIES
# This is probably not the most elegant way to do this, but it works.
# I need to assign the species and all of the "synonyms" for the species
# to accurately assign presence and absence.  So I list all the species 
# and their potential synonyms.  

# Only the species you are looking at should be activated (as opposed to greened out with "#")
# I have a comment "200m resolution" for species with larger territory sizes.
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

# Be sure there are no synonyms or subspecies in the dataset
speciesx<-subset(EB, EB$'SCIENTIFIC NAME'==sppx)
length(which(speciesx$'ALL SPECIES REPORTED'!=1))/nrow(speciesx) # Outputs the recorded presences that are NOT from "complete checklists" (typically ~4%)

# Keep only complete checklists to use as zeros.
# Get rid of imprecise designations that could be the target species:
if (sppx=="Aphelocoma californica") {
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Corvidae sp. (jay sp.)")
} else if (sppx=="Pheucticus melanocephalus") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Pheucticus ludovicianus x melanocephalus")
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Pheucticus ludovicianus/melanocephalus")
} else if (sppx=="Dryobates pubescens") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Dryobates pubescens x nuttallii")
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Dryobates pubescens/villosus")
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Dryobates sp.")
} else if (sppx=="Circus hudsonius") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Accipiter sp.")
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Accipitridae sp. (hawk sp.)")
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Accipitriformes/Falconiformes sp.")
} else if (sppx=="Dryobates nuttallii") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Dryobates nuttallii x scalaris")  # Seems like there are lots of warblers, we aren't going to do this are we?
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Dryobates nuttallii/scalaris")  # Seems like there are lots of warblers, we aren't going to do this are we?
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Dryobates sp.")  # Seems like there are lots of warblers, we aren't going to do this are we?
} else if (sppx=="Baeolophus inornatus") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Baeolophus inornatus/ridgwayi") 
} else if (sppx=="Passerculus sandwichensis") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Passerellidae sp. (sparrow sp.)")
} else if (sppx=="Tachycineta bicolor") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Tachycineta bicolor/thalassina") # Quite a few hybrids
} else if (sppx=="Vireo gilvus") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Vireo sp.")
} else if (sppx=="Sitta carolinensis") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Sitta sp.") 
} else if (sppx=="Sialia mexicana") {  
  putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
  putative_zeros<-subset(putative_zeros, putative_zeros$'SCIENTIFIC NAME'!="Sialia sp.")
} else {
putative_zeros<-subset(EB, EB$'ALL SPECIES REPORTED'==1)
}


# Need to fix the "X" values in the dataset.  Will switch to the median value for the dataset.
# POSSIBILITY FOR REFINEMENT: set X to the mean count, or set X to 1. Make the median count vary by year...
speciesx$counts<-speciesx$`OBSERVATION COUNT`
length(which(speciesx$`OBSERVATION COUNT`=="X"))/nrow(speciesx) # Out puts the fraction of recorded presences that have counts as "X" (typically ~4%)
speciesxnox<-subset(speciesx, speciesx$`OBSERVATION COUNT`!="X") 
speciesx$counts[which(speciesx$`OBSERVATION COUNT`=="X")]<-median(as.numeric(speciesxnox$`OBSERVATION COUNT`)) 

# Need to remove duplicates of speciesx. 
# Where there are duplicates, we average the value reported at the date and location
dupspeciesx<-unique(speciesx$date_locality[which(duplicated(speciesx$date_locality)==TRUE)])
# The following should be have no duplicates 
speciesx_nodups<-speciesx[which(duplicated(speciesx$date_locality)==FALSE),]

# For table with one row for each location-date, we go back and average the counts for location-dates with duplicate counts.
# POSSIBILITY FOR REFINEMENT: take the larger of the two counts, or the count that occurred during the most active part of the day (e.g. morning)...
for (ii in dupspeciesx) {
  speciesx_nodups$counts[which(speciesx_nodups$date_locality==ii)]<-mean(as.numeric(speciesx$counts[which(speciesx$date_locality==ii)]))
}

# Trim all duplicates in ZERO values based on date-locality
putative_zeros<-putative_zeros[which(duplicated(putative_zeros$date_locality)==FALSE),]
# Set all species counts to zero for the checklist data; at this point
# in the script some of these unduplicated zeros may overlap with a 
# location where the species was present
putative_zeros$counts<-rep(0, nrow(putative_zeros))
# Combine species presences - which MUST GO FIRST - with species absences
combo<-rbind(speciesx_nodups, putative_zeros)
# Given that the presences are first in the data table, any zero counts with 
# duplicate location-dates will be removed.
speciesx_all<-combo[which(duplicated(combo$date_locality)==FALSE),]

# Trim to April 1 to July 15
speciesx_apr<-subset(speciesx_all, speciesx_all$month=="04")
speciesx_may<-subset(speciesx_all, speciesx_all$month=="05")
speciesx_jun<-subset(speciesx_all, speciesx_all$month=="06")
speciesx_jul<-speciesx_all[intersect(which(speciesx_all$month=="07"), which(speciesx_all$day<16)),]
speciesx45<-rbind(speciesx_apr, speciesx_may)
speciesx67<-rbind(speciesx_jun, speciesx_jul)
blue<-rbind(speciesx45, speciesx67)

# Rename the data table headings to remove spaces in the names.
names(blue)<-c("GLOBAL_ID", "LAST_EDITED", "TAXMC_ORDER", "CATEGORY", "COM_NAME", "SCI_NAME", "SUBSPECIES_COM", "SUBSPECIES_SCI",  
                "OBS_COUNT", "BB_ATLAS_CODE", "BB_ATLAS_CAT", "AGE_SEX", "COUNTRY", "COUNTRY_CODE", "STATE", "STATE_CODE",              
                "COUNTY", "COUNTY_CODE", "IBA_CODE", "BCR_CODE", "USFWS_CODE", "ATLAS_BLOCK", "LOCALITY", "LOCALITY_ID",                
                "LOCALITY_TYPE", "LATITUDE", "LONGITUDE", "DATE", "TIME_OBS_START", "OBS_ID", "SAMPL_EVENTID", "PROTOCOL", 
                "PROTCL_CODE", "PROJECT_CODE", "DUR_MINUTES", "EFFORT_KM", "EFFORT_HA", "NUM_OBS", "ALL_SPP_YN", "GROUPID",            
                "HAS_MEDIA", "APPROVED", "REVIEWED", "REASON", "TRIP_CMMT", "SPECIES_CMMT", "V47", "yr", "month", "day", 
                "hour", "date_locality", "counts" )

# Remove extra data columns for ease of export. Even shortening the names, some 
# may be further abbreivated when they are exported to a shapefile.
blue2<-blue[,c("GLOBAL_ID", "COM_NAME", "SCI_NAME", "SUBSPECIES_SCI", "OBS_COUNT", "COUNTY_CODE", 
              "LOCALITY_ID", "LOCALITY_TYPE", "LATITUDE", "LONGITUDE", "TIME_OBS_START", "OBS_ID", "SAMPL_EVENTID", 
              "PROTOCOL", "PROTCL_CODE", "DUR_MINUTES", "EFFORT_KM", "EFFORT_HA", "NUM_OBS", "ALL_SPP_YN", "GROUPID",            
              "yr", "month", "day", "hour", "date_locality", "counts")]


################################################
# Make eBird Points into spatial file and clip 
# to region of interest (EBSN). 
# Save the files for later use.
################################################

# To make files we abbreviate the species name to the first two characters of the 
# genus and the species.  For example, Aimophila ruficeps=Airu
spaceis<-regexpr(" ", sppx)[1]
sppname<-paste(substr(sppx, 1,2), substr(sppx, (spaceis+1),(spaceis+2)), sep="") 

# Create a folder to house files:
dir.create(paste(ebirdpath, sppname, "_eBird_data_by_yr", sep=""))

# Upload the reprojected polygon that will be used to trim and assign CRS for points
#eb_poly<-readOGR("G:/EcoHealthEastBay", layer="EHA_dissolved_Landscapeunits_rprj") # Larger landscape units
# For shapefiles, we have to set a second path without the "/" at the end of the folder name.
ebirdpath2<-"G:/EcoHealthEastBay"
eb_poly<-readOGR(ebirdpath2, layer="AgencyBoundary_dissolve_rprj") # Agency boundaries

# Transform to spatial points based on latitude and longitude
for (years in seq(2010, 2020, 1)) { # Separate data by years
  blue<-subset(blue2, as.numeric(blue2$yr)==years)
  coordssppx<-cbind(blue$LONGITUDE, blue$LATITUDE)
  eval(parse(text=paste(sppname, years, " <- SpatialPointsDataFrame(coords = coordssppx, data=blue, # set the coordinates to the latitude and longitude and keep the data.
                                     proj4string = CRS('+proj=longlat +ellps=GRS80 +no_defs'))", sep="")))  # Assign CRS
  eval(parse(text=paste(sppname, years, "_subset <- ", sppname, years, "[eb_poly, ]", sep=""))) # Mask the points by the agency boundaries
  eval(parse(text=paste("writeOGR(", sppname, years, "_subset, '", ebirdpath, sppname, "_eBird_data_by_yr', ",  "layer='", sppname, years, "', driver='ESRI Shapefile')", sep="")))
}
