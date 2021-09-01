# EBSNEcoHealth
The scripts associated with the EBSN Ecohealth assessment.  The objective was to describe trends in presence-absence and abundance across 28 species, using eBird data.

The following Google doc has a narrative explaination of the methods used in the scripts:

https://docs.google.com/document/d/11aJk4GtvEyvmt35_tbMy_7Ou1yBdyCalw1XEjRxXi9k/edit?pli=1

A description of the scripts is as follows:

################################

eBird_allspp_text_to_shapefile.R  

################################

This code takes an eBird text file (an excerpt of the full eBird data) and creates GIS shapefiles for the observations that lie within the EBSN agency boundaries for a specific species in a specific year.  In the resulting shapefiles, repeat locations on a given day are removed but multiple observations within a 100 or 200-m grid cell have not been averaged nor have the grid cells been averaged over the breeding season. 

The folder that is created for the files will be made in this script and named as a four-letter code made up of the first two letters of the genus combined with the first two letters of the species.  For example, rufous-crowned sparrow’s scientific name is Aimophila ruficeps and thus the folder will have the label “Airu” (followed by “_eBird_data_by_yr”).  Within the folder will be files for each year, labeled with the same four-letter code followed by the year. 

Before beginning the script, you need:

1. To download other eBird data.  You will need to "register", sign in, and request access.  The process starts here: https://ebird.org/data/download
1.i. I downloaded data for CA across January 2010 to December 2020.  The resulting file is titled: "ebd_US-CA_201001_202012_relJan-2021.txt"
The title of this file reflects its contents.  The "ebd_US-" means eBird USA, followed by "CA", the state, and the dates 201001 for Jan 2010 and 202012 for Dec 2020
1.ii. This excerpt of eBird data is very large and takes a long time to load. In order to be able to continue loading eBird data directly, you need to make sure that you request a manageable excerpt of data.  Trying to upload the whole eBird data is not possible in R without using the eBird package “auk”. I found this package to be difficult to use and so did the data processing myself.
1.iii. Note in the script there is a command: 
save.image(paste(ebirdpath, "eBird_loaded.RData", sep=""))
to save the text file as an R project after some initial processing has been completed.  This is advisable.  
You can then load the data with the command:
		load(paste(ebirdpath, "eBird_loaded.RData", sep=""))
where “ebirdpath” is the directory to where you want to save the file (e.g. “C:/”).
1.iv. Finally, the eBird data file is too large to upload to the repository or share by email. Here is a link to the eBird data used in this assessment:
XXXXXXXXXXXXXXXXXXXXXX

2. The file "AgencyBoundary_dissolve_rprj" which defines the EBSN agency boundaries. 

3. To choose the species.  See the list of scientific names in the script, starting around line 150. 
4. To install the packages listed at the top with elements “library(PACKAGE NAME)”. To install packages use the command: install.packages("PACKAGE NAME").



################################

eBird_allspp_abundance_rasters.R

################################

This code downloads all the eBird relative abundance rasters within the East Bay Stewardship Network (EBSN) for each week for a specific species and averages them across the breeding season. Relative abundance rasters were created by Cornell to use all of the North American eBird data to model locations where abundances are expected to be highest or lowest for each week across the year (they do not include trends across years).  

These rasters are not necessary for the EBSN Eco-Health analysis of trends.  Instead, these rasters were used as a background check for biases in the EBSN eBird data.  For example, observations made in 2010 could theoretically occur only in locations where Cornell predicted the highest relative abundance whereas observations made in 2020 could occur only in locations with the lowest relative abundance.  If this were the case then the 2010 data would be occurring in locations with high habitat suitability whereas the 2020 data would be in locations with low habitat suitability.  This would introduce significant bias.  

Unfortunately, the eBird relative abundance rasters were not significantly resolved enough to expose biases.  While we saw little bias in the eBird relative abundance values across locations and years, we did see bias when we looked at the land cover types underlying the eBird observations.  Because the land cover background check exposed more bias, we focused on land cover as the primary means to mitigate biased sampling in the EBSN data.

Before beginning the script, you will need:

1. ArcGIS to properly process the relative abundance rasters that are created in this script.  See line ~275 to see a description of the tasks requiring ArcGIS.  I also provide the beginning of R code to do the same work, but this code has not been vetted.  

2. To choose the species.  See the list of scientific names in the script, starting around line 110. 

3. To install the packages listed at the top with elements “library(PACKAGE NAME)”. To install packages use the command: install.packages("PACKAGE NAME"). These packages should be the same packages as were needed for “eBird_allspp_text_to_shapefile.R”



#######################################

eBird_allspp_grid_analyze_and_figures.R

#######################################

This code takes the output of "eBird_text_to_shapefile_allapp.R" and averages the observations within a grid cell on a specific day. Then the resulting daily grids are averaged across the breeding season. The result is one grid per year of average counts (across the grid cells and across days) that is used in trend analyses. As described in the methods above, two analyses are performed on (a) presence-absence data across all eBird observations, and (b) abundance within locations with repeat observations. 

Before beginning the script, you will need:

1. The file "focal_buffer_dissolve" to create boundaries for point density calculations. The script runs without this file, but without the figure you will see error messages for the calculations of point density. 

2. To set the variable ArcGIS to either 0 or 1, depending on whether you have calculated, resampled, and reprojected the eBird relative abundance files (these rasters are created with the script “eBird_allspp_abundance_rasters.R”).

3. The file "focal_rast" OR "foc_rast200m"
the former is for species with max territory size less than 4 ha (which uses a grid cell resolution of 100m) and 
the latter is for all species with max territory size greater than 4 ha (which uses a grid cell resolution of 200 m). 

4. To choose the species.  See the list of scientific names in the script, starting around line 120. 

5. To install the packages listed at the top with elements “library(PACKAGE NAME)”. To install packages use the command: install.packages("PACKAGE NAME"). These packages should be the same packages as were needed for “eBird_allspp_text_to_shapefile.R”

