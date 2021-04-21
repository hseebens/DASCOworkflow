##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# This script constitutes the main file, which can be used to execute the whole 
# workflow in sequence. The user has to provide basic information about the 
# checklist to be used, the paths for storing data and different options for 
# the application of the workflow.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################


graphics.off()
rm(list=ls())


library(data.table) # for clean_GBIF_records, request_GBIF_download
library(rgbif) # for clean_GBIF_records, request_GBIF_download
library(worrms)
library(robis)
library(CoordinateCleaner) # for clean_GBIF_records
library(httr)
library(sf)   # for transform_coords_to_regions



###################################################################################
## load functions #################################################################
source(file.path("R","load_functions.R")) # load all required functions



###################################################################################
### Global variables ##############################################################

### Within this workflow, files will be downloaded and stored in these folders
### Note: All files in that folder will be considered as relevant files
path_to_GBIFdownloads <- "/home/hanno/Storage_large/GBIF/FirstRecords_Mar2021"
path_to_OBISdownloads <- "/home/hanno/Storage_large/OBIS"

## has to be stored in Data/Input/ and has to include a column named 'scientificName'
## for taxon names and 'Location' for region names and 'Taxon' (no authority) for habitat check
# name_of_specieslist <- "SInAS_AlienSpeciesDB_2.3.1_FullTaxaList.csv"
filename_inputData <- "IntroDat_22Mar2021.csv"
column_scientificName <- "scientificName" # taxon name with or without authority; require for GBIF
column_taxonName <- "TaxonName" # taxon name without authority; required for OBIS
column_location <- "Region" # column name of location records
column_eventDate <- "FirstRecord" # column name of year of first record of occurrence

## Name of file with the information of alien species and regions
# name_of_TaxonLoc <- "IntroDat_22Mar2021.csv"

## name of shapefile providing polygons for the new delineation
name_of_shapefile <- "RegionsTerrMarine"

file_name_extension <- "FirstRecords"


###################################################################################
### Check and create folder structure #############################################
create_folders() # creates folders only if they are not existing yet

prepare_dataset(filename_inputData,column_scientificName,column_taxonName,column_location,column_eventDate,file_name_extension)
  

###################################################################################
### Obtaining data ################################################################
### send requests to GBIF #########################################################

## GBIF account details ############
## Note that multiple accounts are required for n_accounts>1.
## The accounts have to numbered x=1...n_accounts, while x is part of
## user name and email address. For example, user name and email should be:
## (ekinhanno1, ekinhanno1@gmail.com), (ekinhanno2, ekinhanno2@gmail.com) and so on.

n_accounts <- 7

## login details for first account (x=1) (the '1' in user name and email
## address will be replaced be account number)
user <- "ekinhanno1"                                  # your gbif.org username
pwd <- "seebenskaplan1234"                                     # your gbif.org password (set the same password for all accounts for convenience)
email <- "ekinhanno1@outlook.com"                 # your email which you will recieve the download link

###################################################################################

## send requests to GBIF 
send_GBIF_request(file_name_extension,n_accounts,user=user,pwd=pwd,email=email)

## get downloads from GBIF (requires running 'send_GBIF_request' first)
get_GBIF_download(path_to_GBIFdownloads)

### extract relevant information from GBIF downloads ##############################
extract_GBIF_columns(path_to_GBIFdownloads,file_name_extension)


###################################################################################
### get OBIS records ##############################################################

get_OBIS_records(path_to_OBISdownloads,file_name_extension)


###################################################################################
### Cleaning data #################################################################

## Thinning of coordinates:
## High numbers of records may cause memory issues. The number of records
## can be reduced by setting thin_records to TRUE. Then, coordinates are
## rounded to the second digit and duplicates are removed. For the single
## remaining record at this site, the original (not the rounded) record
## is kept. Thinning may result in imprecise results when the regions
## considered later in the workflow are small.
thin_records <- T

### clean GBIF records ############################################################

clean_GBIF_records(path_to_GBIFdownloads,file_name_extension,thin_records)

### clean OBIS records ############################################################

thin_records <- F

clean_OBIS_records(path_to_OBISdownloads,file_name_extension,thin_records)
  

###################################################################################
### get alien regions based on coordintates #######################################

## Create shapefile of terrestrial and marine polygons
## loads and combines shapefiles and stores the final shapefile in Data/Input/Shapefiles
terrestrial_polygons <- "RegionsShapefile_200121" # name of terrestrial shapefile
marine_polygons <- "meow_ecos" # name of marine shapefile

create_shapefile(terrestrial_polygons,marine_polygons)


## Assign coordinates to different realms (terrestrial, freshwater, marine)
## depending on geographic location and additional tests
realm_extension <- TRUE 

## assigning country checklists to marine polygons
# checklist_to_marine <- TRUE

# Region shapefile requires a consistent structure for marine and terrestrial polygons !!!!!

coords_to_regions_GBIF(name_of_shapefile,
                       path_to_GBIFdownloads,
                       realm_extension,
                       file_name_extension)

coords_to_regions_OBIS(name_of_shapefile,
                       path_to_OBISdownloads,
                       realm_extension,
                       file_name_extension)


########################################################################
## add first records per region (requires 'eventDate' column) ##########
## and produce final output file containing GBIF and OBIS records ######
dat <- add_first_records(file_name_extension,path_to_GBIFdownloads)
