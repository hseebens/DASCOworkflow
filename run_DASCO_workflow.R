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
# Authors: Hanno Seebens, Ekin Kaplan, 13.06.2021
##################################################################################


graphics.off()
rm(list=ls())


library(data.table) # for clean_GBIF_records, request_GBIF_download
library(httr)
library(sf)   # for transform_coords_to_regions
library(utils)
library(rgbif) # for clean_GBIF_records, request_GBIF_download
library(worrms)
library(robis)
library(CoordinateCleaner) # for clean_GBIF_records
library(rfishbase)

###################################################################################
## load functions #################################################################
source(file.path("R","load_functions.R")) # load all required functions



###################################################################################
### Global variables ##############################################################

### Within this workflow, files will be downloaded and stored in these folders
### Note: All files in that folder will be considered as relevant files; old files should be removed
path_to_GBIFdownloads <- file.path("path","to","GBIF","folder")
path_to_OBISdownloads <- file.path("path","to","OBIS","folder")
# path_to_GBIFdownloads <- "/home/hanno/Storage_large/GBIF/SInASdata/Germany_200522"
# path_to_OBISdownloads <- "/home/hanno/Storage_large/OBIS/SInASdata/Germany_200522"

## has to be stored in Data/Input/ and has to include a column named 'scientificName'
## for taxon names and 'Location' for region names and 'Taxon' (no authority) for habitat check
filename_inputData <- "SInAS_AlienSpeciesDB_2.4.1_Full+Taxonomy.csv"

column_scientificName <- "scientificName" # taxon name with or without authority; require for GBIF
column_taxonName <- "Taxon" # taxon name without authority; required for OBIS
column_location <- "Location" # column name of location records
column_eventDate <- "eventDate" # column name of year of first record of occurrence
column_habitat <- "habitat" # column name of year of first record of occurrence

## Name of file with the information of alien species and regions
# name_of_TaxonLoc <- "IntroDat_22Mar2021.csv"

## name of shapefile providing polygons for the new delineation
name_of_shapefile <- "RegionsTerrMarine_160621"

## term to be added to the names of the output files; can be blank
file_name_extension <- "SInAS_2.4.1"


## check if folders and files exist
if (!dir.exists(path_to_GBIFdownloads)) stop(paste0("Folder '",path_to_GBIFdownloads,"’ does not exist!"))
if (!dir.exists(path_to_OBISdownloads)) stop(paste0("Folder '",path_to_OBISdownloads,"’ does not exist!"))
if (!file.exists(file.path("Data","Input",filename_inputData))) stop(paste0("File '",filename_inputData,"’ could not be found in 'Input' folder!"))
 

###################################################################################
## GBIF account details ###########################################################
## Note that multiple accounts are required for n_accounts>1.
## The accounts have to numbered x=1...n_accounts, while x is part of
## user name and email address. For example, user name and email should be:
## (ekinhanno1, ekinhanno1@gmail.com), (ekinhanno2, ekinhanno2@gmail.com) and so on.

n_accounts <- 1

## login details for first account (x=1) (the '1' in user name and email
## address will be replaced be account number (i.e., 1:n_accounts)
user <- ""                                  # your gbif.org username
pwd <- ""                                     # your gbif.org password (set the same password for all accounts for convenience)
email <- ""                 # your email which you will recieve the download link



###################################################################################
### EXECUTION OF DASCO WORKFLOW ###################################################
###################################################################################

### 1. Check and create folder structure #############################################
create_folders() # creates folders only if they are not existing yet

prepare_dataset(filename_inputData,
                column_scientificName,
                column_taxonName,
                column_location,
                column_eventDate,
                column_habitat,
                file_name_extension)


###################################################################################
### 2. Obtaining data #############################################################

## send requests to GBIF
send_GBIF_request(file_name_extension,
                  path_to_GBIFdownloads,
                  n_accounts,
                  user=user,
                  pwd=pwd,
                  email=email)

## get downloads from GBIF (requires running 'send_GBIF_request' first)
get_GBIF_download(path_to_GBIFdownloads,
                  file_name_extension)

### extract relevant information from GBIF downloads ##############################
extract_GBIF_columns(path_to_GBIFdownloads,
                     file_name_extension)


### get OBIS records ##############################################################

get_OBIS_records(path_to_OBISdownloads,
                 file_name_extension)
## Intermediate download files are stored under Data/Output/Intermediate


###################################################################################
### 3. Cleaning data ##############################################################

## Thinning of coordinates:
## High numbers of records may cause memory issues. The number of records
## can be reduced by setting thin_records to TRUE. Then, coordinates are
## rounded to the second digit and duplicates are removed. For the single
## remaining record at this site, the original (not the rounded) record
## is kept. Thinning may result in imprecise results when the regions
## considered later in the workflow are small.

### clean GBIF records ############################################################

clean_GBIF_records(path_to_GBIFdownloads,
                   file_name_extension,
                   thin_records=FALSE)

### clean OBIS records ############################################################

clean_OBIS_records(path_to_OBISdownloads,
                   file_name_extension,
                   thin_records=FALSE)


###################################################################################
### 4. get alien regions based on coordintates ####################################

## get habitat information for taxa (terrestrial, freshwater, marine, brackish)
get_habitats_DASCO(file_name_extension,path_to_GBIFdownloads,path_to_OBISdownloads)

## Assign coordinates to different realms (terrestrial, freshwater, marine)
## depending on geographic location and additional tests
## realm_extension <- TRUE

# Region shapefile requires a consistent structure for marine and terrestrial polygons !!!!!

coords_to_regions_GBIF(name_of_shapefile,
                       realm_extension=TRUE,
                       file_name_extension)

coords_to_regions_OBIS(name_of_shapefile,
                       realm_extension=TRUE,
                       file_name_extension)


########################################################################
## 5. produce final output of the DASCO workflow #######################
## add first records per region (requires 'eventDate' column) ##########
dat <- final_DASCO_output(file_name_extension)

