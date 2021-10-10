##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# This script loads all required functions of the workflow.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# Authors: Hanno Seebens, Ekin Kaplan, 23.02.2021
##################################################################################



########################################################################
### load scripts #######################################################

## check and create folder structure ###################################
source(file.path("R","create_folders.R")) # 

## prepare data set ####################################################
source(file.path("R","prepare_dataset.R")) # a function to decompress large zip files

## obtain and process GBIF records #####################################
source(file.path("R","send_GBIF_request.R")) # a function to decompress large zip files
source(file.path("R","get_GBIF_download.R")) # downloads requested data from GBIF
source(file.path("R","decompress_file.R")) # a function to decompress large zip files
source(file.path("R","extract_GBIF_columns.R")) # a function to extract from zipped GBIF downloads

## obtain and process OBIS records #####################################
source(file.path("R","get_OBIS_records.R")) # downloads requested data from GBIF

## clean coordinates ####################################################
source(file.path("R","clean_GBIF_records.R")) # a function to clean GBIF records
source(file.path("R","clean_OBIS_records.R")) # create shapefile of marine and terrestrial polygons

## assign coordinates to regions and identify alien populations #########
source(file.path("R","get_WoRMS_habitats.R")) # get habitat information from WoRMS
source(file.path("R","get_habitats.R")) # get habitat information from WoRMS, Fishbase and Sealifebase
source(file.path("R","get_habitats_DASCO.R")) # get habitat information from WoRMS, Fishbase and Sealifebase
source(file.path("R","standardise_location_names.R")) # standardise location names (for matching with shapefile)

source(file.path("R","coords_to_regions_GBIF.R")) # identify region for each coordinate
source(file.path("R","coords_to_regions_OBIS.R")) # identify region for each coordinate

## add first records to final output file ###############################
source(file.path("R","final_DASCO_output.R")) # add first records per species and region (if available)
