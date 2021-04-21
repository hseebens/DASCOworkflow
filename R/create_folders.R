##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# The scripts checks the existence of data folders and generates them if required.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################


create_folders <- function (){
  
  ## create output folder #####
  if (!file.exists("Data")){
    dir.create("Data")
  }
  if (!file.exists(file.path("Data","Output"))){
    dir.create(file.path("Data","Output"))
  }
  if (!file.exists(file.path("Data","Output","Intermediate"))){
    dir.create(file.path("Data","Output","Intermediate"))
  }
  if (!file.exists(file.path("Data","Input"))){
    dir.create(file.path("Data","Input"))
  }
  if (!file.exists(file.path("Data","Input","Shapefiles"))){
    dir.create(file.path("Data","Input","Shapefiles"))
  }
}  