##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# This script prepares the provided data set for further processing. Names of 
# required columns are standardised and a taxa list is created.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



prepare_dataset <- function(
  filename_inputData,
  column_scientificName,
  column_taxonName,
  column_location,
  column_eventDate,
  file_name_extension){

  FullTaxaList <- fread(file.path("Data","Input",filename_inputData))
  
  colnames(FullTaxaList)[names(FullTaxaList) == column_scientificName] <- "scientificName"
  colnames(FullTaxaList)[names(FullTaxaList) == column_taxonName] <- "Taxon"
  colnames(FullTaxaList)[names(FullTaxaList) == column_location] <- "Location"
  colnames(FullTaxaList)[names(FullTaxaList) == column_eventDate] <- "eventDate"
  
  ## standardise location names 
  newLocNames <- standardise_location_names(FullTaxaList$Location,file_name_extension)
  if (nrow(newLocNames)!=nrow(FullTaxaList)){
    stop("\n Standardisation of location names went wrong. Check standardise_location_names.R in coords_to_regions.R \n")
  } 
  FullTaxaList$Location <- newLocNames$Location
  # ind <- which(SpecRegionData$Location!=newLocNames$Location)
  
  ## output full taxon x location data set
  fwrite(FullTaxaList,file.path("Data","Output",paste0("FullDataSet_Standardised_",file_name_extension,".gz")))
  
  ## output taxon list
  col_names <- c("scientificName","Taxon")
  if (all(c("Habitat_marine","Habitat_freshwater","Habitat_terrestrial")%in%colnames(FullTaxaList))){
    col_names <- c(col_names,"Habitat_marine","Habitat_freshwater","Habitat_terrestrial")
  } 
  if (all(c("phylum","class")%in%tolower(colnames(FullTaxaList)))){
    colnames(FullTaxaList)[tolower(colnames(FullTaxaList))=="class"] <- "class"
    colnames(FullTaxaList)[tolower(colnames(FullTaxaList))=="phylum"] <- "phylum"
    col_names <- c(col_names,"phylum","class")
  } else {
    cat("\n Warning: Columns 'phylum' and 'class' needed to distinguish realm! \n")
  }
  fwrite(unique(FullTaxaList[,..col_names]),file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))
}

