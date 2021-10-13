##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# This script adds event dates (i.e., years of first record) to the final data set.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



final_DASCO_output <- function(
  file_name_extension
  ){

  SpecRegionData <-  fread(file.path("Data","Output",paste0("FullDataSet_Standardised_",file_name_extension,".gz")),stringsAsFactors = F,header=T)
  taxlist <- fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))
  
  ### GBIF records ##############################################################
  
  ## load data ###############################

  meow_records <- fread(file.path("Data","Output","Intermediate","MarineRecords_GBIF.gz"))

  if (!"eventDate"%in%colnames(SpecRegionData)){
    cat("\n Column 'eventDate' is missing. Cannot add first records. \n")
  }

  all_records_spec <- fread(file.path("Data","Output",paste0("DASCO_GBIFregions_",file_name_extension,".csv")))

  GBIF_keys <- fread(file.path("Data","Output",paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  GBIF_keys <- unique(GBIF_keys[,c("speciesKey","scientificName")])
  
  SpecRegionData_keys <- merge(SpecRegionData,GBIF_keys,by="scientificName",all.x=T)
  
  ## create empty dummy variables to ensure proper merging; if records exists, these files are overwritten further down
  all_regspec_fr_GBIF <- data.frame(location=character(),scientificName=character(),speciesKey=integer(),Realm=character(),eventDate=integer())
  marine_regspec_fr_GBIF <- data.frame(MEOW=character(),scientificName=character(),speciesKey=integer(),Realm=character(),eventDate=integer())
  terr_regspec_fr_GBIF <- data.frame(location=character(),scientificName=character(),speciesKey=integer(),Realm=character(),eventDate=integer())

  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="location"] <- "MEOW"
  # if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  col_names <- colnames(meow_records)[which(colnames(meow_records)!="scientificName")]
  marine_regs_species_GBIF <- merge(all_records_spec,meow_records[,..col_names],by=c("speciesKey","MEOW"),all.x=T)
  marine_regs_species_GBIF <- subset(marine_regs_species_GBIF,Realm=="marine")
  
  if (nrow(marine_regs_species_GBIF)>0){
    marine_regs_species_GBIF$eventDate[is.na(marine_regs_species_GBIF$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
    marine_regspec_fr_GBIF <- aggregate(eventDate ~ MEOW + scientificName + speciesKey + Realm,data=marine_regs_species_GBIF,FUN=min) # + Source
    marine_regspec_fr_GBIF$eventDate[marine_regspec_fr_GBIF$eventDate==2500] <- NA # remove dummy variable
  }
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "location"
  col_names <- colnames(all_records_spec)[which(colnames(all_records_spec)!="scientificName")]
  terr_regs_species_GBIF <- merge(unique(all_records_spec[,..col_names]),SpecRegionData_keys,by=c("speciesKey","location"),all.y=T)
  terr_regs_species_GBIF <- subset(terr_regs_species_GBIF,Realm!="marine")
  
  if (nrow(terr_regs_species_GBIF)>0){
    terr_regs_species_GBIF$eventDate[is.na(terr_regs_species_GBIF$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
    terr_regspec_fr_GBIF <- aggregate(eventDate ~ location + scientificName + speciesKey + Realm,data=terr_regs_species_GBIF,FUN=min)# + Source
    terr_regspec_fr_GBIF$eventDate[terr_regspec_fr_GBIF$eventDate==2500] <- NA
  }
  
  ## combine terrestrial and marine first records #################################
  colnames(marine_regspec_fr_GBIF)[colnames(marine_regspec_fr_GBIF)=="MEOW"] <- "location"
  all_regspec_fr_GBIF <- rbind(marine_regspec_fr_GBIF,terr_regspec_fr_GBIF)
  all_regspec_fr_GBIF <- unique(all_regspec_fr_GBIF[,c("location","scientificName","Realm","eventDate")])
  
  ## add canonical names to GBIF ###################################################
  col_names <- c("scientificName","taxon")
  taxlist <- unique(taxlist[,..col_names])
  all_regspec_fr_GBIF <- merge(all_regspec_fr_GBIF,taxlist,by="scientificName",all.x=T)
  all_regspec_fr_GBIF <- all_regspec_fr_GBIF[,c("location","scientificName","taxon","Realm","eventDate")]

  
  
  ### OBIS records ##############################################################
  
  ## load data ###############################
  
  meow_records <- fread(file.path("Data","Output","Intermediate","MarineRecords_OBIS.gz"))

  all_records_spec <- fread(file.path("Data","Output",paste0("DASCO_OBISregions_",file_name_extension,".csv")))

  ## create empty dummy variables to ensure proper merging; if records exists, these files are overwritten further down
  all_regspec_fr_OBIS <- data.frame(location=character(),scientificName=character(),Realm=character(),eventDate=integer())
  marine_regspec_fr_OBIS <- data.frame(location=character(),taxon=character(),Realm=character(),eventDate=integer())
  terr_regspec_fr_OBIS <- data.frame(location=character(),taxon=character(),Realm=character(),eventDate=integer())

  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="location"] <- "MEOW"
  # if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  col_names <- colnames(meow_records)[which(colnames(meow_records)!="scientificName")]
  marine_regs_species_OBIS <- merge(all_records_spec,unique(meow_records[,..col_names]),by=c("taxon","MEOW"),all.x=T)
  marine_regs_species_OBIS <- subset(marine_regs_species_OBIS,Realm=="marine")
  
  if (nrow(marine_regs_species_OBIS)>0){
    marine_regs_species_OBIS$eventDate[is.na(marine_regs_species_OBIS$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
    marine_regspec_fr_OBIS <- aggregate(eventDate ~ MEOW + taxon + Realm,data=marine_regs_species_OBIS,FUN=min) # + Source
    marine_regspec_fr_OBIS$eventDate[marine_regspec_fr_OBIS$eventDate==2500] <- NA # remove dummy variable
    colnames(marine_regspec_fr_OBIS)[colnames(marine_regspec_fr_OBIS)=="MEOW"] <- "location"
  }
  
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "location"
  col_names <- colnames(all_records_spec)[which(colnames(all_records_spec)!="scientificName")]
  terr_regs_species_OBIS <- merge(unique(all_records_spec[,..col_names]),SpecRegionData,by=c("taxon","location"))#,all.y=T
  terr_regs_species_OBIS <- subset(terr_regs_species_OBIS,Realm!="marine")
  
  if (nrow(terr_regs_species_OBIS)>0){
    terr_regs_species_OBIS$eventDate[is.na(terr_regs_species_OBIS$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
    terr_regspec_fr_OBIS <- aggregate(eventDate ~ location + taxon + Realm,data=terr_regs_species_OBIS,FUN=min)# + Source
    terr_regspec_fr_OBIS$eventDate[terr_regspec_fr_OBIS$eventDate==2500] <- NA
  }

  ## combine terrestrial and marine first records #################################
  all_regspec_fr_OBIS <- rbind(marine_regspec_fr_OBIS,terr_regspec_fr_OBIS)
  
  ## add scientificName to OBIS ###################################################
  col_names <- c("scientificName","taxon")
  taxlist <- unique(taxlist[,..col_names])
  all_regspec_fr_OBIS <- merge(all_regspec_fr_OBIS,taxlist,by="taxon",all.x=T)
  # all_regspec_fr_OBIS <- all_regspec_fr_OBIS[,-which(colnames(all_regspec_fr_OBIS)=="taxon")]
  all_regspec_fr_OBIS <- all_regspec_fr_OBIS[,c("location","scientificName","taxon","Realm","eventDate")]

  
  ## combine GBIF and OBIS ########################################################
  all_regspec_fr <- rbind(all_regspec_fr_GBIF,all_regspec_fr_OBIS)
  all_regspec_fr <- unique(all_regspec_fr[,c("location","scientificName","taxon","eventDate")])
    
  ## output ###################################################
  fwrite(all_regspec_fr,file.path("Data","Output",paste0("DASCO_AlienRegions_",file_name_extension,".csv")),sep=";",row.names=F)
  # dat_new <- fread(file.path("Data","Output",paste0("DASCO_AlienRegions_",file_name_extension,".csv")))
  # all_regspec_fr <- readRDS(file.path("Data","FirstRecords_TerrMarRegions_min3.rds"))

  return(all_regspec_fr)
}
