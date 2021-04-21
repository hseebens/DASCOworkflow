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



add_first_records <- function(
  file_name_extension,
  path_to_GBIFdownloads
  ){
  
  
  SpecRegionData <-  fread(file.path("Data","Output",paste0("FullDataSet_Standardised_",file_name_extension,".gz")),stringsAsFactors = F,header=T)
  taxlist <- fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))
  
  ### GBIF records ##############################################################
  
  ## load data ###############################

  meow_records <- fread(file.path("Data","Output","Intermediate","MarineRecords_GBIF.gz"))

  if (!"eventDate"%in%colnames(SpecRegionData)){
    cat("\n Column 'eventDate' is missing. Cannot add first records. \n")
  }

  all_records_spec <- fread(file.path("Data","Output",paste0("AlienRegions_GBIF_",file_name_extension,".csv")))

  GBIF_keys <- fread(file.path(path_to_GBIFdownloads,paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  GBIF_keys <- unique(GBIF_keys[,c("speciesKey","scientificName")])
  
  SpecRegionData_keys <- merge(SpecRegionData,GBIF_keys,by="scientificName",all.x=T)

  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="Location"] <- "MEOW"
  # if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  col_names <- colnames(meow_records)[which(colnames(meow_records)!="scientificName")]
  marine_regs_species <- merge(all_records_spec,meow_records[,..col_names],by=c("speciesKey","MEOW"),all.x=T)
  marine_regs_species <- subset(marine_regs_species,Realm=="marine")
  marine_regs_species$eventDate[is.na(marine_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
  marine_regspec_fr <- aggregate(eventDate ~ MEOW + scientificName + speciesKey + Realm,data=marine_regs_species,FUN=min) # + Source
  marine_regspec_fr$eventDate[marine_regspec_fr$eventDate==2500] <- NA # remove dummy variable
  
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Location"
  col_names <- colnames(all_records_spec)[which(colnames(all_records_spec)!="scientificName")]
  terr_regs_species <- merge(unique(all_records_spec[,..col_names]),SpecRegionData_keys,by=c("speciesKey","Location"),all.y=T)
  terr_regs_species <- subset(terr_regs_species,Realm!="marine")
  terr_regs_species$eventDate[is.na(terr_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
  terr_regspec_fr <- aggregate(eventDate ~ Location + scientificName + speciesKey + Realm,data=terr_regs_species,FUN=min)# + Source
  terr_regspec_fr$eventDate[terr_regspec_fr$eventDate==2500] <- NA
  
  ## combine terrestrial and marine first records #################################
  colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Location"
  all_regspec_fr_GBIF <- rbind(marine_regspec_fr,terr_regspec_fr)
  all_regspec_fr_GBIF <- all_regspec_fr_GBIF[,c("Location","scientificName","Realm","eventDate")]
  

  
  ### OBIS records ##############################################################
  
  ## load data ###############################
  
  meow_records <- fread(file.path("Data","Output","Intermediate","MarineRecords_OBIS.gz"))

  all_records_spec <- fread(file.path("Data","Output",paste0("AlienRegions_OBIS_",file_name_extension,".csv")))

  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="Location"] <- "MEOW"
  # if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  col_names <- colnames(meow_records)[which(colnames(meow_records)!="scientificName")]
  marine_regs_species <- merge(all_records_spec,unique(meow_records[,..col_names]),by=c("Taxon","MEOW"),all.x=T)
  marine_regs_species <- subset(marine_regs_species,Realm=="marine")
  marine_regs_species$eventDate[is.na(marine_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
  marine_regspec_fr <- aggregate(eventDate ~ MEOW + Taxon + Realm,data=marine_regs_species,FUN=min) # + Source
  marine_regspec_fr$eventDate[marine_regspec_fr$eventDate==2500] <- NA # remove dummy variable
  
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Location"
  col_names <- colnames(all_records_spec)[which(colnames(all_records_spec)!="scientificName")]
  terr_regs_species <- merge(unique(all_records_spec[,..col_names]),SpecRegionData,by=c("Taxon","Location"))#,all.y=T
  terr_regs_species <- subset(terr_regs_species,Realm!="marine")
  if (nrow(terr_regs_species)>0){
    terr_regs_species$eventDate[is.na(terr_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
    terr_regspec_fr <- aggregate(eventDate ~ Location + Taxon + Realm,data=terr_regs_species,FUN=min)# + Source
    terr_regspec_fr$eventDate[terr_regspec_fr$eventDate==2500] <- NA
  }

  ## combine terrestrial and marine first records #################################
  colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Location"
  all_regspec_fr_OBIS <- rbind(marine_regspec_fr,terr_regspec_fr)
  
  ## add scientificName to OBIS ###################################################
  col_names <- c("scientificName","Taxon")
  taxlist <- unique(taxlist[,..col_names])
  all_regspec_fr_OBIS <- merge(all_regspec_fr_OBIS,taxlist,by="Taxon",all.x=T)
  all_regspec_fr_OBIS <- all_regspec_fr_OBIS[,-which(colnames(all_regspec_fr_OBIS)=="Taxon")]
  all_regspec_fr_OBIS <- all_regspec_fr_OBIS[,c("Location","scientificName","Realm","eventDate")]
  
  ## combine GBIF and OBIS ########################################################
  all_regspec_fr <- rbind(all_regspec_fr_GBIF,all_regspec_fr_OBIS)
  all_regspec_fr <- all_regspec_fr[,c("Location","scientificName","Realm","eventDate")]
    
  ## output ###################################################
  fwrite(all_regspec_fr,file.path("Data","Output",paste0("AlienRegions_FinalDB_",file_name_extension,".csv")),sep=";",row.names=F)
  # all_regspec_fr <- read.table(file.path("Data","Output",paste0("AlienRegionsFirstRecords_",file_name_extension,".csv")),sep=";",header=T,stringsAsFactors = F)
  # all_regspec_fr <- readRDS(file.path("Data","FirstRecords_TerrMarRegions_min3.rds"))

  return(all_regspec_fr)
}
