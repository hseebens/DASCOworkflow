##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using habitat records from Worms, Fishbase, and Sealifebase.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# This script utilizes three different database to determine in which habitats
# do the species reside in with four categories: Terrestrial, Freshwater, 
# Brackish and Marine.
#
# Authors: Hanno Seebens, Ekin Kaplan, 30.08.2021
##################################################################################



get_habitats_DASCO <- function(file_name_extension,path_to_GBIFdownloads,path_to_OBISdownloads,update_sealifebase= F){

  ## load taxa list 
  GBIF_specieskeys <- fread(file.path("Data","Output",paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  colnames(GBIF_specieskeys) <- c("speciesKey","scientificName","canonicalName","matchType","Orig_name")
  # GBIF_names <- GBIF_specieskeys[,c("canonicalName","speciesKey")]
  # colnames(GBIF_names) <- c("Taxon","speciesKey")
  
  OBIS_specieskeys <- fread(file.path("Data","Output",paste0("OBIS_SpeciesKeys_",file_name_extension,".csv")))
  # colnames(OBIS_specieskeys) <- c("taxon","OBIS_speciesKey") # Hanno: Check if correct names!!!

  taxon_names <- unique(data.frame(c(GBIF_specieskeys$canonicalName,OBIS_specieskeys$taxon))) # check OBIS names!!!
  colnames(taxon_names) <- "taxon"

  ## check for potential habitat information in the provided data set
  taxa <- fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))
  if (any(grepl("habitat",colnames(taxa)))){
    col_names <- c("taxon",colnames(taxa)[grepl("habitat",colnames(taxa))])
    taxa_hab <- unique(taxa[,..col_names])
    taxa_hab <- unique(taxa_hab[habitat!="",])
    ind <- duplicated(taxa_hab$taxon)
    taxon_names <- merge(taxon_names,taxa_hab,all.x=T,by="taxon")
  }

  # taxon_names_all <- taxon_names
  # taxon_names <- taxon_names_all[1:100,]

  ## get habitat information from WoRMS, Fishbase and Sealifebase
  habitats <- get_habitats(taxon_names)
  fwrite(habitats,file.path("Data","Output","Intermediate","Habitats_rawdata.csv"))
  # habitats <- fread(file.path("Data","Output","Intermediate","Habitats_rawdata.csv"))

    
  ## check results (correct obvious errors; combine duplicated entries)
  ## correct habitat information    
  SpecNames <-  fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))
  taxa_dupl <- habitats$taxon[duplicated(habitats$taxon)]
  # set birds and mammals to terrestrial only
  if (any(taxa_dupl %in% subset(SpecNames,class%in%c("Aves","Mammalia","Arachnida"))$taxon)) { # No marine alien mammal is known 
    birds <- taxa_dupl[taxa_dupl %in% subset(SpecNames,class%in%c("Aves","Mammalia","Arachnida"))$taxon] # set birds to terrestrial
    habitats$habitat_freshwater[habitats$taxon %in% subset(SpecNames,class%in%c("Aves","Mammalia","Arachnida"))$taxon] <- 0
    habitats$habitat_brackish[habitats$taxon %in% subset(SpecNames,class%in%c("Aves","Mammalia","Arachnida"))$taxon] <- 0
    habitats$habitat_marine[habitats$taxon %in% subset(SpecNames,class%in%c("Aves","Mammalia","Arachnida"))$taxon] <- 0
    habitats$habitat_terrestrial[habitats$taxon %in% subset(SpecNames,class%in%c("Aves","Mammalia","Arachnida"))$taxon] <- 1
    habitats <- unique(habitats)
  }
  # set Insecta to non-marine
  if (any(taxa_dupl %in% subset(SpecNames,class%in%c("Insecta","Amphibia"))$taxon)) {
    habitats$habitat_marine[habitats$taxon %in% subset(SpecNames,class%in%c("Magnoliopsida","Insecta","Amphibia"))$taxon] <- 0
    habitats <- unique(habitats)
  }
  # set land plants to non-marine
  if (any(taxa_dupl %in% subset(SpecNames,phylum%in%c("Tracheophyta","Anthocerotophyta","Bryophyta"))$taxon)) {
    habitats$habitat_marine[habitats$taxon %in% subset(SpecNames,class%in%c("Magnoliopsida","Insecta","Amphibia"))$taxon] <- 0
    habitats <- unique(habitats)
  }
  
  ## combine habitat information for duplicated taxa
  habitats <- unique(habitats)
  if (any(duplicated(habitats$taxon))){
    print("Warning: Duplicated taxa in habitat list! Check file Duplicated_Entries_Habitats.csv. Habitat information has been combined.")
    
    write.csv2(habitats$taxon[duplicated(habitats$taxon)],file=file.path("Data","Output","Intermediate","Duplicated_Entries_Habitats.csv"))

    ## combine habitat information
    habitats_agg_terr <- aggregate(habitat_terrestrial ~ taxon,data=habitats,FUN=max,na.rm=T)
    habitats_agg_fresh <- aggregate(habitat_freshwater ~ taxon,data=habitats,FUN=max,na.rm=T)
    habitats_agg_brack <- aggregate(habitat_brackish ~ taxon,data=habitats,FUN=max,na.rm=T)
    habitats_agg_marine <- aggregate(habitat_marine ~ taxon,data=habitats,FUN=max,na.rm=T)
    habitats_agg <- merge(habitats_agg_terr,habitats_agg_fresh,by="taxon",all=T)
    habitats_agg <- merge(habitats_agg,habitats_agg_brack,by="taxon",all=T)
    habitats_agg <- merge(habitats_agg,habitats_agg_marine,by="taxon",all=T)
  } 
  
  write.csv2(habitats_agg,file = file.path("Data","Output",paste0("DASCO_TaxonHabitats_",file_name_extension,".csv")) ,row.names = F)

}

