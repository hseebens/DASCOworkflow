##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# The script obtains information about the habitat (freshwater, terrestrial,
# marine) of a species from WoRMS. The mostly comprised marine species.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



get_WoRMS_habitats <- function(dat){
  
  uni_species <- unique(dat$taxon)
  dat$Habitat_marine <- NA
  dat$Habitat_freshwater <- NA
  dat$Habitat_terrestrial <- NA
  
  pieces <- c(seq(1,length(uni_species),by=100),length(uni_species))
  for (i in 2:length(pieces)){# download habitat information in pieces as WoRMS does not allow large download requests
    
    cat(paste0("\n Getting habitat information from WoRMS ",i-1,"/",length(pieces)-1))
    
    taxa_sub <- uni_species[pieces[i-1]:pieces[i]]
    entries <- try(wm_records_names(name = taxa_sub,marine_only=F),silent=T) # WoRMS
    
    ind_entry <- which(unlist(lapply(entries,nrow))>0)
    entries <- entries[ind_entry]
    
    if (length(entries)==0) next
    
    for (j in 1:length(entries)){ # loop over individual species and add information to data set
      
      taxon <- taxa_sub[ind_entry[j]]
      
      entry <- entries[[j]]
      
      # if (class(entry)[1]=="try-error") next  # if database did not provide info, jump to next species
      if (all(is.na(entry$scientificname))) next # empty entries....
      
      ind_spec <- dat$taxon==taxon
      
      if ("isMarine"%in%colnames(entry)) dat$Habitat_marine[ind_spec] <- entry$isMarine[1]
      if ("isFreshwater"%in%colnames(entry)) dat$Habitat_freshwater[ind_spec] <- entry$isFreshwater[1]
      if ("isTerrestrial"%in%colnames(entry)) dat$Habitat_terrestrial[ind_spec] <- entry$isTerrestrial[1]
    }
  }
  return(dat)
}