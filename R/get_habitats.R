##################################################################################
# 
# This script has been developed to obtain habitat information (terrestrial, brackish, marine, freshwater)
# using habitat records from Worms, Fishbase, and Sealifebase, and optionally 
# information provided in a column 'habitat'.
#
# Authors: Hanno Seebens, Ekin Kaplan, 30.08.2021
##################################################################################


get_habitats= function(taxon_names){ # data.frame of taxon names and optionally a column 'habitat'
  
  # Get worms habitats ###################################################
  cat(paste0("\n  Get habitat records from WoRMS \n"))
  
  w.habitats= get_WoRMS_habitats(taxon_names)
  w.habitats= data.frame(w.habitats$taxon,
                         w.habitats$Habitat_marine,
                         w.habitats$Habitat_freshwater,
                         w.habitats$Habitat_terrestrial)
  colnames(w.habitats) = c("taxon","w.Marine","w.Freshwater","w.Terrestrial")
  
  # Get Fishbase habitats ###################################################
  cat(paste0("\n  Get habitat records from FishBase \n"))
  
  fb.habitats= species(species_list = taxon_names$taxon)
  fb.habitats= data.frame(fb.habitats$Species, fb.habitats$Fresh,
                          fb.habitats$Brack,fb.habitats$Saltwater)
  colnames(fb.habitats)= c("taxon", "fbFresh", "fbBrack", "fbSalt")
  
  # Get Sealifebase dataset ###################################################
  cat(paste0("\n  Get habitat records from SeaLifeBase \n"))
  
  sl.habitats= species(server = "sealifebase")
  
  sl.habitats.final= data.frame(sl.habitats$Fresh,
                                sl.habitats$Brack,
                                sl.habitats$Saltwater,
                                sl.habitats$Species) # Hanno: does speciesKey exist in SLB? Using 'Taxon' instead?
  colnames(sl.habitats.final) = c("sbFresh", "sbBrack","sbSalt",
                                  "taxon")
  
  
  ## Merge records from all databases
  all_DBs = merge(taxon_names, sl.habitats.final, by='taxon', all.x=T) # names and Sealifebase
  all_DBs = merge(all_DBs, fb.habitats, by='taxon', all.x=T) # add fishbase
  all_DBs = merge(all_DBs,w.habitats, by='taxon', all.x=T) # add WoRMS
  
  all_DBs$habitat_terrestrial = NA
  all_DBs$habitat_freshwater= NA
  all_DBs$habitat_brackish= NA
  all_DBs$habitat_marine= NA
  
  ## Add sealifebase habitat records ##################
  all_DBs$habitat_freshwater= ifelse(all_DBs$sbFresh == "1", 1, 0)
  all_DBs$habitat_brackish= ifelse(all_DBs$sbBrack == "1", 1, 0)
  all_DBs$habitat_marine= ifelse(all_DBs$sbSalt == "1", 1, 0)
  
  ## Add Fishbase habitat records #####################
  all_DBs$habitat_freshwater[all_DBs$fbFresh=="-1"] <- 1
  all_DBs$habitat_brackish[all_DBs$fbBrack=="-1"] <- 1
  all_DBs$habitat_marine[all_DBs$fbSalt=="-1"] <- 1
  all_DBs$habitat_freshwater[all_DBs$fbFresh=="0" & (all_DBs$habitat_freshwater!=1 | is.na(all_DBs$habitat_freshwater))] <- 0
  all_DBs$habitat_brackish[all_DBs$fbBrack=="0"   & (all_DBs$habitat_brackish!=1   | is.na(all_DBs$habitat_brackish))] <- 0
  all_DBs$habitat_marine[all_DBs$fbSalt=="0"      & (all_DBs$habitat_marine!=1     | is.na(all_DBs$habitat_marine))] <- 0
  
  ## Add Worms records #################################
  all_DBs$habitat_freshwater[all_DBs$w.Freshwater=="1"] <- 1
  all_DBs$habitat_terrestrial[all_DBs$w.Terrestrial=="1"] <- 1
  all_DBs$habitat_marine[all_DBs$w.Marine=="1"] <- 1
  all_DBs$habitat_freshwater[all_DBs$w.Freshwater=="0"   & (all_DBs$habitat_freshwater!=1  | is.na(all_DBs$habitat_freshwater))] <- 0
  all_DBs$habitat_terrestrial[all_DBs$w.Terrestrial=="0" & (all_DBs$habitat_terrestrial!=1 | is.na(all_DBs$habitat_terrestrial))] <- 0
  all_DBs$habitat_marine[all_DBs$w.Marine=="0"           & (all_DBs$habitat_marine!=1      | is.na(all_DBs$habitat_marine))] <- 0
  
  
  #If there is pre-existing habitat information (works best with SINAS workflow)
  if("habitat" %in% colnames(taxon_names)){
    # dat2= data.frame(dat$taxon,dat$habitat)
    # colnames(dat2) = c("taxon", "habitat")
    # all_DBs= merge(all_DBs,dat2, by='taxon', all.x = T)
    # all_DBs= unique(all_DBs)
    
    all_DBs$habitat_freshwater= ifelse(grepl("freshwater",all_DBs$habitat), 1, all_DBs$habitat_freshwater)
    all_DBs$habitat_terrestrial= ifelse(grepl("terrestrial",all_DBs$habitat), 1, all_DBs$habitat_terrestrial)
    all_DBs$habitat_marine= ifelse(grepl("marine",all_DBs$habitat), 1, all_DBs$habitat_marine)
    all_DBs$habitat_brackish= ifelse(grepl("brackish",all_DBs$habitat), 1, all_DBs$habitat_brackish)
  }  
  
  
  ### Prepare output ####
  habitats= data.frame(all_DBs$taxon,all_DBs$habitat_terrestrial,
                       all_DBs$habitat_freshwater, all_DBs$habitat_brackish,
                       all_DBs$habitat_marine)
  colnames(habitats) = c("taxon", "habitat_terrestrial","habitat_freshwater","habitat_brackish",
                         "habitat_marine")
  
  return(unique(habitats))
}
