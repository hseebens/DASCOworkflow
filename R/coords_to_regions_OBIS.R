##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# Point-wise occurrences of species's populations obtained from OBIS are matched 
# with terrestrial and marine regions to obtain species occurrences per region. 
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################


coords_to_regions_OBIS <- function(
  name_of_shapefile,
  path_to_OBISdownloads,
  realm_extension,
  file_name_extension
  ){
  
  ### load data ###############################################################
  
  ## get OBIS species keys
  OBIS_specieskeys <- fread(file = file.path(path_to_OBISdownloads,paste0("OBIS_SpeciesKeys_",file_name_extension,".csv")))
  # OBIS_specieskeys$speciesid <- as.numeric(OBIS_specieskeys$speciesid)
  
  
  ## Taxon list
  SpecNames <-  fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))
  
  if (realm_extension  # check if realms should be identified
      & !file.exists(file.path("Data","Output","Intermediate",paste0("Habitats_",file_name_extension)))
      & !all(c("Habitat_marine","Habitat_freshwater","Habitat_terrestrial")%in%colnames(SpecNames))){ # check if required columns exist, if not, download data from WoRMS

      cat("\n 'realm_extension==TRUE' requires information about habitats from WoRMS.")
      cat("\n If not provided in the taxon file, it will be obtained now. \n")
      
      SpecNames <- get_WoRMS_habitats(SpecNames) # get habitats for species in WoRMS
      
      fwrite(SpecNames,file.path("Data","Output","Intermediate",paste0("Habitats_",file_name_extension)))
  }

  ## load taxon list with habitats if existing
  if (realm_extension & file.exists(file.path("Data","Output","Intermediate",paste0("Habitats_",file_name_extension)))){
    SpecNames <- read.table(file.path("Data","Output","Intermediate",paste0("Habitats_",file_name_extension)),stringsAsFactors = F)
  }

  SpecNames <- merge(SpecNames,unique(OBIS_specieskeys[,c("scientificName","Taxon_origName")]),by.x="Taxon",by.y="scientificName")  
  
  ## Taxon x region database 
  SpecRegionData <-  fread(file.path("Data","Output",paste0("FullDataSet_Standardised_",file_name_extension,".gz")))

  # SpecRegionData_keys <- merge(SpecRegionData,OBIS_specieskeys,by.x="Taxon",by.y="Taxon_origName")  
  # uni_spec <- unique(SpecRegionData_keys[,c("scientificName","speciesid")])

  ## Polygon file of marine and terrestrial regions
  # regions2 <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=name_of_shapefile,stringsAsFactors = F)
  regions <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=name_of_shapefile,stringsAsFactors = F)
  # colnames(regions)[2] <- "Ecoregion"
  # regions$Ecoregion[!is.na(regions$Ecoregion)] <- paste(regions$Ecoregion[!is.na(regions$Ecoregion)],"MEOW",sep="_")
  # regions <- regions[is.na(regions$featurecla),] # remove lakes !!!! (no alien distinction available yet)
  
  if (realm_extension){
    # regions$Realm <- NA
    # regions$Realm[!is.na(regions$Ecoregion)] <- "marine"
    # regions$Realm[!is.na(regions$Region)] <- "terrestrial"
    
    # regions$Location <- regions$Region
    # regions$Location[!is.na(regions$Ecoregion)] <- regions$Ecoregion[!is.na(regions$Ecoregion)]
    # regions$Region[!is.na(regions$featurecla)] <- regions$name[!is.na(regions$featurecla)]
    
    ## Terrestrial to marine ecoregion file
    neighbours <- fread(file.path("Data","Input","RegionsMEOW_NeighbourList.csv"))
    neighbours <- subset(neighbours,Action!="remove")
    neighbours$MEOW <- paste(neighbours$MEOW,"MEOW",sep="_")
    colnames(neighbours) <- c("Location","MEOW","Action") ## ADJUST shapefile AND REMOVE!!!!!!!!!!!!!!!!!!!!
  }

  ## Identify region of occurrence for each coordinate entry ######################################
  
  ## All taxon-region pairs for the identification of alien populations
  TaxonRegionPairs <- paste(SpecRegionData$Taxon,SpecRegionData$Location,sep="_")
  
  if (realm_extension){
    
    ## add marine ecoregions to taxon-region pairs 
    marine_terr_recs <- merge(neighbours,SpecRegionData,by="Location",allow.cartesian=T)
    marine_speckeys <- unique(marine_terr_recs[,c("MEOW","Taxon")])
    if (any(colnames(marine_terr_recs)%in%c("FirstRecord","eventDate"))){
      colnames(marine_terr_recs)[colnames(marine_terr_recs)%in%c("FirstRecord")] <- "eventDate"
      meow_records <- unique(marine_terr_recs[,c("MEOW","scientificName","Taxon","eventDate")])#,"Source"
    } else {
      meow_records <- unique(marine_terr_recs[,c("MEOW","scientificName","Taxon")])#,"Source"
    }
    fwrite(meow_records,file.path("Data","Output","Intermediate","MarineRecords_OBIS.gz"))
    
    TaxonRegionPairs <- c(TaxonRegionPairs,paste(marine_speckeys$Taxon,marine_speckeys$MEOW,sep="_"))
  
    ## list species which are clearly non-marine, clearly marine and clearly freshwater ######
    if (all(c("class","phylum")%in%tolower(colnames(SpecNames)))){
      
      non_marine <- subset(SpecNames,
                           class=="Insecta" 
                           | phylum=="Tracheophyta"
                           | phylum=="Anthocerotophyta"
                           | phylum=="Bryophyta"
                           | class=="Arachnida"
                           | class=="Aves"
                           | class=="Amphibia"
                           | class=="Mammalia" # not fully correct, but no marine alien mammal known
                           | Habitat_marine=="0")$Taxon
    } else {
      cat("\n No information of class and phylum is provided. Skip identification of non-marine species based on taxonomic information. \n")
      
      non_marine <- subset(SpecNames,Habitat_marine=="0")$Taxon
    }
    
    marine <- subset(SpecNames,Habitat_marine=="1")$Taxon
    
    freshwater <- unique(subset(SpecNames,Habitat_freshwater=="1" & Habitat_marine=="0" & Habitat_terrestrial=="0")$Taxon)
  }
  
  
  ## Identify in which region coordinates fall ############################
  
  ## check available files in folder 'Intermediate' or 'Output' ##################
  folder <- file.path("Output","Intermediate")
  available_files <- list.files(file.path("Data","Output","Intermediate"))
  available_files <- available_files[grep("OBISrecords_Cleaned_All_",available_files)]
  available_files <- available_files[grep(file_name_extension,available_files)]
  
  if (length(available_files)==0){
    folder <- "Output"
    available_files <- list.files(file.path("Data","Output"))
    available_files <- available_files[grep("OBISrecords_Cleaned_All_",available_files)]
    available_files <- available_files[grep(file_name_extension,available_files)]
  }
  if (length(available_files)==0){
    cat("\n No available files with coordinates found! \n")
  }
  
  nchunks <- length(available_files) # set total number of data files, which contain the GBIF records
  nsteps <- 10
  
  all_counter <- 0
  chunk_out <- all_out <- list()
  for (i in 1:nchunks){ # loop over all chunks of GBIF coordinate data
    
    # Import processed GBIF files
    all_coords <- fread(file.path("Data",folder,available_files[i]))
    colnames(all_coords)[colnames(all_coords)=="scientificName"] <- "Taxon"

    # Transform to sf object
    coords_sf <- st_as_sf(all_coords,coords=c("decimalLongitude","decimalLatitude"),crs=st_crs(regions))

    # Arrange the steps of for loop to avoid large memory consumption
    steps <- ceiling(seq(1,nrow(coords_sf),length.out=nsteps))
    
    # Identify region and alien population
    for (j in 1:(length(steps)-1)){# 
      all_counter <- all_counter + 1
      
      print(paste0(round(all_counter/((nsteps-1)*nchunks)*100,2),"%"))
      
      ## identify region of occurrence
      ptspoly <- st_join(coords_sf[steps[j]:steps[j+1],],regions)
      
      ## identify and keep only alien records
      ptspoly$SpeciesRegion <- paste(ptspoly$Taxon,ptspoly$Location,sep="_")
      ptspoly_alien <- ptspoly[ptspoly$SpeciesRegion%in%TaxonRegionPairs,]
      
      ## export
      coords_mat <- as.data.frame(st_coordinates(ptspoly_alien),stringsAsFactors = F)
      if (realm_extension){ #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        output <- cbind.data.frame(ptspoly_alien$Taxon,ptspoly_alien$Location,ptspoly_alien$Realm,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("Taxon","Location","Realm","Longitude","Latitude")#
      } else { #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        output <- cbind.data.frame(ptspoly_alien$Taxon,ptspoly_alien$Location,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("Taxon","Location","Longitude","Latitude")#
      }
      output <- unique(output)
      
      ## remove non-marine species in marine ecoregions and marine species in terrestrial regions
      if (realm_extension){
        output <- subset(output,!(output$Taxon%in%non_marine & output$Realm=="marine"))
        output <- subset(output,!(output$Taxon%in%marine & output$Realm=="terrestrial"))
      }
      
      ## remove entries with a very low number of records per region (requires coordinates in the file 'outpout')
      region_records <- as.data.frame(table(output$Taxon,output$Location),stringsAsFactors = F)
      region_records <- subset(region_records,Freq>0 & Freq<3)
      remove_taxreg <- paste0(region_records$Var1,"_",region_records$Var2)
      output <- subset(output,!(paste0(output$Taxon,"_",output$Location)%in%remove_taxreg))
      output <- unique(output[,c("Taxon","Location","Realm")])
      
      ## identify species with the majority of records in terrestrial realm and remove
      if (realm_extension){ #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        realm_spec <- as.matrix(table(output$Taxon,output$Realm))
        realm_spec_proc <- round(((realm_spec) / rowSums(realm_spec))*100) # percent records per realm
        marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] > 75]
        output <- subset(output,!(Taxon%in%marinespec & Realm=="terrestrial")) # remove terrestrial records of marine species
        non_marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] <= 75]
        output <- subset(output,!(Taxon%in%non_marinespec & Realm=="marine")) # remove marine records of non-marine species
      }
      
      # ## test 
      # test_dat <- unique(merge(unique(output[,c("Taxon","Location","Realm")]),firstrecords_GBIF[,c("Taxon","Species","Class","Order")],by="Taxon"))
      # graphics.off()
      # x11(width=12,height=12)
      # plot(st_geometry(regions),xlim=c(0,10),ylim=c(50,60))
      # points(output[which(output$Taxon==3189846 & output$Location=="North Sea"),c("Longitude","Latitude")],col="black",pch=16)
      # 
      # subset(firstrecords_GBIF,Taxon==9809222) #>75, no vascular plants, no birds, no insects
      # tab_realm[which(rownames(tab_realm)=="5277297"),]
      
      ## output ###############
      # saveRDS(output,file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,"_",i,"_",j,".rds")))
      
      chunk_out[[j]] <- output
    }
    chunk_records <- unique(do.call("rbind",chunk_out))
    
    all_out[[i]] <- chunk_records
    
    # ## output ###############
    # saveRDS(chunk_records,file.path("Data","Output","Intermediate",paste0("AlienRegions_OBIS_",file_name_extension,"_",i,".rds")))
  }
  all_records_spec <- do.call("rbind",all_out)
  all_records_spec <- unique(all_records_spec)
  
  # all_records_spec <- merge(all_records,uni_spec,by="Taxon",all.x=T)

  ## set realm of freshwater species to freshwater ######################
  all_records_spec$Realm[all_records_spec$Taxon%in%freshwater] <- "freshwater"
  
  # ## output ###############
  fwrite(all_records_spec,file.path("Data","Output",paste0("AlienRegions_OBIS_",file_name_extension,".csv")))
  # all_records_spec <- fread(file.path("Data","Output",paste0("AlienRegions_OBIS_",file_name_extension,".gz")))
  
  # ## remove intermediate files if previous saving was successful
  # if (file.exists(file.path("Data","Output","Intermediate",paste0("AlienRegions_OBIS_",file_name_extension,".rds")))){
  #   for (i in 1:nchunks){ # loop over all chunks of coordinate data
  #     for (j in 1:(length(steps)-1)){# 
  #       file.remove(file.path("Data","Output","Intermediate",paste0("AlienRegions_OBIS_",file_name_extension,"_",i,"_",j,".rds")))
  #     }
  #   }
  # }
  # return(all_records_spec)
}
