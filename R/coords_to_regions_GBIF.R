##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# Point-wise occurrences of species's populations obtained from GBIF are matched 
# with terrestrial and marine regions to obtain species occurrences per region. 
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################


coords_to_regions_GBIF <- function(
  name_of_shapefile,
  realm_extension=TRUE,
  file_name_extension=file_name_extension
  ){
  
  ### load data ###############################################################
  
  ## get GBIF species keys
  GBIF_specieskeys <- fread(file.path("Data","Output",paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  GBIF_specieskeys <- unique(GBIF_specieskeys[,c("scientificName","speciesKey")])
  # colnames(GBIF_specieskeys) <- c("speciesKey","scientificName","canonicalName","matchType","Orig_name")
  
  ## Taxon list
  SpecNames <-  fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))

  SpecNames <- merge(SpecNames,GBIF_specieskeys,by="scientificName")
  
  # if (realm_extension  # check if realms should be identified
  #     & !file.exists(file.path("Data","Output",paste0("Habitats_",file_name_extension,".csv")))
  #     & !all(c("habitat_marine","habitat_freshwater","habitat_terrestrial")%in%colnames(SpecNames))){ # check if required columns exist, if not, download data from WoRMS
  # 
  #     cat("\n 'realm_extension==TRUE' requires information about habitats from WoRMS.")
  #     cat("\n If not provided in the taxon file, it will be obtained now. \n")
  #     
  #     SpecNames <- get_WoRMS_habitats(SpecNames) # get habitats for species in WoRMS
  #     
  #     fwrite(SpecNames,file.path("Data","Output",paste0("Habitats_",file_name_extension,".csv")))
  # }

  ## load taxon list with habitats if existing
  if (realm_extension){
    habitats <- fread(file.path("Data","Output",paste0("DASCO_TaxonHabitats_",file_name_extension,".csv")))
    SpecNames <- merge(SpecNames,habitats,by="taxon",all.x=T)  
  }
  
  ## Taxon x region database 
  SpecRegionData <-  fread(file.path("Data","Output",paste0("FullDataSet_Standardised_",file_name_extension,".gz")))

  SpecRegionData_keys <- merge(SpecRegionData,GBIF_specieskeys[,c("scientificName","speciesKey")],by="scientificName")  
  uni_spec <- unique(SpecRegionData_keys[,c("scientificName","speciesKey")])

  ## Polygon file of marine and terrestrial regions
  regions <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=name_of_shapefile,stringsAsFactors = F)
  colnames(regions)[colnames(regions)=="Location"] <- "location"
  # regions$Ecoregion[!is.na(regions$Ecoregion)] <- paste(regions$Ecoregion[!is.na(regions$Ecoregion)],"MEOW",sep="_")
  # regions <- regions[is.na(regions$featurecla),] # remove lakes !!!! (no alien distinction available yet)
  
  ## standardise location names 
  newLocNames <- standardise_location_names(regions$location,file_name_extension,data_set="Shapefile")
  if (nrow(newLocNames)!=nrow(regions)){
    stop("\n Standardisation of location names went wrong. Check standardise_location_names.R in coords_to_regions.R \n")
  } 
  regions$location <- newLocNames$location

  all_locations <- fread(file.path("Data","Input","AllLocations_DASCO.csv"))
  colnames(all_locations)[colnames(all_locations)=="Location"] <- "location"
  all_locations <- all_locations[,c("locationID","location")]
  regions <- merge(regions,all_locations,by="location",all.x=T)

  if (realm_extension){
    # regions$Realm <- NA
    # regions$Realm[!is.na(regions$Ecoregion)] <- "marine"
    # regions$Realm[!is.na(regions$Region)] <- "terrestrial"

    ## Terrestrial to marine ecoregion file
    neighbours <- read.table(file.path("Data","Input","RegionsMEOW_NeighbourList.csv"),sep=";",header=T)
    neighbours <- subset(neighbours,Action!="remove")
    
    ## standardise location names 
    newLocNames <- standardise_location_names(neighbours$Region,file_name_extension,data_set="Neighbours")
    if (nrow(newLocNames)!=nrow(neighbours)){
      stop("\n Standardisation of location names went wrong. Check standardise_location_names.R in coords_to_regions.R \n")
    } 
    neighbours$Region <- newLocNames$location
    
    neighbours$MEOW <- paste(neighbours$MEOW,"MEOW",sep="_")
    colnames(neighbours) <- c("location","MEOW","Action") ## ADJUST shapefile AND REMOVE!!!!!!!!!!!!!!!!!!!!
  }
  
  
  ## Identify region of occurrence for each coordinate entry ######################################
  
  ## All taxon-region pairs for the identification of alien populations
  TaxonRegionPairs <- paste(SpecRegionData_keys$speciesKey,SpecRegionData_keys$location,sep="_")
  
  if (realm_extension){
    
    ## add marine ecoregions to taxon-region pairs 
    marine_terr_recs <- merge(neighbours,SpecRegionData_keys,by="location")
    marine_speckeys <- unique(marine_terr_recs[,c("MEOW","speciesKey")])
    if (any(colnames(marine_terr_recs)%in%c("FirstRecord","eventDate"))){
      colnames(marine_terr_recs)[colnames(marine_terr_recs)%in%c("FirstRecord","First_Record")] <- "eventDate"
      meow_records <- unique(marine_terr_recs[,c("MEOW","scientificName","speciesKey","eventDate")])#,"Source"
    } else {
      meow_records <- unique(marine_terr_recs[,c("MEOW","scientificName","speciesKey")])#,"Source"
    }
    fwrite(meow_records,file.path("Data","Output","Intermediate","MarineRecords_GBIF.gz"))
    
    TaxonRegionPairs <- c(TaxonRegionPairs,paste(marine_speckeys$speciesKey,marine_speckeys$MEOW,sep="_"))
  
    ## list species which are clearly non-marine, clearly marine and clearly freshwater ######
    if (all(c("class","phylum")%in%tolower(colnames(SpecNames)))){
      colnames(SpecNames)[tolower(colnames(SpecNames))=="class"] <- "class"
      colnames(SpecNames)[tolower(colnames(SpecNames))=="phylum"] <- "phylum"
      
      non_marine <- subset(SpecNames,
                           class=="Insecta" 
                           | phylum=="Tracheophyta"
                           | phylum=="Anthocerotophyta"
                           | phylum=="Bryophyta"
                           | class=="Arachnida"
                           | class=="Aves"
                           | class=="Amphibia"
                           | class=="Mammalia" # not fully correct, but no marine alien mammal known
                           | habitat_marine=="0")$speciesKey
    } else {
      cat("\n No information of class and phylum is provided. Skip identification of non-marine species based on taxonomic information. \n")
      
      non_marine <- subset(SpecNames,habitat_marine=="0")$speciesKey
    }
    
    marine <- subset(SpecNames,habitat_marine=="1")$speciesKey
    non_marine <- non_marine[!non_marine%in%marine] # avoid overlaps
    
    freshwater <- unique(subset(SpecNames,habitat_freshwater=="1" & habitat_marine=="0" & habitat_terrestrial=="0")$speciesKey)
  }
  
  
  ## Identify in which region coordinates fall ############################
  
  ## check available files in folder 'Intermediate' or 'Output' ##################
  folder <- file.path("Output","Intermediate")
  available_files <- list.files(file.path("Data","Output","Intermediate"))
  available_files <- available_files[grep("GBIFrecords_Cleaned_",available_files)]
  available_files <- available_files[grep(file_name_extension,available_files)]
  
  if (length(available_files)==0){
    folder <- "Output"
    available_files <- list.files(file.path("Data","Output"))
    available_files <- available_files[grep("GBIFrecords_Cleaned_All_",available_files)]
    available_files <- available_files[grep(file_name_extension,available_files)]
  }
  if (length(available_files)==0){
    cat("\n No available files with coordinates found! \n")
  }
  
  nchunks <- length(available_files) # set total number of data files, which contain the GBIF records
  nsteps <- 50
  
  all_counter <- 0
  chunk_out <- chunk_out_coords <- all_out <- all_out_coords <- list()
  for (i in 1:nchunks){ # loop over all chunks of GBIF coordinate data
    
    # Import processed GBIF files
    all_coords <- fread(file.path("Data",folder,available_files[i])) #rdsfile

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
      ptspoly$SpeciesRegion <- paste(ptspoly$speciesKey,ptspoly$location,sep="_")
      ptspoly_alien <- ptspoly[ptspoly$SpeciesRegion%in%TaxonRegionPairs,]
      
      ## export
      coords_mat <- as.data.frame(st_coordinates(ptspoly_alien),stringsAsFactors = F)
      if (realm_extension){ 
        output <- cbind.data.frame(ptspoly_alien$speciesKey,ptspoly_alien$location,ptspoly_alien$Realm,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("speciesKey","location","Realm","Longitude","Latitude")#
      } else {
        output <- cbind.data.frame(ptspoly_alien$speciesKey,ptspoly_alien$location,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("speciesKey","location","Longitude","Latitude")#
      }
      output <- unique(output)
      
      ## remove non-marine species in marine ecoregions and marine species in terrestrial regions
      if (realm_extension){
        output <- subset(output,!(output$speciesKey%in%non_marine & output$Realm=="marine"))
        output <- subset(output,!(output$speciesKey%in%marine & output$Realm=="terrestrial"))
      }
      
      ## remove entries with a very low number of records per region (requires coordinates in the file 'outpout')
      region_records <- as.data.frame(table(output$speciesKey,output$location),stringsAsFactors = F)
      region_records <- subset(region_records,Freq>0 & Freq<3)
      remove_taxreg <- paste0(region_records$Var1,"_",region_records$Var2)
      output <- subset(output,!(paste0(output$speciesKey,"_",output$location)%in%remove_taxreg))
      
      ## identify species with the majority of records in terrestrial realm and remove
      if (realm_extension){ 
        realm_spec <- as.matrix(table(output$speciesKey,output$Realm))
        realm_spec_proc <- round(((realm_spec) / rowSums(realm_spec))*100) # percent records per realm
        marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] > 75]
        output <- subset(output,!(speciesKey%in%marinespec & Realm=="terrestrial")) # remove terrestrial records of marine species
        non_marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] <= 75]
        output <- subset(output,!(speciesKey%in%non_marinespec & Realm=="marine")) # remove marine records of non-marine species
      }
      # if (nrow(output)>0){
      #   if (output$location=="Hawaiian Islands" & output$Realm=="marine") {print(paste(i,j)); stop()}
      # }
      
      if (realm_extension){
        output_noCoords <- unique(output[,c("speciesKey","location","Realm")])
        output_coords <- unique(output[,c("speciesKey","location","Realm","Longitude","Latitude")])
      } else {
        output_noCoords <- unique(output[,c("speciesKey","location")])
        output_coords <- unique(output[,c("speciesKey","location","Longitude","Latitude")])
      }

      # ## test 
      # test_dat <- unique(merge(unique(output[,c("speciesKey","location","Realm")]),firstrecords_GBIF[,c("speciesKey","Species","Class","Order")],by="speciesKey"))
      # graphics.off()
      # x11(width=12,height=12)
      # plot(st_geometry(regions),xlim=c(0,10),ylim=c(50,60))
      # points(output[which(output$speciesKey==3189846 & output$location=="North Sea"),c("Longitude","Latitude")],col="black",pch=16)
      # 
      # subset(firstrecords_GBIF,speciesKey==9809222) #>75, no vascular plants, no birds, no insects
      # tab_realm[which(rownames(tab_realm)=="5277297"),]
      
      ## output ###############
      # saveRDS(output,file.path("Data","Output","Intermediate",paste0("DASCO_",file_name_extension,"_",i,"_",j,".rds")))
      
      chunk_out_coords[[j]] <- output_coords
      chunk_out[[j]]        <- output_noCoords
    }
    chunk_records <- unique(do.call("rbind",chunk_out))
    chunk_records_coords <- unique(do.call("rbind",chunk_out_coords))
    
    all_out[[i]] <- chunk_records
    all_out_coords[[i]] <- chunk_records_coords
    
    ## output ###############
    fwrite(chunk_records,file=file.path("Data","Output","Intermediate",paste0("DASCO_GBIFregions_",file_name_extension,"_",i,".csv")))
    
    fwrite(chunk_records_coords,file=file.path("Data","Output","Intermediate",paste0("DASCO_GBIFCoords_",file_name_extension,"_",i,".gz")))
  }
  all_records <- rbindlist(all_out)
  all_records <- unique(all_records)
  
  all_records_spec <- merge(all_records,uni_spec,by="speciesKey",all.x=T)
  
  ## set realm of freshwater species to freshwater ######################
  all_records_spec$Realm[all_records_spec$speciesKey%in%freshwater] <- "freshwater"
  
  # ## output ###############
  fwrite(all_records_spec,file=file.path("Data","Output",paste0("DASCO_GBIFregions_",file_name_extension,".csv")))
  # all_records_spec <- readRDS(file.path("Data","Output",paste0("DASCO_GBIFregions_",file_name_extension,".rds")))
  
  ## remove intermediate files if previous saving was successful
  if (file.exists(file.path("Data","Output",paste0("DASCO_GBIFregions_",file_name_extension,".csv")))){
    for (i in 1:nchunks){ # loop over all chunks of coordinate data
      for (j in 1:(length(steps)-1)){# 
        file.remove(file.path("Data","Output","Intermediate",paste0("DASCO_GBIFregions_",file_name_extension,"_",i,".csv")))
      }
    }
  }
  
  ## with coordinates
  
  # all_out_coords <- list()
  # for (i in 1:nchunks){ # loop over all chunks of coordinate data
  #   all_out_coords[[i]] <- readRDS(file.path("Data","Output","Intermediate",paste0("DASCO_GBIFCoords_",file_name_extension,"_",i,".rds")))
  # }
  
  all_coords <- rbindlist(all_out_coords)
  all_coords <- unique(all_coords)
  
  # all_coords_spec <- merge(all_coords,uni_spec,by="speciesKey",all.x=T)
  all_coords$Realm[all_coords$speciesKey%in%freshwater] <- "freshwater"

  fwrite(all_coords,file=file.path("Data","Output",paste0("DASCO_GBIFCoords_",file_name_extension,".gz")))
  
  ## remove intermediate files if previous saving was successful
  if (file.exists(file.path("Data","Output",paste0("DASCO_GBIFCoords_",file_name_extension,".gz")))){
    for (i in 1:nchunks){ # loop over all chunks of coordinate data
      # for (j in 1:(length(steps)-1)){# 
        file.remove(file.path("Data","Output","Intermediate",paste0("DASCO_GBIFCoords_",file_name_extension,"_",i,".gz")))
      # }
    }
  }
  
  # return(all_records_spec)
}
