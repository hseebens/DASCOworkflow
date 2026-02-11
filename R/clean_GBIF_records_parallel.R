##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# This script cleans coordinates obtained from GBIF using the functionality of the 
# package CoordinateCleaner. The 'outlier' test is memory and time consuming, and
# might be switched off.
#
# The DASCO workflow has been published and has to be cited when used:
# Seebens H, Kaplan E (2022) DASCO: A workflow to downscale alien species 
# checklists using occurrence records and to re-allocate species distributions 
# across realms. NeoBiota 74: 75-91. https://doi.org/10.3897/neobiota.74.81082
#
# Authors: Hanno Seebens with support by Ekin Kaplan, 08.01.2026
##################################################################################



clean_GBIF_records <- function(
  path_to_GBIFdownloads,
  file_name_extension,
  thin_records=TRUE,
  tests_for_cleaning = c("centroids", "equal","gbif","institutions","outliers","zeros")  # remove 'seas' and 'capitals' test from default
  ){

  cat(paste0("\n*** Cleaning GBIF records ***\n") ) # notification for the user
  
  ## identify files to import (i.e., all files within all sub-directories ending with .rds and with 'GBIFrecords_NUMBER_NUMBER' in name)
  allfiles <- list.files(path_to_GBIFdownloads)
  # GBIF_records_files <- allfiles[grepl("\\.gz",allfiles)]
  GBIF_records_files <- allfiles[grepl("\\.rds",allfiles)]
  GBIF_records_files <- GBIF_records_files[grepl(file_name_extension, GBIF_records_files)]
  
    
  dat_all <- list()
  for (i in 1:length(GBIF_records_files)){ #
    
    cat(paste0("\n ",i,": ",GBIF_records_files[i],"\n"))
    
    # load file
    # dat_sub <- fread(file.path(path_to_GBIFdownloads,GBIF_records_files[i]))
    dat_sub <- readRDS(file.path(path_to_GBIFdownloads,GBIF_records_files[i]))
    
    # remove duplicates
    dat_sub <- unique(dat_sub)
    
    # remove non-numeric values
    nonnumeric <- is.na(as.numeric(dat_sub$taxonKey)) | is.na(as.numeric(dat_sub$decimalLatitude)) | is.na(as.numeric(dat_sub$decimalLongitude))
    if (any(nonnumeric)){
      dat_sub <- dat_sub[!nonnumeric,]
    }
    dat_sub$decimalLatitude <- as.numeric(dat_sub$decimalLatitude)
    dat_sub$decimalLongitude <- as.numeric(dat_sub$decimalLongitude)
    
    # remove wrong coordinates
    ind <- (dat_sub$decimalLatitude>90 | dat_sub$decimalLatitude< -90) |  (dat_sub$decimalLongitude>180 | dat_sub$decimalLongitude< -180)
    dat_sub <- dat_sub[!ind,]
    
    # remove empty records
    ind <- is.na(dat_sub$taxonKey) | is.na(dat_sub$decimalLatitude) | is.na(dat_sub$decimalLongitude)
    dat_sub <- dat_sub[!ind,]
    
    # remove inprecise coordinates
    ind <- nchar(sub('[0-9]+\\.', '', dat_sub$decimalLatitude))<2
    dat_sub <- dat_sub[!ind,]
    ind <- nchar(sub('[0-9]+\\.', '', dat_sub$decimalLongitude))<2
    dat_sub <- dat_sub[!ind,]
    
    n_split <- 10^5 # number of records per individual chunks (roughly)
    
    if (nrow(dat_sub)>n_split){
      
      cat(paste0("\n Large data set! Split into smaller pieces.\n"))
      
      if (thin_records){
        cat("\n Record thinning is enabled!  \n\n")
      }
        
      tab_rec <- cumsum(table(dat_sub$taxonKey))
      groups <- (ceiling(tab_rec/n_split))
      group_lvl <- unique(groups)

      cores=detectCores()
      cl <- makeCluster(min(c(cores[1]-1,4))) # avoid overloading your computer
      registerDoParallel(cl)
      
      foreach(j=1:length(group_lvl), .packages=c("data.table","CoordinateCleaner"), .errorhandling = "remove") %dopar% {
      # for (j in 1:length(group_lvl)){
        
        cat(paste0("\n ",i,": ",GBIF_records_files[i]," ",j,"/",length(group_lvl),"\n "))
        
        spec_groups <- names(groups)[groups==group_lvl[j]]
        dat_sub_sub <- subset(dat_sub,taxonKey%in%spec_groups)
        
        # thin records by removing duplicated rounded coordinates
        if (thin_records){
          dat_thinned <- list()

          for (k in 1:length(unique(dat_sub_sub$taxonKey))){
            dat_spec <- subset(dat_sub_sub,taxonKey==unique(dat_sub_sub$taxonKey)[k])
            rounded_lat <- round(dat_spec$decimalLatitude,2) # round coordinates for thinning
            rounded_lon <- round(dat_spec$decimalLongitude,2) # round coordinates for thinning
            ind <- duplicated(cbind(rounded_lon,rounded_lat))
            dat_thinned[[k]] <- dat_spec[!ind,] # select first entry (randomly) and remove remaining
          }
          dat_sub_sub <- rbindlist(dat_thinned)
        }

        # clean records #######################################################################
        # outlier test of clean_coordinate is memory-consuming for species with many records, which need to be separated
        max_records <- 10^5
        if (any(table(dat_sub_sub$taxonKey)>max_records)){ # check for species with many records
          
          spec_manyrecords <- names(which(table(dat_sub_sub$taxonKey)>max_records)) # species with many records
          dat_lessrecords <- subset(dat_sub_sub,!taxonKey%in%spec_manyrecords)
          
          ## clean records for species with less records together
          counter <- 0
          dat_cleaned_sub <- list()
          if (nrow(dat_lessrecords)>0){
            counter <- counter + 1
            dat_cleaned_sub[[counter]] <- clean_coordinates(dat_lessrecords, 
                                                      lon = "decimalLongitude", lat = "decimalLatitude", species = "taxonKey", 
                                                      value ="clean",
                                                      tests = tests_for_cleaning,
                                                      outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'          
          }
          
          cat(paste0("\n Data split into further pieces for ",length(spec_manyrecords)," species.\n"))
          
          for (l in 1:length(spec_manyrecords)){ # loop over species with many records
            
            # cat(paste0("\n Data split into pieces for species ",spec_manyrecords[l],"! \n"))
            
            dat_manyrecords <- subset(dat_sub_sub,taxonKey%in%spec_manyrecords[l])
            pieces <- c(seq(1,nrow(dat_manyrecords),by=max_records),nrow(dat_manyrecords)) # split data into smaller pieces
            for (m in 2:length(pieces)){ # m starts with '2'!
              
              cat(paste0("\n Data split into further pieces for species ",spec_manyrecords[l],": ",m-1,"/",length(pieces)-1,"\n"))
              
              counter <- counter + 1
              ## clean records using subsets of data
              dat_cleaned_sub[[counter]] <- clean_coordinates(dat_manyrecords[pieces[m-1]:pieces[m],], 
                                                        lon = "decimalLongitude", lat = "decimalLatitude", species = "taxonKey", 
                                                        value ="clean",
                                                        tests = tests_for_cleaning,
                                                        outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
            }
          }
          dat_cleaned <- rbindlist(dat_cleaned_sub)

        } else {
          
          dat_cleaned <- clean_coordinates(dat_sub_sub, 
                                           lon = "decimalLongitude", lat = "decimalLatitude", species = "taxonKey", 
                                           value ="clean",
                                           tests = tests_for_cleaning,
                                           outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
        }
         
        # intermediate saving of file (just for safety, files can be removed if everything works)
        fwrite(dat_cleaned,file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,"_",j,".gz")))
        
        return()
      }
      ## stop cluster  
      stopCluster(cl)
      
      # collect data store to disk
      dat_sub_all <- list()
      for (j in 1:length(group_lvl)){
        
        dat_cleaned <- fread(file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,"_",j,".gz")))
        dat_sub_all[[j]] <- dat_cleaned
      }
      dat_cleaned <- rbindlist(dat_sub_all)

      # intermediate saving of file (just for safety, files can be removed if everything works)
      fwrite(dat_cleaned,file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".gz")))
      
      ## remove intermediate files if previous saving was successful
      if (file.exists(file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".gz")))){
        for (j in 1:length(group_lvl)){
          file.remove(file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,"_",j,".gz")))
        }
      }
    } else { # for smaller data sets
      
      # clean records
      dat_cleaned <- clean_coordinates(dat_sub, 
                                       lon = "decimalLongitude", lat = "decimalLatitude", species = "taxonKey", 
                                       value ="clean",
                                       tests = tests_for_cleaning,
                                       outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
      
      # intermediate saving of file (just for safety, files can be removed if everything works)
      fwrite(dat_cleaned,file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".gz")))
    }
    
    dat_all[[i]] <- dat_cleaned
  }
  

  ## in case the workflow was restarted at a later for-loop iteration, all data should be loaded again
  if (any(unlist(lapply(dat_all, length))==0)){
    for (i in 1:length(GBIF_records_files)){
      dat_all[[i]] <- fread(file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".gz")))
    }
  }
  
  # output
  dat_all_df <- rbindlist(dat_all)
  
  fwrite(dat_all_df, file.path("Data","Output",paste0("GBIFrecords_Cleaned_All_",file_name_extension,".gz")))
}

