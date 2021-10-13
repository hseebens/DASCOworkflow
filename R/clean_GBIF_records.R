##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# This cleans coordinates obtained from GBIF using the functionality of the 
# package CoordinateCleaner. The 'outlier' test is memory and time consuming, and
# might be switched off.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



clean_GBIF_records <- function(
  path_to_GBIFdownloads,
  file_name_extension,
  thin_records,
  tests_for_cleaning = c("capitals","centroids", "equal","gbif","institutions","outliers","zeros")  # remove 'seas' test from default
  ){

  ## identify files to import (i.e., all files within all sub-directories ending with .rds and with 'GBIFrecords_NUMBER_NUMBER' in name)
  allfiles <- list.files(path_to_GBIFdownloads)
  # GBIF_records_files <- allfiles[grepl("\\.gz",allfiles)]
  GBIF_records_files <- allfiles[grepl("\\.rds",allfiles)]
  
  dat_all <- list()
  for (i in 1:length(GBIF_records_files)){ # 
    
    cat(paste0("\n ",i,": ",GBIF_records_files[i],"\n"))
    
    # load file
    # dat_sub <- fread(file.path(path_to_GBIFdownloads,GBIF_records_files[i]))
    dat_sub <- readRDS(file.path(path_to_GBIFdownloads,GBIF_records_files[i]))
    
    # remove duplicates
    dat_sub <- unique(dat_sub)
    
    # remove non-numeric values
    nonnumeric <- is.na(as.numeric(dat_sub$speciesKey)) | is.na(as.numeric(dat_sub$decimalLatitude)) | is.na(as.numeric(dat_sub$decimalLongitude))
    if (any(nonnumeric)){
      dat_sub <- dat_sub[!nonnumeric,]
    }
    dat_sub$decimalLatitude <- as.numeric(dat_sub$decimalLatitude)
    dat_sub$decimalLongitude <- as.numeric(dat_sub$decimalLongitude)
    
    # remove wrong coordinates
    ind <- (dat_sub$decimalLatitude>90 | dat_sub$decimalLatitude< -90) |  (dat_sub$decimalLongitude>180 | dat_sub$decimalLongitude< -180)
    dat_sub <- dat_sub[!ind,]
    
    # remove empty records
    ind <- is.na(dat_sub$speciesKey) | is.na(dat_sub$decimalLatitude) | is.na(dat_sub$decimalLongitude)
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
        
      tab_rec <- cumsum(table(dat_sub$speciesKey))
      groups <- (ceiling(tab_rec/n_split))
      group_lvl <- unique(groups)
      
      for (j in 1:length(group_lvl)){
        
        cat(paste0("\n ",i,": ",GBIF_records_files[i]," ",j,"/",length(group_lvl),"\n "))
        
        spec_groups <- names(groups)[groups==group_lvl[j]]
        dat_sub_sub <- subset(dat_sub,speciesKey%in%spec_groups)
        
        # thin records by removing duplicated rounded coordinates
        if (thin_records){
          dat_thinned <- list()

          for (k in 1:length(unique(dat_sub_sub$speciesKey))){
            dat_spec <- subset(dat_sub_sub,speciesKey==unique(dat_sub_sub$speciesKey)[k])
            rounded_lat <- round(dat_spec$decimalLatitude,2) # round coordinates for thinning
            rounded_lon <- round(dat_spec$decimalLongitude,2) # round coordinates for thinning
            ind <- duplicated(cbind(rounded_lon,rounded_lat))
            dat_thinned[[k]] <- dat_spec[!ind,]
          }
          dat_sub_sub <- rbindlist(dat_thinned)
        }

        # clean records #######################################################################
        # outlier test of clean_coordinate is memory-consuming for species with many records, which need to be separated
        max_records <- 10^5
        if (any(table(dat_sub_sub$speciesKey)>max_records)){ # check for species with many records
          
          spec_manyrecords <- names(which(table(dat_sub_sub$speciesKey)>max_records)) # species with many records
          dat_lessrecords <- subset(dat_sub_sub,!speciesKey%in%spec_manyrecords)
          
          ## clean records for species with less records together
          counter <- 0
          dat_cleaned_sub <- list()
          if (nrow(dat_lessrecords)>0){
            counter <- counter + 1
            dat_cleaned_sub[[counter]] <- clean_coordinates(dat_lessrecords, 
                                                      lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                                                      value ="clean",
                                                      tests = tests_for_cleaning,
                                                      outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'          
          }
          
          cat(paste0("\n Data split into further pieces for ",length(spec_manyrecords)," species.\n"))
          
          for (l in 1:length(spec_manyrecords)){ # loop over species with many records
            
            # cat(paste0("\n Data split into pieces for species ",spec_manyrecords[l],"! \n"))
            
            dat_manyrecords <- subset(dat_sub_sub,speciesKey%in%spec_manyrecords[l])
            pieces <- c(seq(1,nrow(dat_manyrecords),by=max_records),nrow(dat_manyrecords)) # split data into smaller pieces
            for (m in 2:length(pieces)){
              
              cat(paste0("\n Data split into further pieces for species ",spec_manyrecords[l],": ",m-1,"/",length(pieces)-1,"\n"))
              
              counter <- counter + 1
              ## clean records using subsets of data
              dat_cleaned_sub[[counter]] <- clean_coordinates(dat_manyrecords[pieces[m-1]:pieces[m],], 
                                                        lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                                                        value ="clean",
                                                        tests = tests_for_cleaning,
                                                        outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
            }
          }
          dat_cleaned <- rbindlist(dat_cleaned_sub)
        } else {
          
          dat_cleaned <- clean_coordinates(dat_sub_sub, 
                                           lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                                           value ="clean",
                                           tests = tests_for_cleaning,
                                           outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
          
        }
         
        # intermediate saving of file (just for safety, files can be removed if everything works)
        fwrite(dat_cleaned,file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,"_",j,".gz")))
      }
      
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
                                       lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                                       value ="clean",
                                       tests = tests_for_cleaning,
                                       outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
      
      # intermediate saving of file (just for safety, files can be removed if everything works)
      fwrite(dat_cleaned,file.path("Data","Output","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".gz")))
    }
    
    dat_all[[i]] <- dat_cleaned
  }
  
  # output
  dat_all_df <- rbindlist(dat_all)
  
  fwrite(dat_all_df, file.path("Data","Output",paste0("GBIFrecords_Cleaned_All_",file_name_extension,".gz")))
}

