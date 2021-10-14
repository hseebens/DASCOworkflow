##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# This cleans coordinates obtained from OBIS using the functionality of the 
# package CoordinateCleaner. The 'outlier' test is memory and time consuming, and
# might be switched off.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



clean_OBIS_records <- function(
  path_to_OBISdownloads,
  file_name_extension,
  thin_records,
  tests_for_cleaning = c("capitals","centroids", "equal","gbif","institutions","outliers","zeros")  # remove 'seas' test from default
  ){
  
  # load file
  dat_sub <- fread(file.path(path_to_OBISdownloads,paste0("OBIS_CompleteDownload_",file_name_extension,".gz")))
  
  # remove duplicates
  dat_sub <- unique(dat_sub[,c("scientificName","decimalLongitude","decimalLatitude","speciesid")])
  
  # remove non-numeric values
  nonnumeric <- is.na(as.numeric(dat_sub$decimalLatitude)) | is.na(as.numeric(dat_sub$decimalLongitude))
  if (any(nonnumeric)){
    dat_sub <- dat_sub[!nonnumeric,]
  }
  
  dat_sub$decimalLatitude <- as.numeric(dat_sub$decimalLatitude)
  dat_sub$decimalLongitude <- as.numeric(dat_sub$decimalLongitude)
  
  # remove wrong coordinates
  ind <- (dat_sub$decimalLatitude>90 | dat_sub$decimalLatitude< -90) |  (dat_sub$decimalLongitude>180 | dat_sub$decimalLongitude< -180)
  dat_sub <- dat_sub[!ind,]
  
  # remove empty records
  ind <- is.na(dat_sub$speciesid) | is.na(dat_sub$decimalLatitude) | is.na(dat_sub$decimalLongitude)
  dat_sub <- dat_sub[!ind,]
  
  # remove inprecise coordinates
  ind <- nchar(sub('[0-9]+\\.', '', dat_sub$decimalLatitude))<2
  dat_sub <- dat_sub[!ind,]
  ind <- nchar(sub('[0-9]+\\.', '', dat_sub$decimalLongitude))<2
  dat_sub <- dat_sub[!ind,]
  
  n_split <- 10^6 # number of records per individual chunks (roughly)
  
  if (nrow(dat_sub)>n_split){
    
    cat(paste0("\n Large data set! Split into smaller pieces.\n"))
    
    tab_rec <- cumsum(table(dat_sub$speciesid))
    groups <- (ceiling(tab_rec/n_split))
    group_lvl <- unique(groups)
    
    dat_sub_all <- list()
    for (j in 1:length(group_lvl)){
      
      cat(paste0("\n Working on chunk ",j,"/",length(group_lvl),"\n "))
      
      spec_groups <- names(groups)[groups==group_lvl[j]]
      dat_sub_sub <- subset(dat_sub,speciesid%in%spec_groups)
      
      # thin records by removing duplicated rounded coordinates
      if (thin_records){
        dat_thinned <- list()
        
        cat("\n Record thinning is enabled!  \n\n")
        
        for (k in 1:length(unique(dat_sub_sub$speciesid))){
          dat_spec <- subset(dat_sub_sub,speciesid==unique(dat_sub_sub$speciesid)[k])
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
      if (any(table(dat_sub_sub$speciesid)>max_records)){ # check for species with many records
        
        spec_manyrecords <- names(which(table(dat_sub_sub$speciesid)>max_records)) # species with many records
        dat_lessrecords <- subset(dat_sub_sub,!speciesid%in%spec_manyrecords)
        
        ## clean records for species with less records together
        counter <- 0
        dat_cleaned_sub <- list()
        if (nrow(dat_lessrecords)>0){
          counter <- counter + 1
          dat_cleaned_sub[[counter]] <- clean_coordinates(dat_lessrecords, 
                                                          lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesid", 
                                                          value ="clean",
                                                          tests = tests_for_cleaning,
                                                          outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'          
        }
        
        cat(paste0("\n Data split into further pieces for ",length(spec_manyrecords)," species.\n"))
        
        for (l in 1:length(spec_manyrecords)){ # loop over species with many records
          
          # cat(paste0("\n Data split into pieces for species ",spec_manyrecords[l],"! \n"))
          
          dat_manyrecords <- subset(dat_sub_sub,speciesid%in%spec_manyrecords[l])
          pieces <- c(seq(1,nrow(dat_manyrecords),by=max_records),nrow(dat_manyrecords)) # split data into smaller pieces
          for (m in 2:length(pieces)){
            
            cat(paste0("\n Data split into further pieces for species ",spec_manyrecords[l],": ",m-1,"/",length(pieces)-1,"\n"))
            
            counter <- counter + 1
            ## clean records using subsets of data
            dat_cleaned_sub[[counter]] <- clean_coordinates(dat_manyrecords[pieces[m-1]:pieces[m],], 
                                                            lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesid", 
                                                            value ="clean",
                                                            tests = tests_for_cleaning,
                                                            outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
          }
        }
        dat_cleaned <- rbindlist(dat_cleaned_sub)
        
      } else {
        
        dat_cleaned <- clean_coordinates(dat_sub_sub, 
                                         lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesid", 
                                         value ="clean",
                                         tests = tests_for_cleaning,
                                         outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
      }
      dat_sub_all[[j]] <- dat_cleaned
    }
    dat_cleaned_all <- do.call("rbind",dat_sub_all)
    
  } else { # for smaller data sets
    
    # clean records
    dat_cleaned_all <- clean_coordinates(dat_sub, 
                                     lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesid", 
                                     value ="clean",
                                     tests = tests_for_cleaning,
                                     outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
  }
  
  fwrite(dat_cleaned_all, file.path("Data","Output",paste0("OBISrecords_Cleaned_All_",file_name_extension,".gz")))
  # dat_cleaned_all <- fread(file.path("Data","Output",paste0("OBISrecords_Cleaned_All_",file_name_extension,".gz")))
}

