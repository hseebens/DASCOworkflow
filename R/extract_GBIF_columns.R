##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# The script open downloaded GBIF files and extracts required information.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



extract_GBIF_columns <- function(path_to_GBIFdownloads,file_name_extension){
  
  ############################################################################################
  ## list available zip files (all zip files in the respective will be considered as relevant)
  
  allfiles <- list.files(path_to_GBIFdownloads)
  zippedfiles <- allfiles[grepl("\\.zip",allfiles)]
  csvfiles <- allfiles[grepl("\\.csv",allfiles)]
  
  
  # extract_files <- c("0034400-200221144449610.zip")
  extract_files <- zippedfiles
  
  
  ############################################################################################
  ## decompress zip files (if they are still zipped) and extract required columns ############
  
  for (i in 1:length(extract_files)){#
    
    cat(paste0("\n ",i," Working on ",extract_files[i],"\n"))
    
    data_key <- strsplit(extract_files[i],"-")[[1]][1]
    
    ## unzip files ###################################
    
    if (!file.exists(paste0(path_to_GBIFdownloads,gsub("\\.zip","\\.csv",extract_files[i])))){ # check if file has been unzipped already 
      try(gbif_raw <- decompress_file(path_to_GBIFdownloads,extract_files[i])) # try to unzip
    } 
    if (class(gbif_raw)=="try-error") next
    
    unzipped <- gsub("\\.zip","\\.csv",extract_files[i])
    
    ### process data files #############################################
    
    ## single file...
    # dat <- fread(file=paste0("FirstRecordsSpec/",unzipped),select=c("speciesKey","decimalLatitude","decimalLongitude","basisOfRecord","eventDate","year","dateIdentified"),quote="")
    dat <- fread(file=file.path(path_to_GBIFdownloads,unzipped),select=c("speciesKey","basisOfRecord","decimalLatitude","decimalLongitude"),quote="")
    
    # ## multiple files...
    # dat <- fread(file=paste0("FirstRecordsSpec/",unzipped),select=ind_columns,quote="")
    # colnames(dat) <- colnames(columns) # only required if original data set was split using 'split' in the terminal
    # if (!is.numeric(dat$decimalLatitude) | !is.numeric(dat$decimalLongitude) | !is.integer(dat$speciesKey) | !is.character(dat$basisOfRecord)){
    #   stop("Columns not correctly defined!")
    # }
    
    dat <- dat[basisOfRecord!="FOSSIL_SPECIMEN"]
    
    # dat_sub <- dat[,c("scientificName","decimalLatitude","decimalLongitude")]
    dat_sub <- dat[,c("speciesKey","decimalLatitude","decimalLongitude")]
    
    #######################################################################
    ### initial cleaning of GBIF records ##################################
    
    # remove duplicates
    ind <- duplicated(dat_sub)
    dat_sub <- dat_sub[!ind]
    
    # remove non-numeric values
    nonnumeric <- is.na(as.numeric(dat_sub$speciesKey)) | is.na(as.numeric(dat_sub$decimalLatitude)) | is.na(as.numeric(dat_sub$decimalLongitude))
    if (any(nonnumeric)){
      dat_sub <- dat_sub[!nonnumeric,]
      
      dat_sub$speciesKey <- as.numeric(dat_sub$speciesKey)
      dat_sub$decimalLatitude <- as.numeric(dat_sub$decimalLatitude)
      dat_sub$decimalLongitude <- as.numeric(dat_sub$decimalLongitude)
    }
    
    # remove wrong coordinates
    ind <- (dat_sub$decimalLatitude>90 | dat_sub$decimalLatitude< -90) |  (dat_sub$decimalLongitude>180 | dat_sub$decimalLongitude< -180)
    dat_sub <- dat_sub[!ind,]
    
    # remove empty records
    ind <- is.na(dat_sub$speciesKey) | is.na(dat_sub$decimalLatitude) | is.na(dat_sub$decimalLongitude)
    dat_sub <- dat_sub[!ind,]

    cat(paste0("\n  ",length(unique(dat_sub$speciesKey))," species in ","GBIFrecords_",file_name_extension,"_",data_key,"-",i,".rds \n"))
    
    ## output ###########
    # fwrite(dat_sub,file = file.path(path_to_GBIFdownloads,paste0("GBIFrecords_",file_name_extension,"_",data_key,"-",i,".gz")))
    saveRDS(dat_sub,file = file.path(path_to_GBIFdownloads,paste0("GBIFrecords_",file_name_extension,"_",data_key,"-",i,".rds"))) # stronger compression compared to gz
   
    ## remove unzipped file to save space ##############
    if (file.exists(file.path(path_to_GBIFdownloads,paste0("GBIFrecords_",file_name_extension,"_",data_key,"-",i,".rds")))){
      file.remove(file.path(path_to_GBIFdownloads,unzipped))
    }
  }
}

