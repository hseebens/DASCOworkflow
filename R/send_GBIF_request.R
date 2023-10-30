##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################




send_GBIF_request <- function(file_name_extension,path_to_GBIFdownloads,n_accounts = 1,user=user,pwd=pwd,email=email){
  
  #######################################################################################
  ### Variables #########################################################################
  
  n_chunks <- 21 # select one without n_accounts if the alternative approache with only one account is used 
  n_chunks <- n_accounts * 3       # number of chunks to divide the species records nearly equally
  # note that GBIF API only allows 3 simultaneous downloads per account
  # so, for example, you need to open 7 accounts to simultaneously request
  # the data in 21 chunks 
  
  
  #######################################################################################
  ## load species list to be downloaded #################################################
  
  SpecNames <-  fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))
  SpecNames <- unique(SpecNames$scientificName)                                      # Get unique species names

  #######################################################################################
  ### get GBIF keys for all species #####################################################
  
  cat("\n Get GBIF keys for taxa \n")
  
  GBIF_speclist <- list()
  x <- 0
  for (i in 1:length(SpecNames)){# loop over all species
    
    specname <- name_backbone(name=SpecNames[i], limit = 10, strict = TRUE)      # overrides the max limit to increase speed
    if (all(colnames(specname)!="species")) next
    
    x <- x + 1
    GBIF_speclist[[x]] <- c(specname$speciesKey,specname$scientificName,specname$canonicalName,specname$matchType,SpecNames[i])
    
    if (x%%1000==0) print(paste(round(x/length(SpecNames),2)*100,"%"))
  }
  GBIF_species <- as.data.frame(do.call("rbind",GBIF_speclist),stringsAsFactors = F)
  colnames(GBIF_species) <- c("speciesKey","scientificName","canonicalName","matchType","Orig_name")

  ## save intermediate output ######
  fwrite(GBIF_species, file.path(path_to_GBIFdownloads,paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  fwrite(GBIF_species, file.path("Data","Output",paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  # GBIF_species <- fread(file.path("Data","Output",paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))

  
  #######################################################################################
  ### Get the number of GBIF records per species ########################################
  ### to split the chunks into roughly equal pieces for download ########################
  
  cat("\n Get number of records per taxon from GBIF \n")
  
  # remove entries with duplicated GBIF keys (e.g. synonyms)
  ind <- !duplicated(GBIF_species$speciesKey)
  GBIF_species <- GBIF_species[ind,]
  
  GBIF_species$nRecords <- 0
  for (i in 1:length(GBIF_species$speciesKey)){
    
    nRecords <- try(occ_count(speciesKey=GBIF_species$speciesKey[i]))
    
    if (class(nRecords)=="try-error") next
    
    GBIF_species$nRecords[i] <- nRecords
    
    if (i%%1000==0) print(paste(round(i/length(GBIF_species$speciesKey),2)*100,"%"))
  }
  
  ## save intermediate output #############################
  fwrite(GBIF_species, file.path("Data","Output","Intermediate",paste0("SpeciesGBIFnRecords_",file_name_extension,".csv")))
  # GBIF_species <- fread(file.path("Data","Output","Intermediate","SpeciesGBIFnRecords.csv"))
  
  
  
  #######################################################################################
  ### Split taxon list into n_chunks of similar size (nRecords) for download ############
  ### (i.e., split sum of records of all species into n_chunks)
  
  
  GBIF_species$cumsum <- 0
  GBIF_species$cumsum[1] <- GBIF_species$nRecords[1]
  GBIF_species$group <- 1
  x <- 1
  for (i in 2:nrow(GBIF_species)){ # loop over all species
    if (GBIF_species$cumsum[i-1] > sum(GBIF_species$nRecords)/ n_chunks){ # if $cumsum is larger than 1/n_chunks fraction of all records, start a new cumulative sum
      x <- x + 1 # counter for the number of groups of species belonging to one chunk
      GBIF_species$cumsum[i] <- GBIF_species$nRecords[i]
    } else { # if not, continue with the former cumulative sum
      GBIF_species$cumsum[i] <- GBIF_species$cumsum[i-1] + GBIF_species$nRecords[i]
    }
    GBIF_species$group[i] <- x # set the group number
  }
  
  #######################################################################################
  ### Prepare the requests for GBIF API #################################################
  ### Two approaches are provided below using either several GBIF accounts or a single
  ### account. The first is less convenient but stable, while the latter is a beta version
  
  ### identify extent of area of interest for downloading records #######################
  regions <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=name_of_shapefile,stringsAsFactors = F)
  # regions <- regions[regions$Location=="Germany",]
  
  ## get spatial extent
  bounding_box <- st_bbox(regions)
  
  ## enlarge bounding box to also cover buffer zones
  bounding_box[1] <- bounding_box[1] -1
  bounding_box[2] <- bounding_box[2] -1
  bounding_box[3] <- bounding_box[3] +1
  bounding_box[4] <- bounding_box[4] +1

  ## define WKT string to define area for GBIF request
  WKT_string <- paste('POLYGON((',
                  bounding_box[1],bounding_box[2],",",
                  bounding_box[3],bounding_box[2],",",
                  bounding_box[3],bounding_box[4],",",
                  bounding_box[1],bounding_box[4],",",
                  bounding_box[1],bounding_box[2],
                  "))",sep=" ")

  ### using various GBIF accounts #######################################################
  user_base <- user
  email_base <- email
  
  counter <- 0
  x <- 1
  file_downloads <- list()
  for (j in unique(GBIF_species$group)) {
    
    counter <- counter + 1                                         # counts the loop
    
    ## gbif account details (note that x is part of user name and email address)
    ## in this case, we opened accounts and emails such as:
    ## (ekinhanno1, ekinhanno1@gmail.com), (ekinhanno2, ekinhanno2@gmail.com) and so on for convenience.
    
    user <- gsub("1",x,user_base)                  # your gbif.org username
    email <- gsub("1",x,email_base)                # your email which you will recieve the download link
    
    if (counter %% 3 == 0){                                        # every time counter can be divided by 3,
      x <- x + 1                                                   # set x + 1 => select new GBIF account below.
    }                                                              # note that GBIF API allows up to
    
    ## send query of each chunk to GBIF #######################################
    
    sub_keys <- subset(GBIF_species,group==j)$speciesKey
    
    ## prepare requests for GBIF download
    file_downloads[[j]] <- occ_download(
      pred_in("taxonKey", sub_keys),
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      pred_within(WKT_string),
      format = "SIMPLE_CSV",
      user=user,pwd=pwd,email=email
    )
  }
  
  save(file_downloads,file=file.path("Data","Output",paste0("GBIF_download_requests_",file_name_extension,".RData")))
  # load(file=file.path("Data","Output",paste0("GBIF_download_requests_",file_name_extension,".RData")))
  
  ### using one GBIF account and different queries #######################################################
  ### the rgbif functions are beta versions and may not work as expected! ################################
  
  ## GBIF account details ##############################################################################
  
  # user <- paste0("ekinhanno",1)                                  # your gbif.org username
  # pwd <- "seebenskaplan1234"                                     # your gbif.org password (set the same password for all accounts for convenience)
  # email <- paste0("ekinhanno",1,"@outlook.com" )                 # your email which you will recieve the download link
  
  queries <- list()
  for (j in unique(GBIF_species$group)) {
    
    ## send query of each chunk to GBIF #######################################
    
    sub_keys <- subset(GBIF_species,group==j)$speciesKey
    
    ## prepare requests for GBIF download (no execution!)
    queries[[j]] <- occ_download_prep(
      pred_in("taxonKey", sub_keys),
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      format = "SIMPLE_CSV",
      user=user,pwd=pwd,email=email
    )
  }
  
  ## execute requests in sequence
  file_downloads <- occ_download_queue(.list = queries, status_ping = 60)
  # file_downloads
  save(file_downloads,file=file.path("Data","Output",paste0("GBIF_download_requests_",file_name_extension,".RData")))
  
}

# dat <- readRDS("/home/hanno/Storage_large/GBIF/SInASdata/GBIFrecords_Cleaned_AllSInAS.rds")
