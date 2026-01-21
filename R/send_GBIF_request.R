##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The script prepares and sends requests for the GBIF API to download occurrence
# records of provided taxa. It is designed to handle large quantatities of data.
#
# Using multiple GBIF accounts is currently disabled.
#
# Authors: Hanno Seebens with support by Ekin Kaplan, 08.01.2026
##################################################################################


send_GBIF_request <- function(file_name_extension,
                              path_to_GBIFdownloads,
                              n_chunks=3,
                              user=user,
                              pwd=pwd,
                              email=email){
  
  #######################################################################################
  ## Number of requests sent to GBIF API for downloading occurrence data
  # note that GBIF API only allows 3 simultaneous downloads per account
  # n_chunks <- 3       # number of chunks to divide the species records nearly equally, will be checked below

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
    GBIF_speclist[[x]] <- c(specname$usageKey, specname$speciesKey, specname$scientificName, specname$canonicalName, specname$matchType, SpecNames[i])
    
    if (x%%1000==0) print(paste(round(x/length(SpecNames),2)*100,"%"))
  }
  GBIF_species <- as.data.frame(do.call("rbind",GBIF_speclist),stringsAsFactors = F)
  colnames(GBIF_species) <- c("taxonKey","speciesKey","scientificName","canonicalName","matchType","Orig_name") # usageKey renamed to taxonKey! usageKey does only exist for name_backbone() while all other rgbif functions used in DASCO required taxonKey and both seems equal. 
  
  ## save intermediate output ######
  fwrite(GBIF_species, file.path(path_to_GBIFdownloads,paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  fwrite(GBIF_species, file.path("Data","Output",paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))
  # GBIF_species <- fread(file.path("Data","Output",paste0("GBIF_SpeciesKeys_",file_name_extension,".csv")))

  
  #######################################################################################
  ### Get the number of GBIF records per species ########################################
  ### to split the chunks into roughly equal pieces for download ########################
  
  cat("\n Get number of records per taxon from GBIF \n")
  
  # remove entries with duplicated GBIF keys (e.g. synonyms)
  ind <- !duplicated(GBIF_species$taxonKey)
  GBIF_species <- GBIF_species[ind,]
  
  GBIF_species$nRecords <- 0
  for (i in 1:length(GBIF_species$taxonKey)){
    
    nRecords <- try(occ_count(speciesKey=GBIF_species$taxonKey[i])) # just a rough proxy for sample size
    
    if (class(nRecords)=="try-error") next
    
    GBIF_species$nRecords[i] <- nRecords
    
    if (i%%1000==0) print(paste(round(i/length(GBIF_species$speciesKey),2)*100,"%"))
  }
  
  ## save intermediate output #############################
  fwrite(GBIF_species, file.path("Data","Output","Intermediate",paste0("SpeciesGBIFnRecords_",file_name_extension,".csv")))
  # GBIF_species <- fread(file.path("Data","Output","Intermediate","SpeciesGBIFnRecords.csv"))
  
  ## re-assess number of chunks regarding handling in the workflow
  if (sum(GBIF_species$nRecords)/n_chunks > 10^7){
    n_chunks_old <- n_chunks
    n_chunks <- ceiling(sum(GBIF_species$nRecords)/(40*10^6))+1 # workflow can handle files of 40 Mio records (+1 as the loop below does not count last iteration)
    cat(paste0("\n Number of records too high for ", n_chunks_old, " chunks. ", n_chunks-1, " chunks used instead.\n\n"))
  }
  
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
  regions <- st_read(dsn=file.path("Data","Input","Shapefiles"), layer=name_of_shapefile,stringsAsFactors = F)
  # regions <- regions[regions$Location=="Germany",]
  
  ## get spatial extent
  bounding_box <- st_bbox(regions)
  
  ## enlarge bounding box to also cover buffer zones
  bounding_box[1] <- max(c(-180, bounding_box[1] -1))
  bounding_box[2] <- max(c(-90,  bounding_box[2] -1))
  bounding_box[3] <- min(c(180,  bounding_box[3] +1))
  bounding_box[4] <- min(c(90,   bounding_box[4] +1))

  ## define WKT string to define area for GBIF request
  WKT_string <- paste('POLYGON((',
                  bounding_box[1],bounding_box[2],",",
                  bounding_box[3],bounding_box[2],",",
                  bounding_box[3],bounding_box[4],",",
                  bounding_box[1],bounding_box[4],",",
                  bounding_box[1],bounding_box[2],
                  "))",sep=" ")

  # ### using various GBIF accounts #######################################################
  # user_base <- user
  # email_base <- email
  # 
  # counter <- 0
  # x <- 1
  # file_downloads <- list()
  # for (j in unique(GBIF_species$group)) {
  # 
  #   counter <- counter + 1                                         # counts the loop
  # 
  #   ## gbif account details (note that x is part of user name and email address)
  #   ## in this case, we opened accounts and emails such as:
  #   ## (ekinhanno1, ekinhanno1@gmail.com), (ekinhanno2, ekinhanno2@gmail.com) and so on for convenience.
  # 
  #   user <- gsub("1",x,user_base)                  # your gbif.org username
  #   email <- gsub("1",x,email_base)                # your email which you will recieve the download link
  # 
  #   if (counter %% 3 == 0){                                        # every time counter can be divided by 3,
  #     x <- x + 1                                                   # set x + 1 => select new GBIF account below.
  #   }                                                              # note that GBIF API allows up to
  # 
  #   ## send query of each chunk to GBIF #######################################
  # 
  #   sub_keys <- subset(GBIF_species,group==j)$speciesKey
  # 
  #   ## prepare requests for GBIF download
  #   file_downloads[[j]] <- occ_download(
  #     pred_in("taxonKey", sub_keys),
  #     pred("hasCoordinate", TRUE),
  #     pred("hasGeospatialIssue", FALSE),
  #     pred_within(WKT_string),
  #     format = "SIMPLE_CSV",
  #     user=user,pwd=pwd,email=email
  #   )
  # }
  # 
  # save(file_downloads,file=file.path("Data","Output",paste0("GBIF_download_requests_",file_name_extension,".RData")))
  # load(file=file.path("Data","Output",paste0("GBIF_download_requests_",file_name_extension,".RData")))
  
  ### using one GBIF account and different queries #######################################################
  ### the rgbif functions are beta versions and may not work as expected! ################################

  ## GBIF account details ##############################################################################

  queries <- list()
  for (j in unique(GBIF_species$group)) {

    ## send query of each chunk to GBIF #######################################

    sub_keys <- subset(GBIF_species,group==j)$taxonKey

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
