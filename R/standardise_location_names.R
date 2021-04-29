##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# Location names are standardised accoring to the information provided in 
# AllLocations_DASCO.csv, which has to be provided in Data/Input.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################


standardise_location_names <- function(dat,file_name_extension,data_set=NULL){

  # dat <- SpecRegionData$Location
  
  dat <- cbind.data.frame(dat,dat,1:length(dat),stringsAsFactors=F)  
  colnames(dat) <- c("Location","Location_orig","order")

  ## load location table #################################################
  regions <- read.table(file.path("Data","Input","AllLocations_DASCO.csv"),sep=";",stringsAsFactors = F,header=T)
  # regions$keywords <- gsub("\\(","\\\\(",regions$keywords)
  # regions$keywords <- gsub("\\)","\\\\)",regions$keywords)
  regions$keywords <- tolower(regions$keywords) # set all to lower case for matching
  regions$keywords[regions$keywords==""] <- NA
  regions$Location_lower <- tolower(regions$Location) # set all to lower case for matching
  
  ## prepare data set ############################################

  # dat <- read.table(file.path("Output","Intermediate",paste0(inputfiles[i])),header=T,stringsAsFactors = F)
  
  dat_match1 <- dat #unique(dat[,c("Location","Location_orig")]) ## use another dat set for region matching to keep the original names
  # dat_match1$order <- 1:nrow(dat_match1)
  dat_match1$Location <- gsub("\\xa0|\\xc2", " ",dat_match1$Location) # replace weird white space with recognised white space
  dat_match1$Location <- gsub("^\\s+|\\s+$", "",dat_match1$Location) # trim leading and trailing whitespace
  dat_match1$Location <- gsub("  ", " ",dat_match1$Location) # turn two spaces into one
  dat_match1$Location <- gsub(" \\(the\\)", "",dat_match1$Location) # remove " (the)" 
  dat_match1$Location <- tolower(dat_match1$Location) # set all to lower case for matching
  colnames(dat_match1) <- c("Location_lower","Location_orig","order")
  
  ## step 1: match names of 'dat' with region names of 'regions'
  dat_match1 <- merge(dat_match1,regions,by.x="Location_lower",by.y="Location_lower",all.x=T)
  
  ## step 3: match names by using keywords in 'regions
  ind_keys <- which(!is.na(regions$keywords))
  for (j in 1:length(ind_keys)){ # loop over rows with multiple keywords
    if (any(grepl("; ",regions$keywords[ind_keys[j]]))){ # check if multiple keywords provided
      keywords <- unlist(strsplit(regions$keywords[ind_keys[j]],"; "))
    } else {
      keywords <- regions$keywords[ind_keys[j]]
    }
    for (k in 1:length(keywords)){
      # ind_match <- grep(keywords[k],dat_match1$Location_lower) 
      ind_match <- which(dat_match1$Location_lower==keywords[k]) 
      if (length(unique(regions$Location[ind_keys[j]]))>1) cat(paste0("    Warning: ",keywords[k],"match multiple location names. Refine keywords!"))
      dat_match1$Location[ind_match] <- regions$Location[ind_keys[j]]
      dat_match1$countryCode[ind_match]             <- regions$countryCode[ind_keys[j]]
    }
  }
  
  ## final merging of both data sets with standardised region names
  dat_match1 <- dat_match1[order(dat_match1$order),]
  # if (!identical(dat_match1$Taxon_orig,dat$Taxon_orig)) stop("Data sets not sorted equally!")
  dat$Location <- dat_match1$Location
  
  dat_regnames <- merge(dat,regions[,c("locationID","Location")],by="Location",all.x=T)
  dat_regnames <- dat_regnames[order(dat_regnames$order),]
  dat_regnames <- dat_regnames[,c("Location","Location_orig","locationID")]
  
  dat_regnames$Location[is.na(dat_regnames$Location)] <- dat_regnames$Location_orig[is.na(dat_regnames$Location)] 
  
  ## output ###############################################################################
  
  missing <- sort(unique(dat_regnames$Location_orig[is.na(dat_regnames$locationID)]))
  
  if (length(missing)>0){ # export missing country names
    
    cat(paste0("\n ",length(missing)," location names did not match. Check file 'Missing_Locations_*' \n"))
    
    # write.table(missing,file.path("Data","Output",paste0("Missing_Locations_",file_name_extension,".csv")),row.names = F,col.names=F)
    write.table(missing,file.path("Data","Output",paste0("Missing_Locations_",file_name_extension,"_",data_set,".csv")),row.names = F,col.names=F)
  }
  
  # write.table(reg_names,file.path("Output","Translated_LocationNames.csv"),row.names=F)
  return(dat_regnames)
}

