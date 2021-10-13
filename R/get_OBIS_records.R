##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# Species occurrence records are downloaded from OBIS.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



get_OBIS_records <- function(path_to_OBISdownloads, file_name_extension,
                      fieldsobis = c("scientificName","scientificNameID",
                                     "decimalLongitude","decimalLatitude", 
                                     "basisOfRecord", "country", "speciesid", 
                                     "marine"),
                      remove.fossils = TRUE,
                      omit.rows = c("speciesid", "decimalLongitude",
                                    "decimalLatitude", "scientificName")
  ){
  
  ### Load the Data
  FullTaxaList <- fread(file.path("Data","Output",paste0("TaxaList_Standardised_",file_name_extension,".csv")))

  SpecList <- sort(unique(FullTaxaList$taxon))

  z <- 0
  out_files <- list()
  for (k in (1:length(SpecList)) ) {#[sample(1:length(SpecList),200)]
    
    if (k%%100==0){
      cat(paste0("\n ",round(k/length(SpecList)*100,2),"% of species (",k,"/",length(SpecList),") processed \n"))
    }
    
    ## download OBIS records
    my_occ <- try(robis::occurrence(scientificname = SpecList[k], fields = fieldsobis,verbose=F),silent=T)#
    
    ## check if downloaded did not work
    if(any(class(my_occ) == "try-error")){
      cat(paste0("  \n OBIS API could not download the species ",SpecList[k],"\n  "))
      next
    }
    
    ## store records if available
    if (nrow(my_occ)>0){ # 
      
      z = z + 1
      my_occ$taxon <- SpecList[k]
      out_files[[z]] <- my_occ
      
      ## generate intermediate output
      save_at <- 200
      if (z%%save_at==0){

        OBIS_Download <- rbindlist(out_files, use.names = TRUE, fill = TRUE)
        filename <- file.path("Data","Output","Intermediate",paste0("OBIS_Data_",z/save_at,".gz"))
        fwrite(OBIS_Download, file = filename )
        
        unlink(file.path("Data","Output","Intermediate",paste0("OBIS_Data_",(z-save_at)/save_at,".gz")))
      }
    }
  }
  OBIS_Download <- rbindlist(out_files, use.names = TRUE, fill = TRUE)
  
  ## remove rows with missing information
  OBIS_Download <- na.omit(OBIS_Download, cols = omit.rows)
  
  ## remove fossil records
  if(remove.fossils == TRUE){
    OBIS_Download <- OBIS_Download[OBIS_Download$basisOfRecord != "FossilSpecimen", ]
  }
  
  OBIS_Download <- unique(OBIS_Download)
  
  # colnames(OBIS_Download)[colnames(OBIS_Download)=="scientificName"] <- "Taxon"
  
  ### Save the occurrence records as one big file
  fwrite(OBIS_Download, file = file.path(path_to_OBISdownloads,paste0("OBIS_CompleteDownload_",file_name_extension,".gz")))
  # OBIS_Download <- fread(file = file.path(path_to_OBISdownloads,paste0("OBIS_CompleteDownload_",file_name_extension,".gz")))
  
  OBIS_Codes_Species <- unique(as.data.frame(cbind(OBIS_Download$scientificName, OBIS_Download$taxon, OBIS_Download$speciesid))) # Hanno: Check if correct column names!!!
  # colnames(OBIS_Codes_Species) <- c("scientificName","Taxon_origName")
  colnames(OBIS_Codes_Species) <- c("taxon","taxon_origName","OBIS_speciesKey")
  
  ### Store the species codes
  fwrite(OBIS_Codes_Species, file = file.path("Data","Output",paste0("OBIS_SpeciesKeys_",file_name_extension,".csv")))
  
  ### Delete intermediate files
  filelist <- list.files(path=file.path("Data","Output","Intermediate"),pattern = "OBIS_Data_*")

  if (file.exists(file.path(path_to_OBISdownloads,paste0("OBIS_CompleteDownload_",file_name_extension,".gz")))){
    file.remove(file.path("Data","Output","Intermediate",filelist))
  }
}

#END