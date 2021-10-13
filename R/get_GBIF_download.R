##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# The scripts downloads GBIF files, which were requested by 'send_GBIF_request.R'
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################


get_GBIF_download <- function(path_to_GBIFdownloads,file_name_extension,overwrite=FALSE){
  
  print("If the function is not working, please check if GBIF API finished processing.")
  
  unlink(paste0(path_to_GBIFdownloads,"/*"))

  # the loaded file is called 'file_downloads'
  file_name <- list.files(file.path("Data","Output"))
  file_name <- file_name[file_name==paste0("GBIF_download_requests_",file_name_extension,".RData")]
  load(file=file.path("Data","Output",file_name)) # loads file named 'file_downloads'
  
  for (i in 1:length(file_downloads)){
    occ_download_get(file_downloads[[i]],overwrite=overwrite,path=path_to_GBIFdownloads)
  }
  
  ## alternative if previous does not work, execute the following and copy-paste output to command line
  # for (i in 1:length(file_downloads)){
  #   cat(paste0("wget https://api.gbif.org/v1/occurrence/download/request/",as.vector(file_downloads[[i]]),".zip \n"))
  # }
}
  