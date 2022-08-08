## DASCO: A workflow to down-scale alien species checklists using occurrence records and to re-allocate species distributions across realms
Authors: Hanno Seebens & Ekin Kaplan

Please cite the following paper when using the workflow or parts of it: 
Seebens, H., & Kaplan, E. (2022). DASCO: A workflow to downscale alien species checklists using occurrence records and to re-allocate species distributions across realms. NeoBiota, 74, 75–91. https://doi.org/10.3897/neobiota.74.81082

The workflow "Downscaling Alien Species Checklists using Occurrence data" (DASCO) has been developed to assign alien species records available on regional checklists (i.e. lists of taxa recorded as being alien in a certain geographic region) to individual populations at local scale. The workflow can therefore be used to down-scale information available in alien species checklists. 

The DASCO workflow takes advantage of the coordinate-based occurrences of  taxa provided by the Global Biodiversity Information Facility (GBIF) and Ocen Biodiversity Information System (OBIS). For each taxon listed in the provided checklist, coordinate-based occurences are obtained from GBIF and OBIS. Obtained data were cleaned and standardised. These local data can then be used for analyses and can be aggregate to any delineation provided by the user. Thus, the workflows allows to re-allocate the occurrence of alien species to deviating spatial representations. In this way, information obtained from checklists can be made available for spatial analyses or regional assessments with different geographic boundaries.

An additional component of the workflow represents the determination of the habitat type (terrestrial, marine, freshwater and brackish) of the taxa. This further allows to assign taxa to regions deviating from the original checklist. For example, many checklists include taxa from any habitat type. By integrating the information obtained from checklists, occurrence portals and habitat sources, the taxa on the checklists can be assigned to marine regions as well. By including a series of checks, the DASCO workflow allows the assignment of taxa to marine ecoregions and terrestrial regions with occurrences confirmed through GBIF and OBIS.

All parts of the workflow are open source and can be modified by the user, but it is strongly recommend to report all modifications of the provided codes and tables together with the versions of the individual databases in the publication of the respective results to ensure transparency and reproducibility.

The workflow has been created in the R software, version 4.0.0 (R Core Team 2019). 
