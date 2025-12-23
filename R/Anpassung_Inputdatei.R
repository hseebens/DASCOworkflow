########### Anpassung der ´SInAS-Inputdatei für den DASCO-Workflow #############

# Die zwei Dateien, aus denen die finale Inputdatei erstellt wird, kommen von
# https://zenodo.org/records/17727120 Version 3.1.1
#
# all_taxa: Enthält alle Arten unabhängig vom Invasivitätsstatus (?) mit 
# taxonomischer Einordnung (genus, family...), aber keine Spalte mit location
# species_by_region: Jede Spezies ist den Regionen zugeordnet, in denen sie gebietsfremd ist



# Damit der Workflow funktioniert, wird die Location-Spalte aus species_by_regions
# und die Spalten scientificName, class und phylum aus all_taxa verwendet.
# Beide Dateien besitzen die Spalte taxonID, mithilfe derer gemerged werden kann.
#
# Die Spalte taxon sollte für den merge nicht verwendet werden, auch wenn diese 
# ebenfalls in beiden Dateien vorkommen, da die Gefahr von Synonymen besteht.
# Auch enthält die Spalte NAs, welche im Merge zusammengebracht werden, 
# was keinen Sinn ergibt.

# setwd("C:/Users/karak/Downloads")# bisher habe ich die Dateien im Downloadordner

all_taxa <- read.csv(file.path("Data","Input","SInAS_3.1.1_FullTaxaList.csv"), header = T, sep = " ") 
names(all_taxa)

# Entfernung von Duplikaten in der Spalte taxonID und überprüfen auf NAs
all_taxa <- all_taxa[!duplicated(all_taxa$taxonID), ]
#x <- any(is.na(all_taxa$scientificName)) # Keine NAs

# Nur die notwendigen Spalten werden behalten (muss ggf. noch angepasst werden)
all_taxa <- all_taxa[, c("scientificName", "taxonID", "class", "phylum")]
# all_taxa$scientificName
# Die Spalten aus all_taxa sollen mit species_by_region gemerged werden.
species_by_region <- read.csv(file.path("Data","Input","SInAS_3.1.1.csv"), header = T, sep = " ")
fully_merged <- merge(species_by_region, all_taxa, by = "taxonID", all.x = T)
# NAs in scentific Name mi taxon auffüllen

# Für den Testdurchlauf will ich nur die ersten 100 Zeilen.
# final_test_input <- fully_merged[1:100,]

# Speichern der Datei
# setwd("C:/Users/karak/Documents/Masterarbeit/Orientierung/DASCOworkflow-master")
# write.csv(final_test_input, file.path("Data", "Input","SinAs_3.1.1_mini.csv"))

# fwrite(fully_merged, file.path("Data","Input","SInAS_3.1.1_DASCOinput.csv"))
