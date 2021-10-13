
graphics.off()
rm(list=ls())


library(data.table) # for clean_GBIF_records, request_GBIF_download

alien_records <- fread(file.path("Data","Output","AlienRegions_FinalDB_SInAS_2.3.2.csv"))

## habitats
habitats <- fread(file.path("Data","Output","DASCO_TaxonHabitats_SInAS_2.3.2.csv"))

## canonical names 
alien_habitats <- merge(alien_records,habitats,by="taxon")


## some stats
table(apply(habitats[,2:5],1,sum,na.rm=T))


time_bin <- 5
start_year <- 1700

## marine species
marine_aliens <- subset(alien_habitats,habitat_marine=="1")
marine_aliens$eventDate <- round(marine_aliens$eventDate/time_bin)*time_bin
marine_ts <- as.data.frame(table(marine_aliens$eventDate),stringsAsFactors = F)
colnames(marine_ts) <- c("Year","nSpec")

## brackish species
brackish_aliens <- subset(alien_habitats,habitat_marine=="0" & habitat_brackish=="1")
brackish_aliens$eventDate <- round(brackish_aliens$eventDate/time_bin)*time_bin
brackish_ts <- as.data.frame(table(brackish_aliens$eventDate),stringsAsFactors = F)
colnames(brackish_ts) <- c("Year","nSpec")


## terrestrial species
terr_aliens <- subset(alien_habitats,habitat_terrestrial=="1" & habitat_freshwater=="0")
terr_aliens$eventDate <- round(terr_aliens$eventDate/time_bin)*time_bin
terr_ts <- as.data.frame(table(terr_aliens$eventDate),stringsAsFactors = F)
colnames(terr_ts) <- c("Year","nSpec")


## freshwater species
fresh_aliens <- subset(alien_habitats,habitat_freshwater=="1" & habitat_terrestrial=="0" )
fresh_aliens$eventDate <- round(fresh_aliens$eventDate/time_bin)*time_bin
fresh_ts <- as.data.frame(table(fresh_aliens$eventDate),stringsAsFactors = F)
colnames(fresh_ts) <- c("Year","nSpec")

x11()
op <- par(pch=16)
# plot(terr_ts,xlim=c(start_year,2020))
# plot(marine_ts,xlim=c(start_year,2020))
# plot(fresh_ts,xlim=c(start_year,2020))
plot(terr_ts,xlim=c(start_year,2020),type="b")
points(marine_ts$Year,marine_ts$nSpec,xlim=c(start_year,2020),col="blue",type="b")
points(fresh_ts$Year,fresh_ts$nSpec,xlim=c(start_year,2020),col=3,type="b")
points(brackish_ts$Year,brackish_ts$nSpec,xlim=c(start_year,2020),col=5,type="b")
abline(v=2005,col="gray")
par(op)

x11()
op <- par(pch=16)
plot(terr_ts$Year,terr_ts$nSpec/max(terr_ts$nSpec),xlim=c(start_year,2020),type="b")
points(marine_ts$Year,marine_ts$nSpec/max(marine_ts$nSpec),col="blue",type="b")
points(fresh_ts$Year,fresh_ts$nSpec/max(fresh_ts$nSpec),col="green",type="b")
abline(v=2005,col="gray")
par(op)