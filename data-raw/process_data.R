library(dplyr)
library(tidyr)

# load in mlst profile data
date_ran <- "20201103"

profiles_file <- file.path("data-raw", paste0("profiles_", date_ran, ".txt"))
if (0) {
  source_url <- "https://pubmlst.org/bigsdb?db=pubmlst_campylobacter_seqdef&page=downloadProfiles&scheme_id=1"
  profiles <- read.table(source_url, header=TRUE, sep="\t")
  write.table(profiles, profiles_file, row.names=FALSE, sep="\t")
}
profiles <- read.table(profiles_file, header=TRUE, sep="\t")

# load in raw isolate data and sum up to determine coli status
# the isolates file can be downloaded from:
# https://pubmlst.org/bigsdb?db=pubmlst_campylobacter_isolates&page=plugin&name=Export
# where you select:
# id, isolate, source, species and then Typing->MLST scheme
# deselect include all fields (we don't need ST and CC)
isolates_file <- file.path("data-raw", paste0("isolates_", date_ran, ".txt"))
isolates <- read.table(isolates_file, header=T, sep="\t", comment.char="", quote='')

cols_iso <- c("aspA", "glnA", "gltA", "glyA", "pgm", "tkt", "uncA")

for (i in cols_iso)
  isolates[,i] <- suppressWarnings(as.numeric(as.character(isolates[,i])))

# cleanup the species name column
isolates$species <- factor(isolates$species)
levels(isolates$species) <- gsub("Campylobacter (.*)", "\\1", levels(isolates$species))
levels(isolates$species) <- gsub(" ", ".", levels(isolates$species))

# join to the profiles table
joined <- profiles %>% left_join(isolates)

# aggregate up by reshaping
types <- joined %>% filter(!is.na(species)) %>% group_by(ST, species) %>% summarise(Count = n()) %>% tidyr::spread(species, Count, fill=0)

# and join back to the original database
profiles <- profiles %>% left_join(types)
profiles <- profiles %>% mutate(Coli = ifelse(coli + jejuni == 0, NA, coli > jejuni))

# cleanup the clonal complex
wch <- grepl("^ST-([0-9]+) complex", profiles$clonal_complex)
profiles$CC <- NA
profiles$CC[wch] <- as.numeric(gsub("^ST-([0-9]+) complex", "\\1", profiles$clonal_complex)[wch])
profiles$clonal_complex <- NULL

# grab the columns we want
pubmlst <- profiles %>% select(ST, ASP=aspA, GLN=glnA, GLT=gltA, GLY=glyA, PGM=pgm, TKT=tkt, UNC=uncA, CC, Coli)
usethis::use_data(pubmlst, overwrite=TRUE)

# map the ST to row number by adding empty rows (more efficient)
per_row <- data.frame(ST=1:max(profiles$ST))
per_row <- per_row %>% left_join(pubmlst) %>% mutate(ST = ifelse(is.na(ASP), NA, ST))
.pubmlst_flat <- as.matrix(per_row)
#devtools::use_data(.pubmlst_flat, internal=TRUE, overwrite=TRUE)

# more efficient again is to use a map
.pubmlst_map <- list()
for (i in 1:7) {
  .pubmlst_map[[i]] <- list()
  for (j in 1:max(.pubmlst_flat[, i+1], na.rm=T)) {
    .pubmlst_map[[i]][[j]] <- which(.pubmlst_flat[, i+1] == j)
  }
}
usethis::use_data(.pubmlst_flat, .pubmlst_map, internal=TRUE, overwrite=TRUE)
