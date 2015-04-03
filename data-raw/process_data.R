library(dplyr)
library(reshape2)

# load in mlst profile data
profiles_file <- "profiles_20150403.txt"
#profiles_file <- "http://pubmlst.org/data/profiles/campylobacter.txt"
profiles <- read.table(profiles_file, header=T, sep="\t")

# load in raw isolate data and sum up to determine coli status
isolates_file <- "isolates_20150219.txt"
isolates <- read.table(isolates_file, header=T, sep="\t", comment.char="")

cols_iso <- c("aspA", "glnA", "gltA", "glyA", "pgm", "tkt", "uncA")

for (i in cols_iso)
  isolates[,i] <- suppressWarnings(as.numeric(as.character(isolates[,i])))

# cleanup the species name column
levels(isolates$species) <- gsub("Campylobacter (.*)", "\\1", levels(isolates$species))
levels(isolates$species) <- gsub(" ", ".", levels(isolates$species))

# join to the profiles table
joined <- profiles %>% left_join(isolates)

# aggregate up by reshaping
types <- joined %>% mutate(value=1) %>% dcast(ST ~ species, sum)

# drop the useless "NA" type
types <- types %>% mutate("NA"=NULL)

# and join back to the original database
profiles <- profiles %>% inner_join(types)
profiles <- profiles %>% mutate(Coli = ifelse(coli + jejuni == 0, NA, coli > jejuni))

# cleanup the clonal complex
wch <- grepl("^ST-([0-9]+) complex", profiles$clonal_complex)
profiles$CC <- NA
profiles$CC[wch] <- as.numeric(gsub("^ST-([0-9]+) complex", "\\1", profiles$clonal_complex)[wch])
profiles$clonal_complex <- NULL

# grab the columns we want
pubmlst <- profiles %>% select(ST, ASP=aspA, GLN=glnA, GLT=gltA, GLY=glyA, PGM=pgm, TKT=tkt, UNC=uncA, CC, Coli)
devtools::use_data(pubmlst) 

# map the ST to row number by adding empty rows (more efficient)
per_row <- data.frame(ST=1:max(profiles$ST))
per_row <- per_row %>% left_join(pubmlst) %>% mutate(ST = ifelse(is.na(ASP), NA, ST))
pubmlst_flat <- as.matrix(per_row)
devtools::use_data(pubmlst_flat, internal=TRUE)

#TODO: Construct the more efficient version for matching stuff up