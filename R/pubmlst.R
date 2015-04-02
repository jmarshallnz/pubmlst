
#' Download the latest PubMLST Campylobacter jejuni and coli sequences.
#' 
#' @param url the PubMLST URL to use.
#' @return A data frame of sequence types, their allelic profile, and clonal complex.
#' @seealso
#' \url{http://pubmlst.org}
#' @examples
#' mlst <- download_latest_sequences()
#' head(mlst)
download_latest_sequences <- function(url = "http://pubmlst.org/data/profiles/campylobacter.txt") {

  mlst <- read.table(url, header=T, sep="\t")

  # Format up the clonal complex column
  wch <- grepl("^ST-([0-9]+) complex", mlst$clonal_complex)
  mlst$CC <- NA
  mlst$CC[wch] <- as.numeric(gsub("^ST-([0-9]+) complex", "\\1", mlst$clonal_complex)[wch])
  mlst$clonal_complex <- NULL

  # Fixup the names
  names(mlst)[2:8] <- c("ASP", "GLN", "GLT", "GLY", "PGM", "TKT", "UNC")

  mlst
}

#' Find the matching profiles in the PubMLST Campylobacter jejuni and coli database.
#' 
#' Imputes missing alleles as appropriate.
#' 
#' @param profile the allelic profile
#' @param pubmlst the database to match against
#' @return A boolean vector of which rows in pubmlst match
#' @examples
#' mlst    <- download_latest_sequences()
#' matches <- find_matching_profiles(c(2, 1, 54, NA, 4, 1, 5), mlst)
#' mlst[matches,]
find_matching_profiles <- function(profile, pubmlst) {
  possibles <- rep(T, nrow(pubmlst))
  for (j in 1:length(profile)) {
    if (!is.na(profile[j])) {
      possibles <- possibles & pubmlst[,j+1] == profile[j]
    }
  }
  possibles
}
