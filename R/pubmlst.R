
#' Download the latest PubMLST Campylobacter jejuni and coli sequences.
#' 
#' @param url the PubMLST URL to use.
#' @return A data frame of sequence types, their allelic profile, and clonal complex.
#' @seealso pubmlst
#' \url{http://pubmlst.org}
#' @examples
#' # not run
#' # mlst <- download_latest_sequences()
#' # head(mlst)
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
#' Profiles with missing alleles return all matching profiles.
#' 
#' @param profile the allelic profile
#' @return A boolean vector of which rows in pubmlst match
#' @examples
#' matches <- find_matching_profiles(c(2, 1, 54, NA, 4, 1, 5))
#' pubmlst[matches,]
find_matching_profiles <- function(profile) {
  for (j in seq_along(profile)) {
    if (!is.na(profile[j])) {
      # out of range
      if (profile[j] > length(.pubmlst_map[[j]]))
        return(numeric(0))

      # first allele found
      possibles <- .pubmlst_map[[ j ]][[ profile[j] ]]
      for (k in seq_len(length(profile) - j) + j) {
        if (!is.na(profile[k])) {
          # out of range
          if (profile[k] > length(.pubmlst_map[[k]]))
            return(numeric(0))

          possibles <- intersect(possibles, .pubmlst_map[[k]][[ profile[k] ]])
        }
      }
      return(possibles)
    }
  }
  return(1:length(.pubmlst_flat))
}

.find_matching_profiles <- function(profile, profiles) {
  possibles <- rep(T, nrow(profiles))
  for (j in 1:length(profile)) {
    if (!is.na(profile[j])) {
      possibles <- possibles & profiles[,j] == profile[j]
    }
  }
  which(possibles)
}

.impute_sequences <- function(sequences) {

  # output
  out_seq <- matrix(NA, nrow(sequences), 7)
  out_imp <- matrix(NA, nrow(sequences), 2) # easier to store separately

  # storage for new MLST sequences we find
  newmlst <- matrix(NA, 0, 7)

  for (i in 1:nrow(sequences)) {
    st <- sequences[i,]

    # find hits in pubmlst
    possibles <- find_matching_profiles(st)

    if (length(possibles) == 0 && !any(is.na(st))) {
      # we may have found this new profile before
      possibles <- .find_matching_profiles(st, newmlst)
      if (length(possibles) == 0) {
        # new ST
        newmlst <- rbind(newmlst, st)
        out_seq[i,] <- st
        out_imp[i,] <- c(nrow(newmlst)+10000, 0)
      } else if (length(possibles) == 1) {
        out_seq[i,] <- newmlst[possibles,]
        out_imp[i,] <- c(possibles, 0)
      } else {
        stop("Bug - can't possibly have duplicate new_mlst matches")
      }
    } else if (length(possibles) == 1) {
      # found - fill in the gaps from PubMLST
      out_seq[i,] <- .pubmlst_flat[possibles,2:8]
      out_imp[i,] <- c(possibles, sum(is.na(st)))
    } else {
      # multiple hits -> leave as unimputed
      out_seq[i,] <- st
    }
  }

  cbind(out_seq, out_imp)
}

#' Impute allelic profiles by matching profiles in the PubMLST Campylobacter jejuni and coli database.
#' 
#' Imputes missing alleles if there is a unique match in PubMLST.
#' 
#' @param mlst the allelic profiles to match as a numeric matrix or data.frame of 7 columns.
#' @return A data frame of 9 columns, representing the imputed profile, the ST and the number of imputed alleles.
#' @examples
#' profiles <- matrix(c(2, 1, 54, NA, 4, 1,  5,
#'                      2, 4,  1,  2, 2, 1, NA), nrow=2, byrow=TRUE)
#' impute_sequences(profiles)
impute_sequences <- function(mlst) {

  # speed up by processing the unique allele combinations
  seqs <- factor(apply(mlst, 1, paste, collapse="_"))

  sequences <- t(sapply(levels(seqs), function(x) { suppressWarnings(as.numeric(unlist(strsplit(x, split="_")))) }))

  # impute the sequences
  output <- .impute_sequences(sequences)

  # and paste into the appropriate rows
  result <- output[as.numeric(seqs),]
  rownames(result) <- 1:nrow(result)
  colnames(result) <- c("ASP", "GLN", "GLT", "GLY", "PGM", "TKT", "UNC", "ST", "Imputed")

  as.data.frame(result)
}

#' Impute MLST allelic profiles by matching profiles in the PubMLST Campylobacter jejuni and coli database.
#' 
#' Imputes missing alleles if there is a unique match in PubMLST, and adds sequence type (ST), clonal complex (CC)
#' and c.coli (Coli, logical) columns to the given database or data.frame.
#' 
#' Alleles marked as "NEW" are interpreted as new alleles, and given a large allele number.
#' 
#' @param db the database or data.frame.
#' @param cols_mlst the columns in the data representing the seven loci.
#' @return db, with additional columns ST, CC, and Coli, where cols_mlst have been imputed.
#' @examples
#' df <- data.frame(id=c("A", "B"), ASP=c(2,2), GLN=c(1,4), GLT=c(54,1), 
#'                  GLY=c(NA, 2), PGM=c(4,2), TKT=c(1,1), UNC=c(5,NA))
#' impute_mlst_in_data(df)
impute_mlst_in_data <- function(db, cols_mlst = c("ASP", "GLN", "GLT", "GLY", "PGM", "TKT", "UNC")) {

  # Convert "NEW" to a special value
  for (i in cols_mlst) {
    newbies <- db[,i] == "NEW"
    if (sum(newbies, na.rm=T) > 0)
      db[newbies,i] <- 1000 - 1:sum(newbies)
  }

  # Convert to numeric
  for (i in cols_mlst)
    db[,i] <- suppressWarnings(as.numeric(db[,i]))

  # Impute sequences
  results <- impute_sequences(db[,cols_mlst])

  # store results
  db[,cols_mlst] <- results[,cols_mlst]
  db$ST          <- results$ST
  db$Imputed     <- results$Imputed
  db$CC          <- .pubmlst_flat[ifelse(results$ST < 10000, results$ST, NA), "CC"]
  db$Coli        <- .pubmlst_flat[ifelse(results$ST < 10000, results$ST, NA), "Coli"]

  # return
  db
}