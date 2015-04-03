\name{pubmlst-package}
\alias{pubmlst-package}
\docType{package}
\title{
pubmlst
}
\description{
Infer Sequence Types Based On Allelic Profiles
}
\details{
Utilises the database from PubMLST (http://pubmlst.org) to infer sequence types and species given an allelic profile for Campylobacter jejuni and coli.  Imputes missing alleles if present.
}
\author{
Jonathan Marshall
}
\references{
http://pubmlst.org
}
\keyword{ package }
\examples{
# impute missing alleles, sequence type etc.
df <- data.frame(id=c("A", "B"), ASP=c(2,2), GLN=c(1,4), GLT=c(54,1), 
                 GLY=c(NA, 2), PGM=c(4,2), TKT=c(1,1), UNC=c(5,NA))
impute_mlst_in_data(df)

# take a look at the PubMLST data
head(pubmlst)
}
