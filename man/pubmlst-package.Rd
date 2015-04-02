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
mlst <- download_latest_sequences()
head(mlst)
}
