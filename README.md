# pubmlst

The pubmlst package makes it easy to derive the sequence type and species for _Campylobacter jejuni_ and _Campylobacter coli_ where the Multi-locus sequence typing profile (allelic profile) is available, by looking up the data using PubMLST information.

The package encapsulates the data from PubMLST and allows downloading of the latest profiles directly within R.

Given a data.frame (or database) containing 7 columns for the 7 housekeeping genes, the package will allow determination of the
sequence type, clonal complex, and species.

It allows imputation where one or more of the loci have missing alleles, if there is a unique match within the PubMLST database.

## Installation

Pubmlst is not currently available from CRAN, but you can install it from github with:

```R
# install.packages("devtools")
devtools::install_github("jmarshallnz/pubmlst")
```

## Usage

```R
library(pubmlst)

# assemble some data
df <- data.frame(id=c("A", "B"), ASP=c(2,2), GLN=c(1,4), GLT=c(54,1), 
                 GLY=c(NA, 2), PGM=c(4,2), TKT=c(1,1), UNC=c(5,NA))

# impute the MLST sequence type
impute_mlst_in_data(df)
```
