% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gl_match_rules}
\alias{gl_match_rules}
\title{Genotype List allele names rules (used for common community inconsistencies in usage of IMGT allele names)}
\format{
A list object with five observations three variables:
\itemize{
  \item V1: version matching rules in R
  \item V2: R definitions of possible strings to append and retry matching to official IMGT allele name strings
  \item rules: R definitions of possible strings to append and retry matching to official IMGT allele name strings
}
}
\source{
\url{https://github.com/hughsalamon/disambiguateR/blob/c75bddc8dea7a4926897f3be0fc9229faa0a2b78/data/gl_match_rules.RData}
}
\description{
A list object containing IMGT version information and rules for appending allele names to attempt matching missed with official allele names
}
