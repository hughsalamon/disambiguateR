#' HLA frequencies by accession
#'
#' A list object containing supporting information and a data frame that maps IMGT version 3.9.0 allele names to accession numbers and regional frequency data
#' @format A list object with three items:
#' \itemize{
#'   \item $data: a data frame with 115 observations and 8 variables:
#'   \itemize{
#'     \item region: a factor with 11 levels characterizing geographic regions, "AUS", "EUR", "NAF", "NAM", "NEA", "OCE", "OTH", "SAM", "SEA", "SSA", or "SWA"
#'     \item accession: a factor with 2031 levels charactrizing HLA accession numbers
#'     \item frequency: a numeric giving the frequency of the HLA allele in the geographic region
#'     \item cwdstatus: a factor with three levels, "Neither reg-CWD nor reg-Rare", "reg-CWD", or "reg-Rare"
#'     \item cwdp: a factor with three levels, "", "C", "WD"
#'     \item cwdg: a factor with three levels, "", "C", "WD"
#'     \item mrwaf.allele: a factor with 840 levels
#'     \item X390.allele: a factor with 2031 levels
#'   }
#'   \item $version.date: version date
#'   \item $error: value "". When using updateHLAdata(), a list object of the same format as package_HLA_frequencies_by_accession will be created and named gl_HLA_frequencies_by_accession. In the latter case, errors in creating the object upon attempted download will be recorded in this variable.
#' }
#' @docType data
#' @source \url{http://igdawg.org/pubs/HLA_frequencies_by_accession_and_region.txt}
#' @name package_HLA_frequencies_by_accession
NULL

#' HLA allele history
#'
#' A list object containing supporting information and a data frame that maps multiple IMGT version allele names to accession numbers
#' @format A list object with three items:
#' \itemize{
#'   \item $data: a data frame with 17558 observations and 76 variables:
#'   \itemize{
#'     \item HLA_ID: a factor with 17558 levels characterizing HLA accession numbers
#'     \item All other variables: a factor with levels characterizing allele name strings
#'   }
#'   \item $version.date: version date
#'   \item $error: value "". When using updateHLAdata(), a list object of the same format as package_HLA_allele_history will be created and named gl_HLA_allele_history. In the latter case, errors in creating the object upon attempted download will be recorded in this variable.
#' }
#' @docType data
#' @source \url{ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist_history.txt}
#' @name package_HLA_allele_history
NULL

#' HLA deleted allele records (limited to those where allele names were replaced)
#'
#' A list object containing supporting information and a data frame that documents allele name deprecation and replacements
#' @format A list object with three items:
#' \itemize{
#'   \item $data: a data frame with 115 observations and 6 variables:
#'   \itemize{
#'     \item V1: a factor with 49 levels characterizing HLA loci
#'     \item V2: a factor with 13718 (largely unused) levels characterizing allele names at the specified loci that were replaced
#'     \item V3: an integer with date information characterizing the date the allele in $V2 was assigned
#'     \item V4: an integer with date information characterizing the date the allele in $V5 was assigned (the date of change)
#'     \item V5: a factor with 112 levels characterizing allele names at the specified loci that replaced the names in $V2
#'     \item V6: a factor with 6 levels providing the reason for the replacement
#'   }
#'   \item $version.date: version date
#'   \item $error: value "". When using updateHLAdata(), a list object of the same format as package_HLA_deleted_alleles will be created and named gl_HLA_deleted_alleles. In the latter case, errors in creating the object upon attempted download will be recorded in this variable.
#' }
#' @docType data
#' @source \url{https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt}
#' @name package_HLA_deleted_alleles
NULL

#' Genotype List allele names rules (used for common community inconsistencies in usage of IMGT allele names)
#'
#' A list object containing IMGT version information and rules for appending allele names to attempt matching missed with official allele names
#' @format A list object with five observations three variables:
#' \itemize{
#'   \item V1: version matching rules in R
#'   \item V2: R definitions of possible strings to append and retry matching to official IMGT allele name strings
#'   \item rules: R definitions of possible strings to append and retry matching to official IMGT allele name strings
#' }
#' @docType data
#' @source \url{https://github.com/hughsalamon/disambiguateR/blob/c75bddc8dea7a4926897f3be0fc9229faa0a2b78/data/gl_match_rules.RData}
#' @name gl_match_rules
NULL
