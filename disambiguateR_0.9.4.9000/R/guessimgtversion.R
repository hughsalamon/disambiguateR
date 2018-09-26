#' IMGT/HLA Version Guessing for Genotype List Disambiguation Using Human Leukocyte Antigen (HLA) Allele Frequency Data
#' 
#' This function finds the column name in the data specified by the \code{allelehistory} parameter that matches the most allele names found in the GL string data specified by the \code{glstrings} parameter. When more than one column (typically IMGT version) matches equivalently, the column name that sorts highest (most recent IMGT version) is returned.
#' @param glstrings A vector of glstrings to be disambiguated.
#' @param allelehistory Should be the name of an IMGT-versioned allele name table object. If not provided and \code{alleleload} is set to \code{"default"}, as long as the load or download and load is successful, this parameter will be given the value \code{gl_HLA_allele_history} which can be created by this package (see \code{updateHLAdata}). Any data frame with the appropriate structure can be passed as this argument, useful for providing updated IMGT/HLA version data. For example, \code{package_HLA_allele_history} is such a data frame that is provided with this package. If the parameter is not provided and \code{alleleload} is set to \code{"noload"}, this parameter will be given the value \code{package_HLA_allele_history}.
#' @param alleleload Can be \code{"default"}, in which case if \code{"allelehistory"} is not defined or does not have a \code{$version} variable, the function will attempt to load current HLA data tables, and offer to download them if they can't be loaded. If set to \code{"noload"}, \code{allelehistory} will be set to  \code{package_HLA_allele_history}. Other values for this parameter will result in an error.
#' @return A string containing the latest version of the IMGT/HLA with the maximal number of allele names matching, NA if no matching alleles found, -1 if there was an error found for the alleleload parameter, -2 if allelehistory could not be loaded.
#' @details In this version, Cw and C locus prefixes are treated equivalently.
#' @export guessimgtversion
#' @examples
#' gls <- c("A*1101/A*11012+A*3303^B*40:02:01+B*44:03:01/B*44:03:02^Cw*03041/Cw*0308/Cw*0309/Cw*0310+Cw*07011/Cw*07012/Cw*0705/Cw*0706|Cw*0305+Cw*07011/Cw*07012/Cw*0706")
#' guessimgtversion(glstrings=gls,allelehistory=package_HLA_allele_history)
#' gls <- c("A*1101/A*11012+A*3303^B*40:02:01+B*44:03:01/B*44:03:02^C*03:04:01:01/C*03:08/C*03:09/C*03:10+C*07:01:01/C*07:01:02/C*07:05/C*07:06|C*03:05+C*07:01:01/C*07:01:02/C*07:06")
#' guessimgtversion(glstrings=gls,allelehistory=package_HLA_allele_history)
#'
guessimgtversion <- function(glstrings,allelehistory=NULL,alleleload="default") {
    if(alleleload == "default") {
        if(is.null(get0("allelehistory")$version)) {
            # need an allele history object
            if(!(try(loadHLAdata()) == "try-error")) {
                # seems good to go
                allelehistory <- gl_HLA_allele_history
            } else if(.uquery() == "y" && !(try(updateHLAdata()) == "try-error")) {
                # downloading and loading
                allelehistory <- gl_HLA_allele_history
            } else {
                # calling download failed, exit gracefully
                cat("guessimgtversion: data tables load failed, attempt to download and load data also failed, Exiting...\n")
                return(-2)
            }
        }
    } else if (alleleload == "noload") {
        allelehistory <- package_HLA_allele_history
    } else {
        cat("guessimgtversion: alleleload must be \"default\" or \"noload\", Exiting...\n")
        return(-1)
    }
    allelehistory <- allelehistory$data
    # Evaluate arguments to this function, return with errors as needed
    # remove HLA- throughout glstrings
    glstrings <- gsub("HLA-","",glstrings)
    # convert all Cw to C in glstrings, deleted alleles, and allele history
    glstrings <- gsub("Cw","C",glstrings)
    allelehistory <- apply(allelehistory,2,function(x) {gsub("Cw","C",x)})
    # split glstrings into alleles
    allalleles <- unique(sort(unlist(sapply(1:length(glstrings),function (i) {unlist(strsplit(glstrings[i],"\\^|\\||\\+|\\/"))}))))
    amatches <- sapply(1:(dim(allelehistory)[2]),function (i) {length(intersect(allalleles,allelehistory[,i]))})
    names(amatches) <- colnames(allelehistory)
    if(max(amatches) > 0) {
        ret <- max(colnames(allelehistory)[amatches==max(amatches)])
        if(substr(ret,1,1) == "X") {
            ret <- substring(ret,2)
        }
        return(ret)
    } else {
        return(NA)
    }
}
.uquery <- function () {
    a <- "z"
    while(tolower(a) != "y" && tolower(a) != "n") {
        a <- readline(prompt="HLA data tables couldn't be loaded. Update HLA data tables and try loading again (requires internet access)? Type y or n. ")
    }
    return(a)
}
