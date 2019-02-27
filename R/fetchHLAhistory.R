#' fetch HLA history data function
#' 
#' This function fetches HLA allele name history data and is called by updateHLAhistory
#' @param url A url where the HLA history text data are found.
#' @return A list object with $data storing the data matrix, $version.date providing a date class variable for the version, and $error with a message or NULL upon successful completion
#' @details Called by updateHLAhistory
#' @export
#' @examples hret <- fetchHLAhistory("http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist_history.txt")
#' dim(hret$data)
#' print(hret$version)
#' print(hret$error)
#'

fetchHLAhistory <- function (url) {
    hdat.raw <- try(RCurl::getURLContent(url))
    if(class(hdat.raw) == "try-error") 
    {
        err <- paste("There was an error attempting to retrive data from ",url," : ",fdat.raw[[1]],sep="")
        return(hdat <- list(data = NULL,version = NULL, error = err))
    } else {
        parseret <- try ( {
            hdat.data <- read.table(textConnection(hdat.raw),sep=",",header=T)
            hdat.header <- unlist(strsplit(hdat.raw,"\n"))[1:6]
            n <- grep("date:",hdat.header)
            hdat.version.date <- unlist(strsplit(hdat.header[n],": "))[2]
            n <- grep("version:",hdat.header)
            hdat.version <- unlist(strsplit(hdat.header[n],": "))[2]
            hdat.version <- gsub("/","_",hdat.version)
            hdat.version <- paste(hdat.version,hdat.version.date)
        } )
        if(class(parseret) == "try-error") {
            return(hdat <- NULL, error = parseret[[1]])
        } else {
            return(hdat <- list(data = hdat.data,version = hdat.version, error = ""))
        }
    }
}
