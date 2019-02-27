#' fetch HLA frequency data function
#' 
#' This function fetches HLA frequency data and is called by updateHLAfrequencies
#' @param url A url where the HLA frequency text data are found.
#' @return A list object with $data storing the data matrix, $version.date providing a date class variable for the version, and $error with a message or empty string upon successful completion
#' @details Called by updateHLAfrequencies
#' @export
#' @examples fret <- fetchHLAfrequencies("http://igdawg.org/pubs/HLA_frequencies_by_accession_and_region.txt")
#' dim(fret$data)
#' print(fret$version.date)
#' print(fret$error)
#'

fetchHLAfrequencies <- function (url) {
    fdat.raw <- try(RCurl::getURLContent(url))
    if(class(fdat.raw) == "try-error") 
    {
        err <- paste("There was an error attempting to retrive data from ",url," : ",fdat.raw[[1]],sep="")
        return(fdat <- list(data = NULL,version.date = NULL,error = err))
    } else {
        parseret <- try ( {
            fdat.data <- read.table(textConnection(fdat.raw),sep="\t",header=T)
            fdat.header <- unlist(strsplit(fdat.raw,"\n"))[1:4]
            n <- grep("date:",fdat.header)
            fdat.version.date <- as.Date(unlist(strsplit(fdat.header[n],": "))[2])
        } )
        if(class(parseret) == "try-error") {
            return(fdat <- NULL,version.date = NULL, error = parseret[[1]])
        } else {
            return(fdat <- list(data = fdat.data,version.date = fdat.version.date, error = ""))
        }
    }
}
