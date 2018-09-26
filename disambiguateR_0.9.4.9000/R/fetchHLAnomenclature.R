#' fetch HLA nomenclature data function
#' 
#' This function fetches HLA nomenclature data and is called by updateHLAnomenclature
#' @param url A url where the HLA nomenclature text data are found.
#' @return A list object with $data storing the data matrix, $version.date providing a date class variable for the version, and $error with a message or empty string upon successful completion
#' @details Called by updateHLAnomenclature
#' @export
#' @examples fret <- fetchHLAnomenclature("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt")
#' dim(fret$data)
#' print(fret$version.date)
#' print(fret$error)
#'

fetchHLAnomenclature <- function (url) {
    ddat.raw <- try(RCurl::getURLContent(url))
    if(class(ddat.raw) == "try-error") {
        err <- paste("There was an error attempting to retrive data from ", url," : ",ddat.raw[[1]],sep="")
        return(ddat <- list(data = NULL,version.date = NULL,error = err))
    } else {
        ddat.data <- read.table(textConnection(ddat.raw),sep=";")
        ddat.data <- ddat.data[!(ddat.data[,5]==""),] 
        ddat.header <- unlist(strsplit(ddat.raw,"\n"))[1:6]
        n <- grep("date:",ddat.header)
        ddat.version.date <- as.Date(unlist(strsplit(ddat.header[n],": "))[2])
        return(dat <- list(data = ddat.data,version.date = ddat.version.date, error = ""))
    }
}
