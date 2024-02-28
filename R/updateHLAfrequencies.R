#' HLA allele frequencies data updating function
#' 
#' This function updates HLA frequency data used for disambiguation.
#' @param force Forces updating even if data is older or smaller than the current data.
#' @param current A list object containing the current HLA frequency data and version information.
#' @param url A url where the HLA frequency text data are found.
#' @return A list object with $status 0 when an update is successful and 1 otherwise. See $msg for clarifications.
#' @details Calls fetchHLAfrequencies
#' @export
#' @examples 
#' loadHLAdata()
#' uret <- updateHLAfrequencies()
#' uret$status
#' uret$msg
#'

updateHLAfrequencies <- function (force=FALSE,current=gl_HLA_frequencies_by_accession,url="https://raw.githubusercontent.com/hughsalamon/disambiguateR/master/inst/HLA_frequencies_by_accession_and_region.txt") {
    fredat.fetched <- fetchHLAfrequencies(url)
    if(!(fredat.fetched$error) == "") {
        return(list(status=1,msg=paste("Error fetching HLA allele frequency data. ",fredat.fetched$error,sep="")))
    }
    datadir <- try(paste(path.package("disambiguateR"), "/data/", sep = ""))
    if(class(datadir) == "try-error") {
        return(list(status=1,msg=paste("Error getting package path. ",datadir$message,sep="")))
    } else if(fredat.fetched$version.date > current$version.date) {
        if(!(all(dim(fredat.fetched$data) >= dim(current$data))) && force == FALSE) {
            msg <- "Fetched HLA frequencies data file that was smaller than the current data, aborting update. Use updateHLAfrequencies(force = TRUE) to override."
            msg <- paste(msg,"\n","Current data dimensions: ",dim(current$data)," New data dimensions: ", dim(fredat.fetched$data, sep=""))
            return(list(status=1,msg=msg))
        }  else {
            # update frequencies data
            saveret <- try ( {
                gl_HLA_frequencies_by_accession <- fredat.fetched
                currentHLAfrequencies <- paste(datadir,"HLAfrequencies_",gl_HLA_frequencies_by_accession$version.date,".Rdata",sep="")
                save(gl_HLA_frequencies_by_accession,file=currentHLAfrequencies)
                save(currentHLAfrequencies,file=paste(datadir,"currentHLAfrequencies.info.RData",sep=""))
            } )
            if(class(saveret) == "try-error") {
                return(list(status=1,msg=paste("Error saving HLA frequencies: ",saveret$message,sep="")))
            } else {
                return(list(status=0,msg=paste("Updated HLA frequencies to: ",fredat.fetched$version.date,sep="")))
            }
        }
    } else {
        if(force == FALSE) {
            msg <- paste("Fetched HLA frequencies data with version ",fredat.fetched$version.date,". Current version is ",current$version.date,". HLA frequencies not updated.",sep="")
            return(list(status=1,msg=msg))
        } else if (force == TRUE) {
            # update frequencies data
            saveret <- try ( {
                gl_HLA_frequencies_by_accession <- fredat.fetched
                currentHLAfrequencies <- paste(datadir,"HLAfrequencies_",gl_HLA_frequencies_by_accession$version.date,".Rdata",sep="")
                save(gl_HLA_frequencies_by_accession,file=currentHLAfrequencies)
                save(currentHLAfrequencies,file=paste(datadir,"currentHLAfrequencies.info.RData",sep=""))
            } )
            if(class(saveret) == "try-error") {
                return(list(status=1,msg=paste("Error saving HLA frequencies (force = TRUE): ",saveret$message,sep="")))
            } else {
                return(list(status=0,msg=paste("Forced update of HLA frequencies to: ",fredat.fetched$version.date,sep="")))
            }
        } else {
            msg <- "force must be TRUE or FALSE. Not sure what do with downloaded frequencies data, aborting."
            return(list(status=1,msg=msg))
        }
    }
}
