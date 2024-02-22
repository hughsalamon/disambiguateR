#' HLA history data updating function
#' 
#' This function updates HLA history data used for disambiguation.
#' @param force Forces updating even if data is older or smaller than the current data.
#' @param quickcheck A logical, defaults to FALSE. If TRUE, performs a quick check on HLA history data at \code{url} for the most recent version therein and compares to currently downloaded data. Assumes that the table header of the online data gives the most recent IMGT version string at row 1 column 2. \code{force} takes precedence and if TRUE, \code{quickcheck} will be ignored.
#' @param current A list object containing the current HLA history data and version information.
#' @param url A url where the HLA history text data are found.
#' @return A list object with $status 0 when an update is successful and 1 otherwise. See $msg for clarifications.
#' @details Calls fetchHLAhistory. \code{quickcheck} is useful if download times for HLA history data are inconvenient.
#' @export
#' @examples 
#' #uret <- updateHLAhistory()
#' #uret$status
#' #uret$msg
#'

updateHLAhistory <- function (force=FALSE,quickcheck=FALSE,current=gl_HLA_allele_history,url="http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist_history.txt") {
    if(quickcheck==TRUE && force==FALSE && !is.null(gl_HLA_allele_history$version)) {
    	tmp <- gl_HLA_allele_history$version
        check.header <- read.table(file="http://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Allelelist_history.txt",nrows = 6,comment.char="")
        n <- grep("date:",check.header)
        check.version.date <- unlist(strsplit(check.header[n],": "))[2]
        n <- grep("version:",check.header)
        check.version <- unlist(strsplit(check.header[n],": "))[2]
        check.version <- gsub("/","_",check.version)
        check.version <- paste(check.version,check.version.date)
    	if(check.version <= tmp) {
            msg <- "A quick check suggests that HLA history is not newer than the current data. Not fetching, aborting update. Use updateHLAhistory(force = TRUE) to override."
            return(list(status=1,msg=msg))
        } else {
            msg <- "A quick check suggests that HLA history is newer than the current data. Attemping update..."
        }
    }
    hisdat.fetched <- fetchHLAhistory(url)
    if(!(hisdat.fetched$error) == "") {
        return(list(status=1,msg=paste("Error fetching HLA history data. ",hisdat.fetched$error,sep="")))
    }
    datadir <- try(paste(path.package("disambiguateR"), "/data/", sep = ""))
    if(class(datadir) == "try-error") {
        return(list(status=1,msg=paste("Error getting package path. ",datadir$message,sep="")))
    } else if(hisdat.fetched$version > current$version) {
        if(!(all(dim(hisdat.fetched$data) >= dim(current$data))) && force == FALSE) {
            msg <- "Fetched HLA history that was smaller than the current data, aborting update. Use updateHLAhistory(force = TRUE) to override."
            msg <- paste(msg,"\n","Current data dimensions: ",dim(current$data)," New data dimensions: ", dim(hisdat.fetched$data, sep=""))
            return(list(status=1,msg=msg))
        }  else {
            # update history data
            saveret <- try ( {
                gl_HLA_allele_history <- hisdat.fetched
                currentHLAhistory <- paste(datadir,"HLAhistory_",gl_HLA_allele_history$version,".Rdata",sep="")
                save(gl_HLA_allele_history,file=currentHLAhistory)
                save(currentHLAhistory,file=paste(datadir,"currentHLAhistory.info.RData",sep=""))
            } )
            if(class(saveret) == "try-error") {
                return(list(status=1,msg=paste("Error saving HLA history: ",saveret$message,sep="")))
            } else {
                return(list(status=0,msg=paste("Updated HLA history to: ",hisdat.fetched$version,sep="")))
            }
        }
    } else {
        if(force == FALSE) {
            msg <- paste("Fetched HLA history data with version ",hisdat.fetched$version,". Current version is ",current$version,". HLA history data not updated.",sep="")
            return(list(status=1,msg=msg))
        } else if (force == TRUE) {
            # update history data
            saveret <- try ( {
                gl_HLA_allele_history <- hisdat.fetched
                currentHLAhistory <- paste(datadir,"HLAhistory_",gl_HLA_allele_history$version,".Rdata",sep="")
                save(gl_HLA_allele_history,file=currentHLAhistory)
                save(currentHLAhistory,file=paste(datadir,"currentHLAhistory.info.RData",sep=""))
            } )
            if(class(saveret) == "try-error") {
                return(list(status=1,msg=paste("Error saving HLA history (force = TRUE): ",saveret$message,sep="")))
            } else {
                return(list(status=0,msg=paste("Forced update of HLA history to: ",hisdat.fetched$version,sep="")))
            }
        } else {
            msg <- "force must be TRUE or FALSE. Not sure what do with downloaded history data, aborting."
            return(list(status=1,msg=msg))
        }
    }
}
