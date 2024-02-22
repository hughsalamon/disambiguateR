#' HLA nomenclature data updating function
#' 
#' This function updates HLA nomenclature data used for disambiguation.
#' @param force Forces updating even if data is older or smaller than the current data.
#' @param current A list object containing the current HLA nomenclature data and version information.
#' @param url A url where the HLA nomenclature text data are found.
#' @return A list object with $status 0 when an update is successful and 1 otherwise. See $msg for clarifications.
#' @details Calls fetchHLAnomenclature
#' @export
#' @examples uret <- updateHLAnomenclature()
#' uret$status
#' uret$msg
#'

updateHLAnomenclature <- function (force=FALSE,current=gl_HLA_deleted_alleles,url="https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt") {
    nomdat.fetched <- fetchHLAnomenclature(url)
    if(!(nomdat.fetched$error == "")) {
        return(list(status=1,msg=paste("Error fetching HLA nomenclature data. ",nomdat.fetched$error,sep="")))
    }
    datadir <- try(paste(path.package("disambiguateR"), "/data/", sep = ""))
    if(class(datadir) == "try-error") {
        return(list(status=1,msg=paste("Error getting package path. ",datadir$message,sep="")))
    } else if(nomdat.fetched$version.date > current$version.date) {
        if(!(all(dim(nomdat.fetched$data) >= dim(current$data))) && force == FALSE) {
            msg <- "Fetched HLA nomenclature that was smaller than the current data, aborting update. Use updateHLAnomenclature(force = TRUE) to override."
            msg <- paste(msg,"\n","Current data dimensions: ",dim(current$data)," New data dimensions: ", dim(nomdat.fetched$data, sep=""))
            return(list(status=1,msg=msg))
        }  else {
            # update nomenclature data
            saveret <- try ( {
                gl_HLA_deleted_alleles <- nomdat.fetched
                currentHLAnomenclature <- paste(datadir,"HLAnomenclature_",gl_HLA_deleted_alleles$version.date,".Rdata",sep="")
                save(gl_HLA_deleted_alleles,file=currentHLAnomenclature)
                save(currentHLAnomenclature,file=paste(datadir,"currentHLAnomenclature.info.RData",sep=""))
            } )
            if(class(saveret) == "try-error") {
                return(list(status=1,msg=paste("Error saving HLA nomenclature: ",saveret$message,sep="")))
            } else {
                return(list(status=0,msg=paste("Updated HLA nomenclature to: ",nomdat.fetched$version.date,sep="")))
            }
        }
    } else {
        if(force == FALSE) {
            msg <- paste("Fetched HLA nomenclature data with version ",nomdat.fetched$version.date,". Current version is ",current$version.date,". HLA nomenclature not updated.",sep="")
            return(list(status=1,msg=msg))
        } else if (force == TRUE) {
            # update nomenclature data
            saveret <- try ( {
                gl_HLA_deleted_alleles <- nomdat.fetched
                currentHLAnomenclature <- paste(datadir,"HLAnomenclature_",gl_HLA_deleted_alleles$version.date,".Rdata",sep="")
                save(gl_HLA_deleted_alleles,file=currentHLAnomenclature)
                save(currentHLAnomenclature,file=paste(datadir,"currentHLAnomenclature.info.RData",sep=""))
            } )
            if(class(saveret) == "try-error") {
                return(list(status=1,msg=paste("Error saving HLA nomenclature (force = TRUE): ",saveret$message,sep="")))
            } else {
                return(list(status=0,msg=paste("Forced update of HLA nomenclature to: ",nomdat.fetched$version.date,sep="")))
            }
        } else {
            msg <- "force must be TRUE or FALSE. Not sure what do with downloaded nomenclature data, aborting."
            return(list(status=1,msg=msg))
        }
    }
}

