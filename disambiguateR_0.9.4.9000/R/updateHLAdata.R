#' HLA data updating function
#' 
#' This function updates HLA data used for disambiguation.
#' @param force FALSE by default. TRUE to force updating independent of loaded object version records.
#' @param quickcheck A logical, defaults to FALSE. If TRUE, performs a quick check on HLA history data at \code{url} for the most recent version therein and compares to currently downloaded data. Assumes that the table header of the online data gives the most recent IMGT version string at row 1 column 2. \code{force} takes precedence and if TRUE, quickcheck will be ignored. See \code{updateHLAhistory}.
#' @return NULL. See printed messages for information on update successes or failures.
#' @details When upgrading the package, if you do not start a clean R session you may have stale objects \code{gl_HLA_nomenclature}, \code{gl_HLA_frequencies_by_accession}, and \code{gl_HLA_history} in the workspace that were reloaded when you started your session. In such cases, \code{updateHLAdata} may fail. You may choose to start a clean R session without restoring the workspace or simply run \code{updateHLAdata(force=TRUE)}.
#' @export
#' @examples updateHLAdata()
#'

updateHLAdata <- function(force=FALSE, quickcheck=FALSE) {
    nulldat <- list(data=NULL,version.date=as.Date("0000-01-01"),version="0000",error=NULL)
    cat("Attempting to load any downloaded data currently installed with the package.\n")
    retpre <- loadHLAdata()
    if(retpre != 0) {
        cat("Did not load any previously downloaded data. Proceeding with update.\n")
        gl_HLA_deleted_alleles <<- nulldat
        gl_HLA_frequencies_by_accession <<- nulldat
        gl_HLA_allele_history <<- nulldat
    } else {
        cat("Previously downloaded data found. Versions will be checked. Proceeding with update.\n")
    }
    if(!exists("gl_HLA_deleted_alleles")) {
        gl_HLA_deleted_alleles <<- nulldat
    }
    nomret <- updateHLAnomenclature(force=force)
    cat(nomret$msg,"\n")
    if(!exists("gl_HLA_frequencies_by_accession")) {
        gl_HLA_frequencies_by_accession <<- nulldat
    }
    freret <- updateHLAfrequencies(force=force)
    cat(freret$msg,"\n")
    if(!exists("gl_HLA_allele_history")) {
        gl_HLA_allele_history <<- nulldat
    }
    hisret <- updateHLAhistory(force=force,quickcheck=quickcheck)
    cat(hisret$msg,"\n")
    if(nomret$status + freret$status + hisret$status == 0) {
        cat("HLA data updates were completed\n")
        uret <- 0
    } else {
        cat("HLA data updates were not completed. See messages above.\n")
        uret <- 1
    }
    ret <- loadHLAdata()
    if(ret != 0) {
        cat("loadHLAdata() failed, see errors above.\n")
        uret <- uret + 2
    }
    return(uret)
}
