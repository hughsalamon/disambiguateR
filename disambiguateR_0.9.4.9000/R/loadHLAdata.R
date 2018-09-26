#' HLA data loading function
#' 
#' This function loads the most current HLA data that has been saved using updateHLAdata()
#' @return 0 when an load is successful 1 otherwise. See printed messages for clarifications.
#' @export
#' @examples loadHLAdata()
#'

loadHLAdata <- function() {
	datadir <- try(paste(path.package("disambiguateR"), "/data/", sep = ""))
	if(class(datadir) == "try-error") {
		cat("Error getting package path. ",datadir$message,"\n")
		return(1)
	} 
	if(length(f <- list.files(datadir))) {
		if(!("currentHLAnomenclature.info.RData" %in% f) & !("currentHLAfrequencies.info.RData" %in% f) & !("currentHLAhistory.info.RData" %in% f)) {
			cat("Looks like there is no information on downloaded data, nothing loaded.\n")
			return(1)
		}
	}
	nomret <- try( {
		load(paste(datadir,"currentHLAnomenclature.info.RData",sep=""),.GlobalEnv)
		load(currentHLAnomenclature,.GlobalEnv)
	} )
	if(class(nomret) == "try-error") {
		cat("There was an error loading nomenclature data used for tracking deleted allele names.","\n")
		return(1)
	}
	freret <- try( {
		load(paste(datadir,"currentHLAfrequencies.info.RData",sep=""),.GlobalEnv)
		load(currentHLAfrequencies,.GlobalEnv)
	} )
	if(class(freret) == "try-error") {
		cat("There was an error loading allele frequency data.","\n")
		return(1)
	}
	hisret <- try( {
		load(paste(datadir,"currentHLAhistory.info.RData",sep=""),.GlobalEnv)
		load(currentHLAhistory,.GlobalEnv)
	} )
	if(class(hisret) == "try-error") {
		cat("There was an error loading allele name history data.","\n")
		return(1)
	}
	return(0)
}
