#' Check GL Strings syntax for locus names appearing in multiple locus fields
#'
#' This function parses a GL string sufficiently to check if locus names appear in multiple locus fields. This function is used by canonicalize().
#' @param glstring A character string containing valid GL string information on alleles at one or more loci.
#' @return NULL if no locus identifiers are found in multiple locus fields, otherwise returns a locus name that appears in multiple locus fields.
#' @details Any string containing the GL string delimiters ^, +, /, or | and locus names separated from allele names by * will be parsed and analyzed.
#' @export
#'
checkglstring <- function(glstring) {
    loci <- sort(unlist(strsplit(glstring,"\\^")))
    acheck <- list()
	for(locusn in 1:(length(loci))) {
		alleles <- unique(sort(unlist(strsplit(loci[locusn],"[\\|\\+\\/]"))))
        for(acount in 1:length(alleles)) {
            alocus <- unlist(strsplit(alleles[acount],"[\\*]"))[1]
            if(!is.null(acheck[[alocus]]) && acheck[[alocus]] < locusn) {
                return(alocus)
            }
            acheck[[alocus]] <- locusn
        }
    }
    return(NULL)
}
