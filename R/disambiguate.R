#' Genotype List Disambiguation Using Human Leukocyte Antigen (HLA) Allele Frequency Data
#' 
#' This function disambiguates Genotype List (GL) string data based on human leukocyte antigen (HLA) allele frequency data, deleted allele data, geographic region, and IMGT version. Arguements must also be defined to define the behavior of the function with respect to desired output, including complete or partial disambiguation, whether allele names should be standardized to IMGT version data or HLA accession numbers should be output in GL string format.
#' @param glstrings A vector of GL strings to be disambiguated.
#' @param freqbyacc Frequency-by-accession data frame. If not provided and if \code{loadHLAdata()} can be called succesfully, the parameter defaults to \code{gl_HLA_frequencies_by_accession}, an object that can be created using \code{updateHLAdata()}. Also, a data object in this package, \code{package_HLA_frequencies_by_accession}, contains frequency data. Any data frame with the latter's structure can be passed as this argument and should be useful for providing updated or alternative allele-frequency-by-HLA-accession data.
#' @param deletedalleles Deleted alleles data frame. If not provided and if \code{loadHLAdata()} can be called succesfully, the parameter defaults to \code{gl_HLA_deleted_alleles}, an object that can be created using \code{updateHLAdata()}. Also, a data object in this package, \code{package_HLA_deleted_alleles}, contains deleted allele data. Any data frame with this the latter's structure can be passed as this argument, useful for providing updated deleted allele information.
#' @param allelehistory IMGT-versioned allele name table. If not provided and if \code{loadHLAdata()} can be called succesfully, the parameter defaults to \code{gl_HLA_allele_history}, an object that can be created using \code{updateHLAdata()}. Also, a data object in this package, \code{package_HLA_allele_history}, contains IMGT allele history data. Any data frame with the latter's structure can be passed as this argument, useful for providing updated IMGT/HLA version data.
#' @param probratio Probability ratio <=1 and >= 0 used to maintain ambiguity for possible genotypes at a locus with predicted frequency >= \code{probratio} * maximum genotype frequency calculated for that locus. Note that setting \code{dislevel} to 1 will override \code{probratio} < 1.0. Defaults to 0.5.
#' @param regionstring Geographic region designation from the frequency by accession data frame. This is the region for which HLA allele frequency data will be used in genotype frequency calculations. For example, after running \code{updateHLAdata()} successfully, \code{tmp <- disambiguate(region="")} will show the available region designations on the console.
#' @param imgtversion Defaults to \code{"guess"}. This parameter is required to specify the version of IMGT/HLA allele namespace as defined in data passed with parameter \code{allelehistory}. For example, after running \code{updateHLAdata()} successfully, \code{tmp <- disambiguate(region="AUS",probratio=1,imgtversion="")} will show the available IMGT versions on the console. If this parameter is set to \code{"guess"}, \code{guessimgtversion} will be passed the GL string(s) to attempt to determine the IMGT/HLA allele version.
#' @param ruledat Disambiguation rules for matching HLA allele names that do not comply perfectly with IMGT version names, defaults to \code{gl_match_rules}, provided by this package.
#' @param dislevel Disambiguation level, where 1 indicates complete disambiguation to single genotypes at each locus, and 0 permits partial disambiguation.
#' @param mode Disambugation mode, where "normal" returns GL strings with alleles as named in input with only "HLA-" stripped off (the prefix can be replaced using the \code{hlaprefix} parameter), "version" returns GL strings with alleles standardized to the names matching the specified version, and "accession" returns GL strings containing accession numbers.
#' @param hlaprefix Defaulting to FALSE, if TRUE this parameter will force "HLA-" prefix onto all reported allele names, with the exception of accession numbers.
#' @param log Logical to determine if verbose logging of disambiguation calculations is returned with disambiguated GL strings. Defaults to \code{FALSE}.
#' @param allelelog Logical to determine if a table of allele name mappings to accession numbers should be returned for each GL string. Defaults to \code{TRUE}.
#' @param verbose Default value \code{"default"} provides some console feedback on each disambiguated GL string, \code{"quiet"} silences all messaging by the package with the exception of inconsistency errors, and \code{"verbose"} provides information regarding allele name mapping in addition to default messaging on the console
#' @param dataload Controls usage of HLA data tables and by default attempts to load or update and load the tables required for disambiguation. This parameter has the default value \code{"default"}, in which case when versions cannot be determined for \code{allelehistory}, \code{freqbyacc}, or \code{deletedalleles}, the function attempts to load or download and load all the HLA data tables. If the parameter has the value \code{"package"}, then the tables delivered with the package install will be used. If the parameter has the value \code{"argsonly"}, then only when versions can be determined for \code{allelehistory}, \code{freqbyacc}, and \code{deletedalleles} will the function proceed.
#'
#' @return A list with the components:
#'  \item{data}{A data frame with two variables:}
#'    \item{$glstrings}{A vector of disambiguated GL strings, one for each input GL string, in the same order as input. Returns NA upon encountering inconsistency processing arguments.}
#'    \item{$glscores}{A vector of scores, one for each input GL string, in the same order as input.}
#'  
#'  \item{alleletable}{A data frame containing five variables:}
#'    \item{$allele}{The allele for which a match was sought (character).}
#'    \item{$matched}{The allele string that matched, possibly with special rules applied, or NA (character).}
#'    \item{$accession}{The accession number assigned the allele, or NA (character).}
#'    \item{$status}{A factor with possible levels "deleted" (allele was found to have been deleted and replaced), "imgt" (allele was found in the specified IMGT version), "accession-imgt" (allele was identified as an HLA accession number), "rules" (allele was found in the specified IMGT version after applying HLA GL string matching rules), or "unidentified" if an accession number could not be assigned.}
#'
#'  \item{inconsistency.accession}{TRUE if failures to identify accession numbers occurred, FALSE otherwise. }
#'  \item{inconsistency.parse}{TRUE if failures to parse input GL strings or process them were encountered, FALSE otherwise.}
#'  \item{inconsistency.call}{A character string describing an inconsistency. NA if GL strings are returned.}
#'  \item{log}{Verbose details on GL String disambiguation. NA if argument \code{log} is FALSE.}
#'
#' @details In this version, "Cw" and "C" are both permitted in glstrings, even if inappropriate to the IMGT version. Use \code{mode="version"} to see the correct allele names for the version. "HLA-" preceding locus designations will be stripped off into order to avoid incorrect allele pairs calculations that could occur in the case that the prefix is used inconsistently within a locus. Haplotype information is not parsed or sorted in this version, so the behavior of disambiguateR for such data is not well-defined. Because disambiguateR calls canonicalize before returning GL strings, the GL strings returned will be reproducible and all input GL strings that describe logically equivalent information will be reported with the same output GL string. If the mode="accession" parameter is set and disambiguation levels are set identically (dislevel,probratio), two GL strings that use different variations of allele names or different IMGT version namespaces and are determined to be semantically equivalent will result in the same output GL string. The glscores are the means of the maximum allele pair frequencies across all loci. When allele pair frequencies cannot be calculated, the maximum allele pair frequency value for a locus is set to zero. The glscores are not affected by complete or partial disambiguation (for example, the \code{dislevel} and \code{probratio} paramaters do not affect the calculation of glscores).
#' @export
#' @examples
#' gls <- c("A*11011/A*11012+A*3303^C*03041/C*0308/C*0309+C*07011/C*07012/C*0706|C*03041/C*0308/C*0309+C*0705|C*0305+C*07011/C*07012/C*0706|C*0310+C*07011/C*07012/C*0706|C*0310+C*0705^B*4002+B*44031/B*44032","HLA-A*0203+HLA-A*2402101/HLA-A*2402101L/HLA-A*24022/HLA-A*2409N/HLA-A*2411N/HLA-A*2420/HLA-A*2426^HLA-B*40011/HLA-B*40012+HLA-B*5401")
#' # dislevel=0 will perform partial disambiguation. Note extensive log information to stderr.
#' gls.out <- disambiguate(probratio=0.5,freqbyacc=package_HLA_frequencies_by_accession,deletedalleles=package_HLA_deleted_alleles,allelehistory=package_HLA_allele_history,region="NAF",imgtversion = "1050",dislevel=0,mode="normal",glstrings=gls)
#' str(gls.out)
#' gls.out$data
#' # dislevel=1 will perform complete disambiguation, irrespective of probratio value.
#' gls.out <- disambiguate(probratio=0.5,freqbyacc=package_HLA_frequencies_by_accession,deletedalleles=package_HLA_deleted_alleles,allelehistory=package_HLA_allele_history,region="NAF",imgtversion = "1050",dislevel=1,mode="normal",glstrings=gls,log=TRUE)
#' # Top 10 rows of log
#' cat(paste(strsplit(gls.out$log,"\n")[[1]][1:10],collapse="\n"))
#' gls.out$data
#'
#' # not run
#' # Normally, one will want the latest HLA data tables:
#' # updateHLAdata()
#' # The update may take a little time to download.
#' # Once updated data are updated (or later simply loaded with load_HLA_data), there is no need to specify freqbyacc, deletedalleles, or allelehistory for typical usage:
#' # gls.out <- disambiguate(probratio=0.5,region="NAF",imgtversion = "1050",dislevel=0,mode="normal",glstrings=gls,log=TRUE)
#'
#' # For a vector or list of GL Strings, gls, use guessimgtversion to calculate a different imgt version for each GL String, and create a data frame of disambiguated results
#' # require(dplyr)
#' # ret <- lapply(1:length(gls), function(X) {print(X);imgt <- guessimgtversion(gls[X]); ret<-cbind(disambiguate(gls[X],regionstring="EUR",dislevel=0,mode="normal",probratio=0.8,imgtversion=imgt,log=TRUE,allelelog=FALSE)$data,imgt); return(ret)}) %>% bind_rows
#' # ret[1:2,]
#' # write the verbose log for the second glstring to a file for inspection
#' # write(ret$log[2],file="HLA_disambiguateR_tmplog.txt")
#' #
#' # When there is doubt about the geographic region appropriate for disambiguating GL strings, one can investigate the effect of region on the glscores:
#' # sapply(unique(gl_HLA_frequencies_by_accession$data$region), function(X) {d <- disambiguate(test.glstrings,region=X,verbose="quiet",dislevel=1);return(sum(d$data$glscores))})
#'
disambiguate <- function(glstrings,freqbyacc=NULL,deletedalleles=NULL,allelehistory=NULL,regionstring="SSA",imgtversion="guess",probratio=0.5,dislevel=0,ruledat=gl_match_rules,mode="normal",hlaprefix=FALSE,log=FALSE,allelelog=TRUE,verbose="default",dataload="default") {
    # At the end of this code is a brief outline of the algorithm used to disambiguate GL strings. 
    # Much of the code is concerned with processing the function call, handling table loading, and processing the GL string(s).
    # Future code might be easier to understand were the algorithm itself in another function.
    version <- "1.1.2"
    if(log == FALSE) {
        dlog <- NA
    } else {
        dlog <- paste("disambiguate ", version, "\n", sep="")
        dlog <- paste(dlog,paste("freqbyacc version date:",freqbyacc$version.date),sep="\n")
        dlog <- paste(dlog,paste("deletedalleles version date:",deletedalleles$version.date),sep="\n")
        dlog <- paste(dlog,paste("allelehistory version:",allelehistory$version),sep="\n")
    }
    alleletable <- data.frame(allele=NA,matched=NA,accession=NA,status=NA)
    if(!(verbose %in% c("default","quiet","verbose"))) {
        cat(paste("disambiguate",version,"Error: verbose  must be \"default\", \"quiet\", or \"verbose\".\n Exiting..."))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: verbose  must be \"default\", \"quiet\", or \"verbose\".\n Exiting..."),sep="\n")
        }
        data <- data.frame(glstrings=NA,glscores=NA)
        ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="verbose parameter value not recognized",log=dlog)
        return(ret)
    }
    if(verbose == "default" || verbose == "verbose") {
        cat("disambiguateR version ",version,"\n")
    }
    if(!(dataload %in% c("default","package","argonly"))) {
        cat(paste("disambiguate",version,"Error: dataload must be \"default\", \"package\", or \"argonly\".\n Exiting..."))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: dataload  must be \"default\", \"package\", or \"argonly\".\n Exiting..."),sep="\n")
        }
        data <- data.frame(glstrings=NA,glscores=NA)
        ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="dataload parameter value not recognized",log=dlog)
        return(ret)
    }
    # check HLA data loading parameters, load or download and load data if needed, depending on dataload parameter
    if(is.null(get0("freqbyacc")$version.date) || is.null(get0("allelehistory")$version) || is.null(get0("deletedalleles")$version.date)) {
        if(dataload == "default") {
            # try loading
            if(verbose == "verbose") {
                cat("disambiguateR freqbyacc, allelehistory, and deletedalleles must all be defined and point to valid objects or default data will be loaded\n")
                cat("disambiguateR attempting to load HLA data tables...\n")
            }
            if(!((ltest <- try(loadHLAdata())) == "try-error") && ltest == 0) {
                # seems good to go
                freqbyacc <- gl_HLA_frequencies_by_accession
                deletedalleles <- gl_HLA_deleted_alleles
                allelehistory <- gl_HLA_allele_history
                if(verbose == "verbose") {
                    cat("disambiguateR is using gl_HLA_frequencies_by_accession, gl_HLA_deleted_alleles, and gl_HLA_allele_history\n")
                }
            } else if(!(try(updateHLAdata()) == "try-error")) {
                # downloading and loading
                cat("disambiguateR freqbyacc, allelehistory, and deletedalleles must all be defined and point to valid objects or default data will be loaded...which failed\n")
                cat("disambiguateR could not load HLA data tables and attempted to update (see updateHLAdata).\nPlease read console messages and note any errors above. Trying downloaded data...\n")
                freqbyacc <- gl_HLA_frequencies_by_accession
                deletedalleles <- gl_HLA_deleted_alleles
                allelehistory <- gl_HLA_allele_history
                cat("disambiguateR is using gl_HLA_frequencies_by_accession, gl_HLA_deleted_alleles, and gl_HLA_allele_history\n")
            } else {
                # calling download failed, exit gracefully
                cat("disambiguateR data tables load failed, attempt to download and load data also failed, Exiting...\n")
                data <- data.frame(glstrings=NA,glscores=NA)
                ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="freqbyacc, allelehistory, and deletedalleles parameters values were not all valid and attempts to load and download failed",log=dlog)
                return(ret)
            }
        } else if(dataload == "package") {
            # define tables as package-shipped data
            freqbyacc <- package_HLA_frequencies_by_accession
            deletedalleles <- package_HLA_deleted_alleles
            allelehistory <- package_HLA_allele_history
            if(verbose == "verbose") {
                cat("disambiguateR is using package_HLA_frequencies_by_accession, package_HLA_deleted_alleles, and package_HLA_allele_history\n")
            }
        } else if(dataload == "argsonly") {
            # warn regarding missing data and exit gracefully
            cat("disambiguateR freqbyacc, allelehistory, and deletedalleles must all be defined and point to valid HLA data objects, Exiting...\n")
            data <- data.frame(glstrings=NA,glscores=NA)
            ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="freqbyacc, allelehistory, and deletedalleles parameters values were not all valid",log=dlog)
            return(ret)
        }
    }
    deletedalleles <- cbind(paste(deletedalleles$data[,1],deletedalleles$data[,2],sep=""),paste(deletedalleles$data[,1],deletedalleles$data[,5],sep=""))
    allelehistory <- allelehistory$data

    # Evaluate remaining arguments to this function, return with inconsistencies as needed
    if(!(regionstring %in% unique(freqbyacc$data$region))) {
        cat(paste("disambiguate",version,"Error: region",regionstring,"not recognized. Possible regions are:",paste(as.list(unique(freqbyacc$data$region)),collapse=", "),". Exiting...\n"))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: region",regionstring,"not recognized. Possible regions are:",paste(as.list(unique(freqbyacc$data$region)),collapse=", "),". Exiting..."),sep="\n")
        }
        data <-  data.frame(glstrings=NA,glscores=NA)
        ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="region not recognized",log=dlog)
        return(ret)
    }
    if(verbose == "default" || verbose == "verbose") {
        cat(paste("...region requested:",regionstring),"\n")
    }
    if(log == TRUE) {
        dlog <- paste(dlog,paste("...region requested:",regionstring),sep="\n")
    }
    if(!(dislevel == 0 | dislevel == 1)) {
        cat(paste("disambiguate",version,"Error: disambiguation level must be 1 (complete disambiguation) or 0 (partial disambiguation). Exiting...\n"))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: disambiguation level must be 1 (complete disambiguation) or 0 (partial disambiguation). Exiting..."),sep="\n")
        }
        data <- data.frame(glstrings=NA,glscores=NA)
        ret<- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="disambiguation level must be 0 or 1",log=dlog)
        return(ret)
    }
    if(verbose == "default" || verbose == "verbose") {
        cat(paste("...disambiguation level:",dislevel),"\n")
    }
    if(log == TRUE) {
        dlog <- paste(dlog,paste("...disambiguation level:",dislevel),sep="\n")
    }
    if(!(mode %in% c("normal","accession","version"))) {
        cat(paste("disambiguate",version,"Error: mode must be \"normal\" (returns GL strings with alleles as named in input, with only minor standardization), \"accession\" (returns GL strings with accession numbers), or \"version\" (returns GL strings with alleles standardized to the names matching the specified version)\n Exiting..."))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: mode must be \"normal\" (returns GL strings with alleles as named in input, with only minor standardization), \"accession\" (returns GL strings with accession numbers), or \"version\" (returns GL strings with alleles standardized to the names matching the specified version)\n Exiting..."),sep="\n")
        }
        data <- data.frame(glstrings=NA,glscores=NA)
        ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="mode not recognized",log=dlog)
        return(ret)
    }
    if(verbose == "default" || verbose == "verbose") {
        cat(paste("...mode:",mode),"\n")
    }
    if(log == TRUE) {
        dlog <- paste(dlog,paste("...mode:",mode),sep="\n")
    }
    if(!is.logical(hlaprefix)) {
        cat(paste("disambiguate",version,"Error: hlaprefix must be a logical.\n Exiting..."))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: hlaprefix must be a logical\n Exiting..."),sep="\n")
        }
        data <- data.frame(glstrings=NA,glscores=NA)
        ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="hlaprefix not a logical",log=dlog)
        return(ret)
    }
    if(verbose == "default" || verbose == "verbose") {
        cat(paste("...hlaprefix:",hlaprefix),"\n")
    }
    if(log == TRUE) {
        dlog <- paste(dlog,paste("...hlaprefix:",hlaprefix),sep="\n")
    }
    if(!(probratio <= 1.0 & probratio >= 0)) {
        cat(paste("disambiguate",version,"Error: probratio",probratio,"not <= 1 and >= 0. Exiting...\n"))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: probratio",probratio,"not <= 1 and >= 0. Exiting..."),sep="\n")
        }
        data <- data.frame(glstrings=NA,glscores=NA)
        ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="probratio not <= 1 and >= 0",log=dlog)
        return(ret)
    } else if(dislevel == 1 && probratio < 1.0) {
        if(verbose == "default" || verbose == "verbose") {
            cat("...disambiguate ",version," Warning: dislevel set to 1, probratio will be reset to 1.0.\n")
        }
        if(log == TRUE) {
            dlog <- paste(dlog,"Warning: dislevel set to 1, probratio will be reset to 1.0.",sep="\n")
        }
        probratio <- 1.0
    }
    if(verbose == "default" || verbose == "verbose") {
        cat(paste("...probratio:", probratio),"\n")
    }
    if(log == TRUE) {
        dlog <- paste(dlog,paste("...probratio:", probratio),sep="\n")
    }
    if(imgtversion == "guess") {
        imgtversion <- try(guessimgtversion(glstrings,allelehistory=allelehistory))
        if(class(imgtversion) == "try-error") {
            cat(paste("disambiguate",version,"Error: guessimgtversion failed. Exiting...\n"))
            if(log == TRUE) {
                dlog <- paste(dlog,paste("disambiguate",version,"Error: guessimgtversion failed. Exiting..."),sep="\n")
            }
            data <- data.frame(glstrings=NA,glscores=NA)
            ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="guessimgtversion failed",log=dlog)
            return(ret)
        }
    } else if(!(paste("X",imgtversion,sep="") %in% names(allelehistory))) {
        cat(paste("disambiguate",version,"Error: IMGT version ",imgtversion,"not recognized. Exiting...\n"))
        cat("...Allowed IMGT versions:\n")
        cat(paste(c(gsub("X","",names(allelehistory)[-1]),"\n"),collapse="\n"))
        if(log == TRUE) {
            dlog <- paste(dlog,paste("disambiguate",version,"Error: IMGT version ",imgtversion,"not recognized. Exiting..."),sep="\n")
            dlog <- paste(dlog,"...Allowed IMGT versions:\n",sep="\n")
            dlog <- paste(dlog,paste(c(gsub("X","",names(allelehistory)[-1]),"\n"),collapse="\n"),sep="\n")
        }
        data <- data.frame(glstrings=NA,glscores=NA)
        ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=FALSE,inconsistency.call="IMGT version not recognized",log=dlog)
        return(ret)
    }
    if(verbose == "default" || verbose == "verbose") {
        cat(paste("...IMGT version:", imgtversion),"\n")
    }
    if(log == TRUE) {
        dlog <- paste(dlog,paste("...IMGT version:", imgtversion),sep="\n")
    }

    Sys.setlocale("LC_COLLATE", "C") # turn off locale-specific sorting, usually, but not on all platforms
    inconsistency.accession = FALSE
    inconsistency.parse = FALSE

    name_version <- paste("v",imgtversion,sep="")

    # remove HLA- throughout glstrings
    glstrings <- gsub("HLA-","",glstrings)

    # split glstrings into alleles
    allalleles <- unique(sort(unlist(sapply(1:length(glstrings),function (i) {unlist(strsplit(glstrings[i],"\\^|\\||\\+|\\/"))}))))

    # track Cw in alleles
    Cw <- grepl("Cw",allalleles)
    if(!(name_version < "v3") & any(Cw) == TRUE) {
        if(verbose == "default" || verbose == "verbose") {
            cat(paste("disambiguate",version,"Warning: Cw locus designations are not valid for IMGT versions > 2"),"\n")
        }
        if(log == TRUE) {
            dlog <- paste(dlog,cat(paste("disambiguate",version,"Warning: Cw locus designations are not valid for IMGT versions > 2")),sep="\n")
        }
    }
    # track C in alleles
    C <- grepl("C",allalleles)
    if(name_version < "v3" & all(Cw == C) == FALSE) {
        if(verbose == "default" || verbose == "verbose") {
            cat(paste("disambiguate",version,"Warning: C locus designations are not valid for IMGT versions < 3"),"\n")
        }
        if(log == TRUE) {
            dlog <- paste(dlog,cat(paste("disambiguate",version,"Warning: C locus designations are not valid for IMGT versions < 3")),sep="\n")
        }
    }
    if(any(Cw) == TRUE & all(Cw == C) == FALSE) {
        # C and Cw designations are both found in this collection of alleles
        if(verbose == "default" || verbose == "verbose") {
            cat(paste("disambiguate",version,"Warning: C and Cw locus designations are not both valid any IMGT version, yet both are found in GL Strings"),"\n")
        }
        if(log == TRUE) {
            dlog <- paste(dlog,cat(paste("disambiguate",version,"Warning: C and Cw locus designations are not both valid any IMGT version, yet both were found in GL Strings.")),sep="\n")
        }
    }

    # convert all Cw to C in glstrings, deleted alleles, and allele history and keep in separate objects
    m.allalleles <- allalleles
    m.allalleles <- gsub("Cw","C",m.allalleles)
    m.deletedalleles <- deletedalleles
    m.deletedalleles[,1] <- gsub("Cw","C",m.deletedalleles[,1])
    m.deletedalleles[,2] <- gsub("Cw","C",m.deletedalleles[,2])
    ahistory <- allelehistory[c("HLA_ID",paste("X",imgtversion,sep=""))]
    m.ahistory <- ahistory
    m.ahistory[,2] <- gsub("Cw","C",m.ahistory[,2])

    if(log == TRUE) {
        dlog <- paste(dlog,paste(length(allalleles),"unique allele strings found in glstrings"),sep="\n")
    }

    allele2acc <- vector()
    alleleout <- vector()
    if(allelelog == TRUE) {
        alleletable <- matrix(rep(NA,length(allalleles)*4),length(allalleles),4)
    }
    for(i in 1:length(allalleles)) {
        found <- 0
        if (length((res <- t(ahistory[which(m.ahistory[,2] == m.allalleles[i]),]))) > 0) {
            # this is an allele name in the designated IMGT version
            allele2acc[allalleles[i]] <- res[1]
            if(log == TRUE) {
                dlog <- paste(dlog,paste(allalleles[i]," was found in IMGT version ", imgtversion," alleles and assigned accession ",res[1],".",sep=""),sep="\n")
            }
            if(allelelog == TRUE) {
                alleletable[i,1] <- allalleles[i]
                alleletable[i,2] <- res[2]
                alleletable[i,3] <- res[1]
                alleletable[i,4] <- "imgt"
            }
            found <- 1
            if(mode == "accession") {
                alleleout[allalleles[i]] <- res[1]
            } else if (mode == "version"){
                alleleout[allalleles[i]] <- res[2]
            } else if (mode == "normal") {
                alleleout[allalleles[i]] <- allalleles[i]
            }
        } else if(length((res <- t(ahistory[which(m.ahistory[,1] == m.allalleles[i]),]))) > 0) {
            # this is an HLA accession number
            allele2acc[allalleles[i]] <- res[1]
            if(log == TRUE) {
                dlog <- paste(dlog,paste(allalleles[i]," was identified as an HLA accession number",sep=""),sep="\n")
            }
            if(allelelog == TRUE) {
                alleletable[i,1] <- allalleles[i]
                alleletable[i,2] <- res[2]
                alleletable[i,3] <- res[1]
                alleletable[i,4] <- "accession-imgt"
            }
            found <- 1
            if(mode == "accession") {
                alleleout[allalleles[i]] <- res[1]
            } else if (mode == "version"){
                alleleout[allalleles[i]] <- res[2]
            } else if (mode == "normal") {
                alleleout[allalleles[i]] <- allalleles[i]
            }
        } else if(length((res <- t(deletedalleles[which(m.deletedalleles[,1] == m.allalleles[i]),]))) > 0 & length((res2 <- as.character(ahistory[which(m.ahistory[,2] == res[2]),1]))) > 0) {
            # this is a deleted allele name with known mapping to an allele with an accession number
            allele2acc[allalleles[i]] <- res2
            if(log == TRUE) {
                dlog <- paste(dlog,paste(allalleles[i]," was found in deleted alleles and assigned accession ",res2,".",sep=""),sep="\n")
            }
            if(allelelog == TRUE) {
                alleletable[i,1] <- allalleles[i]
                alleletable[i,2] <- res[2]
                alleletable[i,3] <- res2
                alleletable[i,4] <- "deleted"
            }
            found <- 1
            if(mode == "accession") {
                alleleout[allalleles[i]] <- res2
            } else if(mode == "version") {
                alleleout[allalleles[i]] <- res[2]
            } else if(mode == "normal") {
                alleleout[allalleles[i]] <- allalleles[i]
            }
        } else {
            # try special rules to modify the allele name and check against the IMGT version
            name_numlength <- nchar(strsplit(m.allalleles[i],"\\*")[[1]][2])
            for(j in 1:(dim(ruledat)[1])) {
                if(eval(parse(text=ruledat[j,1]))) {
                    for(k in 1:length(ruledat$rule[j][[1]])) {
                        append <- ""
                        eval(parse(text=ruledat$rule[j][[1]][k]))
                        tryallelename <- paste(m.allalleles[i],append,sep="")
                        if(length((res <- t(ahistory[which(m.ahistory[,2] == tryallelename),]))) > 0) {
                            # we have a match with a modified name
                            allele2acc[allalleles[i]] <- res[1]
                            if(log == TRUE) {
                                dlog <- paste(dlog,paste(allalleles[i]," was found in IMGT version ", imgtversion," alleles and assigned accession ",res[1]," using special rules. The modified string that matched was ",tryallelename,".",sep=""),sep="\n")
                            }
                            if(allelelog == TRUE) {
                                alleletable[i,1] <- allalleles[i]
                                alleletable[i,2] <- res[2]
                                alleletable[i,3] <- res[1]
                                alleletable[i,4] <- "rules"
                            }
                            found <- 1
                            if(mode == "accession") {
                                alleleout[allalleles[i]] <- res[1]
                            } else if(mode == "version") {
                                alleleout[allalleles[i]] <- res[2]
                            } else if(mode=="normal") {
                                alleleout[allalleles[i]] <- allalleles[i]
                            }
                            break
                        }
                    }
                }
            }
        }
        if(found == 0) {
            # failed to assign an accession number to this allele string
            if(verbose == "default" || verbose == "verbose") {
                cat(paste("disambigluate ",version," Warning: ",allalleles[i]," could not be assigned an accession number from deleted alleles, IMGT version " , imgtversion,", or with special rules to match version alleles.",sep=""),"\n")
            }
            if(log == TRUE) {
                dlog <- paste(dlog,paste("disambiguate ",version," Warning: ",allalleles[i]," could not be assigned an accession number from deleted alleles, IMGT version ", imgtversion,", or with special rules to match version alleles.",sep=""),sep="\n")
            }
            if(allelelog == TRUE) {
                alleletable[i,1] <- allalleles[i]
                alleletable[i,2] <- NA
                alleletable[i,3] <- NA
                alleletable[i,4] <- "unidentified"
            }
            inconsistency.accession = TRUE
            if(mode == "accession" | mode == "version") {
                alleleout[allalleles[i]] <- paste(allalleles[i],"Unidentified",sep="")
            } else {
                alleleout[allalleles[i]] <- allalleles[i]
            }
        }
    }
    alleletable <- data.frame(allele=as.character(alleletable[,1]),matched=as.character(alleletable[,2]),accession=as.character(alleletable[,3]),status=as.factor(alleletable[,4]))

    if(hlaprefix==TRUE) {
	    alleleout <- gsub("^","HLA-",alleleout)
    }

    glout <- vector()
    glscores <- vector()

    # Iterate through the glstrings and disambiguate each one. An outline of the algorithm is at the end of this code.
    for(ng in 1:length(glstrings)) {
        glstring = glstrings[ng]
        if(log == TRUE) {
            dlog <- paste(dlog,paste("\ndisambiguate",version,": Disambiguating GL String (cleaned):\n    ",glstring,"\n",sep=""),sep="\n")
        }
        loci <- sort(unlist(strsplit(glstring,"\\^")))
        lociout <- vector()
        lociscores <- vector()
        for(locusn in 1:(length(loci))) {
            scoretmp <- 0
            if(log == TRUE) {
                dlog <- paste(dlog,paste("locus =",loci[locusn]),sep="\n")
            }
            genotypelist = sort(unlist(strsplit(loci[locusn],"\\|")))
            locusout <- ""
            locusfreq <- -1
            pairs <- data.frame(matrix(NA,nrow=1,ncol=4))
            names(pairs) <- c("a1","a2","prob","onefreq")
            npairs <- 0
            foundallelepair <- vector()
            for(ngl in 1:length(genotypelist)) {
                if(log == TRUE) {
                    dlog <- paste(dlog,paste("genotypelist =",genotypelist[ngl]),sep="\n")
                }
                allelelist <- sort(unlist(strsplit(genotypelist[ngl],"\\+")))
                if(length(allelelist) > 2) {
                    cat(paste("disambiguate ",version," Error: in GL String:\n    ",  glstring, "\n    \n    genotype:\n    ", genotypelist[ngl], "\n    has more than than two ambiguous alleles, because multiple + delimiters were found. Exiting...\n",sep=""))
                    if(log == TRUE) {
                        dlog <- paste(dlog,paste("disambiguate ",version," Error: in GL String:\n    ",  glstring, "\n    \n    genotype:\n    ", genotypelist[ngl], "\n    has more than than two ambiguous alleles, because multiple + delimiters were found. Exiting...\n",sep=""),sep="\n")
                    }
                    inconsistency.parse = TRUE
                    data <- data.frame(glstrings=NA,glscores=NA)
                    ret <- list(data=data,alleletable=alleletable,inconsistency.accession=FALSE,inconsistency.parse=inconsistency.parse,inconsistency.call="multiple + delimiters found",log=dlog)
                    return(ret)
                }
                alleles1 <- sort(unlist(strsplit(allelelist[1],"\\/")))
                alleles2 <- sort(unlist(strsplit(allelelist[2],"\\/")))
                if(length(alleles1) >= 1 & length(alleles2) >= 1) {
                    for(i in 1:length(alleles1)) {
                        res <- subset(freqbyacc$data, region == regionstring & accession == allele2acc[alleles1[i]], select = "frequency")
                        if(length(res[[1]]) == 0) {
                            freq1 <- NA
                        }
                        else {
                            freq1 <- res[1][[1]]
                        }
                        for(j in 1:length(alleles2)) {
                            res <- subset(freqbyacc$data, region == regionstring & accession == allele2acc[alleles2[j]], select = "frequency")
                            if(length(res[[1]]) == 0) {
                                freq2 <- NA
                            }
                            else {
                                freq2 <- res[1][[1]]
                            }
                            # at this point we should have a pair of alleles and their frequencies in the designated geographic region
                            a1 <- alleleout[alleles1[i]]
                            a2 <- alleleout[alleles2[j]]
                            if(a2 < a1) {
                                atmp <- a1
                                a1 <- a2
                                a2 <- atmp
                            }
                            if(is.na(foundallelepair[paste(a1,a2)])) {
                                foundallelepair[paste(a1,a2)] <- 1
                                if(is.na(sort(c(freq1,freq2),na.last=TRUE)[1])) {
                                    prob <- NA
                                    onefreq <- NA
                                } else if (is.na(freq1*freq2)) {
                                    prob <- NA
                                    onefreq <- sort(c(freq1,freq2),na.last=TRUE)[1]
                                } else {
                                    # This is the allele pair frequency product (APFP). See notes on the algorithm at the end of this code.
                                    prob <- freq1 * freq2
                                    onefreq <- NA
                                }

                                npairs <- npairs + 1
                                if(npairs == 1) {
                                    pairs$a1[1] <- a1
                                    pairs$a2[1] <- a2
                                    pairs$prob[1] <- prob
                                    pairs$onefreq[1] <- onefreq
                                } else {
                                    pairs <- rbind(pairs,c(a1,a2,prob,onefreq))
                                }
                            }
                        }
                    }
                } else if (length(alleles1) >= 1) {
                    # we have only one allele, so for purposes of disambiguation will assume homozygosity for this allele
                    if(log == TRUE) {
                        dlog <- paste(dlog,paste("Assuming genotype data implies homozygosity in GL String:\n  ",  glstring, "\n       \n      genotype:\n     ", genotypelist[ngl],"\n", sep=""),sep="\n")
                    }
                    for(i in 1:length(alleles1)) {
                        res <- subset(freqbyacc$data, region == regionstring & accession == allele2acc[alleles1[i]], select = "frequency")
                        if(length(res[[1]]) == 0) {
                            freq1 <- NA
                        }
                        else {
                            freq1 <- res[1][[1]]
                        }
                        freq2 <- freq1
                        a1 <- alleleout[alleles1[i]]
                        a2 <- a1
                        if(is.na(foundallelepair[paste(a1,a2)])) {
                            foundallelepair[paste(a1,a2)] <- 1
                            if(is.na(sort(c(freq1,freq2),na.last=TRUE)[1])) {
                                prob <- NA
                                onefreq <- NA
                            } else {
                                prob <- freq1*freq2
                                onefreq <- NA
                            }
                            npairs <- npairs + 1
                            if(npairs == 1) {
                                pairs$a1[1] <- a1
                                pairs$a2[1] <- a1
                                pairs$prob[1] <- prob
                                pairs$onefreq[1] <- onefreq
                            } else {
                                pairs <- rbind(pairs,c(a1,a2,prob,onefreq))
                            }
                        }
                    }
                }
            }
            if(dim(pairs)[1] >= 1) {
                pairs$prob <- as.numeric(pairs$prob)
                pairs$onefreq <- as.numeric(pairs$onefreq)
            }
            if(log == TRUE) {
                dlog <- paste(dlog,paste(capture.output(write.table(as.matrix(pairs), quote=FALSE, row.names=FALSE)),sep="\n",collapse="\n"),sep="\n")
            }

            # Sort pairs by a1, a2
            pairs <- pairs[order(pairs$a1,pairs$a2),]

            if(!all(is.na(pairs$prob)) & any(pairs$prob >= 0)) {
                # This is the maximum allele pair frequency product (MAPFP). See notes on the algorithm at the end of this code.
                locusfreq <- max(pairs$prob,na.rm=TRUE)
                scoretmp <- locusfreq
                locusfreq <- locusfreq * probratio
                pairstmp <- subset(pairs, prob >= locusfreq)
                pairstmp <- pairstmp[order(pairstmp$a1,pairstmp$a2),]
                if(dislevel == 1) {
                    locusout <- paste(pairstmp$a1[1],"+",pairstmp$a2[1],sep="")
                } else {
                    locusout <- paste(pairstmp$a1[1],"+",pairstmp$a2[1],sep="")
                    if(dim(pairstmp)[1] > 1) {
                        for(k in 2:(dim(pairstmp)[1])) {
                            locusout <- paste(locusout,"|",pairstmp$a1[k],"+",pairstmp$a2[k],sep="")
                        }
                    }
                }
            } else if(!all(is.na(pairs$onefreq)) & any(pairs$onefreq >= 0)) {
                if(log == TRUE) {
                    dlog <- paste(dlog,paste("disambiguate ",version," Warning: for GL String:\n    ",  glstring, "\n    locus data:\n    ", loci[locusn], "\n    no frequency pair data were defined that could distinguish possible allele pairs, so frequencies of single alleles were used when one allele frequency in a pair was undefined\n",sep=""),sep="\n")
                }
                locusfreq <- max(pairs$onefreq,na.rm=TRUE)
                locusfreq <- locusfreq * probratio
                pairstmp <- subset(pairs, onefreq >= locusfreq)
                pairstmp <- pairstmp[order(pairstmp$a1,pairstmp$a2),]
                if(dislevel == 1) {
                    locusout <- paste(pairstmp$a1[1],"+",pairstmp$a2[1],sep="")
                } else {
                    locusout <- paste(pairstmp$a1[1],"+",pairstmp$a2[1],sep="")
                    if(dim(pairstmp)[1] > 1) {
                        for(k in 2:(dim(pairstmp)[1])) {
                            locusout <- paste(locusout,"|",pairstmp$a1[k],"+",pairstmp$a2[k],sep="")
                        }
                    }
                }
            } else {
                if(log == TRUE) {
                    dlog <- paste(dlog,paste("disambiguate ",version," Warning: for GL String:\n    ",  glstring, "\n    locus data:\n    ", loci[locusn], "\n    no frequency data were defined that could distinguish possible allele pairs\n",sep=""),sep="\n")
                }
                pairstmp <- pairs
                pairstmp <- pairstmp[order(pairstmp$a1,pairstmp$a2),]
                if(dislevel == 1) {
                    locusout <- paste(pairstmp$a1[1],"+",pairstmp$a2[1],sep="")
                } else {
                    locusout <- paste(pairstmp$a1[1],"+",pairstmp$a2[1],sep="")
                    if(dim(pairstmp)[1] > 1) {
                        for(k in 2:(dim(pairstmp)[1])) {
                            locusout <- paste(locusout,"|",pairstmp$a1[k],"+",pairstmp$a2[k],sep="")
                        }
                    }
                }
                
            }
            lociout[locusn] <- locusout
            lociscores[locusn] <- scoretmp
        }
        glscores[ng] <- mean(lociscores)
        gltmp <- canonicalize(paste(lociout,collapse="^"))
        glout[ng] <- gltmp$glstring
        if(gltmp$inconsistency == "") {
            if(log == TRUE) {
                dlog <- paste(dlog,"Disambiguated GL String:",glout[ng],sep="\n")
            }
        } else {
                inconsistency.parse = TRUE
                dlog <- paste(dlog,paste("disambiguate ",version," Error: processing GL String:\n    ",  glstring, "\n    canonicalize failed for:\n    ", paste(lociout,collapse="^"), "\n    with inconsistency: ", gltmp$inconsistency,"\n",sep=""),sep="\n")
            
        }
    }
    data <- data.frame(glstrings=I(as.character(glout)),glscores)
    ret <- list(data=data,alleletable=alleletable,inconsistency.accession=inconsistency.accession,inconsistency.parse=inconsistency.parse,inconsistency.call=NA,log=dlog)
    return(ret)

    # Notes on the overall algorithm
    # for each allele
    # identify deleted alleles
        # deleted
        # not deleted
            # search for accession, applying special rules as needed
            # report alleles for which accession cannot be identified
    # if dislevel == 0 & probratio < 1 the algorithm is
    # for each locus
        # create all possible allele pairs consistent with GL String
        # calculate allele pair frequency products (APFPs)
            # Maximum allele pair frequency product (MAPFP) > 0
                # keep allele pairs for which APFP > R*MAPFP, for probratio R
            # Maximum allele pair frequency product (MAPFP) == 0
                # all allele frequncies = 0
                # yes - keep all allele pairs
                # no - keep allele pairs for which one allele has a frqeuency in the region >= max frequency for pairs with only one allele
}
