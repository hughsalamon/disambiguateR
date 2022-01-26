#' Canonicalize Genotype List (GL) Strings
#'
#' This function parses a GL string, represents information at each locus as an adjacency matrix, and returns a sorted, canonicalized GL string. Two GL strings that contain the same information regarding genotypes, not matter how different, both will be converted to an identical, compact GL string.
#' @param glstring A character string containing valid GL string information on alleles at one or more loci.
#' @return A list with the components:
#'  \item{glstring}{A character string transformed according to the sorting and simplifying algorithm in the function or an empty string in the case of an inconsistency error.}
#'  \item{inconsistency}{An inconsistency error message string or empty string.}
#' @details Any string containing the GL string delimiters ^, +, /, or | will be parsed, analyzed, and (possibly) transformed.  Strings not containing these symbols should be returned unchanged. No validation of human leukocyte antigen (HLA) or killer cell immunoglobulin-like receptor (KIR) GL strings is performed.
#' @export
#' @examples
#' gls <- c("HLA-A*11011/HLA-A*11012+HLA-A*2410^HLA-C*0702/HLA-C*0710+HLA-C*0704/HLA-C*0711/HLA-C*0712|HLA-C*0703+HLA-C*0704/HLA-C*0711/HLA-C*0712^HLA-B*1532/HLA-B*1535+HLA-B*1803","HLA-A*02011/HLA-A*02012/HLA-A*02014/HLA-A*0209/HLA-A*0230/HLA-A*0231+HLA-A*0207/HLA-A*0215N/HLA-A*0218|HLA-A*0207/HLA-A*0215N/HLA-A*0218+HLA-A*0207/HLA-A*0215N/HLA-A*0218^HLA-C*03041/HLA-C*0308/HLA-C*0309+HLA-C*0702/HLA-C*0710|HLA-C*0310+HLA-C*0702/HLA-C*0710^HLA-B*4601+HLA-B*4601")
#' sapply(gls, canonicalize)
#'
#' # the following represent identical genotypes
#' gls <- c("HLA-C*0308/HLA-C*0309+HLA-C*0702/HLA-C*0710|HLA-C*0310+HLA-C*0702/HLA-C*0710","HLA-C*0710+HLA-C*0310/HLA-C*0308/HLA-C*0309|HLA-C*0702+HLA-C*0310/HLA-C*0308/HLA-C*0309")
#' x <- as.data.frame(t(sapply(gls,canonicalize)))
#' x$glstring[[1]] == x$glstring[[2]]
#'
canonicalize <- function (glstring) {
    checkval <- checkglstring(glstring)
    if(!is.null(checkval)) {
        return(list(glstring="",inconsistency=paste("Error: a locus cannot appear in multiple gene fields! ",checkval," appears in multiple '^'-delimited fields.",sep="")))
    }
    loci <- sort(unlist(strsplit(glstring,"\\^")))
    lociout <- vector()
	for(locusn in 1:(length(loci))) {
		alleles <- unique(sort(unlist(strsplit(loci[locusn],"[\\|\\+\\/]"))))
		m <- matrix(0,length(alleles),length(alleles))
		rownames(m) <- alleles
		colnames(m) <- rownames(m)
		genotypesout <- vector()
		ng <- 0
		genotypes <- sort(unlist(strsplit(loci[locusn],"[\\|]")))
		for(genotypen in 1:(length(genotypes))) {
			copies <- sort(unlist(strsplit(genotypes[genotypen],"[\\+]")))
			if(length(copies) == 2) {
				# This is the most typical case, pairs of alleles, possibly with ambiguity.  
				# We will populate the adjacency matrix with information representing the pairs of alleles.
				alleles1 <- sort(unlist(strsplit(copies[1],"\\/")))
				alleles2 <- sort(unlist(strsplit(copies[2],"\\/")))
				for(j in 1:length(alleles1)) {
					for(k in 1:length(alleles2)) {
						m[alleles1[j],alleles2[k]] <- 1
						m[alleles2[k],alleles1[j]] <- 1
					}
				}
			} else if(length(copies) == 1) {
				ng <- ng + 1
				genotypesout[ng] <- paste(sort(unlist(strsplit(copies[1],"\\/"))),collapse="/")
			} else {
				for(j in 1:length(copies)) {
					copies[j] <- paste(sort(unlist(strsplit(copies[j],"\\/"))),collapse="/")
				}
				ng <- ng + 1
				genotypesout[ng] <- paste(sort(unlist(copies)),collapse="+")
			}
		}
		# Now we analyze the adjacency matrix
		while(sum(m) > 0) {
			cs <- colSums(m)
			mm <- m
			if((dim(mm)[2] > 1)) {
				mm <- mm[,order(apply(mm,2,paste,collapse=""))]
				# find clusters
				last <- 1
				counter <- 1
				countern <- 1
				cluster <- list()
				cluster[[counter]] <- vector()
				cluster[[counter]][countern] <- colnames(mm)[1]
				for(i in 2:(dim(mm)[2])) {
					if(all(mm[,i]==mm[,last])) {
						countern <- countern + 1
						cluster[[counter]][countern] <- colnames(mm)[i]
					} else {
						last <- i
						counter <- counter + 1
						countern <- 1
						cluster[[counter]] <- vector()
						cluster[[counter]][countern] <- colnames(mm)[i]
					}
				}
				clustersizes <- unlist(lapply(cluster,length))
				clusterweights <- clustersizes * unlist(lapply(cluster,function (X) {sum(m[,X[1]])}))
				clusterlength <- unlist(lapply(cluster,function (X) {sum(m[,X[1]])}))
				if(sum(clusterweights == max(clusterweights)) > 1) {
					# need to use tiebreaker
					ties <- clusterlength[which(clusterweights == max(clusterweights))]
					names(ties) <- which(clusterweights == max(clusterweights))
					minties <- ties[ties == min(ties)]
					if(length(minties) == 1) {
						# we have broken the tie
						# create string for the cluster with the max weight, and remove cluster members from matrix
						thiscluster <- cluster[[as.numeric(names(minties))]]
						set1 <- paste(sort(thiscluster),collapse="/")
						set2 <- paste(sort(rownames(m)[which(m[,thiscluster[1]] == 1)]),collapse="/")
						ng <- ng + 1
						genotypesout[ng] <- paste(sort(c(set1,set2)),collapse="+")
						m <- m[,!(colnames(m) %in% thiscluster),drop=FALSE]
						m <- m[!(rownames(m) %in% thiscluster),,drop=FALSE]
					
					} else {
						# ties remain, use alphabetic sort of cluster allele names
						tmpstring <- unlist(lapply(names(minties),function (X) {thiscluster <- cluster[[as.numeric(X)]];paste(sort(c(thiscluster,rownames(m)[which(m[,thiscluster[1]] == 1)])),collapse="/")}))
						# we have broken the tie
						# create string for the cluster with the alphabetically first max weight, and remove cluster members from matrix
						thiscluster = cluster[[as.numeric(names(minties[which(order(tmpstring)==1)]))]]
						set1 <- paste(sort(thiscluster),collapse="/")
						set2 <- paste(sort(rownames(m)[which(m[,thiscluster[1]] == 1)]),collapse="/")
						ng <- ng + 1
						genotypesout[ng] <- paste(sort(c(set1,set2)),collapse="+")
						m <- m[,!(colnames(m) %in% thiscluster),drop=FALSE]
						m <- m[!(rownames(m) %in% thiscluster),,drop=FALSE]
					}
				} else {
					# create string for the cluster with the max weight, and remove cluster members from matrix
					thiscluster <- cluster[[which(clusterweights == max(clusterweights))]]
					set1 <- paste(sort(thiscluster),collapse="/")
					set2 <- paste(sort(rownames(m)[which(m[,thiscluster[1]] == 1)]),collapse="/")
					ng <- ng + 1
					genotypesout[ng] <- paste(sort(c(set1,set2)),collapse="+")
					m <- m[,!(colnames(m) %in% thiscluster),drop=FALSE]
					m <- m[!(rownames(m) %in% thiscluster),,drop=FALSE]
				}
			}
			else {
				# we have only one pattern
				# create string for the cluster with the max number of pairs, and remove cluster members from matrix
				thiscluster <- colnames(mm)
				set1 <- paste(sort(thiscluster),collapse="/")
				set2 <- paste(sort(rownames(m)[which(m[,thiscluster[1]] == 1)]),collapse="/")
				ng <- ng + 1
				genotypesout[ng] <- paste(sort(c(set1,set2)),collapse="+")
				m <- m[,!(colnames(m) %in% thiscluster),drop=FALSE]
				m <- m[!(rownames(m) %in% thiscluster),,drop=FALSE]
			}
		}
		lociout[locusn] <- paste(sort(genotypesout),collapse="|")
		rm(m)
	}
	return(list(glstring=paste(sort(lociout),collapse="^"),inconsistency=""))
}
