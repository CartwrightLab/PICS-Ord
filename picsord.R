#!/usr/bin/Rscript --vanilla
################################################################################
# Copyright (c) 2010 Reed A. Cartwright, PhD <reed@scit.us>
# Released under GNU General Public License, version 2
# http://www.gnu.org/licenses/gpl-2.0.html
#
# $Revision: 1735 $
#
# Usage: ./picsord.R sample.fasta
#
################################################################################
## Configurable Variables

# the path to ngila
ngila <- "ngila"

# the arguments to ngila
ngilacmds <- "-m zeta -o dist-c:-"

# do not change unless you know what you are doing
ngilacmds2 <- "--pairs=all" 

# number of bins range: 2-64
nbins <- 10

# remove constant columns
rmconst <- TRUE

# quantization method
quant <- "default"

# maximum length of sequence names
seqlenmax <- 100

################################################################################
## Program
options(warn = -1)

cargs <- commandArgs(trailingOnly = TRUE)
# read the name of the fasta file
aln <- cargs[1]

# print usage if needed
if(is.na(aln)) {
	cat("
PICS-Ord: Pairwise Identity and Cost Scores Ordination
Copyright (c) 2010 Reed A. Cartwright, PhD <reed@scit.us>
with Robert Luecking, Brendan Hodkinson, and Alexandros Stamatakis

Usage: Rscript picsord.R input.fasta > output.phy
To change parameters, simply edit the file.

")
	quit()
}

# character encoding information
b64 <- function(x) {
	v <- c(as.character(0:9),
		LETTERS[1:26], letters[1:26], "%","@")
	v[x%%64+1]
}

# run ngila and fetch output
tab <- scan(pipe(paste(ngila, ngilacmds, ngilacmds2, aln)), quiet=TRUE,
	what=character(0))
# massage output into an R matrix
n <- 0.5*(sqrt(1+4*length(tab))-1)
names <- tab[1:n]
tab <- matrix(as.numeric(tab[-(1:n)]), n, n)

# md scaling
res <- cmdscale(tab, k=n-1,eig=TRUE)
mm <- res$points[,which(res$eig > 0)]

# quantization method
if(quant == "kmeans") {
	# build clusters using kmeans
	k <- kmeans(c(mm),nbins,nstart=10)
	o <- order(order(k$centers))
	m <- matrix(o[k$cluster]-1,nrow=n)
} else {
	# rescale resulting points
	m <- t(mm)
	smin <- apply(m,1,min)
	smax <- apply(m,1,max)
	srm <- max(smax-smin)
	m <- t((m-smin)/srm)

	# create characters
	#m <- round((m*0.999 - 0.04995)*nbins)
	m <- trunc(nbins*m*(1-.Machine$double.eps/2))
}

if(rmconst) {
	smin <- apply(m,2,min)
	smax <- apply(m,2,max)
	m <- m[,smin != smax,drop=FALSE]
}
m <- apply(m,2,b64)

# find and format otu names
seqlenmin <- 10
names <- sub("^\\s+", "", names)
names <- sub("\\s+$", "", names)
names <- sub("\\s+", "_", names)
name.width <- max(nchar(names),seqlenmin-1)
name.width <- min(name.width,seqlenmax-1)
name.width <- sprintf("%% -%d.%ds", name.width, name.width)
names <- sprintf(name.width, names)

# print character matrix in phylip format
cat(sprintf("    %d    %d\n", n, ncol(m)))
mm <- apply(m,1,function(x) paste(x,collapse=""))
cat(paste(names,mm,collapse="\n"))
cat("\n")
