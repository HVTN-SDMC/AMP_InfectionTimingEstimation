# I, Paul Edlefsen, obtained this July 21, 2020 directly from Elena Giorgi at lanl. I then merged in my edits to the previous version that I had worked on. The newer version accepts a fourth argument as the sample name, so I moved my acceptance of the prior name as the fourth argument to the fifth.
# You should be able to re-obtain the latest version by doing what I did the first time: I went to
# http://www.hiv.lanl.gov/content/sequence/POISSON_FITTER/pfitter.html,
# selected the [sample input] link, which put text in the textbox
# labeled "Paste your alignment(s) here", and then I hit submit
# without changing anything from its default value.  On the resulting
# page, which is called
# http://www.hiv.lanl.gov/cgi-bin/POISSON_FITTER/v2/pfitter.cgi (but
# you can't just go there circumventing the form, or you will just see
# an error message), I followed link 'Download R Script'. 

#########################################################################
##########       PURPOSE: fit Poisson, calculate lambdas,      ##########
##########    U-stat standard deviation, and goodness of fit   ##########
##########      written by EEG, last modified on 5/26/15       ##########
##########          send questions to egiorgi@lanl.gov         ##########

# INPUT: 
#  pairwise hamming distances file, mutation rate, length of a sequence, sample_name
#  distances file: tab-delimited 3-column file. seqname1(1st col), 
#                  seqname2(2nd) and distance between seq1 and seq2(3rd).
#                  based on large-scale formatted sequnce input, which 
#                  means every unique sequence is represented only once 
#                  and a seqname should end with _nnn wher nnn is the 
#                  the multiplicity of such sequence in the alignment.

#  example:  R CMD BATCH '--vanilla --args sample.dist  2.16e-05 2517 1018' this_script 

# OUTPUT:
#  2 files, one with lambdas, maxhd, jackknife standard deviation 
#  and estimated days, goodness of fit p-values, and the other with the 
#  star-phylogeny estimated an dobserved numbers (if they coincide, you have a star)
## NEW BY PAUL: if a fifth argument is present and is "gamma00" or "gamma11", output Bayesian credible interval too using the given prior eg Gamma( 0, 0 ).
#########################################################################


library(tools)
args = commandArgs(trailingOnly=TRUE)

infile <- args[1]
epsilon <- c(as.numeric(args[2]))
nbases <- c(as.numeric(args[3]))

sample  <- args[4]
prefix  <- gsub("[^0-9A-Za-z]","_" , sample ,ignore.case = TRUE)  # replace all chars other than alphanumerics. outfile prefix

## NEW BY PAUL: fifth argument is the prior spec; if present, output Bayesian credible interval too.
if( length( args > 4 ) ) {
    prior <- args[5];
} else {
    prior <- NULL;
}

prior.gamma.pseudocounts <-
    list(
         "gamma00" = c( shape = 0, rate = 0 ),
         "gamma11" = c( shape = 1, rate = 1 )
        );


dir <- paste(dirname(infile), '/', sep='')
outfile <- paste(dir, prefix,  "_LOG_LIKELIHOOD.txt", sep="")
outfile2 <- paste(dir, prefix,  "_CONVOLUTION.txt", sep="")

## ADDED BY PAUL:
if( !is.null( prior ) ) {
    if( prior %in% names( prior.gamma.pseudocounts ) ) {
        # good.
        gamma.prior.shape <- prior.gamma.pseudocounts[[ prior ]][ "shape" ];
        gamma.prior.rate <- prior.gamma.pseudocounts[[ prior ]][ "rate" ];
        print( paste( "Including Bayesian credibility interval using Gamma prior: dgamma( ., shape = ", gamma.prior.shape, ", rate = ", gamma.prior.rate, ")", sep = "" ) );
    } else {
        stop( paste( "ERROR: UNRECOGNIZED PRIOR `", prior, "'", sep = "" ) )
    }
    outfile3 <- paste(dir, prefix, "_BAYESIAN_", prior, ".txt", sep="")
}

### FUNCTIONS ###
iseven <- function(c) {
	c1 <- c/2-as.integer(c/2)
	if(c1==0){ ev <- TRUE } else { ev <- FALSE }
	return(ev)
}
		
phi <- sqrt(1+4/3) # \sqrt{ 1 + \frac{8}{R_0} }
#gens <- function(l,nb,epsilon) (l/(nb*epsilon)-(1-phi)/(phi^2))*((phi)/(1+phi))
days <- function(l,nb,epsilon) 1.5*((phi)/(1+phi))*(l/(epsilon*nb) - (1-phi)/(phi^2))
### end FUNCTIONS ###

write(paste("Sample", "Lambda", "St.Dev", "NSeq", "NBases", "LengthYvec0", "SumYvec0", "MeanHD", "MaxHD","Days(CI)", "Chi2","DF","Goodness_of_pval", sep="\t"), file=outfile, append=FALSE)
	
dlist <- read.table(file=infile, sep="\t", stringsAsFactors=F, quote="",comment.char = "")   


### calc HD with consensus
d0 <- dlist[which(dlist[,1]==dlist[1,1]),]
mult0 <- as.numeric(sub('.+_(\\d+)$', '\\1', d0[,2]))
nseq <- sum(mult0)
yvec0 <- rep(0, (1+max(d0[,3])))
for(i in 1:(1+max(d0[,3]))){ yvec0[i] <- sum(mult0[which(d0[,3]==(i-1))]) }

nl0 <- length(yvec0)
clambda <- sum((1:(nl0-1))*yvec0[-1])/sum(yvec0) #### THIS IS THE LAMBDA THAT FITS THE CONSENSUS ONLY DISTRIBUTION
sumyvec0 <- sum(yvec0)

### calc intersequence HD
d1 <- dlist[-which(dlist[,1]==dlist[1,1]),]	
yvec <- rep(0, (1+max(d1[,3])))
## PTE July 2020 to fix bug in which all entries in the deduplicated fasta file that is given to the hamming_dist.pl script have > 1 multiplicity, or are not in descending order by multiplicity; the previous code ignored the last sequence name in computing the yvec[0] count so undercounted if it had multiplicity >1. This fixes that bug that if the last entry has > 1 multiplicity it was ignored, and it is also now robust to the order of the sequences in the input fasta file.
seqnames <- unique(c( d1[,1], d1[,2]))
# PTE July 2020 seqnames <- unique(d1[,1])
for(i in 1:length(seqnames)){
	m0 <- as.numeric(sub('.+_(\\d+)$', '\\1', seqnames[i]))
	yvec[1] <- yvec[1] + 0.5*m0*(m0-1) ## 0 bin
        # PTE July 2020
        if( !any(d1[,1]==seqnames[i]) ) {
            next;
        }
	tmp <- d1[which(d1[,1]==seqnames[i]),]
	for(j in 1:dim(tmp)[1]){
		m1 <- as.numeric(sub('.+_(\\d+)$', '\\1', tmp[j,2]))
		val <- tmp[j,3]
		yvec[val+1] <- yvec[val+1] + m0*m1
	}
}

### Fitting

nl <- length(yvec)
lambda <- sum((1:(nl-1))*yvec[-1])/sum(yvec)
estdays <- days(lambda, nbases, epsilon)
		
print(paste("Estimated Lambda", format(lambda, digits=4), sep=" "))
#print(paste("Estimated Lambda based on consensus HD", format(clambda, digits=4), sep=" "))



#### U STAT ESTIMATE OF ST DEV
#### FORMULAE
#### Var(HD) = (N(N-1)/2)^(-1) (2(N-2)sigma1^2 + sigma2^2)
#### sigma1^2 = (N(N-1)(N-2)/3 -1)^(-1) sum_{i<j<l} ((Dij-mu)(Dil-mu)+(Dij-mu)(Djl-mu))
#### sigma2^2 = (N(N-1)/2-1)^(-1) sum_{i<j} (Dij-mu)^2

### construct a matrix of Dij's
### number of unique sequences
nuni <- dim(d0)[1]
TX <- matrix(rep(0,nuni^2), ncol=nuni)

for(i in 1:(dim(d0)[1]-1)){
	useq <- d0[i,2]
	TX[((i+1):dim(TX)[1]),i] <- d1[which(d1[,1]==useq),3]
}

sigma1 <- 0
sigma2 <- 0
muhat <- 0
denmu <- (nseq*(nseq-1)/2)^(-1)
den1 <- 12*(nseq*(nseq-1)*(nseq-2)*(nseq-3))^(-1)  
den2 <- den1/4   


for(n in 1:(nuni-1)){
	for(m in (n+1):nuni){
		muhat <- muhat + mult0[n]*mult0[m]*denmu*TX[m,n]
	}
}

for(n in 1:nuni){
	dnn <- 0
	sigma1 <- sigma1 + choose(mult0[n],3)*den1*2*(dnn-muhat)^2 
	sigma2 <- sigma2 + choose(mult0[n],2)*den2*(dnn-muhat)^2
	if(n != nuni){
		for(m in (n+1):nuni){
			dnm <- TX[m,n]
			dmm <- 0
			sigma2 <- sigma2 + mult0[n]*mult0[m]*(dnm - muhat)^2
			sigma1 <- sigma1 + (2/3)*choose(mult0[n],2)*mult0[m]*(dnm-muhat)*(dnm+2*dnn-3*muhat)
			sigma1 <- sigma1 + (2/3)*mult0[n]*choose(mult0[m],2)*(dnm-muhat)*(dnm+2*dmm-3*muhat)
			if(m != nuni){
				for(l in (m+1):nuni){
					dnl <- TX[l,n]
					dlm <- TX[l,m]
					sigma1 <- sigma1 + (2/3)*mult0[n]*mult0[m]*mult0[l]*((dnm-muhat)*(dnl-muhat)+(dnm-muhat)*(dlm-muhat)+(dnl-muhat)*(dlm-muhat)) 
				}
			}
		}
 	}
}

## varhd <- sqrt(denmu*(2*(nseq-2)*sigma1 + sigma2))
A <- 8/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
B <- 4/(nseq*(nseq-1)*(nseq-2)*(nseq-3))
newvarhd <- sqrt(A*sigma1 + B*sigma2)
	
upplim <- days(lambda + 1.96*newvarhd, nbases, epsilon)
lowlim <- days(lambda - 1.96*newvarhd, nbases, epsilon)
uppdays <- round(upplim)
lowdays <- round(lowlim)
	
formatteddays <- paste(round(estdays), " (", lowdays, ", ", uppdays, ")", sep="") 

### output figures
dvec1 <- 0
for(i in 1:length(yvec)){ dvec1 <- c(dvec1, rep((i-1),yvec[i])) }
dvec1 <- dvec1[-1]

histdatafile <- paste(dir, prefix, "_hd.txt",sep="")
write('HDIST', file=histdatafile)
write.table(dvec1, file=histdatafile, append=T, row.names=F, col.names=F)

meanhd <- mean(dvec1)
maxhd <- max(dvec1)

figure <- paste(dir, prefix, "_hd_freq_dist.jpg", sep="")
jpeg(file=figure, pointsize=14, width=800, height=800)
par(mar=c(5,5,4,2))

layout(matrix(c(1,2,3), nr=1, ncol=3), widths=c(5), heights=c(1), FALSE)
h <- hist(dvec1, breaks=seq(-1,max(dvec1),1), plot=FALSE)
hist(dvec1, breaks=seq(-0.5,max(dvec1)+0.5,1), freq=TRUE, xlab=sample, ylim = c(0,15+max(h$counts)), labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Hamming Distance Frequency Distribution", sep=" "), cex.lab=2, cex.axis=2,cex.main=2)
lines(seq(0,max(dvec1)+1,1), 0.5*nseq*(nseq-1)*dpois(seq(0,1+max(dvec1),1), lambda=lambda),col="red", lwd=2)

dev.off()	
	
figure <- paste(dir, prefix, "_hd_freq_dist.ps", sep="")
postscript(file=figure, width=400, height=300)
par(mar=c(5,5,4,3))

h <- hist(dvec1, breaks=seq(-1,max(dvec1),1), plot=FALSE)
hist(dvec1, breaks=seq(-0.5,max(dvec1)+0.5,1), freq=TRUE, xlab=sample, ylim = c(0,15+max(h$counts)),
labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Hamming Distance Frequency Distribution", sep=" "), cex.lab=2, cex.axis=2,cex.main=2)
lines(seq(0,max(dvec1)+1,1), 0.5*nseq*(nseq-1)*dpois(seq(0,1+max(dvec1),1), lambda=lambda),col="red", lwd=2)

dev.off()

figurelinedata <- paste(dir, prefix, "_hd_freq_dist.txt",sep="")
m<-cbind(seq(0,max(dvec1)+1,1), 0.5*nseq*(nseq-1)*dpois(seq(0,1+max(dvec1),1), lambda=lambda))
d<-as.data.frame(m)
names(d) <- c('HDIST', 'FREQ')
write.table(d, file=figurelinedata, row.names=F,  quote=F, sep="\t")

#### FIT THE CONSENSUS ONLY HD DISTRIBUTION
	
if (lambda!=0) {
			
	xvec1 <- c(yvec0, rep(0,nl0))
	yvec1 <- rep(0,2*nl0)
	yvec1[1] <- 1/2*yvec0[1]*(yvec0[1]-1)  ### freq at zero 
	mvals <- seq(2,2*nl0,2)

	for(m in mvals) {
		delta <- rep(0,m)
		delta[1+m/2] <- 1
		for(hj in 1:m){			
			yvec1[m] <- yvec1[m] + 1/2*xvec1[hj]*xvec1[m-hj+1]
                        ## PAUL NOTES: This appears to be a bug; I believe that "nl" should be "nl0" here. The former is the length of the intersequence HD histogram vector (yvec), the latter is the corresponding length for the distances-to-consensus histogram vector (yvec0).  The "nl" two lines down should also be "nl0" I believe.
			if(m<2*nl) { yvec1[m+1] <- yvec1[m+1] + 1/2*xvec1[hj]*(xvec1[m-hj+2] - delta[hj]) } 
		}
		if(m<2*nl) { yvec1[m+1] <- yvec1[m+1] + 1/2*xvec1[m+1]*xvec1[1] }
	}
	
	dvec2 <- rep(0, yvec1[1])
	w <- which(yvec1>0)	
	for(hk in w[-1]) { dvec2 <- c(dvec2, rep((hk-1), yvec1[hk])) }
	
	mmax <- 1.5*(max(c(yvec0, yvec1)))			

	cfigure <- paste(dir, prefix, "_conv_plot.jpg", sep="")
	jpeg(file=cfigure, pointsize=14, width=800, height=800)
	par(mar=c(5,5,4,2))

	layout(matrix(c(1,2,3), nr=1, ncol=3), widths=c(5), heights=c(1), FALSE)
	hh1 <- hist(dvec1, breaks=seq(-0.5,max(dvec1+0.5),1), freq=TRUE, xlab=sample, labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Convolution Plot", sep=" "), col="blue",cex.lab=2, cex.axis=2,cex.main=2)
	lines(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], col="red", lwd=2)
	points(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], pch=23, col="red", lwd=2)
	legend("topright", legend=c("CONV","OBS"), fill=c("red","blue"), text.width=0.8)

	dev.off()
			
	figure <- paste(dir, prefix, "_conv_plot.ps", sep="")
	postscript(file=figure, width=400, height=300)
	par(mar=c(5,5,4,3))
	
	hh1 <- hist(dvec1, breaks=seq(-0.5,max(dvec1+0.5),1), freq=TRUE, xlab=sample, labels=ifelse(max(dvec1)<20,TRUE,FALSE), main=paste("Convolution Plot", sep=" "), col="blue",cex.lab=2, cex.axis=2,cex.main=2)
	lines(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], col="red", lwd=2)
	points(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))], pch=23, col="red", lwd=2)
	legend("topright", legend=c("CONV","OBS"), fill=c("red","blue"), text.width=0.8)

	dev.off()
	
    figurelinedata <- paste(dir, prefix, "_conv_plot.txt",sep="")
    m<-cbind(seq(0,(2*nl0-1),1), yvec1[1:length(seq(0,(2*nl0-1),1))])
    d<-as.data.frame(m)
    names(d) <- c('HDIST', 'FREQ')
    write.table(d, file=figurelinedata, row.names=F,  quote=F, sep="\t")

	check <- 0 

	write(sample, file=outfile2, append=TRUE)
	write(paste("HD", "OBS", "CONV", sep="\t"), file=outfile2, append=TRUE)
	for(jj in 1:length(yvec1)){ 
		write(paste(jj-1, hh1$counts[jj], yvec1[jj], sep="\t"), file=outfile2, append=TRUE) 
		hey <- abs(hh1$counts[jj]-yvec1[jj])
		if(!is.na(hey)) { check <- check + hey }
	}

	check <- check/sum(yvec1, na.rm=T)
	ifclause <- ifelse(check <= 0.1, "FOLLOWS", "DOES NOT FOLLOW")
	astring <- paste(sample, ifclause, "A STAR-PHYLOGENY", sep=" ")

	write(astring, file=outfile2, append=TRUE)
	write(" ", file=outfile2, append=TRUE)	
}
		

### CONSTRUCT SIGMA_ij MATRIX THEN INVERT IT
pk <- function(x) ((nseq^2)*(2^x)*exp(-2*clambda)*(clambda^x))/factorial(x)
mui <- function(x) nseq*dpois(x, lambda=clambda)
eyvec <- 0.5*pk(0:(2*nl0-1))

if (lambda!=0) {
	
	sigmaij <- matrix(nrow=(2*nl0), ncol=(2*nl0))
	coeff <- (nseq^3)*exp(-3*clambda)
	
	#### RICORDATI!!!! EYVEC[K] == E(Y_{K-1}) !!!!!!
	
	for(k in 0:(2*nl0-1)){  
  
		for(l in 0:(2*nl0-1)){   
			
			if(k>=l){ 
				c1 <- ((clambda^k)/factorial(k))*sum(choose(k,l:0)*((clambda^(0:l))/factorial(0:l)))
				c2 <- ((clambda^l)/factorial(l))*sum(choose(l,l:0)*((clambda^((k-l):k))/factorial((k-l):k))) 
			}
			if(k<l){
				c1 <- ((clambda^l)/factorial(l))*sum(choose(l,k:0)*((clambda^(0:k))/factorial(0:k)))
				c2 <- ((clambda^k)/factorial(k))*sum(choose(k,k:0)*((clambda^((l-k):l))/factorial((l-k):l))) 
			}
						
			sigmaij[k+1,l+1] <- 0.5*coeff*(c1+c2)
			
			if(k==l){ sigmaij[k+1,l+1] <- sigmaij[k+1,l+1] + (0.5)*pk(k) }
			if((k==l)&(iseven(k))){ sigmaij[k+1,l+1] <- sigmaij[k+1,l+1] - (0.25)*mui(k/2) }
		}
	}
									
	sdec <- La.svd(sigmaij)
	diag <- ifelse(sdec$d>1e-4,sdec$d,0)
	diagmat <- matrix(rep(0,(2*nl0)^2), ncol=2*nl0)
	for(ii in 1:(2*nl0)){diagmat[ii,ii]<-ifelse(diag[ii]==0,0,1/diag[ii])}
	sigmainv <- sdec$u%*%diagmat%*%sdec$vt
	
	h <- hist(dvec1, breaks=seq(-1,max(dvec1),1), plot=FALSE)
	xvec <- h$breaks
	yvec <- h$counts
	nl1 <- length(yvec)
	aplambda <- sum((1:(nl1-1))*yvec[-1])/sum(yvec)
			
	pesce <- 0.5*nseq*(nseq-1)*dpois(0:(nl1-1), lambda=aplambda)

	if (length(yvec)<2*nl0) { 
		ccvv <- 2*nl0 - length(yvec) 
		yvec <- c(yvec, rep(0,ccvv)) 
	}
	
	chisq <- t(abs(yvec-eyvec))%*%sigmainv%*%(abs(yvec-eyvec))
	pval <- ifelse(chisq<0,2e-16,1-pchisq(chisq,df=nl0-1))		
	if(pval==0){ pval <- 2e-16 }
	if(chisq<0){ chisq <- NA }
} else { 
	chisq <- NA
	nl <- NA
	pval <- 0
}
       
write(paste(sample, format(lambda, digits=4), format(newvarhd, digits=4), nseq, nbases, nl0, sumyvec0, format(meanhd, digits=2), maxhd, formatteddays, format(as.numeric(chisq), digits=4), nl0-1, format(as.numeric(pval), digits=4), sep="\t"), file=outfile, append=TRUE)

## New Bayesian code by PAUL    
if( !is.null( prior ) ) {

    ## This file has to be sourced inline, here. It depends on eg "phi" defined in this script.
    source( "PFitter_Bayesian_functions.R" )

    if(!file.exists(outfile3)){
	write(paste("Sample", "CLambda", "Interval_width", "NSeq", "NBases", "k", "n","Days(CI)", sep="\t"), file=outfile3, append=FALSE)
    }

    dvec0 <- 0
    for(i in 1:length(yvec0)){ dvec0 <- c(dvec0, rep((i-1),yvec0[i])) }
    dvec0 <- dvec0[-1]

    # These are the total number of cells with a mutation, and the total number of cells.
    k <- sum( dvec0 );
    n <- nbases * length( dvec0 );

    meanhd0 <- mean(dvec0)
    maxhd0 <- max(dvec0)

    credible.interval.and.posterior.median <-
        compute.poisson.desired.quantiles.results.recursively( k, n, prior.gamma.shape = prior.gamma.pseudocounts[[ prior ]][ "shape" ], prior.gamma.rate = prior.gamma.pseudocounts[[ prior ]][ "rate" ] );

    ## This is the posterior median.
    posterior.estdays0 <- credible.interval.and.posterior.median[ "0.5" ] 

    posterior.upplim0 <- credible.interval.and.posterior.median[ "0.975" ]
    posterior.lowlim0 <- credible.interval.and.posterior.median[ "0.025" ]
    posterior.interval.width <- unname( posterior.upplim0 - posterior.lowlim0 );
    posterior.uppdays0 <- round(posterior.upplim0)
    posterior.lowdays0 <- round(posterior.lowlim0)
    	
    posterior.formatteddays0 <- paste(round(posterior.estdays0), " (", posterior.lowdays0, ", ", posterior.uppdays0, ")", sep="") 
    
    write(paste(sample, format(clambda, digits=4), format(posterior.interval.width, digits=4), nseq, nbases, k, n, posterior.formatteddays0, sep="\t"), file=outfile3, append=TRUE)
    
} # End if( !is.null( prior ) )







