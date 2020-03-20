circularize <- function(x, len) {
	if (length(x) < len) {
		x <- rep(x, length.out=len)
	}
	x
}

#plot densities from columns of a matrix (or variables)
plot.densities <- function(m, m2, m3, m4, lty=rep(1:6, each=length(palette())), lwd=1, cex=1, legends=TRUE, legend.pos="topleft", legend.cols=1, legend.cex=1, horiz=FALSE, xlim=if(horiz) c(mny,mxy) else c(mnx,mxx),ylim=if(!horiz) c(mny,mxy) else c(mnx,mxx), pch=NULL, type="l", npch=20, col=1:length(d), normalize=FALSE, q=NULL, overall=length(d) > 4, ylab="Density", xlab="", main="", ...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- as.vector(m)
		if (!missing(m2)) { m <- cbind(m,as.vector(m2)); xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m <- cbind(m,as.vector(m3)); xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m <- cbind(m,as.vector(m4)); xn <- c(xn,deparse(substitute(m4))) }
		colnames(m) <- xn
	}
	else if (!is.list(m)){
		m <- as.matrix(m)
	}
	if (is.list(m)) d <- lapply(m, function(x) { suppressWarnings(density(x[is.finite(x)],...)) })
	else d <- apply(m, 2, function(x) { suppressWarnings(density(x[is.finite(x)],...)) })
	lty <- circularize(lty, length(d))
	lwd <- circularize(lwd, length(d))
	cex <- circularize(cex, length(d))
	if (!is.null(q)) {
		if (is.list(m)) qq <- lapply(m, function(x) { suppressWarnings(quantile(x[is.finite(x)], q)) })
		else qq <- apply(m, 2, function(x) { suppressWarnings(quantile(x[is.finite(x)], q)) })
	}
	if (overall) {
		o <- unlist(m)
		n <- names(d)
		d[[length(d)+1]] <- density(o[is.finite(o)], ...)
		names(d) <- c(n, "OVERALL")
		lty[length(d)] <- 1
		lwd[length(d)] <- 3
		col[length(d)] <- 1
		cex[length(d)] <- 1
		if (!is.null(q)) qq[length(d)] <- suppressWarnings(quantile(o[is.finite(o)], q))
	}
	if (normalize) {
		for (i in 1:length(d)) {
			d[[i]]$y <- d[[i]]$y/max(d[[i]]$y)
		}
	}
	mxy <- max(unlist(lapply(d, function(x) max(x$y))))
	mny <- min(unlist(lapply(d, function(x) min(x$y))))
	mxx <- max(unlist(lapply(d, function(x) max(x$x))))
	mnx <- min(unlist(lapply(d, function(x) min(x$x))))
	lty <- rep(lty, length(d), length.out=length(d))
	pch <- ifelse(type=="l",rep(0,length(d)),rep(pch, length(d), length.out=length(d)))
	suppressWarnings(plot(0,0,type="n",ylim=ylim,xlim=xlim, main=main, xlab=xlab, ylab=ylab, ...))
	cols <- if (horiz) c("y","x") else c("x","y")
	for (i in length(d):1) {
		suppressWarnings(lines(d[[i]][cols], col=col[i], lty=lty[i], pch=pch[i], type=type, lwd=lwd[i], cex=cex[i]))
		if (!is.null(q)) {
			for (k in 1:length(q)) {
				if (horiz) abline(h=qq[[i]][k], lty=k, col=col[i])
				else abline(v=qq[[i]][k], lty=k, col=col[i])
			}
		}
	}
	if (type != "l") {
		pch <- rep(pch, length(d), length.out=length(d))
		if (horiz) for (i in 1:length(d)) suppressWarnings(points(d[[i]]$y[seq(1,length(d[[i]]$y), length.out=npch)], d[[i]]$x[seq(1,length(d[[i]]$y), length.out=npch)], col=col[i], pch=pch[i], ...))
		else	   for (i in 1:length(d)) suppressWarnings(points(d[[i]]$x[[seq(1,length(d[[i]]$y), length.out=npch)]], d[[i]]$y[seq(1,length(d[[i]]$y), length.out=npch)], col=col[i], pch=pch[i], ...))
	}
	if (length(legends) > 1 || (!is.character(legends) && legends)) {
		if ((length(legends) == length(d)) && (length(legends) > 0)) xl <- legends
		else {
			xl <- if (is.list(m)) names(m) else colnames(m)
			if (overall && !is.null(xl)) xl <- c(xl, "OVERALL")
		}
		if (!is.null(xl)) {
			plot.legend(legend.pos, xl, col=col[1:length(d)], lty=lty, ncol=legend.cols, cex=legend.cex, pch=if (is.null(pch[1])) NULL else pch)
		}
	}
	invisible(d)
}

my.quantile <- function(v, probs, na.rm=FALSE) {
	s <- sort(v, na.last=ifelse(na.rm,NA,TRUE))
	l <- length(s)
	s[pmax(pmin(l,probs*l),1)]
}

#center - median centering before quantile (almost no effect, just create an offset)
#qmin - scale the distribution forcing the minimum to be qmin shriniking or stretching the distribution
#qmax - scale the distribution forcing the maximum to be qmax shriniking or stretching the distribution
#qmid - scale the centered qmid part of the distribution to be between qmin and qmax shriniking or stretching the distribution
#qshrink - scale every element to have the centered qshrink part of the distribution equal before estimating the averages (helps to smooth eventual wild distributions)
quantile.normalization <- function(l.mx, verbose=FALSE, averages=FALSE, center=FALSE, q.func=my.quantile, m.func=mean, qmin=min(avg), qmax=max(avg), qmid=1, qshrink=NULL) {
	#q.func=quantile	# if you want to use the embed slower R quantile (slightly different results)
	
	#convert to list
	retfn <- function(x) x
	if (is.data.frame(l.mx) || is.matrix(l.mx)) { 
		.rn <- rownames(l.mx)
		l.mx <- lapply(as.list(data.frame(l.mx)), function(x) { names(x) <- .rn; x; })
		if (is.data.frame(l.mx))
			retfn <- function(x) { x <- data.frame(x); rownames(x) <- .rn; x; }
		else
			retfn <- function(x) { x <- data.matrix(data.frame(x)); rownames(x) <- .rn; x; }
	}

	#scale values to min~max range, or to part of the median/center quantile (qmid)
	q.scale <- function(v, qmin=min(v), qmax=max(v), qmid=1) {
		qs <- quantile(v, p=c((1-qmid)/2,(qmid+1)/2))
		vmin <- qs[1]
		vmax <- qs[2]
		#3-rule eq: (v - vmin) / (vmax - vmin) = (q - qmin) / (qmax - qmin)
		((v - vmin) / (vmax - vmin)) * (qmax - qmin) + qmin
	}

	#median centering
	if (center) {
		if (verbose) cat("Centering...\n")
		# center to 
		md <- unlist(lapply(l.mx, median, na.rm = TRUE))
		mmd <- median(md)
		l.mx <- lapply(l.mx, function(l) l-(median(l,na.rm=TRUE) - mmd))
	}

	#shrinking distributions to median/center quantile range
	if (!is.null(qshrink)) {
		if (verbose) cat("Shrinking...\n")
		minq <- median(unlist(lapply(l.mx, quantile, p=(1-qshrink)/2, na.rm = TRUE)))
		maxq <- median(unlist(lapply(l.mx, quantile, p=(qshrink+1)/2, na.rm = TRUE)))
		l.mx <- lapply(l.mx, q.scale, qmin=minq, qmax=maxq, qmid=qshrink)
	}

	#maximum number of values
	m <- 0 
	for (i in 1:length(l.mx)) m <- max(m,sum(is.finite(l.mx[[i]])))
	if (verbose) cat("Max-Quantiles:",m,"\n")

	#quantiling
	q.mx <- matrix(NA,nrow=m, ncol=length(l.mx))
	for (i in 1:length(l.mx)) { q.mx[,i] <- q.func(l.mx[[i]],probs=(1:m)/m,na.rm=TRUE); if (verbose) cat("Quantiled:",i,"\n"); }

	#averaging
	avg <- numeric(length=m)
	if (verbose) cat("Averaging...\n")
	for (i in 1:m) avg[i] <- m.func(q.mx[i,], na.rm=TRUE) 

	#scale averages?
	if (min(avg) != qmin  ||  max(avg) != qmax  ||  qmid < 1) avg <- q.scale(avg, qmin=qmin, qmax=qmax, qmid=qmid)

	#normalizing
	q.mx <- list()
	if (verbose) cat("Normalizing...\n")
	for (i in 1:length(l.mx)) { 
		if (verbose) cat("Computing ecdf",i,"...\n");
		avg.cdf <- ecdf(l.mx[[i]][is.finite(l.mx[[i]])])
		r <- pmin(pmax(1,avg.cdf(l.mx[[i]])*m),m)
		xq <- avg[r]
		names(xq) <- names(l.mx[[i]])
		q.mx[[i]] <- xq
	}
	names(q.mx) <- names(l.mx)

	#return
	if (verbose) cat("Done!\n");
	if (averages) avg else retfn(q.mx)
}

#differential expression testER - de.test
#ref: selecting genes by test statistics Chen Hua J Biomed Biotech 2005
#ref:Comparison and evaluation of methods for generating differentially expressed gene lists from microarray data - Jeffrey Culhane - BMC Bioinf 2006
#ref: Comparison of various statistical methods for identifying differential gene expression in replicated microarray data Kim Lee Stat Met Med Res 2006
#x - matrix or data frame
#classes
de.test <- function(x, classes, test=c("all", "ttest", "kolmogorov", "kruskal", "ftest", "welch", "cochran", "wilcoxon", "copa", "s2n", "osum", "sam", "infcont", "yeohchi2", "pca", "affinity", "slide","bimodal","brown","perm1"), nperms=100, ret.all = FALSE, debug=FALSE, except=c("affinity","bimodal","brown"), ...) {
	
	if (!is.factor(classes)) classes <- factor(classes)	
	tx <- t(x)
	nl <- nlevels(classes)
	l <- levels(classes)
	m <- matrix(Inf, nrow = nrow(x), ncol=nl)
	colnames(m) <- l	
	p <- list()
	if (("ttest" %in% test  ||  test == "all")  &&  ! "ttest" %in% except) {
		if (debug) cat("ttest...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) t.test(xwi[,k], xwj[,k])$p.value)
		}
		p[["ttest"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("bimodal" %in% test  ||  test == "all")  &&  ! "bimodal" %in% except) {
		if (debug) cat("bimodal...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) p.bimodal(xwi[,k]))
		}
		m <- cbind(m, apply(t(x), 2, p.bimodal))
		p[["bimodal"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("infcont" %in% test  ||  test == "all") &&  ! "infcont" %in% except) {
		if (debug) cat("infcont...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			m[,i] <- p.information.content(x, classes, l[i], large=FALSE)
		}
		p[["infcont"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("sliding" %in% test  ||  test == "all")  &&  ! "sliding" %in% except) {
		if (debug) cat("sliding...\n")
		uclasses <- unclass(classes)
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			.s <- sum(uclasses == i)
			.v <- numeric(nrow(x))
			for (j in 1:nrow(x)) {
				a <- order(x[j,])
				ac <- uclasses[a]
				.v[j] <- min(sapply(1:(ncol(x)-.s+1), function(k) {
					p.overlap(ncol(x), .s, .s, sum(ac[1:.s+k-1] == i))
				}))
			}
			m[,i] <- .v
		}
		p[["sliding"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("yeohchi2" %in% test  ||  test == "all")  &&  ! "yeohchi2" %in% except) {
		if (debug) cat("yeohchi2...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			m[,i] <- p.yeohchi2(x, classes, l[i])
		}
		p[["yeohchi2"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("kolmogorov" %in% test  ||  test == "all")  &&  ! "kolmogorov" %in% except) {
		if (debug) cat("kolmogorov...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) ks.test(xwi[,k], xwj[,k])$p.value)
		}
		p[["kolmogorov"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("kruskal" %in% test  ||  test == "all")  &&  ! "kruskal" %in% except) {
		if (debug) cat("kruskal...\n")
		p[["kruskal"]] <- apply(tx, 2, function(x) kruskal.test(x, classes)$p.value)
	}
	if ("ftest" %in% test  ||  test == "all") {
		if (debug) cat("ftest...\n")
		p[["ftest"]] <- apply(tx, 2, function(x) {
			sm <- summary(lm(x~classes))$fstatistic
			pf(sm[1], sm[2], sm[3], lower.tail=FALSE)
		} )
	#	if (debug) cat("ftest...\n")
	#	p[["ftest"]] <- apply(tx, 2, function(x) { tdf <- data.frame(x=x, g=classes); oneway.test(x ~ g, tdf, var.equal=TRUE)
	}
	if (("wilcoxon" %in% test  ||  test == "all")  &&  ! "wilcoxon" %in% except) {
		if (debug) cat("wilcoxon...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) wilcox.test(xwi[,k], xwj[,k])$p.value)
		}
		p[["wilcoxon"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("s2n" %in% test  ||  test == "all")  &&  ! "s2n" %in% except) {
		if (debug) cat("s2n...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) (mean(xwi[,k],na.rm=TRUE)-mean(xwj[,k],na.rm=TRUE)) / (sd(xwi[,k],na.rm=TRUE)+sd(xwj[,k],na.rm=TRUE)))
		}
		m <- glog(abs(m))
		m <- apply(m, 2, function(x) 1 - x / max(x,na.rm=TRUE))  # transform to a p-value like
		p[["s2n"]] <- if (ret.all) m else apply(m,1,min,na.rm=TRUE)
	}
	if (("welch" %in% test  ||  "oneway" %in% test  ||  test == "all")   &&  ! "welch" %in% except) {
		if (debug) cat("welch...\n")
		#xdf <- as.data.frame(tx)
		#xdf[,"testcluster"] <- classes
		#print(colnames(xdf))
		tdf <- data.frame(x=tx[,1], g=classes)
		p[["welch"]] <- apply(tx, 2, function(x) { tdf$x <- x; oneway.test(x ~ g, tdf, var.equal=FALSE)$p.value })
	}
	if (("brown" %in% test  ||  "brown-forsythe" %in% test  ||  
		  "levene"  %in% test  || test == "all")   
		  &&  ! "brown" %in% except) {
		if (debug) cat("brown-forsythe/levene...\n")
		#xdf <- as.data.frame(tx)
		#xdf[,"testcluster"] <- classes
		#print(colnames(xdf))
		library(lawstat)
		p[["brown"]] <- apply(tx, 2, function(x) { levene.test(x, g, location="median")$p.value })
	}
	if (("osum" %in% test  ||  test == "all")  &&  ! "osum" %in% except) {
		if (debug) cat("osum...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			m[,i] <- outlier.sum(x, group=classes, ref=l[i])
		}
		m <- glog(abs(m))
		m <- apply(m, 2, function(x) 1 - x / max(x))  # transform to a p-value like
		p[["osum"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("copa" %in% test  ||  test == "all")  &&  ! "copa" %in% except) {
		if (debug) cat("copa...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			#m[,i] <- pmax(copa(x[,classes == l[i]], c(.75,.90,.95)), copa(x[,classes == l[i]], .90), copa(x[,classes == l[i]], .95))
			m[,i] <- apply(copa(x[,classes == l[i]], c(.75,.90,.95)), 2, max)
		}
		m <- glog(abs(m))
		m <- apply(m, 2, function(x) 1 - x / max(x))  # transform to a p-value like
		p[["copa"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("sam" %in% test  ||  test == "all")  &&  ! "sam" %in% except) {
		if (debug) cat("sam...\n")
		library(samr)
		data <- list(x=x, y=unclass(classes), geneid=as.character(1:nrow(x)),genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)		
		samr.obj <- samr(data,  resp.type=ifelse(nl == 2, "Two class unpaired", "Multiclass"), nperms=nperms)
		p[["sam"]] <- if(ret.all) samr.obj else samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
	}
	if (("pca" %in% test  ||  test == "all")  &&  ! "pca" %in% except) {
		if (debug) cat("pca...\n")
		if (any(is.na(x))) {
			x <- na.impute(x)
		}
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			xsd <- apply(t(x[,classes == l[i]]),2,sd,na.rm=TRUE)
			wsd <- which(xsd == 0 | xsd < 0.0000000001)
			if (length(wsd) > 0) {
				# all values are the same, sd = 0, problems
				pc <- prcomp(t(x[-wsd,classes == l[i]]), scale=TRUE)
				m[-wsd,i] <- apply(abs(t(pc$rotation[,1:4])), 2, max)
				m[wsd,i] <- -Inf
			} else {
				pc <- prcomp(t(x[,classes == l[i]]), scale=TRUE)
				m[,i] <- apply(abs(t(pc$rotation[,1:4])), 2, max)
			}
		}
		pc <- apply(t(m), 2, max)
		p[["pca"]] <- if (ret.all) m else 1 - pc / max(pc)  # transform to a p-value like
	}
	if (("affinity" %in% test  ||  test == "all")  &&  ! "affinity" %in% except) {
		### this is not run in "all" mode because it takes long
		if (debug) cat("affinity propagation...\n")
		m[] <- 0
		npart <- trunc(nrow(x) / 4000 + 1)
		ngenes <- round(nrow(x) / npart,0)
		for (j in 0:(npart-1)) {
			g <- 1:ngenes + j * ngenes
			if (debug) cat("partition ",j,", genes ",min(g),":",max(g),"...\n")
			ptx <- tx[,min(ncol(tx),g)]
			for (i in 1:nl) {
				if (debug) cat(l[i],"\n")
				kpk <- affinity.propagation.clustering(data=tx[classes==l[i],g], debug=debug, max.iter=1000)
				m[g[kpk$exemplars],i] <- 1
			}
		}
		pc <- apply(t(m), 2, sum)
		pc[pc > 0] <- 1 / pc[pc > 0]
		p[["affinity"]] <- if (ret.all) m else 1 - pc / max(pc)  # transform to a p-value like
	}
	
	if (length(grep("perm1:",test)) > 0  &&  ! "perm1" %in% except) {
		xtest <- sub("perm1:","",test[grep("perm1:",test)])		
		thepis <- de.test(x, classes, test=xtest, debug=debug)
		minperm = 30
		if (nperms < minperm) nperms = minperm
		xstats <- list()
		xcounts <- c()
		#BREAKS ARE for p-values
		xbrks <- c(-1e300, 0, 10^-(100:20), 10^-(399:1/20), 1, +1e300)
		for (i in 1:nperms) {
			if (debug  ||  TRUE) {
				if ((i-1) %% 20 == 0) cat("\nPermutation ", i,"") else cat(".")
			}
			xp <- de.test(x, sample(classes), test=xtest, debug=debug)
			if (i <= minperm) {
				xstats[[i]] <- xp
			}
			if (i == minperm) {
				#xr <- range(unlist(lapply(xstats, range)))
				#xd <- abs(xr[2]-xr[1])				
				#xbrks <- c(-Inf, seq(xr[1]-xd,xr[2]+xd,len=100000), +Inf)
				xcounts <- numeric(length(xbrks))
				for (j in 1:minperm) {
					xh <- hist(xstats[[j]], breaks=xbrks, plot=FALSE)$counts
					xcounts <- xcounts + xh
				}
				xstats <- NULL
			} else {
				xh <- hist(xp, breaks=xbrks, plot=FALSE)$counts
				xcounts <- xcounts + xh
			}
		}
		if (debug) cat("\n")
		xsums <- cumsum(xcounts)
		thesum <- sum(xcounts)
		#wisna <- which(is.na(thepis))
		#if (length(wisna) > 0) thepis[wisna] <- 1
		xc <- cut(thepis, breaks=xbrks, labels=FALSE)
		newp <- (xsums[xc]+1)/thesum
		#if (length(wisna) > 0) {
		#	newp[wisna] <- NA
		#	thepis[wisna] <- NA
		#}
		p[["perm1"]] <- if (ret.all) data.frame(p.raw=thepis, perm.counts=xcounts, estimated.p=newp) else newp
	}

	#if ("ftest2" %in% test  ||  test == "all") {
	#	if (debug) cat("ftest2...\n")
	#	tdf <- data.frame(x=tx[,1], g=classes)
	#	p[["ftest2"]] <- apply(tx, 2, function(x) { tdf$x <- x; oneway.test(x ~ g, tdf, var.equal=TRUE)$p.value })
	#}
	#if ("cochran" %in% test  ||  test == "all") {
	#	### this seems to be the same than ftest
	#	if (debug) cat("cochran...\n")
	#	tdf <- data.frame(x=tx[,1], g=classes)
	#	p[["cochran"]] <- apply(tx, 2, function(x) { tdf$x <- x; oneway.test(x ~ g, tdf, var.equal=TRUE)$p.value })
	#}
	# Brown-Forsythe seems to be f-test with unequal variance so the same than welch or very similar
	#ebayes - limma
	#rank product
	#template matching
	#bga - between groups analysis
	#maxT
	#ROC
	#logistic regression
	if (length(p) == 1) p[[1]] else p
}
