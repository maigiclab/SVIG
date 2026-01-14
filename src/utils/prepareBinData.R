    prepareBinData <- function(maxBins=Inf, # how many segments to downsample to
	bedList=c(), # vector of paths to bed files, with names of variables
	bpList=NULL,
	cnDf=NULL,
	repTimingBed=NULL,
    genes.high.expr=NULL,
    genes.low.expr=NULL, addSummary=TRUE, gc_wig_path=NULL,     exprm=NULL, expr.tissue=NULL, gene.table.fn='../data/genes.table.csv', binSize=NULL) {
	
    if (is.null(binSize)) {
        stop("binSize needs to be provided")
    }
        
	# design the bins: 500kb bins 
    binList <- list()
    for (c in c(as.character(1:22), 'X')) {       
        chr.length <- seqlengths(Hsapiens)[[paste('chr',c,sep='')]]       
        binStarts <- seq(from=1, to=chr.length-binSize, by=binSize ) # ER: edit so that bin size is not hard-coded in "to"
        binEnds <- binStarts + binSize-1
        binList[[c]] <- data.frame(chr=c, chromStart=binStarts, chromEnd=binEnds)
    }
    allBins <- do.call('rbind', binList)

    gr.allBins <- GRanges(seqnames=Rle(paste0('chr',allBins$chr)),
                          ranges=IRanges(allBins$chromStart, allBins$chromEnd),
                          strand=rep(c("*"), nrow(allBins)))
    	# down-sample the bins
    
    print('bins are ready')
    
	if (maxBins!=Inf) {
		si <- sample(nrow(allBins))
		allBins <- allBins[si,]
		gr.allBins <- gr.allBins[si,]
	}


	if (!is.null(genes.high.expr) & !is.null(genes.low.expr)) {
	    gr.highly.expr.genes <- GRanges(seqnames=Rle(paste0('chr', genes.high.expr$chr)),
	                         ranges=IRanges(genes.high.expr$chromStart+1, genes.high.expr$chromEnd),
	                         strand=rep(c("*"), nrow(genes.high.expr)) ,
	                             seqlengths=seqlengths(Hsapiens)                        
	                     )
	    gr.lowly.expr.genes <- GRanges(seqnames=Rle(paste0('chr', genes.low.expr$chr)),
	                                    ranges=IRanges(genes.low.expr$chromStart+1, genes.low.expr$chromEnd),
	                                    strand=rep(c("*"), nrow(genes.low.expr)) ,
	                                    seqlengths=seqlengths(Hsapiens))                        

	    coverage.expr.genes <- coverage(gr.highly.expr.genes)
    	coverage.low.expr.genes <- coverage( gr.lowly.expr.genes)
        allBins$highExpGenes <-  rep(NA, nrow(allBins))
    	allBins$lowExpGenes <-  rep(NA, nrow(allBins))
	}

    # replication timing: medianRepTime
    if (!is.null(repTimingBed)) {
    	rep.timing.df <- read.csv(repTimingBed, sep='\t', header=FALSE)
		colnames(rep.timing.df) <- c('chr', 'chromStart', 'chromEnd', 'timing')
		gr.timing <- GRanges(seqnames=Rle(paste0(rep.timing.df$chr)),
        	             ranges=IRanges(rep.timing.df$chromStart+1, rep.timing.df$chromEnd),
            	         strand=rep(c("*"), nrow(rep.timing.df )),
                	     timing=rep.timing.df$timing,
                         seqinfo=seqlevels(gr.allBins)
                     	)
		overlaps.rep.domains <- as.data.frame(findOverlaps(gr.allBins, gr.timing ))
	    allBins$medianRepTime <- rep(NA, nrow(allBins))
	}	

	# extract sequences of the bins, in preparation to count the non-mapping bases
    binSequences <-     as.character(getSeq(Hsapiens, paste0('chr',allBins$chr),
                                            start=allBins$chromStart,
                                            end=allBins$chromEnd))

    # copy number data
    if (!is.null(cnDf)) {

    	gr.ascat <-  GRanges(seqnames=Rle(paste0(cnDf$Chromosome)),
                      ranges=IRanges(cnDf$chromStart, cnDf$chromEnd),
                      strand=rep(c("*"), nrow(cnDf)),
                     seqlengths=seqlengths(Hsapiens),
                     totalCn=cnDf$total.copy.number.inTumour
                     )
	    #bin.midpoints <- rowMeans(allBins[,c('chromStart', 'chromEnd')])
    	#gr.bin.midpoints <- GRanges(seqnames=Rle(paste0('chr',allBins$chr)),
        #                  ranges=IRanges(bin.midpoints , bin.midpoints ),
        #                  strand=rep(c("*"), nrow(allBins)))
        ascat.cov <- coverage(gr.ascat, weight=cnDf$major.copy.number.inTumour)


    	allBins$meanCn <- NA
    }


	covList <- list()
    grList <- list()
	# all the other bed files
    if (length(bedList)>0) {
        
        for (bi in 1:length(bedList)) {
            print(names(bedList)[bi])
            aBed <- read.table(bedList[bi])
            colnames(aBed) <- c('chr', 'chromStart', 'chromEnd')
            aBed <- subset(aBed, chr %in% seqnames(gr.allBins))
            gr.bed <- GRanges(seqnames=Rle(paste0(aBed[,1])),
                        ranges=IRanges(pmax(aBed[,2],1), aBed[,3]),
                        strand=rep(c("*"), nrow(aBed)),
                        seqlengths=seqlengths(Hsapiens)[seqlevels(gr.allBins)]
                        )
            grList[[names(bedList)[bi]]] <- gr.bed
            coverage.bed <-coverage(gr.bed)
            covList[[bi]] <- coverage.bed

            allBins[,names(bedList)[bi]] <- NA
            # added April 24
            ba<- binnedAverage(gr.allBins, coverage.bed, "mean_cvg")
            allBins[,names(bedList)[bi]] <- ba$mean_cvg    
            

        }
    }

	# counting the breakpoints
	if (!is.null(bpList)) {
		# for all breakpoint lists
		for (bps.i in 1:length(bpList)) {
			allBins[,names(bpList)[bps.i]] <- NA
			gr.bps <- GRanges(seqnames=Rle(paste0('chr',bpList[[bps.i ]]$chr)),
                    ranges=IRanges(bpList[[bps.i ]]$position, bpList[[bps.i ]]$position),
                        strand=rep(c("*"), nrow(bpList[[bps.i ]])))
			allBins[,names(bpList)[bps.i]] <-  countOverlaps(gr.allBins,gr.bps)
		}
	}

    
	# summarize all bins
    if (addSummary) {
        
       
        print('loop over bed files')
        for (bi in 1:length(bedList)) {
            featureName <- names(bedList)[bi]
            print(featureName)
            aGr <- grList[[featureName]]
            olps <- as.data.frame(findOverlaps(gr.allBins, aGr))
            olps$overlapEnd <- pmin(end(gr.allBins[olps$queryHits]), end(aGr[olps$subjectHits]))
            olps$overlapStart <- pmax(start(gr.allBins[olps$queryHits]), start(aGr[olps$subjectHits]))
            olps$overlapSize <-  olps$overlapEnd - olps$overlapStart
            
            allBins[,featureName] <- 0
            allBins[olps$queryHits,featureName] <- olps$overlapSize
        }
        

        
        # replication domains
        print('replication origins')
        rep.cov <- coverage(gr.timing, weight=gr.timing$timing)
        a<- binnedAverage(gr.allBins, rep.cov, "timing")
        allBins$medianRepTime <- a$timing

        # copy number
        print('copy number')   
        for (chr in c(1:22, 'X', 'Y')) {
            print(chr)
            is.bin.chr<- allBins$chr==chr
            allBins$meanCn[is.bin.chr] <- ascat.cov[[paste0('chr', chr)]][round(rowMeans(allBins[is.bin.chr,c('chromStart', 'chromEnd')]))]/560
        }
        
        #olps_cnv <- as.data.frame(findOverlaps(gr.allBins, gr.ascat))
        #sampling_x <- 10
        #olps_cnv <- olps_cnv[sample(1:nrow(olps_cnv), round(nrow(olps_cnv)/sampling_x)),]
        #olps_cnv$totalCn <- gr.ascat[olps_cnv$subjectHits,]$totalCn
        #cnv.aggr <- aggregate(olps_cnv$totalCn, by=list(Category=olps_cnv$queryHits), FUN=mean)
        #allBins$meanCn[cnv.aggr$Category] <- cnv.aggr$x
        
        # GC content
        # this page https://research.stowers.org/cws/CompGenomics/Tutorial/GRanges/guide.html
        print('GC content')
        gc.gr<-import(gc_wig_path) 
        gc.gr <- gc.gr[seqnames(gc.gr) %in% paste0('chr', c(1:22,'X'))]
        seqlevels(gc.gr) <- seqlevels(gr.allBins)
        gr.cov <- coverage(gc.gr, weight=as.numeric(gc.gr$name))
        a <- binnedAverage(gr.allBins, gr.cov, "gc")
        allBins$gc <- a$gc
        
            
        #    exprm=NULL, expr.tissue=NULL
        if (!is.null(exprm)) {
            load(exprm)
            # expr.m, disease.m

            genes.table <- as.data.frame(read.csv(gene.table.fn, row.names=NULL)[,2:7])
            disease.df <- as.data.frame(disease.m)
            disease.df$ensid <- gsub('\\..*', '', rownames(disease.m))
            disease.m2 <- merge(disease.df, genes.table, all.x=TRUE,  by.x='ensid', by.y='ensembl_gene_id')
            disease.m2 <- subset(disease.m2, !is.na(hgnc_symbol) & (chr %in% c(1:22, 'X')))
            
            gr.genes <- GRanges(seqnames=Rle(paste0('chr', disease.m2$chr)),
                     ranges=IRanges(disease.m2$chromStart, disease.m2$chromEnd),
                     strand=rep(c("*"), nrow(disease.m2)),
                    seqinfo=seqlevels(gr.allBins)
                 )
            gr.genes$expression.sel <- disease.m2[,expr.tissue]
            gr.genes$length <- disease.m2$chromEnd - disease.m2$chromStart
            
            expr.coverage <- coverage(gr.genes, weight=log(gr.genes$expression.sel+1) )
            ba<- binnedAverage(gr.allBins, expr.coverage, "mean_cvg")
            allBins$expr <- ba$mean_cvg
            
            long.thresh <- quantile(gr.genes$length, 0.8)
            is.long <- gr.genes$length>long.thresh
            expr.coverage <- coverage(gr.genes[is.long], weight=log(gr.genes$expression.sel[is.long]+1) )
            ba<- binnedAverage(gr.allBins, expr.coverage, "mean_cvg")
            allBins$exprLong <- ba$mean_cvg         
            
        }
        
    }

    return(allBins)

}