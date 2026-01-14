evalExprOverlap <- function(test.bedpe, segments.gr, disease.m2.sel, expr.tissue, no.bins=100, plot.tile=exp.name) {

    query_samples <- unique(test.bedpe$sample)
    
    v.list <- list()
    #print(length(query_samples))
    for (si in query_samples) {
        segments.sample.gr <- segments.gr[segments.gr$sample==si]

        genes.gr <- GRanges(seqnames=Rle(paste0('chr',disease.m2.sel$chr)),
                      ranges=IRanges(disease.m2.sel$chromStart,disease.m2.sel$chromEnd), seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        #print(head(as.numeric(disease.m2.sel[,expr.tissue])))
        genes.gr$expr <- as.numeric(disease.m2.sel[,expr.tissue])
        genes.gr$id <- disease.m2.sel$hgnc_symbol

        gaps.gr <- gaps(genes.gr)
        gaps.gr <- gaps.gr[strand(gaps.gr)=='*']
        gaps.gr$expr <- 0
        gaps.gr$id <- 'gap'
        #genes.gaps.gr <- c(genes.gr, gaps.gr)
        genes.gaps.gr <- genes.gr
        #print(head(genes.gaps.gr))
        #print(length(genes.gaps.gr))
        
        o.df <- as.data.frame(findOverlaps(segments.sample.gr, genes.gaps.gr))
        if (nrow(o.df)>0) {
            o.df$rs.midpoint <-segments.sample.gr$midpoint[o.df$queryHits]
            o.df$distGeneStartToMidpoint <- start(genes.gaps.gr)[o.df$subjectHits] - o.df$rs.midpoint 
            o.df$distGeneEndToMidpoint <- end(genes.gaps.gr)[o.df$subjectHits] - o.df$rs.midpoint 
            o.df$rs.start1 <- segments.sample.gr$rearr.start[o.df$queryHits]
            o.df$rs.start2 <- segments.sample.gr$rearr.end[o.df$queryHits]
            o.df$rs.length <- o.df$rs.start2 - o.df$rs.start1
            o.df$rs.class <- segments.sample.gr$rs.class[o.df$queryHits]
            o.df$gene.start <- start(genes.gaps.gr)[o.df$subjectHits]
            o.df$gene.end <- end(genes.gaps.gr)[o.df$subjectHits]
            o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint/o.df$rs.length 
            o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint/o.df$rs.length   
            o.df$id <- NA
            o.df$id <- genes.gaps.gr$id[o.df$subjectHits]
            o.df$expr <- NA
            #print(head(o.df))
            
            o.df$expr <- genes.gaps.gr$expr[o.df$subjectHits]

            
            o.df$distGeneStartToMidpointScaledInt <- round(o.df$distGeneStartToMidpointScaled,2)*no.bins
            o.df$distGeneEndToMidpointScaledInt <- round(o.df$distGeneEndToMidpointScaled,2) * no.bins
            minBin <- max(3.5*no.bins + 1, abs(min(o.df$distGeneStartToMidpointScaledInt)) + 1)
            #print(paste('min', min(o.df$distGeneStartToMidpointScaledInt)))
            o.df$distGeneStartToMidpointScaledIntPos <- o.df$distGeneStartToMidpointScaledInt + minBin
            o.df$distGeneEndToMidpointScaledIntPos <- o.df$distGeneEndToMidpointScaledInt + minBin
            maxBinPos <- max(o.df$distGeneEndToMidpointScaledIntPos)
            o.df$geneBinFootprint <- o.df$distGeneEndToMidpointScaledInt - o.df$distGeneStartToMidpointScaledInt + 1

            #print(head(o.df))
            cvg.genes = coverage(IRanges(start=o.df$distGeneStartToMidpointScaledIntPos, 
                                   end=o.df$distGeneEndToMidpointScaledIntPos), 
                         #weight=log(o.df$expr+1)/o.df$rs.length)
                         weight=log(o.df$expr+1)/(o.df$rs.length))
                         #weight=log(o.df$expr+1)/(o.df$rs.length *  o.df$geneBinFootprint))
            if (length(cvg.genes)<round(3.5*no.bins + minBin)) {
                cvg.genes <- c(cvg.genes, rep(0,3.5*no.bins + minBin - length(cvg.genes) +1))
            }
                
            
            #plot(cvg.genes, type='l', xlim=c(-3.5*no.bins + minBin, 3.5*no.bins + minBin), main=si)
            #abline(v=-3.5*no.bins + minBin)
            #abline(v=3.5*no.bins + minBin)
            #abline(v=minBin)
            #abline(v=-0.5*no.bins + minBin)
            #abline(v=0.5*no.bins + minBin)
            
            #print(paste('no.bins', no.bins))
            #print(paste('minBin', minBin))
            #print(paste('minDist',min(o.df$distGeneStartToMidpointScaledInt)))  
            #print(round(-3.5*no.bins + minBin))
            #print(round(3.5*no.bins + minBin))
            #print(paste('length',length(cvg.genes)))
            
            v <- cvg.genes[round(-3.5*no.bins + minBin):round(3.5*no.bins + minBin)]
            v.list[[si]] <- as.vector(v)

        }

    }
    
    v.m<-do.call(rbind,v.list)
    v <- colSums(v.m)/nrow(test.bedpe)
    plot(v, type='l', main=plot.tile, ylab='sum expression log', ylim=c(0, max(v)))
    abline(v=3*no.bins)
    abline(v=4*no.bins)
}