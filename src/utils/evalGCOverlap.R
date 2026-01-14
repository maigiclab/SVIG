evalGCOverlap <- function(test.bedpe, segments.gr, gc.gr) {

    query_samples <- unique(test.bedpe$sample)
    
    v.list <- list()
    for (si in query_samples) {
        segments.sample.gr <- segments.gr[segments.gr$sample==si]

        o.df <- as.data.frame(findOverlaps(segments.sample.gr, gc.gr))
        if (nrow(o.df)>0) {
            o.df$rs.midpoint <-segments.sample.gr$midpoint[o.df$queryHits]
            o.df$distGeneStartToMidpoint <- start(gc.gr)[o.df$subjectHits] - o.df$rs.midpoint 
            o.df$distGeneEndToMidpoint <- end(gc.gr)[o.df$subjectHits] - o.df$rs.midpoint 
            o.df$rs.start1 <- segments.sample.gr$rearr.start[o.df$queryHits]
            o.df$rs.start2 <- segments.sample.gr$rearr.end[o.df$queryHits]
            o.df$rs.length <- o.df$rs.start2 - o.df$rs.start1
            o.df$rs.class <- segments.sample.gr$rs.class[o.df$queryHits]
            o.df$gene.start <- start(gc.gr)[o.df$subjectHits]
            o.df$gene.end <- end(gc.gr)[o.df$subjectHits]
            o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint/o.df$rs.length 
            o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint/o.df$rs.length   
            o.df$id <- NA
            
            o.df$gc <- as.numeric(gc.gr$name[o.df$subjectHits])
            
            no.bins <- 100
            o.df$distGeneStartToMidpointScaledInt <- round(o.df$distGeneStartToMidpointScaled,2)*no.bins
            o.df$distGeneEndToMidpointScaledInt <- round(o.df$distGeneEndToMidpointScaled,2) * no.bins
            minBin <- abs(min(o.df$distGeneStartToMidpointScaledInt)) + 1
            o.df$distGeneStartToMidpointScaledIntPos <- o.df$distGeneStartToMidpointScaledInt + minBin
            o.df$distGeneEndToMidpointScaledIntPos <- o.df$distGeneEndToMidpointScaledInt + minBin
            maxBinPos <- max(o.df$distGeneEndToMidpointScaledIntPos)
            o.df$geneBinFootprint <- o.df$distGeneEndToMidpointScaledInt - o.df$distGeneStartToMidpointScaledInt + 1

            o.df <- subset(o.df, !is.na(gc))
            cvg.genes = coverage(IRanges(start=as.integer(o.df$distGeneStartToMidpointScaledIntPos), 
                                   end=as.integer(o.df$distGeneEndToMidpointScaledIntPos)), 
                         weight=o.df$gc)
            cvg.overlaps = coverage(IRanges(start=as.integer(o.df$distGeneStartToMidpointScaledIntPos), 
                                   end=as.integer(o.df$distGeneEndToMidpointScaledIntPos)))
            
            v <- cvg.genes[round(-3.5*no.bins + minBin):round(3.5*no.bins + minBin)]
            v2 <- cvg.overlaps[round(-3.5*no.bins + minBin):round(3.5*no.bins + minBin)]
            v.list[[si]] <- as.vector(v)/as.vector(v2)

        }

    }
    
    v.m<-do.call(rbind,v.list)
    v <- colMeans(v.m, na.rm=TRUE)
    plot(v, type='l', main=exp.name, ylab='average GC content')
    abline(v=3*no.bins)
    abline(v=4*no.bins)

}