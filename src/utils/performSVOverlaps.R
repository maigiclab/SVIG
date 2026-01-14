performSVOverlaps.R <- function(muts.gr, test.bedpe, simulate, margin.size) {
# find and characterize overlaps between point mutations and SVs.
# should be run per sample
# muts.gr - a genomic ranges object with point mutations
# test.bedpe - a data frame with SVs: chrom1, start1, start2, svclass
# simulate - whether perturb the observed SVs
# margin.size - margin size around SVs [bps]
    
    if (nrow(test.bedpe)>0) {
        if (simulate) {
            test.bedpe$start1 <- test.bedpe$start2
            # fake duplicaitons!
            test.bedpe$start2 <- test.bedpe$start1+ test.bedpe$length         
        }
        # prepare a GR object with SVs
        rs.gr <- GRanges(seqnames=Rle(paste0('chr', test.bedpe$chrom1)),
                      ranges=IRanges(test.bedpe$start1, test.bedpe$start2))
        rs.gr$rs.class = test.bedpe$svclass
    
        
        # prepare the GR object with margin regons     
        margin.left.gr <- GRanges(seqnames=Rle(paste0('chr', test.bedpe$chrom1)),
                              ranges=IRanges(test.bedpe$start1-margin.size, test.bedpe$start1-1),seqinfo= seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        margin.left.gr$dup.rel.pos <- -1
        margin.left.gr$dup.start <- test.bedpe$start1
        margin.left.gr$dup.end <- test.bedpe$end2

        # overlaps of SVs with margins
        margin.right.gr <- GRanges(seqnames=Rle(paste0('chr', test.bedpe$chrom1)),
                              ranges=IRanges(test.bedpe$end2+1, test.bedpe$end2+margin.size),seqinfo= seqinfo(BSgenome.Hsapiens.UCSC.hg19))
        margin.right.gr$dup.rel.pos <- 1
        margin.right.gr$dup.start <- test.bedpe$start1
        margin.right.gr$dup.end <- test.bedpe$start2
        
         # Combine GRanges objects: left and right margins
        all.margins.gr <- c(margin.left.gr, margin.right.gr)
        # Identify overlaps
        # Identify self-overlapping ranges and drop the "subject" member of each overlap pair
        # (with drop.redundant=TRUE, this typically drops the later index in each overlapping pair)
        overlaps <- findOverlaps(all.margins.gr, drop.self = TRUE, drop.redundant = TRUE)
        keep <- setdiff(seq_along(all.margins.gr), subjectHits(overlaps))
        all.margins.gr <- all.margins.gr[keep]
        
        # find the overlaps between mutations and SVs
        overlap.df <- as.data.frame(findOverlaps(muts.gr, rs.gr))
        # find the overlaps between mutations and SV margins
        overlap.margin.df <- as.data.frame(findOverlaps(muts.gr, all.margins.gr))
            
        # deal with clustered SVs by creating a separate GRanges object    
        clustered.bedpe <- subset(do.call('rbind',sample.rearrs[sample_id]),(is.clustered==TRUE) & (svclass %in% c('deletion', 'duplication', 'inversion')))   
        if (nrow(clustered.bedpe)>0) {
            clustered.gr <- GRanges(seqnames=Rle(paste0('chr', clustered.bedpe$chrom1)),
                      ranges=IRanges(clustered.bedpe$start1, clustered.bedpe$start2))                
        } else {
            clustered.gr <- GRanges()
        }
        clustered.overlap.df <- as.data.frame(findOverlaps(muts.gr, clustered.gr))
        
        # annotate the mutations according to the overlaps
        muts_df$is.within.dup <- FALSE
        muts_df$is.within.dup[overlap.df$queryHits] <- TRUE
        muts_df$is.within.cluster <- FALSE
        muts_df$is.within.cluster[clustered.overlap.df$queryHits] <- TRUE
        muts_df$is.within.margin <- FALSE
        muts_df$is.within.margin[overlap.margin.df$queryHits] <- TRUE    
        muts_df$overlapsMultipleDups <- FALSE
        muts_df$overlapsMultipleDups[overlap.df$queryHits[duplicated(overlap.df$queryHits)]] <- TRUE

        muts_df$dup.ID <- NA
        muts_df$dup.ID[overlap.df$queryHits] <- paste(test.bedpe$sample[overlap.df$subjectHits], test.bedpe$sv_id[overlap.df$subjectHits])
        muts_df$dup.start <- NA
        muts_df$dup.start[overlap.df$queryHits] <- start(rs.gr[overlap.df$subjectHits])
        muts_df$dup.end <- NA
        muts_df$dup.end[overlap.df$queryHits] <- end(rs.gr[overlap.df$subjectHits])
        muts_df$dup.length <- muts_df$dup.end - muts_df$dup.start
        muts_df$position.witin.dup <- NA
        muts_df$position.witin.dup[overlap.df$queryHits] <- (muts_df$starts[overlap.df$queryHits] - muts_df$dup.start[overlap.df$queryHits]) / (muts_df$dup.length[overlap.df$queryHits])

        muts_df$margin.pos <- NA
        muts_df$margin.pos[overlap.margin.df$queryHits] <- all.margins.gr[overlap.margin.df$subjectHits]$dup.rel.pos
        
        muts_df$position.witin.margin <- NA
        
        muts_df$is.within.dup.label <- mapvalues(muts_df$is.within.dup, 
          from=c(TRUE,FALSE), 
          to=c("within a duplication","outside"))    
    
        # table of verlap frequesncies for muts
        t <- table(muts_df$is.within.dup.label, 
                   muts_df$max_sig)
        t.norm <- t/rowSums(t)  
    
        # now check the apobec muts specificaly
        sample_muts_apobec <- subset(muts_df, max_sig %in% c('SBS2', 'SBS13') & MutCN<1.5)
        if (nrow(sample_muts_apobec)>0) {
            sample_muts_apobec <- sample_muts_apobec[order(sample_muts_apobec$chroms, sample_muts_apobec$starts),]
            sample_muts_apobec$distPrev <- sample_muts_apobec$starts - c(NA, sample_muts_apobec$starts[1:(nrow(sample_muts_apobec)-1)])
            sample_muts_apobec$distPrev[sample_muts_apobec$chroms != c('-1',sample_muts_apobec$chroms[1:(nrow(sample_muts_apobec)-1)])] <- NA
    
            sample_muts_apobec$isC <- sample_muts_apobec$wt=='C'
            sample_muts_apobec$position.witin.dup.cut <- NA
            if (sum(!is.na(sample_muts_apobec$position.witin.dup))) {
                sample_muts_apobec$position.witin.dup.cut <- cut(sample_muts_apobec$position.witin.dup, breaks=10 )
                t <- table(sample_muts_apobec$position.witin.dup.cut, sample_muts_apobec$isC)
            }
        sample_muts_apobec$sample <- sample_id
        sample_muts_apobec$dcc_project_code <- dcc_project_code

        sv.footprint <- sum(sum(coverage(reduce(rs.gr))))
        margin.footprint <- sum(sum(coverage(all.margins.gr)))
        bg.footprint <- sum(sum(coverage(mapability.gr))) - sum(sum(coverage(reduce(rs.gr)))) - sum(sum(coverage(all.margins.gr)))

        # create summary statistics for the sample
        sample.summary <- data.frame(sample=sample_id,
                                               apobec.rate.dups=sum(sample_muts_apobec$is.within.dup)/sv.footprint,
                                               apobec.rate.margin=sum(sample_muts_apobec$is.within.margin)/margin.footprint,
                                               apobec.rate.background=sum(!sample_muts_apobec$is.within.dup & !sample_muts_apobec$is.within.margin)/ bg.footprint,
                                               apobec.number.dups=sum(sample_muts_apobec$is.within.dup),
                                               apobec.number.margin=sum(sample_muts_apobec$is.within.margin),
                                               apobec.number.background=sum(!sample_muts_apobec$is.within.dup & !sample_muts_apobec$is.within.margin), # it will be background in the same samples
                                               margin.footprint=margin.footprint,
                                               sv.footprint=sv.footprint,
                                               bg.footprint=bg.footprint,
                                               no.dups=nrow(test.bedpe)
                                               )
        r <- list()
        r[['sample.summary']] <- sample.summary
        r[['sample_muts_apobec']] <- sample_muts_apobec
        r[['muts_df']] <- muts_df
        return(r)
        }
    }
}