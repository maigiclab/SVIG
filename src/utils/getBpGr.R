getBpGr <- function(bedpe, addMidpoint=FALSE, addMargin=FALSE) {
# bedpe fields required
# chrom1, start1,
# chrom2, start2
# sv_id_global

    if (addMidpoint) {

        bedpe_nontransloc <- subset(bedpe, chrom1==chrom2)
        bps.df <- rbind(data.frame(chr=bedpe$chrom1, position=bedpe$start1),
                       data.frame(chr=bedpe$chrom2, position=bedpe$start2),
                       data.frame(chr=bedpe_nontransloc$chrom2, position=rowMeans(bedpe_nontransloc[,c('start1','start2')])) 
        )
        if (addMargin) {
            bps.df <- rbind(data.frame(chr=bedpe$chrom1, position=bedpe$start1-3*bedpe$length),
                       data.frame(chr=bedpe$chrom2, position=bedpe$start2 + 3*bedpe$length),
                       data.frame(chr=bedpe_nontransloc$chrom2, position=rowMeans(bedpe_nontransloc[,c('start1','start2')])) 
            )
        }
        if ('sv_id_global' %in% colnames(bedpe)) {
        bps.df$sv_id_global <- c(
                        paste0(bedpe$sv_id_global, '_bp1'),
                        paste0(bedpe$sv_id_global, '_bp2'),
                        paste0(bedpe_nontransloc$sv_id_global, '_midpoint')
                        )

        }
        bps.df$type <- c(rep('bp1', nrow(bedpe)),
                 rep('bp2', nrow(bedpe)),
                 rep('midpoint', nrow(bedpe_nontransloc))
                 )
    
    } else if (('start3' %in% colnames(bedpe)) && ('start4' %in% colnames(bedpe))) {
        bps.df <- rbind(data.frame(chr=bedpe$chrom1, position=bedpe$start1),
                       data.frame(chr=bedpe$chrom2, position=bedpe$start2),
                       data.frame(chr=bedpe$chrom1, position=bedpe$start3),
                        data.frame(chr=bedpe$chrom1, position=bedpe$start4) 
        )
        if ('sv_id_global' %in% colnames(bedpe)) {
            bps.df$sv_id_global <- c(
                        paste0(bedpe$sv_id_global, '_bp1'),
                        paste0(bedpe$sv_id_global, '_bp2'),
                        paste0(bedpe$sv_id_global, '_bp3'),
                        paste0(bedpe$sv_id_global, '_bp4')
                        )
        }
        bps.df <- subset(bps.df, !is.na(position))
    } else {
        bps.df <- rbind(data.frame(chr=bedpe$chrom1, position=bedpe$start1),
                       data.frame(chr=bedpe$chrom2, position=bedpe$start2)
        )
        if ('sv_id_global' %in% colnames(bedpe)) {
        bps.df$sv_id_global <- c(
                        paste0(bedpe$sv_id_global, '_bp1'),
                        paste0(bedpe$sv_id_global, '_bp2')
                        )
        }
    }

    

    
    if (grepl('chr', bps.df$chr[1])) {
    bps.gr <- GRanges(seqnames=Rle(paste0( bps.df$chr)),
                  ranges=IRanges(bps.df$position, 
                                 bps.df$position))            
    } else {
    bps.gr <- GRanges(seqnames=Rle(paste0('chr', bps.df$chr)),
                  ranges=IRanges(bps.df$position, 
                                 bps.df$position))        
    }

    if ('sv_id' %in% colnames(bedpe)) {
        bps.gr$sv_id_global <- bps.df$sv_id_global
    }
    if (addMidpoint) {
        bps.gr$type <- bps.df$type
    }
    return(bps.gr)
    
}