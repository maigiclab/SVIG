#marginmode <- '3len'
#marginmode <- '1Mb'

# define bins
# no.bins bins per Mb or SV length
# when working with scaled SVs, gene counts overlaping SV will be scaled inversely to SV length
# that's not true when working with fixed 1Mb borders

rs_gene_overlaps <- function(genes_df, rs.breast.gr, rs.groups=NULL, marginmode='3len', 
                             plot_tile='',ylab='expressed genes norm. wrt. gene and SV size', 
                             ymax=NA,   
                             no.bins =50, # per SV unit length, or per Mb
                            add.newplot = TRUE,
                             lineCol='black'
                            ) {
    # SV rescaling is additionally scaling by inverse of SV length on the length rescaled analysis
    
    
    genes.high.gr <- GRanges(seqnames=Rle(paste0('chr',genes_df$chr)),
                  ranges=IRanges(genes_df$chromStart,genes_df$chromEnd))
    genes.high.gr$hgnc_symbol <- genes_df$hgnc_symbol
    
    
    # query, subject
    o.df <- as.data.frame(findOverlaps(rs.breast.gr, genes.high.gr))
    o.df$hgnc_symbol <- genes_df$hgnc_symbol[o.df$subjectHits]
   
    # gene start and end to midpoint
    o.df$rs.midpoint <-rs.breast.gr$midpoint[o.df$queryHits]
    o.df$distGeneStartToMidpoint <- genes_df$chromStart[o.df$subjectHits] - o.df$rs.midpoint 
    o.df$distGeneEndToMidpoint <- genes_df$chromEnd[o.df$subjectHits] - o.df$rs.midpoint 
    o.df$geneStrand <- genes_df$strand[o.df$subjectHits]

    o.df$rs.start1 <- rs.breast.gr$start1[o.df$queryHits]
    o.df$rs.start2 <- rs.breast.gr$start2[o.df$queryHits]
    o.df$rs.length <- o.df$rs.start2 - o.df$rs.start1
    o.df$rs.class <- rs.breast.gr$rs.class[o.df$queryHits]


    if (marginmode=='3len') {
        o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint/o.df$rs.length 
        o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint/o.df$rs.length        
    } else if (marginmode == '1Mb') {
       o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint/1e6
       o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint/1e6
    }

    o.df$distGeneStartToMidpointScaledInt <- round(o.df$distGeneStartToMidpointScaled,2)*no.bins
    o.df$distGeneEndToMidpointScaledInt <- round(o.df$distGeneEndToMidpointScaled,2) * no.bins
    o.df$geneBinFootprint <- o.df$distGeneEndToMidpointScaledInt - o.df$distGeneStartToMidpointScaledInt + 1

    minBin <- abs(min(o.df$distGeneStartToMidpointScaledInt)) + 1

    
    o.df$distGeneStartToMidpointScaledIntPos <- o.df$distGeneStartToMidpointScaledInt + minBin
    o.df$distGeneEndToMidpointScaledIntPos <- o.df$distGeneEndToMidpointScaledInt + minBin
    maxBinPos <- max(o.df$distGeneEndToMidpointScaledIntPos)
    

    if (!is.null(rs.groups)) {
        options(repr.plot.width=28, repr.plot.height=5)
        par(mfrow=c(1,6))
        for (rs in c('RS1', 'RS2', 'RS3', 'RS4', 'RS5', 'RS6')) {

            median_len <- median(rs.breast.gr[rs.breast.gr$rs.class==rs]$length)
            no.rearrs <- sum(rs.breast.gr$rs.class==rs)
            
            # coverage of all genes
            xvals <- (1:length(cvg.rs1) - minBin) / no.bins
            yvals <- (cvg.rs1 / no.rearrs) * 1e6 * no.bins / 7
            keep <- xvals >= -3.5 & xvals <= 3.5
            xvals <- xvals[keep]
            yvals <- yvals[keep]
     
            if (marginmode=='3len') {
                cvg.rs1 = coverage(IRanges(start=subset(o.df, rs.class==rs)$distGeneStartToMidpointScaledIntPos, 
                                       end=subset(o.df, rs.class==rs)$distGeneEndToMidpointScaledIntPos), 
                               weight=1/((subset(o.df, rs.class==rs)$rs.length)))# * (subset(o.df, rs.class==rs)$geneBinFootprint)))
                
                        # (h$counts/(nrow(sample.rearrs.all.conf.rs))) * (no.bins * 1e6) / (7*median.size.rs)
                plot(xvals, yvals, 
                  cex.lab=2,cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlim=c(-3.5, +3.5), type='l', lwd=2,
                 xlab='distance to rearrangement midpoint [1Mb]',
                ylab=ylab,
                 main=rs,xaxs="i", ylim=c(0,ymax)
                )
            abline(v=0, col='gray')
            abline(v=-3.5, col='gray')
            abline(v=+3.5, col='gray')
            abline(v=-0.5, col='gray')
            abline(v=+0.5, col='gray')
                
            } else if (marginmode == '1Mb') {
                cvg.rs1 = coverage(IRanges(start=subset(o.df, rs.class==rs)$distGeneStartToMidpointScaledIntPos, 
                                       end=subset(o.df, rs.class==rs)$distGeneEndToMidpointScaledIntPos))
                
                plot((1:length(cvg.rs1)-minBin)/no.bins, 
                 (cvg.rs1 / no.rearrs),  
                  cex.lab=2,cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlim=c(-1, +1), type='l', lwd=2,
                 xlab='scaled distance to rearrangement midpoint',
                ylab='normalized count of expressed genes',
                 main=rs,xaxs="i", ylim=c(0,.5)
                )
                abline(v=0, col='gray')
            }
        }        
    } else {
        #options(repr.plot.width=7, repr.plot.height=7)
        #par(mfrow=c(1,1))
        
        no.rearrs <- length(rs.breast.gr)     
        if (marginmode=='3len') {
            
            cvg.rs1 = coverage(IRanges(start=o.df$distGeneStartToMidpointScaledIntPos, 
                                       end=o.df$distGeneEndToMidpointScaledIntPos)) 
                             #weight=1/(o.df$rs.length))
            #ys <- (cvg.rs1 / no.rearrs) * 1e6 * no.bins / ( 7 )
            yvals <- (cvg.rs1 / no.rearrs) #* 1e6 
            if (is.na(ymax)) {
                ymax=max(yvals)
            }

            # coverage of all genes
            xvals <- (1:length(cvg.rs1) - minBin) / no.bins # no.bins per SV length
            keep <- xvals >= -3.5 & xvals <= 3.5
            xvals <- xvals[keep]
            yvals <- yvals[keep]

            
            
            if (add.newplot) {
                plot(xvals, 
                 yvals,  
                 cex.lab=2, cex.axis=2, cex.main=2, cex.sub=1.5, 
                 xlim=c(-3.5, 3.5), 
                 ylim=c(0, ymax),
                 type='l', 
                 xlab='Dist. to SV midpoint [SV lengths]',
                 ylab=ylab,
                 main=plot_tile,
                 lwd=6,
                 xaxt = "n", xaxs = "i", 
                 yaxt = "s",  # keep default y-axis
                 bty = "l",
                 col = lineCol)   # draw only left and bottom box lines

                            # Add your vertical lines
                abline(v=-.5, lwd=6, col="#4EA72E")
                abline(v=0.5, lwd=6 ,col="#4EA72E")
                abline(v=0, lty=2)
                
                # Add custom x-axis only
                xticks <- c(-0.5, 0, 0.5)
                axis(side=1, at=xticks, labels=c('-.5', '0', '+.5'), cex.axis=2)
                 
                     
            } else {
                 lines(xvals, 
                 yvals,
                     lwd=6, col=lineCol)           
                }
            

            
        } else if (marginmode=='1Mb') {
            cvg.rs1 = coverage(IRanges(start=o.df$distGeneStartToMidpointScaledIntPos, 
                                       end=o.df$distGeneEndToMidpointScaledIntPos)) 
                             #weight=1/(o.df$rs.length))
            
            if (add.newplot) {
                plot((1:length(cvg.rs1)-minBin)/no.bins, 
                     (cvg.rs1 / no.rearrs),
                     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=1.5,
                     type='l', lwd=2,
                     xlab=' distance to rearrangement midpoint [1Mb]',
                    ylab=ylab,
                     main=plot_tile,xaxs="i",
                     col=lineCol,
                     
                     xlim=c(-1, 1)
                    )
                points((1:length(cvg.rs1)-minBin)/no.bins, (cvg.rs1 / no.rearrs), pch=19)
                abline(v=0, col='gray')                  
            } else {
                lines((1:length(cvg.rs1)-minBin)/no.bins, 
                     (cvg.rs1 / no.rearrs),
                     lwd=2, col=lineCol)
            }
     
        }
    }

    r <- list()
    r[['o.df']] <- o.df
    r[['max_y']] <- max((cvg.rs1 / no.rearrs))
    return(r)
}