plotSignature <- function(svCounts,newSubtypes, plotTitle='', col='gray', cis=NULL, fixCategs=TRUE, ylim=NULL) {
    options(repr.plot.width=24, repr.plot.height=10)  
    par(mar=c(12, 4, 4, 2))  # Increase bottom margin (first value)

    allCategs <- getCategs(newSubtypes)
    
    cols <- rep(col, length(svCounts))
    #cols[grepl('non-clustered_tds',allCategs)] <- 'green'
    if (fixCategs) {
        bp <- barplot(svCounts[allCategs], las=2, cex.names=0.7, border=NA, col=cols, main=paste(plotTitle, sum(svCounts),'SVs'))
    } else {
        allCategs <- allCategs[allCategs %in% names(svCounts)]
        if (is.null(ylim)) {
            bp <- barplot(svCounts[allCategs], las=2, cex.names=0.7, border=NA, col=cols, main=paste(plotTitle, sum(svCounts),'SVs'))
        } else {
            bp <- barplot(svCounts[allCategs], las=2, cex.names=0.7, border=NA, col=cols, main=paste(plotTitle, sum(svCounts),'SVs'), ylim=ylim)
        }
    }
    
    if (!is.null(cis)) {
        y <- svCounts[allCategs]    
        arrows(x0 = bp, y0 = pmax(0, y - cis),  # lower end (not below 0)
           x1 = bp, y1 = y + cis,           # upper end
           angle = 90, code = 3, length = 0.05, lwd = 1.2, col='black')
    }
    if (fixCategs) {
        for (i in seq(from=length(newSubtypes),to = length(allCategs), by =length(newSubtypes))) {
            abline(v=bp[i]+0.5, col='darkgray')
        }
    } else {
        for (i in seq(from=length(newSubtypes),to = length(svCounts), by =length(newSubtypes))) {
            abline(v=bp[i]+0.5, col='darkgray')    
        }
    }

    
}
