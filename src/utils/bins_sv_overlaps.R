bins_sv_overlaps <- function(common.origins.1Mb.gr, rs.gr, bin.size=1e3) {
    o.to.segment.df <- as.data.frame(findOverlaps(common.origins.1Mb.gr, rs.gr))

    
    o.to.segment.df$origin.midpoint <- core.origins$midopoint[o.to.segment.df$queryHits]
    o.to.segment.df$rs.segment.start <- start(rs.gr)[o.to.segment.df$subjectHits]
    o.to.segment.df$rs.segment.end <- end(rs.gr)[o.to.segment.df$subjectHits]

    
    o.to.segment.df$rs.segment.start.rel.unbound <- (o.to.segment.df$rs.segment.start - o.to.segment.df$origin.midpoint + 1e6)/bin.size + 1
    o.to.segment.df$rs.segment.end.rel.unbound <- (o.to.segment.df$rs.segment.end - o.to.segment.df$origin.midpoint + 1e6)/bin.size + 1

    o.to.segment.df$rs.segment.start.rel <- (pmax(o.to.segment.df$rs.segment.start - o.to.segment.df$origin.midpoint, -1e6) + 1e6)/bin.size + 1
    o.to.segment.df$rs.segment.end.rel <- (pmin(o.to.segment.df$rs.segment.end - o.to.segment.df$origin.midpoint, 1e6) + 1e6)/bin.size + 1
    # select only rearrangements RS1, RS3, RS5
    o.to.segment.df$length <- o.to.segment.df$rs.segment.end -  o.to.segment.df$rs.segment.start
    return(o.to.segment.df)
      
}