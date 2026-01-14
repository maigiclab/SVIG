loadRepeats <- function(repeat_fn) {

    # Column names per UCSC rmsk schema
    colnames_rmsk <- c(
      "bin","swScore","milliDiv","milliDel","milliIns",
      "genoName","genoStart","genoEnd","genoLeft","strand",
      "repName","repClass","repFamily",
      "repStart","repEnd","repLeft","id"
    )
    # Read with base R
    con <- gzfile(repeat_fn, "rt")
    rmsk <- read.table(
      con, header = FALSE, sep = "\t", quote = "", comment.char = "",
      col.names = colnames_rmsk,
      colClasses = c(
        "integer","integer","numeric","numeric","numeric",
        "character","integer","integer","integer","character",
        "character","character","character",
        "integer","integer","integer","integer"
      )
    )
    close(con)
    # (Optional) keep only canonical chromosomes
    keep_chr <- rmsk$genoName %in% c(paste0("chr", 1:22), "chrX", "chrY", "chrM", "chrMT")
    rmsk <- rmsk[keep_chr, ]
    # Build GRanges
    # Note: UCSC coords are 0-based start, 1-based end → add 1 to start for GRanges
    rmsk_gr <- GRanges(
      seqnames = rmsk$genoName,
      ranges   = IRanges(start = rmsk$genoStart + 1L, end = rmsk$genoEnd),
      strand   = rmsk$strand,
      repName  = rmsk$repName,
      repClass = rmsk$repClass,
      repFamily= rmsk$repFamily
    )


    return(rmsk_gr)
}