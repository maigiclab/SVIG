# replication timing data
replis.gr <- toGRanges('~/projects/rsignatures/data//external/BASIS/Repliseq_data//RepliTime/MCF7_data//MCF7_RepliSeq.bedGraph')
replis.gr <- replis.gr[lengths(replis.gr)<=1001]

#rfd_fn <- '~/park_dglodzik/OK-Seq/HeLa_S3/Hela.EdC.Combined_OkaSeq.RFD.bw'
#rfd.gr <- rtracklayer::import(rfd_fn, format = "bigWig")
rfd_fn <- '~/park_dglodzik/OK-Seq/HeLa_MRL2/Hela.EdC.Combined_OkaSeq.RFD.bw'
rfd.gr <- rtracklayer::import(rfd_fn, format = "bigWig")
rfd.gr <- rfd.gr[rfd.gr$score!=-2]

# initiation and termination zones
iz_fn <- '~/park_dglodzik/OK-Seq/HeLa_MRL2/HeLa_hmm_HMMsegments_IZ.bed'
iz.gr <- toGRanges(iz_fn)
tz_fn <- '~/park_dglodzik/OK-Seq/HeLa_MRL2/HeLa_hmm_HMMsegments_TZ.bed'
tz.gr <- toGRanges(tz_fn)
total_footprint_iz <- sum(width(iz.gr))
total_footprint_tz <- sum(width(tz.gr))

load('~/projects//rsignatures/data/processed/disease.RData')
dim(disease.m)


gene.table.fn <- '~/park_dglodzik/repos/hotspots/data/genes.table.csv'
genes.table <- as.data.frame(read.csv(gene.table.fn, row.names=NULL)[,2:7])
disease.df <- as.data.frame(disease.m)
disease.df$ensid <- gsub('\\..*', '', rownames(disease.m))
disease.m2 <- merge(disease.df, genes.table, all.x=TRUE,  by.x='ensid', by.y='ensembl_gene_id')
disease.m2 <- subset(disease.m2, !is.na(hgnc_symbol) & (chr %in% c(1:22, 'X')))
genes.gr <- GRanges(seqnames=Rle(paste0('chr',disease.m2$chr)),
              ranges=IRanges(disease.m2$chromStart,disease.m2$chromEnd))
genes.gr$hgnc_symbol <- disease.m2$hgnc_symbol
genes.gr$ensid <- disease.m2$ensid
genes.gr$expr_fpkm <- disease.m2[,expr_tissue]