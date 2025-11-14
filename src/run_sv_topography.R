if (interactive()) {
    # I did not get it to work for RS2, RS6a
    #exp.name <- 'non_clustered'
    #filter_str <- " (is.clustered==FALSE)"
    exp.name <- 'RS1'
    filter_str <- "max.Ref.Sig=='Ref.Sig.R1'  & (is.clustered==FALSE)"
    #exp.name <- 'RS2'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R2'"
    #exp.name <- 'RS6a'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R6a'"
    #exp.name <- 'RS6b'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R6b'"
    
    #exp.name <- 'RS8'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R8'  & (is.clustered==FALSE)"
    #exp.name <- 'RS3'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R3'  & (is.clustered==FALSE) "
    #exp.name <- 'R7'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R7'  & (is.clustered==FALSE) "
    #exp.name <- 'cdk12_dups'
    #filter_str <- "sample %in% curated_cdk12 & label2=='_tds' & is.clustered==FALSE  & length < 5e6 & length >1e5"
    #exp.name <- 'ccne1_dups'
    #filter_str <- "(sample %in% ccne.samples.df$aliquot_id) & label2=='_tds' & is.clustered==FALSE  & length < 5e6 & length >1e5"
    #exp.name <- 'RS11'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R11'  & (is.clustered==FALSE)"    
    #exp.name <- 'RS15'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R15'  & (is.clustered==FALSE)"    
    #exp.name <- 'RS1_esophagus'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R1' & label2=='_tds'& V15=='Esophagus'"   
    #exp.name <- 'RS1_bone'
    #filter_str <- "max.Ref.Sig=='Ref.Sig.R1' & label2=='_tds'& V15=='Bone/SoftTissue'"
    
    #exp.name <- 'cdk12_dups_shorter'
    #filter_str <- "sample %in% curated_cdk12 & label2=='_tds' & is.clustered==FALSE  & length < 1e6 & length >1e5
    
    shuffle <- FALSE
    shuffleMode <- 'simple'
    doRegression <- FALSE
} else {
    args <- commandArgs(trailing = T)
    exp.name <-  args[1]
    filter_str <- args[2]
    if (length(args)>2) {
        shuffle_str <- args[3]
        if (shuffle_str=='shuffleFlip') {
            shuffle<-TRUE
            shuffleMode <- 'flip'
        } else if (shuffle_str=='shuffleSimple'){
            shuffle<-TRUE
            shuffleMode <- 'simple'
        } else {
            shuffle<-FALSE 
            shuffleMode <- 'simple'
        }
    } else {
        shuffle<-FALSE
        shuffleMode <- 'simple'
    }
    
    doRegression <- FALSE
}
print(exp.name)
print(filter_str)
print(paste('shuffle:',shuffle, shuffleMode))

makePDFs <- TRUE
# create plots of replication origin density around the rearrangement midpoints
marginmode <- '3len'


chr.margin <- 1e6

if (shuffle==TRUE) {
    exp.name <- paste0(exp.name, '_shuffle_', shuffleMode)
}

result_folder <- paste0('/home/dg204/park_dglodzik/svig/', exp.name, '/')
if (!dir.exists(result_folder)) {
  dir.create(result_folder)
}
suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg19))
suppressWarnings(library("GenomicRanges"))
suppressWarnings(library(regioneR))
suppressWarnings(library("plotrix"))
suppressWarnings(library(dplyr))
suppressWarnings(library(tidyverse))
suppressWarnings(library(Rsamtools))
suppressWarnings(library(plyr))
library(signature.tools.lib)
library(lmtest)
library(MASS)
library(stringr)
library(rhdf5)
source('~/repos/hotspots/R/utils/prepareBinData.R')
source('~/projects/rsignatures/src/utils/evalExprOverlap.R')
source('~/projects/rsignatures/src/utils/evalGCOverlap.R')
source('~/projects/rsignatures/./src/utils/qcut.R')
source('~/projects/rsignatures//src/utils/bins_sv_overlaps.R')
source('~/projects/rsignatures//src/utils/rs_gene_overlaps.R')
source('~/projects/rsignatures//src/utils/normalizeVector.R')
source('~/projects/rsignatures//src/utils/getBpGr.R')
source('~/projects/rsignatures//src/utils/performOverlaps.R')
source('~/projects/rsignatures//src/utils/getBpGr.R')
source('~/projects/rsignatures//src/utils/loadRepeats.R')
source('~/projects/rsignatures/src/utils/plotSignature.R')
source('~/projects/rsignatures/src/utils/getCategs.R')
source('~/repos/hotspots/R/utils/plotCoefficients.R')
# this is where the results will be stored
resultList <- list()
sample_metadata = read.csv('~/park_dglodzik/data_repo/PanCan/WGS.metadata.txt', sep='\t')
sample_metadata <- subset(sample_metadata, !duplicated(icgc_specimen_id))
rownames(sample_metadata) <- sample_metadata$aliquot_id
specimen_hist <- read.table('~/park_dglodzik/data_repo/PanCan//pcawg_specimen_histology_August2016_v7.tsv', sep='\t', header=FALSE)
sample_metadata.m <- merge(sample_metadata, specimen_hist, by.x='icgc_sample_id', by.y='V6')
rownames(sample_metadata.m) <- sample_metadata.m$aliquot_id
# load the reference data
replis.gr <- toGRanges('~/projects/rsignatures/data//external/BASIS/Repliseq_data//RepliTime/MCF7_data//MCF7_RepliSeq.bedGraph')
replis.gr <- replis.gr[lengths(replis.gr)<=1001]
replis.gr$chr.len<- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[as.character(seqnames(replis.gr))]
replis.gr<-replis.gr[(start(replis.gr)>chr.margin) & ((replis.gr$chr.len - end(replis.gr))>chr.margin)]

replis.gr$quantile <- qcut(replis.gr$V4, 10)
# Prepare replication origins, SNS-seq

core.origins <- read.table('~/projects/rsignatures/data/external/origins/akerman_core_origins_hg19.bed', sep='\t')
core.origins$midopoint <- rowMeans(core.origins[,c('V2', 'V3')])
common.origins.gr <- GRanges(seqnames=Rle(core.origins$V1),
                  ranges=IRanges(core.origins$midopoint, core.origins$midopoint))

common.origins.width.gr <- GRanges(seqnames=Rle(core.origins$V1),
                  ranges=IRanges(core.origins$V2, core.origins$V3))

core.origins$distPrev <- core.origins$V2 - c(NA,core.origins$V3[1:(nrow(core.origins)-1)])

core.origins$leftMargin <- core.origins$midopoint - 1e6
core.origins$rightMargin <- core.origins$midopoint + 1e6

common.origins.1Mb.gr <- GRanges(seqnames=Rle(core.origins$V1),
                  ranges=IRanges(core.origins$leftMargin, core.origins$rightMargin))

# repeats
repeat_fn <- '/home/dg204/park_dglodzik/ref/repeats/hg19/rmsk.txt.gz'
#repeats.gr <- loadRepeats(repeat_fn)
#save(repeats.gr, file='~/projects/rsignatures/data/processed/repeats.gr.RData')
load('~/projects/rsignatures/data/processed/repeats.gr.RData')

# sample groups of interest
ccne.samples.df <- read.csv('~/projects/rsignatures/data/processed/CCNE1.amp.sample.csv')
curated_cdk12 <- c('0009b464-b376-4fbc-8a56-da538269a02f',
                  '84ca6ab0-9edc-4636-9d27-55cdba334d7d',
                  'b243adb4-b3e7-4e0e-bc0d-625aa8dbb1be',
                   '89dad92e-5b3f-479a-a6da-a94ee7df7f8a',
                   'bc0dee07-de20-44d6-be65-05af7e63ac96',
                   '36d1a85e-a09b-4537-86e0-eaf1eb03aed8',
                   '0bfd1043-816e-e3e4-e050-11ac0c4860c5' # this one has two hits
                  )
hrd_brca_loh <- read.table('~/projects/rsignatures/data/processed/BRCA_and_LOH_PCAWG.tsv', header=TRUE)


# end of references
# load the sample RS signature attributions
# load the bedpe file 
load('~/projects/rsignatures//data/processed/sample.rearrs.RData')
# the the notebooks: rs.ipynb and rearr_catalogues.ipynb
length(sample.rearrs)
all.rearrs <- do.call('rbind',sample.rearrs)
all.rearrs.m <- merge(all.rearrs, sample_metadata.m, by.x='sample', by.y='aliquot_id', all.x=TRUE)
eval(parse(text = paste('test.bedpe <- subset(all.rearrs.m,', filter_str,')')))

# samples with most SVs
st <- sort(table(test.bedpe$sample), decreasing=TRUE)
top_df <- data.frame(no_svs=st)
colnames(top_df) <- c('sample', 'no_svs')
top_df$tissue <- sample_metadata.m[names(st),'V13']
write.csv(top_df, file=paste0(result_folder,exp.name,'_top_samples.csv'), row.names=FALSE, quote=FALSE)



test.bedpe$dcc_project_code <- sample_metadata[test.bedpe$sample,'dcc_project_code']
# filter rearrangements

#  pos.1, pos.2, Chromosome.1, Chromosome.2, Signature.SV1, Signature.SV6
test.bedpe$mid <- rowMeans(test.bedpe[,c('start1', 'start2')])
test.bedpe$length <- test.bedpe$start2 - test.bedpe$start1
class.maxlen <- quantile( test.bedpe$length, 0.80)

test.bedpe$chr1.len<- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[paste0('chr', test.bedpe$chrom1) ]
test.bedpe$chr2.len<- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[paste0('chr', test.bedpe$chrom2) ]
test.bedpe$is.len.common <- test.bedpe$length < class.maxlen

# shuffling of the SVs
if (shuffle) {
    test.bedpe$chr1.len<- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[paste0('chr', test.bedpe$chrom1) ]
    test.bedpe$chr2.len<- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[paste0('chr', test.bedpe$chrom2) ]
    test.bedpe$length <- test.bedpe$start2 - test.bedpe$start1
    test.bedpe$repliseq.new <- NA
    
    transloc <- test.bedpe$chrom1!= test.bedpe$chrom2
    if (sum(transloc)>0) {
        test.bedpe$start1[transloc] <- chr.margin + sapply(test.bedpe$chr1.len[transloc]-2*chr.margin, sample, 1) 
        test.bedpe$start2[transloc] <- chr.margin + sapply(test.bedpe$chr2.len[transloc]-2*chr.margin, sample, 1) 
    }
        
    not_transloc <- test.bedpe$chrom1==test.bedpe$chrom2

    if (shuffleMode == 'simple') { 
        test.bedpe$start1[!transloc] <- chr.margin + sapply(test.bedpe$chr1.len[!transloc]-2*chr.margin, sample, 1) 
        test.bedpe$start2[!transloc] <- test.bedpe$start1[!transloc] + test.bedpe$length[!transloc]
    } else if (shuffleMode == 'flip'){     
        sides <- sample(c("left", "right"), size = nrow(test.bedpe), replace = TRUE)
        leftFlag <- !transloc & sides=='left'
        test.bedpe$start2[leftFlag] <- test.bedpe$start1[leftFlag] 
        test.bedpe$start1[leftFlag] <- test.bedpe$start1[leftFlag] - test.bedpe$length[leftFlag]
        rightFlag <- !transloc & sides=='right'
        test.bedpe$start1[rightFlag] <- test.bedpe$start2[rightFlag] 
        test.bedpe$start2[rightFlag] <- test.bedpe$start2[rightFlag] + test.bedpe$length[rightFlag]

        
    } else {
        # obsolete: shuffling while preserving repliseq distribution
        #shuffle_right <- sample(c(TRUE, FALSE), size=nrow(test.bedpe), replace=TRUE)
        shuffle_right <- abs(test.bedpe$repliseq.sim.right - test.bedpe$repliseq.actual ) < abs(test.bedpe$repliseq.sim.left - test.bedpe$repliseq.actual )
        shuffle_right[is.na(shuffle_right)] <- sample(c(TRUE, FALSE), size=sum(is.na(shuffle_right)), replace=TRUE)
        
        test.bedpe$start1[not_transloc&shuffle_right] <- test.bedpe$start2[not_transloc&shuffle_right]
        test.bedpe$start2[not_transloc&shuffle_right] <- test.bedpe$start1[not_transloc&shuffle_right] + 
        test.bedpe$length[not_transloc&shuffle_right]
        test.bedpe$repliseq.new[not_transloc&shuffle_right] <- test.bedpe$repliseq.sim.right[not_transloc&shuffle_right]
    
        test.bedpe$start2[not_transloc&!shuffle_right] <- test.bedpe$start1[not_transloc&!shuffle_right]
        test.bedpe$start1[not_transloc&!shuffle_right] <- test.bedpe$start2[not_transloc&!shuffle_right] - 
        test.bedpe$length[not_transloc&!shuffle_right]
        test.bedpe$repliseq.new[not_transloc&!shuffle_right] <- test.bedpe$repliseq.sim.left[not_transloc&!shuffle_right]
    } 


        
    test.bedpe <- subset(test.bedpe, start2 < (chr2.len-2*chr.margin))
    test.bedpe <- test.bedpe[order(test.bedpe$repliseq.new),]
    # hist(subset(test.bedpe,chrom1==1)$start1)
    # maybe uneven repliseq profile results from chromosome distribution?
}

test.bedpe$is.proper <- (test.bedpe$start1>chr.margin) & (test.bedpe$start2 <(test.bedpe$chr1.len - chr.margin) ) & test.bedpe$is.len.common
test.bedpe.sel <- subset(test.bedpe, (chrom1==chrom2) & (is.proper==TRUE))
test.bedpe2 <- test.bedpe
test.bedpe2$sample <- 'X'
test.bedpe2$svclass[test.bedpe2$svclass=='duplication'] <- 'tandem-duplication'

catalogue <- bedpeToRearrCatalogue(test.bedpe2)
m <- as.matrix(rowSums(catalogue$rearr_catalogue))
colnames(m) <- 'all samples'


if (makePDFs) {
    pdf(paste0(result_folder,exp.name,'_rearr_catalogue_all.pdf'), width=10, height=3)
}
plotRearrSignatures(m) 
if (makePDFs) {
    dev.off()
}

r_dens <- density(log10(subset(test.bedpe, chrom1==chrom2)$length))

if (makePDFs) {
    pdf(paste0(result_folder,exp.name,'_length_density.pdf'), width=8, height=8)
}


# Kernel density plot
plot(r_dens, lwd = 2, xlab='rearrangement length', col='dark gray', xlim=c(2.5,8), main=paste(exp.name, 'length density'))

if (makePDFs) {
    dev.off()
}


# this used to be the MH section


# these are the segments that denote SV footprint
# add 3xSV margin
segments.gr <- GRanges(seqnames=Rle(paste0('chr', test.bedpe.sel$chrom1)),
          ranges=IRanges(test.bedpe.sel$start1 - 3*test.bedpe.sel$length, 
                         test.bedpe.sel$start2 + 3*test.bedpe.sel$length), seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19))
segments.gr <- trim(segments.gr)
segments.gr$rs.class <- test.bedpe.sel$max.Ref.Sig
segments.gr$length <- test.bedpe.sel$start2 - test.bedpe.sel$start1
segments.gr$midpoint = rowMeans(test.bedpe.sel[,c('start2', 'start1')])
segments.gr$rearr.start <- test.bedpe.sel$start1
segments.gr$rearr.end <- test.bedpe.sel$start2
segments.gr$sample <- test.bedpe.sel$sample
segments.no.margin.gr <- segments.gr
start(segments.no.margin.gr) <- test.bedpe.sel$start1
end(segments.no.margin.gr) <- test.bedpe.sel$start2
# SV segments + 1Mb "margim":
margim <- 1e6
segments.1Mb.gr <- GRanges(seqnames=Rle(paste0('chr', test.bedpe.sel$chrom1)),
                  ranges=IRanges(test.bedpe.sel$start1-margim, 
                                 test.bedpe.sel$start2+margim), seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19))
segments.1Mb.gr <- trim(segments.1Mb.gr)
segments.1Mb.gr$midpoint = rowMeans(test.bedpe.sel[,c('start2', 'start1')])
segments.1Mb.gr$rearr.start <- test.bedpe.sel$start1
segments.1Mb.gr$rearr.end <- test.bedpe.sel$start2
segments.1Mb.gr$rs.length <- segments.1Mb.gr$rearr.end-segments.1Mb.gr$rearr.start
segments.1Mb.gr$rs.class <- test.bedpe.sel$max.Ref.Sig
segments.1Mb.gr$rs.length <- segments.1Mb.gr$length
segments.1Mb.gr$start1 <- segments.1Mb.gr$rearr.start
segments.1Mb.gr$start2 <-  segments.1Mb.gr$rearr.end

# overlap of SVs (1Mb margin) with replication origins
o.breast.df.1Mb <- as.data.frame(findOverlaps(segments.1Mb.gr, common.origins.gr))

o.breast.df.1Mb$rs.midpoint <- segments.1Mb.gr$midpoint[o.breast.df.1Mb$queryHits]
o.breast.df.1Mb$origin <- start(common.origins.gr)[o.breast.df.1Mb$subjectHits]
o.breast.df.1Mb$rs.length <- segments.1Mb.gr$length[o.breast.df.1Mb$queryHits]
o.breast.df.1Mb$distanceToOrigin <- o.breast.df.1Mb$origin - o.breast.df.1Mb$rs.midpoint
o.breast.df.1Mb$isWithinSV <-  (start(common.origins.gr)[o.breast.df.1Mb$subjectHits] > segments.1Mb.gr$rearr.start[o.breast.df.1Mb$queryHits]) &
                                (end(common.origins.gr)[o.breast.df.1Mb$subjectHits] < segments.1Mb.gr$rearr.end[o.breast.df.1Mb$queryHits])
origin.v <- c(sum(o.breast.df.1Mb$isWithinSV==TRUE), sum(o.breast.df.1Mb$isWithinSV==FALSE))
trials.v <-  c(sum(segments.1Mb.gr$rearr.end - segments.1Mb.gr$rearr.start ),(length(segments.1Mb.gr) * 2 * margim))
comp.data <- matrix(c(origin.v, 
               trials.v), nrow = 2, byrow = TRUE)
# enrichment of replication origins compared to the margins
v <- origin.v/trials.v
RR <- v[1]/v[2]
pt <- prop.test(origin.v, trials.v)
resultList[['originEnrichmentES']] <- RR
resultList[['originEnrichmentP']] <- pt$p.value

SE_logRR <- sqrt((1-v[1])/origin.v[1] + (1-v[2])/origin.v[2])
CI_logRR <- log(RR) + c(-1,1)*1.96*SE_logRR
CI_RR <- exp(CI_logRR)
resultList[['originEnrichmentES_lower95']] <- CI_RR[1]
resultList[['originEnrichmentES_higher95']] <- CI_RR[2]

# Step 2: Perform Fisher's Exact Test
#result <- fisher.test(data)
# overlap of SVs (1Mb margin) with repli-seq
o.df.1Mb <- as.data.frame(findOverlaps(segments.1Mb.gr, replis.gr))
o.df.1Mb$rs.midpoint <- segments.1Mb.gr$midpoint[o.df.1Mb$queryHits]
o.df.1Mb$rearr.start <- segments.1Mb.gr$rearr.start[o.df.1Mb$queryHits]
o.df.1Mb$rearr.end <- segments.1Mb.gr$rearr.end[o.df.1Mb$queryHits]
o.df.1Mb$rs.length <- o.df.1Mb$rearr.end - o.df.1Mb$rearr.start
o.df.1Mb$repli_measure <- replis.gr$V4[o.df.1Mb$subjectHits]
o.df.1Mb$distGeneStartToMidpoint <-start(replis.gr)[o.df.1Mb$subjectHits] - o.df.1Mb$rs.midpoint 
o.df.1Mb$distGeneEndToMidpoint <- end(replis.gr)[o.df.1Mb$subjectHits] - o.df.1Mb$rs.midpoint 
o.df.1Mb$distGeneStartToMidpointScaled <- o.df.1Mb$distGeneStartToMidpoint
o.df.1Mb$distGeneEndToMidpointScaled <- o.df.1Mb$distGeneEndToMidpoint       
o.df.1Mb$distGeneStartToMidpointScaledInt <- round(o.df.1Mb$distGeneStartToMidpointScaled/1e4,0)
o.df.1Mb$distGeneEndToMidpointScaledInt <- round(o.df.1Mb$distGeneEndToMidpointScaled/1e4,0) 


minBin <- abs(min(o.df.1Mb$distGeneStartToMidpointScaledInt, na.rm=TRUE)) + 1
maxBin <- abs(max(o.df.1Mb$distGeneStartToMidpointScaledInt, na.rm=TRUE)) + 1

o.df.1Mb$distGeneStartToMidpointScaledIntPos <- o.df.1Mb$distGeneStartToMidpointScaledInt + minBin
o.df.1Mb$distGeneEndToMidpointScaledIntPos <- o.df.1Mb$distGeneEndToMidpointScaledInt + minBin
maxBinPos <- max(o.df.1Mb$distGeneEndToMidpointScaledIntPos)

cvg = coverage(IRanges(start=o.df.1Mb$distGeneStartToMidpointScaledIntPos, 
                                   end=o.df.1Mb$distGeneEndToMidpointScaledIntPos))

measure = coverage(IRanges(start=o.df.1Mb$distGeneStartToMidpointScaledIntPos, 
                                   end=o.df.1Mb$distGeneEndToMidpointScaledIntPos), 
                           weight=as.integer(o.df.1Mb$repli_measure))

options(repr.plot.width=8, repr.plot.height=8)
if (makePDFs) {
    pdf(file=paste0(result_folder, 'replication_association_',exp.name,'_1Mb.pdf'), width=8, height=8)
}    
par(mar = c(6, 6, 6, 6))
    h <- hist(o.breast.df.1Mb$distanceToOrigin, breaks=100, plot=FALSE)
    # Scale the counts
    h$counts <- h$counts/nrow(test.bedpe)
    # Generate pretty ticks based on new count values
    yticks <- pretty(h$counts)
    xticks <- c(-1e6, 0, 1e6)
    # Plot the histogram with improved y-axis
    plot(h, 
         main=paste(exp.name, round(median(test.bedpe.sel$length/1000),0), 'kb'), col='darkgray', border=NA, 
         xlab='Dist. to SV midpoint [Mb]', 
         ylab='Repl. origins per bin per SV', 
         cex.lab=2, cex.axis=2, cex.main=2, cex.sub=1.5,
         xlim=c(-1e6, 1e6),
         ylim=range(yticks),
     xaxs = "i", xaxt='n', yaxt='n')  # suppress both axes
    # Add custom y-axis
    axis(side=2, at=yticks, labels=yticks, cex.axis=2)
    # Add custom x-axis (no scientific notation)
    axis(side=1, at=xticks, labels=c('-1', '0', '+1'), cex.axis=2)
# Add vertical line at zero
abline(v=0, lty=2)
# repliseq line
par(new = TRUE)                            
 plot((1:length(cvg)-minBin)*1e4, measure/cvg,pch=19, ylim=c(30,60), 
    xlim=c(-1e6, 1e6), 
      xaxs="i", col='#7570b3',  
    cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlab='', ylab='', main='',axes=FALSE)
axis(side = 4, 
     at = c(30, 40, 50, 60), 
     cex.lab = 2, 
     cex.axis = 2, 
     col.axis = "#7570b3",  # tick label color
     col = "#7570b3")
mtext("Average Repli-seq signal", side = 4, line = 3,cex=2, col='#7570b3')     
if (makePDFs) {
    dev.off()
}

# the overlap data frame
# left is query, right is subject
# segments: SVs + 3 SV lengths to the left and to the right
o.breast.df <- as.data.frame(findOverlaps(segments.gr, common.origins.gr))
o.breast.df$rs.midpoint <- segments.gr$midpoint[o.breast.df$queryHits]
# note: for origins, we use the start position
o.breast.df$origin <- start(common.origins.gr)[o.breast.df$subjectHits]
o.breast.df$rs.length <- segments.gr$length[o.breast.df$queryHits]
o.breast.df$distanceToOrigin <- o.breast.df$origin - o.breast.df$rs.midpoint

# o.breast.3widths.df$distanceToOrigin/o.breast.3widths.df$rs.length
o.breast.df$distanceToOrigin.scaled <- o.breast.df$distanceToOrigin/o.breast.df$rs.length
o.breast.df$len.inv.weight <- 1/o.breast.df$rs.length 

o.breast.df$rs.start1 <- segments.gr$rearr.start[o.breast.df$queryHits]
o.breast.df$rs.start2 <- segments.gr$rearr.end[o.breast.df$queryHits]

if (marginmode=='3len') {
    o.breast.df$distanceToOrigin.scaled <- o.breast.df$distanceToOrigin/o.breast.df$rs.length
} else if (marginmode=='1Mb') {
    # really, this is unscaled
    o.breast.df$distanceToOrigin.scaled <- o.breast.df$distanceToOrigin
}
o.breast.df$rs.class <- segments.gr$rs.class[o.breast.df$queryHits]
o.breast.df$len.inv.weight <- 1/o.breast.df$rs.length

# overlap with Repli-seq (margin flexible, depending on SV length)
o.df <- as.data.frame(findOverlaps(segments.gr, replis.gr))
o.df$rs.midpoint <- segments.gr$midpoint[o.df$queryHits]
o.df$rearr.start <- segments.gr$rearr.start[o.df$queryHits]
o.df$rearr.end <- segments.gr$rearr.end[o.df$queryHits]
o.df$rs.length <- o.df$rearr.end - o.df$rearr.start
o.df$rs.class <- segments.gr$rs.class[o.df$queryHits]
o.df$repli_measure <- replis.gr$V4[o.df$subjectHits]
no.bins <- 100 # per one unit: 1Mb or SV length
if (marginmode=='3len') {
    o.df$distGeneStartToMidpoint <-start(replis.gr)[o.df$subjectHits] - o.df$rs.midpoint 
    o.df$distGeneEndToMidpoint <- end(replis.gr)[o.df$subjectHits] - o.df$rs.midpoint 
    o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint/o.df$rs.length 
    o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint/o.df$rs.length
    o.df$distGeneStartToMidpointScaledInt <- round(o.df$distGeneStartToMidpointScaled,2)*no.bins
    o.df$distGeneEndToMidpointScaledInt <- round(o.df$distGeneEndToMidpointScaled,2) * no.bins
} else if (marginmode=='1Mb') {
    o.df$distGeneStartToMidpoint <-start(replis.gr)[o.df$subjectHits] - o.df$rs.midpoint 
    o.df$distGeneEndToMidpoint <- end(replis.gr)[o.df$subjectHits] - o.df$rs.midpoint 
    o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint
    o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint       
    o.df$distGeneStartToMidpointScaledInt <- round(o.df$distGeneStartToMidpointScaled/1e4,0)
    o.df$distGeneEndToMidpointScaledInt <- round(o.df$distGeneEndToMidpointScaled/1e4,0) 
} else if (marginmode=='1Mb.left') {
     o.df$distGeneStartToMidpoint <-start(replis.gr)[o.df$subjectHits] - o.df$rs.start1 
    o.df$distGeneEndToMidpoint <- end(replis.gr)[o.df$subjectHits] - o.df$rs.start1 
    o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint
    o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint       
    o.df$distGeneStartToMidpointScaledInt <- round(o.df$distGeneStartToMidpointScaled/1e4,0)
    o.df$distGeneEndToMidpointScaledInt <- round(o.df$distGeneEndToMidpointScaled/1e4,0) 
} else if (marginmode=='1Mb.right') {
     o.df$distGeneStartToMidpoint <-start(replis.gr)[o.df$subjectHits] - o.df$rs.start2 
    o.df$distGeneEndToMidpoint <- end(replis.gr)[o.df$subjectHits] - o.df$rs.start2 
    o.df$distGeneStartToMidpointScaled <- o.df$distGeneStartToMidpoint
    o.df$distGeneEndToMidpointScaled <- o.df$distGeneEndToMidpoint       
    o.df$distGeneStartToMidpointScaledInt <- round(o.df$distGeneStartToMidpointScaled/1e4,0)
    o.df$distGeneEndToMidpointScaledInt <- round(o.df$distGeneEndToMidpointScaled/1e4,0) 
}
minBin <- abs(min(o.df$distGeneStartToMidpointScaledInt, na.rm=TRUE)) + 1
maxBin <- abs(max(o.df$distGeneStartToMidpointScaledInt, na.rm=TRUE)) + 1
o.df$distGeneStartToMidpointScaledIntPos <- o.df$distGeneStartToMidpointScaledInt + minBin
o.df$distGeneEndToMidpointScaledIntPos <- o.df$distGeneEndToMidpointScaledInt + minBin
maxBinPos <- max(o.df$distGeneEndToMidpointScaledIntPos)

 # counts of repliseq segments overlapping with each bin
cvg = coverage(IRanges(start=o.df$distGeneStartToMidpointScaledIntPos, 
                                   end=o.df$distGeneEndToMidpointScaledIntPos)
                           )
# sum of repliseq at a given position
measure = coverage(IRanges(start=o.df$distGeneStartToMidpointScaledIntPos, 
                                   end=o.df$distGeneEndToMidpointScaledIntPos), 
                           weight=as.integer(o.df$repli_measure))
# measure/cvg gives average repliseq signal at a given position
resultList[['measure']] <- measure
resultList[['cvg']] <- cvg
resultList[['minBin']] <- minBin
resultList[['repliseq_o.df']] <- o.df

# a histogram of replication origins
# where bins correspond to rescaled distance from SV midpoint to replication origin
options(repr.plot.width=8, repr.plot.height=8)
if (makePDFs) {
    pdf(file=paste0(result_folder, 'replication_association_',exp.name,'.pdf'), width=8, height=8)
}    
# number of bins in the plotted regions will depend on how much is not shown because they fall out of the -3.5, 3.5 range
par(mar = c(6, 6, 6, 6))
    h <- hist(o.breast.df$distanceToOrigin.scaled, breaks=40, plot=FALSE)
    # Scale the counts
    h$counts <- h$counts/nrow(test.bedpe)
    yticks <- pretty(h$counts)
    xticks <- c(-0.5, 0, 0.5)
    plot(h, 
         col = 'darkgray', border = NA,
         xlab = 'Dist. to SV midpoint [SV lengths]',
         ylab = 'Repl. origins per bin per SV',
         cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 1.5,
         xlim = c(-3.5, 3.5),
         main = paste(exp.name, round(median(test.bedpe.sel$length / 1000), 0), 'kb'),
         xaxs = "i",
         xaxt = "n")  # suppress x-axis

    axis(side=2, at=yticks, labels=yticks, cex.axis=2)
    axis(side=1, at=xticks, labels=c('-.5', '0', '+.5'), cex.axis=2)
    abline(v=-.5, lwd=6, col="#4EA72E")
    abline(v=.5, lwd=6 ,col="#4EA72E")
    abline(v=0, lty=2)
    par(new = TRUE)                             # Add new plot
       plot((1:length(cvg)-minBin)/100, measure/cvg,  pch=19, ylim=c(30,60), 
        xlim=c(-3.5, 3.5), 
        xaxs="i",  ylab='',
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xlab='' ,axes=FALSE, col='#7570b3', main='')
    
    axis(side = 4, 
     at = c(30, 40, 50, 60), 
     cex.lab = 2, 
     cex.axis = 2, 
     col.axis = "#7570b3",  # tick label color
     col = "#7570b3")
mtext("Average Repli-seq signal", side = 4, line = 3,cex=2, col='#7570b3')     
if (makePDFs) {
    dev.off()
}


# gene expression analysis
# expr.m, disease.m
load('~/projects//rsignatures/data/processed/disease.RData')
cancer_types_tab <- table(gsub('-.*', '', test.bedpe$dcc_project_code))
cancer_types_tab <- cancer_types_tab[names(cancer_types_tab) %in% colnames(disease.m)]
if (length(cancer_types_tab)>0) {
    expr.tissue <- names(cancer_types_tab)[which.max(cancer_types_tab)]
} else {
    print('No matching tissue found')
    expr.tissue <- 'BRCA' 
}

# load the genes
print(paste('Tissue for expression analysis:', expr.tissue))
gene.table.fn <- '~/repos/hotspots/data/genes.table.csv'
genes.table <- as.data.frame(read.csv(gene.table.fn, row.names=NULL)[,2:7])
disease.df <- as.data.frame(disease.m)
disease.df$ensid <- gsub('\\..*', '', rownames(disease.m))
disease.m2 <- merge(disease.df, genes.table, all.x=TRUE,  by.x='ensid', by.y='ensembl_gene_id')
disease.m2 <- subset(disease.m2, !is.na(hgnc_symbol) & (chr %in% c(1:22, 'X')))
disease.m2$footprint <- disease.m2$chromEnd - disease.m2$chromStart
disease.m2.bed <- disease.m2[,c('chr', 'chromStart', 'chromEnd', 'strand', 'ensid', 'OV', 'PRAD', 'BRCA', 'hgnc_symbol')]
disease.m2.bed$chr <- paste0('chr',disease.m2.bed$chr )
genes_high_disease <- subset(disease.m2, disease.m2[,expr.tissue] > quantile(disease.m2[,expr.tissue], 0.5))
genes_low_disease <- subset(disease.m2, disease.m2[,expr.tissue] < quantile(disease.m2[,expr.tissue], 0.5))
segments.gr$rs.length <- segments.gr$length
segments.gr$start1 <- segments.gr$rearr.start
segments.gr$start2 <-  segments.gr$rearr.end

# all the gene plots (3x3)
if (makePDFs) {
    pdf(file=paste0(result_folder, 'genes_bps_',exp.name,'.pdf'), width=24, height=24)
    par(mfrow = c(3, 3))

} else {
    options(repr.plot.width = 12, repr.plot.height = 12)
    par(mfrow = c(3, 3))
}
# add SVs + 1Mb margin 
# perform overlaps
# for each gene, calculate distance from the gene to the midpoint
# come up with new units
# calculate coverage of genes in this new frame of reference, weigh by SV size
r1<- rs_gene_overlaps(genes_high_disease, segments.1Mb.gr, rs.groups=NULL, marginmode='1Mb',  ylab='highly expressed genes', no.bins=100,  plot_tile=paste('highly expressed: in ', expr.tissue))
r2 <- rs_gene_overlaps(genes_low_disease, segments.gr, rs.groups=NULL, marginmode=marginmode,  ylab='low expressed genes', no.bins=100,  plot_tile=paste('low expressed: in ', expr.tissue),
                      add.newplot = TRUE)

r3 <- rs_gene_overlaps(genes_high_disease,  segments.gr, rs.groups=NULL, marginmode=marginmode,  ylab='highly expressed genes', no.bins=100,  plot_tile=paste('high expressed in:', expr.tissue)) #ymax=2.5,
r4 <- rs_gene_overlaps(subset(genes_high_disease, strand==-1), segments.gr, rs.groups=NULL, marginmode=marginmode,  ylab='highly expressed genes', no.bins=100,  plot_tile=paste('high expressed - in:', expr.tissue)) #ymax=2.5,
r5 <- rs_gene_overlaps(subset(genes_high_disease, strand==1), segments.gr, rs.groups=NULL, marginmode=marginmode,  ylab='highly expressed genes', no.bins=100,  plot_tile=paste('high expressed + in:', expr.tissue)) #ymax=2.5,
# the above counts genes
# the below adds their expression
evalExprOverlap(test.bedpe, segments.gr, disease.m2, expr.tissue=expr.tissue, plot.tile=paste(exp.name, 'expression around SVs'))
evalExprOverlap(test.bedpe, segments.gr, subset(disease.m2, strand==-1), expr.tissue, plot.tile=paste(exp.name, 'expression around SVs (-)'))
evalExprOverlap(test.bedpe, segments.gr, subset(disease.m2, strand==1), expr.tissue, plot.tile=paste(exp.name, 'expression around SVs (+)'))
disease.m2$footprint <- disease.m2$chromEnd - disease.m2$chromStart
disease.m2.long <- subset(disease.m2,  footprint > quantile(disease.m2$footprint, 0.75) )
evalExprOverlap(test.bedpe, segments.gr, disease.m2.long, expr.tissue, plot.tile=paste(exp.name, 'expression around SVs (long)'))
if (makePDFs) {
    dev.off()
}

# gene plots (2,1)
if (makePDFs) {
    pdf(file=paste0(result_folder, 'genes_high_expr_bps_',exp.name,'.pdf'), width=8, height=8)
    par(mar = c(6, 6, 6, 6))
    #par(mfrow = c(2, 1))
}
#r1<- rs_gene_overlaps(genes_high_disease, segments.1Mb.gr, rs.groups=NULL, marginmode='1Mb',  
#                      ylab=paste0('highly expressed genes per SV, ',expr.tissue), 
#                      no.bins=100,  
#                      plot_tile=paste(nrow(test.bedpe), 'SVs'))

ymax <- max(r3[['max_y']], r2[['max_y']])+ 0.01
ymax <- max(0.4,ymax)

r2 <- rs_gene_overlaps(genes_low_disease, segments.gr, rs.groups=NULL, marginmode=marginmode,  ylab=paste0('High expr. genes (',expr.tissue,') per bin per SV'), no.bins=100,  plot_tile=exp.name,
                      add.newplot = TRUE, ymax=ymax, lineCol='gray')

r3 <- rs_gene_overlaps(genes_high_disease,  
                       segments.gr, 
                       rs.groups=NULL, 
                       marginmode=marginmode,  
                       plot_tile=exp.name,
                       ylab=paste0('High expr. genes (',expr.tissue,') per bin per SV'), no.bins=100, ymax=ymax,  add.newplot = FALSE)

resultList[['expression_profile_high']] <- r2
resultList[['expression_profile_low']] <- r3
resultList[['segments.gr']] <- segments.gr

if (makePDFs) {
    dev.off()
}    

# repliseq at bp1, midpoint, bp2
# replication time at breakpoints
test.bedpe$sv_id_global <- paste(test.bedpe$sample, test.bedpe$sv_id)
rownames(test.bedpe) <- test.bedpe$sv_id_global 
bps.temp.gr <- getBpGr(test.bedpe, addMidpoint=TRUE)
bps_repliseq.bps.df <- as.data.frame(findOverlaps(bps.temp.gr, replis.gr))
bps_repliseq.bps.df$sv_id_global <- bps.temp.gr$sv_id_global[bps_repliseq.bps.df$queryHits]
bps_repliseq.bps.df$sv_id_global <- gsub('_.*', '', bps_repliseq.bps.df$sv_id_global)
bps_repliseq.bps.df$repliseq <- replis.gr$V4[bps_repliseq.bps.df$subjectHits]
# type can be bp1, bp2, midpoint
bps_repliseq.bps.df$type <- bps.temp.gr$type[bps_repliseq.bps.df$queryHits]
# Aggregate using mean (or another function like median, min, max, etc.)
agg_df <- aggregate(repliseq ~ sv_id_global + type, data = bps_repliseq.bps.df, FUN = mean)
# Now reshape (safe since there are no duplicates)
wide_df <- reshape(
  agg_df,
  timevar = "type",
  idvar = "sv_id_global",
  direction = "wide"
)
wbp1 <- wilcox.test(wide_df$repliseq.bp1, wide_df$repliseq.midpoint)
resultList[['repliseq_midpoint_bp1_wilcoxonP']] <- wbp1$p.value
repliseq_midpoint_bp1_mean_diff <- mean(wide_df$repliseq.midpoint - wide_df$repliseq.bp1, na.rm=TRUE)
resultList[['repliseq_midpoint_bp1_mean_diff']] <- repliseq_midpoint_bp1_mean_diff
wbp2 <- wilcox.test(wide_df$repliseq.bp2, wide_df$repliseq.midpoint)
resultList[['repliseq_midpoint_bp2_wilcoxonP']] <- wbp2$p.value
repliseq_midpoint_bp2_mean_diff <-  mean(wide_df$repliseq.midpoint - wide_df$repliseq.bp2, na.rm=TRUE)
resultList[['repliseq_midpoint_bp2_mean_diff']] <- repliseq_midpoint_bp2_mean_diff
wbps <- wilcox.test(c(wide_df$repliseq.bp1, wide_df$repliseq.bp2), wide_df$repliseq.midpoint)
resultList[['repliseq_midpoint_both_bp2s_wilcoxonP']] <- wbps$p.value

# gene stats
# checking genes overlapping with breakpoints only
bps.temp.gr <- getBpGr(test.bedpe, addMidpoint=TRUE)
genes.high.gr <- GRanges(seqnames=Rle(paste0('chr',genes_high_disease$chr)),
              ranges=IRanges(genes_high_disease$chromStart,genes_high_disease$chromEnd))
genes.high.gr$strand <- genes_high_disease$strand
bps_genes.df <- as.data.frame(findOverlaps(bps.temp.gr, genes.high.gr))
bps_genes.df$type <- bps.temp.gr$type[bps_genes.df$queryHits]
bps_genes.df$strand <- genes.high.gr$strand[bps_genes.df$subjectHits]
bps_genes.df <- subset(bps_genes.df, type!='midpoint')
bp_strand_table <- table(bps_genes.df$strand, bps_genes.df$type)
strand_t <- fisher.test(bp_strand_table)
# effect size between whether a gene overlaps bp1, bp2 and gene strand
resultList[['gene_strand_bp_P']] <- strand_t$p.value
resultList[['gene_strand_bp_estimate']] <- strand_t$estimate[['odds ratio']]
resultList[['gene_strand_bp_lower']] <- strand_t$conf.int[1]
resultList[['gene_strand_bp_higher']] <- strand_t$conf.int[2]
r3[['o.df']]$orientationVsMidpoint <- NA
r3[['o.df']]$orientationVsMidpoint[r3[['o.df']]$distGeneEndToMidpoint<0] <- 'geneLeftOfMidpoint'
r3[['o.df']]$orientationVsMidpoint[r3[['o.df']]$distGeneStartToMidpointScaled>0] <- 'geneRightOfMidpoint'
midpoint_strand_table <- table(r3[['o.df']]$orientationVsMidpoint, r3[['o.df']]$geneStrand)
strand_midpoint_gene <- fisher.test(midpoint_strand_table)
resultList[['gene_strand_sv_P']] <- strand_midpoint_gene$p.value
resultList[['gene_strand_sv_estimate']] <- strand_midpoint_gene$estimate[['odds ratio']]
resultList[['gene_strand_sv_lower']] <- strand_midpoint_gene$conf.int[1]
resultList[['gene_strand_sv_higher']] <- strand_midpoint_gene$conf.int[2]
resultList[['median_length']] <- median(test.bedpe.sel$length/1000)

test.bedpe$sv_id_global <- paste0(test.bedpe$sample, '_', test.bedpe$id)
rownames(test.bedpe) <- test.bedpe$sv_id_global
bps.temp.gr <- getBpGr(test.bedpe)
nrow(test.bedpe)

expr_tissue <- 'BRCA'
source('~/projects/rsignatures//src/utils/loadGRs.R')

bps.temp.gr <- getBpGr(test.bedpe)
ovl_r <- performOverlaps(test.bedpe, bps.temp.gr, replis.gr, rfd.gr, genes.gr )




# repeat elements
test.bedpe$bp1_repeat_class <- 'no'
test.bedpe$bp2_repeat_class <- 'no'
test.bedpe$bp1_repID <-  NA
test.bedpe$bp2_repID <- NA

repeat_overlaps.df <- as.data.frame(findOverlaps(bps.temp.gr, repeats.gr))
repeat_overlaps.df$sv_id_global <- bps.temp.gr$sv_id_global[repeat_overlaps.df$queryHits]
repeat_overlaps.df$repFamily <- repeats.gr$repFamily[repeat_overlaps.df$subjectHits]
#sort(table(repeat_overlaps.df$repFamily )/(length(bps.temp.gr)))
repeat_overlaps.df_bp1 <- subset(repeat_overlaps.df, grepl('_bp1', sv_id_global))
repeat_overlaps.df_bp2 <- subset(repeat_overlaps.df, grepl('_bp2',sv_id_global))
test.bedpe[gsub('_bp.*', '',repeat_overlaps.df_bp1$sv_id_global),'bp1_repeat_class'] <- repeat_overlaps.df_bp1$repFamily
test.bedpe[gsub('_bp.*', '',repeat_overlaps.df_bp2$sv_id_global),'bp2_repeat_class'] <- repeat_overlaps.df_bp2$repFamily
test.bedpe[gsub('_bp.*', '',repeat_overlaps.df_bp1$sv_id_global),'bp1_repID'] <- repeat_overlaps.df_bp1$subjectHits
test.bedpe[gsub('_bp.*', '',repeat_overlaps.df_bp2$sv_id_global),'bp2_repID'] <- repeat_overlaps.df_bp2$subjectHits

test.bedpe$bp1_isL1 <- test.bedpe$bp1_repeat_class=='L1'
test.bedpe$bp2_isL1 <- test.bedpe$bp2_repeat_class=='L1'
test.bedpe_reps <-subset(test.bedpe, bp1_repeat_class!='no' & bp2_repeat_class!='no' & bp1_repID!=bp2_repID) # ensure that this is not a single repeat at both breakpoints
m_rep <- table(test.bedpe$bp1_repeat_class, test.bedpe$bp2_repeat_class)
m_l1_rep <- table(test.bedpe$bp1_isL1, test.bedpe$bp2_isL1)
rep_fisher <- fisher.test(m_l1_rep)
resultList[['repeatOverlapBpProp']] <- nrow(repeat_overlaps.df)/(2*nrow(test.bedpe))
#sort(table(repeat_overlaps.df$repFamily )/(length(bps.temp.gr)))
resultList[['repeatOverlapCoincidence']] <- sum(test.bedpe_reps$bp1_repeat_class==test.bedpe_reps$bp2_repeat_class)/nrow(test.bedpe_reps)
resultList[['repeatL1FisherEstimate']] <- rep_fisher$estimate
resultList[['repeatL1FisherEstimate_lower']] <- rep_fisher$conf.int[1]
resultList[['repeatL1FisherEstimate_higher']] <- rep_fisher$conf.int[2]
resultList[['repeatL1FisherPvalue']] <- rep_fisher$p.value

if ('t_rfd_2d' %in%  names(ovl_r)) {
    t_rfd_2d <- ovl_r[['t_rfd_2d']]
    fn <- paste0(result_folder, exp.name,'_RFD.csv')
    write.csv(t_rfd_2d, file=fn)
}

if ('t_repli_2d' %in% names(ovl_r)) {
    t_repli_2d <- ovl_r[['t_repli_2d']]
    fn <- paste0(result_folder,exp.name,'_repliTiming.csv')
    write.csv(t_repli_2d, file=fn)
}

if (('t_expr_2d' %in% names(ovl_r)) && ('t_rfd_2d' %in% names(ovl_r))) {
    t_expr_2d <- ovl_r[['t_expr_2d']]
    fn <- paste0(result_folder, exp.name,'_expression.csv')
    if (makePDFs) {
        pdf(file=paste0(result_folder, 'catalogue_rfd_',exp.name,'.pdf'), width=20, height=8)
    }    
    plotSignature(colSums(ovl_r[['t_rfd_2d']]),c('L_L', 'L_R', 'R_L', 'R_R')) 
    df <- as.data.frame(colSums(ovl_r[['t_rfd_2d']]))
    colnames(df) <- c( 'count')
    write.csv(df, file=paste0(result_folder, 'rfd_counts_',exp.name,'.csv'))
    if (makePDFs) {
        dev.off()
    }
}
resultList[['ovl_r']] <- ovl_r


# repeat elements


if (doRegression) {
    bps.df <- rbind(data.frame(chr=test.bedpe$chrom1, position=test.bedpe$start1, position=test.bedpe$start1),
                   data.frame(chr=test.bedpe$chrom2, position=test.bedpe$start2, position=test.bedpe$start1)
               )
    
    load('~/repos/hotspots//data/bins.regression.1000_OV.RData')
    # in future runs, model.variables and 
    #model.variables <- colnames(allBins)[ 4:ncol(allBins)]
    model.variables <- model.variables[-which(model.variables=='meanCn' | model.variables=='gaps'  )]
    
    bins.regression$exprLog <- log(bins.regression$expr + 1)
    bins.regression$exprLog.norm <- normaliseVector(bins.regression$exprLog)
    
    bpList=list()
    bpList[[1]] <- bps.df
    names(bpList)[1] <- 'test.bps'

    allBins2 <- prepareBinData(
        bedList=c(), 
        bpList=bpList, 
        cnDf=NULL , 
        repTimingBed=NULL,
        addSummary=FALSE,
        binSize=1e3
        )
    
    bins.regression.m <- merge(bins.regression, allBins2, by=c('chr', 'chromStart', 'chromEnd'))

    formula.test <- as.formula(paste('test.bps ~ ' , paste(paste0(model.variables, '.norm'), collapse=' + ')))

    bins.regression.mat <- as.matrix(bins.regression.m[,4:ncol(bins.regression.m)])
    
    print('Fitting regression model...')
    fm_nbin.test <- glm.nb(formula.test , data = bins.regression.m ) # negative binomial fit
    
    save(fm_nbin.test, model.variable.names, 
         file=paste0(result_folder,'/regression_', exp.name, '.RData'))
    
    if (makePDFs) {
    pdf(paste0(result_folder,exp.name,'regression_topography.pdf'), width=8, height=8)
    }
    options(repr.plot.width=8, repr.plot.height=8)
    par(mar = c(5, 10, 5, 5))
    c <- plotCoefficients(fm_nbin.test, nbfit.ref=NULL, exp.name, model.variable.names, reorder=TRUE, 
                          myxlim=c(0,2.5))
    if (makePDFs) {
        dev.off()
    }

}

save(resultList, file=paste0(result_folder, 'stats_',exp.name,'.RData'))
resultList[['ovl_r']] <- NULL
resultList[['expression_profile_low']] <- NULL
resultList[['expression_profile_high']] <- NULL
resultList[['cvg']] <- NULL
resultList[['measure']] <- NULL
resultList[['repliseq_o.df']] <- NULL
resultList[['expression_profile_high']] <- NULL
resultList[['expression_profile_low']] <- NULL
resultList[['segments.gr']] <- NULL
save(resultList, file=paste0(result_folder, 'stats_light_',exp.name,'.RData'))

# mh analysis
# microhomology analysis based on Hartwig data
# 
load('~/projects/rsignatures//data/processed/sample.rearrs.hartwig.RData')
all.rearrs.m <- merge(sample.rearrs.df, sample_metadata.m, by.x='sample', by.y='aliquot_id', all.x=TRUE)
eval(parse(text = paste('test.bedpe <- subset(all.rearrs.m,', filter_str,')')))

t_mh <- table(test.bedpe$sample, test.bedpe$label, test.bedpe$HOMLEN_lbl)
t_mh_2d <- collapseTable3D2D(t_mh)
order1 <- c('clustered', 'non-clustered')
order2 <- c('del',  'tds','inv', 'trans')
order3 <- c('1-10Kb', '10-100Kb', '100Kb-1Mb', '1Mb-10Mb','>10Mb')
order4 <- c('0','1', '2', '3', '4', '>=5' )
totalCats <- length(order1) * length(order2) * length(order3) * length(order4)   
categs <- paste0( rep(order1, each=totalCats/length(order1)), '_', 
                 rep(rep(order2, each=totalCats/length(order1)/length(order2)), length(order1)), '_',
                rep(rep(order3, each=totalCats/length(order1)/length(order2)/length(order3)), length(order1) * length(order2)), '_',
                rep(order4,length(order1) * length(order2) * length(order3))
                )
categs[grepl( 'trans_.*b_',categs)] <- gsub('trans_.*b_','trans_',categs[grepl( 'trans_.*b_',categs)])
categs <- categs[!duplicated(categs)]
categs <- categs[categs %in% colnames(t_mh_2d)]
t_mh_2d <- t_mh_2d[,categs]
cols <- rep('gray', ncol(t_mh_2d))
cols[grepl('tds',colnames(t_mh_2d))] <- 'darkgreen'
cols[grepl('del',colnames(t_mh_2d))] <- 'darkred'
cols[grepl('inv',colnames(t_mh_2d))] <- 'darkblue'
cols[grepl('trans',colnames(t_mh_2d))] <- 'purple'
if (makePDFs) {
    pdf(paste0(result_folder,exp.name,'_mh_catalogue.pdf'), width=20, height=8)
    par(mar = c(12.1, 4.1, 4.1, 2.1))
}
    bp <- barplot(colSums(as.matrix(t_mh_2d)) , las=2, cex.names=0.7, border=NA, col=cols, ylab='total number of SVs', main=paste(exp.name, 'microhomology'))
    for (i in seq(from=length(order4),to = ncol(t_mh_2d), by =length(order4))) {
        abline(v=bp[i]+0.5)
    }
if (makePDFs) {
    dev.off()
}

if (nrow(top_df)>0) {
    if (makePDFs) {
        pdf(paste0(result_folder,exp.name,'_mh_catalogue_top10.pdf'), width=20, height=8)
        par(mar = c(12.1, 4.1, 4.1, 2.1))
    }   
    for (i in 1:min(nrow(top_df), 10)) {
        if (top_df$sample[i] %in% rownames(t_mh_2d)) {
            bp <- barplot(as.matrix(t_mh_2d)[as.character(top_df$sample[i]),] , las=2, cex.names=0.7, border=NA, col=cols, ylab='total number of SVs', main=paste(top_df$sample[i], top_df$tissue[i],'microhomology'))
            for (i in seq(from=length(order4),to = ncol(t_mh_2d), by =length(order4))) {
                abline(v=bp[i]+0.5)
            }
        }
    }
    if (makePDFs) {
        dev.off()
    }

    arr <- unclass(t_mh)
    fn <- paste0(result_folder,exp.name,'_mh.h5')
    if (file.exists(fn)) file.remove(fn)
    h5createFile( fn)
    h5write(arr, file=fn, name="t_mh")

    # micro-homology vector
    v_mh <- as.vector(test.bedpe$HOMLEN)
    fn <- paste0(result_folder,exp.name,'_mh_v.h5')
    if (file.exists(fn)) file.remove(fn)
    h5createFile( fn)
    h5write(v_mh, file=fn, name="v_mh")
    
}

t_mh_1d <- table(test.bedpe$HOMLEN_lbl)
if (makePDFs) {
    pdf(paste0(result_folder,exp.name,'_mh_catalogue_1D.pdf'), width=3, height=4)
    par(mar = c(3.1, 3.1, 3.1, 2.1))
}
t_mh_1d <- t_mh_1d[order4]
t_mh_1d[is.na(t_mh_1d)] <- 0
barplot(t_mh_1d[order4], horiz=TRUE, border=NA, main=exp.name, las=2)
if (makePDFs) {
    dev.off()
}
write.csv(t_mh_1d, file=paste0(result_folder,exp.name,'_mh_catalogue_1D.csv'),row.names=FALSE)
# mh analysis end




