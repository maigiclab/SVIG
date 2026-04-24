
suppressPackageStartupMessages({
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library("GenomicRanges")
library(regioneR)
library(stringr)
library(lmtest)
library(MASS)
library(signature.tools.lib)
library(plyr)
library(readxl)
require(MutationTimeR) #
library(reshape2)
})    

source('utils/plotCoefficients.R')
source('utils/prepareBinData.R')
source('utils/performSVOverlaps.R') # the main worker function


if (interactive()) {
    # default values as an example and for debugging
    is.clustered <- FALSE
    svclass <- "duplication"
    size_min <- -1
    size_max <- 3e5
    max_sig <- 'Ref.Sig.R1'
    sample_subset <- 'CDK12'
    output_folder <- '../data/processed/'
} else {

    option_list <- list(
      make_option(c( "--svclass"), type = "character", help = "SV class"),
      make_option(c("--size_min"), type = "integer", help = "Minimum size"),
      make_option(c("--size_max"), type = "integer", help = "Maximum size"),
      make_option(c( "--output"), type = "character", help = "Output folder"),
      make_option(c( "--max_sig"), type = "numeric", help = "RS signature (eg. Ref.Sig.R1)"),
      make_option(c("--subset"), type = "character", help = "Sample subset, eg 'CCNE1', 'CDK12', 'HRD'; can be NULL")
    )
    # parse the command arguments
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    svclass <- opt$svclass
    size_min <- opt$size_min
    size_max <- opt$size_max
    output_folder <- opt$output
    max_sig <- if (opt$max_sig == "null") NULL else opt$max_sig
    sample_subset <- opt$subset
    # by default, only analyze non-clustered SVs
    is.clustered <- FALSE
}

# key settings
simulate <- FALSE # the naive generation of synthetic data is obsolete
margin.size <- 1e5 # size of the region around SVs
loadTimeR <- TRUE

# inputs:
pcawg_muts_folder <- '../data/interim/ICGC/SNV_indel_calls/final_consensus_12oct_passonly/snv_mnv/'
sample_metadata_fp <- '../data/interim/WGS.metadata.txt'
sv_data_fp <- '../data/interim/sample.rearrs.RData'
# TimeR result path
# this is missing the context field
# we will be loading sample files in /data/*.RData folder, which come from TimeR
# see projects/brca_timing/src/script_hrd_timer_pcawg.sh
muts_timer_path <- "../data/interim/TimeR_bb_all_sigs/"
sbs_signatures_pcawg_fp <- '../data/interim/SigProfilier_PCAWG_WGS_probabilities_SBS.csv'
sbs_signatures_pcawg_exposures_fp <- '../data/interim/PCAWG_sigProfiler_SBS_signatures_in_samples.csv'
sv_sigs_fp <- '../data/interim/43018_2020_27_MOESM3_ESM.xlsx'
# sample lists
ccne1_sample_list_fp <- '../data/interim/CCNE1.amp.sample.csv'
hrd_sample_list_fp <- '../data/interim/BRCA_and_LOH_PCAWG.tsv'
# others
mappability_fp <- '../data/interim/hg19.CRC.100mer.bed'


# load the metadata
sample_metadata = read.csv(sample_metadata_fp, sep='\t')
sample_metadata <- subset(sample_metadata, !duplicated(icgc_specimen_id))

# make a name for the experiment
exp.name <- paste0('PCAWG_', svclass,'_',size_min, '_', as.integer(size_max), '_', max_sig) 
if (!is.null(sample_subset)) {
    exp.name <- paste0(exp.name, '_', sample_subset)
}
test.samples.df <- sample_metadata
# this is used if maxSig is null



#notebooks/signature calc/rearr_catalogue
load(sv_data_fp)

# SBS signatures and exposures
pcawg_attributions <- read.csv(sbs_signatures_pcawg_fp)
pcawg_attributions$context <- paste0(substr(pcawg_attributions$Mutation.Subtype,1,1), '[', pcawg_attributions$Mutation.Type, ']', substr(pcawg_attributions$Mutation.Subtype,3,3))
pcawg_attributions$max_sig<-colnames(pcawg_attributions[,5:(ncol(pcawg_attributions)-1)])[apply(pcawg_attributions[,5:(ncol(pcawg_attributions)-1)],1,which.max)]
pcawg_attributions$max_prob <- apply(pcawg_attributions[,5:(ncol(pcawg_attributions)-2)],1,max)
pcawg_exposures <- read.csv(sbs_signatures_pcawg_exposures_fp)
rs_exposures <- as.data.frame(read_excel(sv_sigs_fp, sheet='S7'))
rownames(rs_exposures) <- rs_exposures[,1]
colnames(rs_exposures)[1] <- 'sample'
rs_exposures$sample <- NULL
# merge the sample metadata table with SBS and RS exposures
test.samples.df.sinatures <- merge(test.samples.df, pcawg_exposures, by.x = 'icgc_specimen_id', by.y='Sample.Names')
test.samples.df.sinatures.abobec <- merge(test.samples.df.sinatures, rs_exposures, by.x='aliquot_id', by.y=0)

# sample groups of interest
if (!is.null(sample_subset)) {
    if (sample_subset=='CCNE1') {
        ccne.samples.df <- read.csv(ccne1_sample_list_fp)
        test.samples.df.sinatures.abobec <- subset(test.samples.df.sinatures.abobec, aliquot_id %in% ccne.samples.df$aliquot_id)
    } else if (sample_subset=='CDK12') {
        curated_cdk12 <- c('0009b464-b376-4fbc-8a56-da538269a02f',
                  '84ca6ab0-9edc-4636-9d27-55cdba334d7d',
                  'b243adb4-b3e7-4e0e-bc0d-625aa8dbb1be',
                   '89dad92e-5b3f-479a-a6da-a94ee7df7f8a',
                   'bc0dee07-de20-44d6-be65-05af7e63ac96',
                   '36d1a85e-a09b-4537-86e0-eaf1eb03aed8',
                   '0bfd1043-816e-e3e4-e050-11ac0c4860c5' # this one has two hits
                  )
        test.samples.df.sinatures.abobec <- subset(test.samples.df.sinatures.abobec, aliquot_id %in% curated_cdk12)
    } else if (sample_subset=='HRD'){
        hrd_brca_loh <- read.table(hrd_sample_list_fp, header=TRUE)
        test.samples.df.sinatures.abobec <- subset(test.samples.df.sinatures.abobec, aliquot_id %in% hrd_brca_loh$aliquot_id)
    }    
}
print(paste(nrow(test.samples.df.sinatures.abobec), 'samples after  filtering'))

# mappability of the reference genome
mapability <- read.table(mappability_fp)
mapability.gr <-GRanges(seqnames=Rle(mapability$V1),
                                  ranges=IRanges(mapability$V2, mapability$V3),seqinfo= seqinfo(BSgenome.Hsapiens.UCSC.hg19))
mapability.gr <- trim(reduce(mapability.gr))



apobec.mut.list <- list()
mut.list<- list()
sample.summary.list <- list()

# counter of samples to be loaded
lf <- 0

# looping over samples
# samples that were filtered according to the criteria
for (si in 1:nrow(test.samples.df.sinatures.abobec)) {
    sample_id <- test.samples.df.sinatures.abobec$aliquot_id[si]
    dcc_project_code <- test.samples.df.sinatures.abobec$dcc_project_code[si]
    # these attributions are on 3-nucleotide level
    pcawg_attributions_si <- subset(pcawg_attributions, Sample==test.samples.df.sinatures.abobec$icgc_specimen_id[si])  
    rownames(pcawg_attributions_si) <- pcawg_attributions_si$context
    
    muts.gr <- NULL
    # loading point mutations processed with MutationTimeR, which additionally have information on mutation multiplicity
    if (loadTimeR) {
        fn <- paste(muts_timer_path, "/data/", sample_id, ".RData", sep="")
        # prepare muts.gr
        if (file.exists(fn)) {
            lf <- lf + 1
            load(fn) # loads the vcf variable
            muts.gr <- rowRanges(vcf)
            mut.info.df <- as.data.frame(info(vcf))
            # context is missing
            muts.gr$context <- mut.info.df$context
            muts.gr$MutCN <- mut.info.df$MutCN
            muts.gr$pMutCN <- mut.info.df$pMutCN
            muts.gr$MajCN <- mut.info.df$MajCN
            muts.gr$MinCN <- mut.info.df$MinCN
            muts.gr$max_sig <- pcawg_attributions_si[muts.gr$context,'max_sig']
            
            seqlevels(muts.gr) <- paste0('chr', seqlevels(muts.gr))
            muts_df <- as.data.frame(mut.info.df)
            muts_df$chroms <- as.character(seqnames(muts.gr))
            muts_df$starts <- start(muts.gr)
            muts_df$wt <- as.character(muts.gr$REF)
            muts_df$Validation_status <- NULL
        } else {
            #print('file not found')
        }
    } else {
        # starting with a default VCF file
        sample_muts <- vcfToSNVcatalogue(paste0(pcawg_muts_folder,sample_id,'.consensus.20160830.somatic.snv_mnv.vcf.gz'))
        muts_df <- sample_muts$muts
        muts_df$max_sig <- pcawg_attributions_si[muts_df$context,'max_sig']
        muts_df$max_prob <- pcawg_attributions_si[muts_df$context,'max_prob']        
        
        muts.gr <- GRanges(seqnames=Rle(paste0('chr', muts_df$chroms)),
                  ranges=IRanges(muts_df$starts, muts_df$ends))
        muts.gr$max_sig <- muts_df$max_sig
    }

    # if file was loaded
    if (!is.null(muts.gr)) {

            if (sample_id %in% names(sample.rearrs)) {
                if (is.null(max_sig)) {
                    filter.string <- paste0('(is.clustered==',is.clustered,') & (svclass=="',svclass,'")' )
                    eval(parse(text=paste("test.bedpe <- subset(do.call('rbind',sample.rearrs[sample_id]),", filter.string, ")" )))
                } else {
                    test.bedpe <- subset(do.call('rbind',sample.rearrs[sample_id]),(max.Ref.Sig==max_sig) & (chrom1==chrom2 ))
                }
                
                if (!is.null(size_min) & size_min>0) {
                    test.bedpe <- subset(test.bedpe, length>size_min)
                }
                if (!is.null(size_max) & size_max>0) {
                    test.bedpe <- subset(test.bedpe, length<size_max)
                }
                
                
                if (nrow(test.bedpe)>0) {
                    if (simulate) {
                        # this was an attempt to simulate SVs by shifting them to the right, preserving length (not used in the manuscript
                        test.bedpe$start1 <- test.bedpe$start2
                        # fake duplicaitons!
                        test.bedpe$start2 <- test.bedpe$start1+ test.bedpe$length         
                    }
                   
                    r <- performSVOverlaps.R(muts.gr, test.bedpe, simulate, margin.size)
                    sample.summary.list[[si]] <-  r[['sample.summary']]
                    apobec.mut.list[[si]] <- r[['sample_muts_apobec']]
                    mut.list[[si]] <- r[['muts_df']]
                                      
                }
            }
    }
}

print(paste(lf, 'files loaded'))
# save sample-level statistics
save(sample.summary.list, apobec.mut.list, mut.list, file=paste0(output_folder, exp.name, '.RData'))



