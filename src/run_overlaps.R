if (interactive()) {
    is.clustered <- FALSE
    svclass <- "duplication"
    size_min <- -1
    size_max <- 10000000
    outputFolder <- '/home/dg204/park_dglodzik/APOBEC_overlaps/sigs_all/'
    max_sig <- NULL
    sample_subset <- 'CCNE1'
    output_folder <- '/home/dg204/park_dglodzik/APOBEC_overlaps/sigs_all/'
} else {

    library(optparse)
    option_list <- list(
      make_option(c( "--svclass"), type = "character", help = "SV class"),
      make_option(c("--size_min"), type = "integer", help = "Minimum size"),
      make_option(c("--size_max"), type = "integer", help = "Maximum size"),
      make_option(c( "--output"), type = "character", help = "Output folder"),
      make_option(c( "--max_sig"), type = "numeric", help = "RS signature (eg. Ref.Sig.R15)"),
      make_option(c("--subset"), type = "character", help = "Sample subset (CCNE1, CDK12, HRD)")
    )
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)
    print(opt)
    svclass <- opt$svclass
    size_min <- opt$size_min
    size_max <- opt$size_max
    output_folder <- opt$output
    max_sig <- if (opt$max_sig == "null") NULL else opt$max_sig
    sample_subset <- opt$subset
    
    is.clustered <- FALSE

}

sample_metadata = read.csv('~/park_dglodzik/data_repo//PanCan/WGS.metadata.txt', sep='\t')
sample_metadata <- subset(sample_metadata, !duplicated(icgc_specimen_id))

exp.name <- paste0('PCAWG_', svclass,'_',size_min, '_', as.integer(size_max), '_', max_sig) 
if (!is.null(sample_subset)) {
    exp.name <- paste0(exp.name, '_', sample_subset)
}
test.samples.df <- sample_metadata
# this is used if maxSig is null


# key settings
simulate <- FALSE
margin.size <- 1e5
loadTimeR <- TRUE
# TimeR result path
# this is missing the context field
# see projects/brca_timing/src/script_hrd_timer_pcawg.sh
result_path <- "/home/dg204/park_dglodzik/TimeR_bb_all_sigs/"

# to be used with CCNE1
# to be used with CDK12
# this is INCOMPLETE!
#result_path <- "~/projects/rsignatures/data/processed/TimeR_bb"


suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg19))
suppressWarnings(library("GenomicRanges"))
suppressWarnings(library(regioneR))
library(stringr)
library(lmtest)
library(MASS)
source('~/repos/hotspots/R/utils/plotCoefficients.R')
source('~/repos/hotspots/R/utils/prepareBinData.R')
library(signature.tools.lib)
library(plyr)
library(readxl)
require(MutationTimeR) #
library(reshape2)
options(repr.matrix.max.cols=200, repr.matrix.max.rows=100)
source('~/projects/rsignatures/src/utils/performSVOverlaps.R')

#notebooks/signature calc/rearr_catalogue
load('~/projects/rsignatures//data/processed/sample.rearrs.RData')

# SBS signatures and exposures
pcawg_attributions <- read.csv('~/park_dglodzik/data_repo/PanCan/repertoire of signatures/attributions/SigProfilier_PCAWG_WGS_probabilities_SBS.csv')
pcawg_attributions$context <- paste0(substr(pcawg_attributions$Mutation.Subtype,1,1), '[', pcawg_attributions$Mutation.Type, ']', substr(pcawg_attributions$Mutation.Subtype,3,3))
pcawg_attributions$max_sig<-colnames(pcawg_attributions[,5:(ncol(pcawg_attributions)-1)])[apply(pcawg_attributions[,5:(ncol(pcawg_attributions)-1)],1,which.max)]
pcawg_attributions$max_prob <- apply(pcawg_attributions[,5:(ncol(pcawg_attributions)-2)],1,max)
pcawg_exposures <- read.csv('~/park_dglodzik/data_repo/PanCan/repertoire of signatures/signatures_in_samples/PCAWG_sigProfiler_SBS_signatures_in_samples.csv')
rs_exposures <- as.data.frame(read_excel('~/park_dglodzik/data_repo//PanCan//PCAWG Andrea//43018_2020_27_MOESM3_ESM.xlsx', sheet='S7'))
rownames(rs_exposures) <- rs_exposures[,1]
colnames(rs_exposures)[1] <- 'sample'
rs_exposures$sample <- NULL
# merge the sample metadata table with SBS and RS exposures
test.samples.df.sinatures <- merge(test.samples.df, pcawg_exposures, by.x = 'icgc_specimen_id', by.y='Sample.Names')
test.samples.df.sinatures.abobec <- merge(test.samples.df.sinatures, rs_exposures, by.x='aliquot_id', by.y=0)
if (!is.null(sample_subset)) {
    if (sample_subset=='CCNE1') {
        ccne.samples.df <- read.csv('~/projects/rsignatures/data/processed/CCNE1.amp.sample.csv')
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
        hrd_brca_loh <- read.table('~/projects/rsignatures/data/processed/BRCA_and_LOH_PCAWG.tsv', header=TRUE)
        test.samples.df.sinatures.abobec <- subset(test.samples.df.sinatures.abobec, aliquot_id %in% hrd_brca_loh$aliquot_id)
    }    
}
print(paste(nrow(test.samples.df.sinatures.abobec), 'samples after  filtering'))


# sample groups of interest



mapability <- read.table('~/repos/hotspots//data/breastData/hg19CRG.100bp/hg19.CRC.100mer.bed')
mapability.gr <- trim(reduce(GRanges(seqnames=Rle(mapability$V1),
                                  ranges=IRanges(mapability$V2, mapability$V3),seqinfo= seqinfo(BSgenome.Hsapiens.UCSC.hg19))))


options(error = recover)

apobec.mut.list <- list()
mut.list<- list()
sample.summary.list <- list()

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
    if (loadTimeR) {

        fn <- paste(result_path, "/data/", sample_id, ".RData", sep="")
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
        sample_muts <- vcfToSNVcatalogue(paste0('/home/dg204/park_data/ICGC/SNV_indel_calls/final_consensus_12oct_passonly/snv_mnv/',sample_id,'.consensus.20160830.somatic.snv_mnv.vcf.gz'))
        muts_df <- sample_muts$muts
        muts_df$max_sig <- pcawg_attributions_si[muts_df$context,'max_sig']
        muts_df$max_prob <- pcawg_attributions_si[muts_df$context,'max_prob']        
        
        muts.gr <- GRanges(seqnames=Rle(paste0('chr', muts_df$chroms)),
                  ranges=IRanges(muts_df$starts, muts_df$ends))
        muts.gr$max_sig <- muts_df$max_sig
    }

    
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
                        test.bedpe$start1 <- test.bedpe$start2
                        # fake duplicaitons!
                        test.bedpe$start2 <- test.bedpe$start1+ test.bedpe$length         
                    }
                   
                    r <- performSVOverlaps.R(muts.gr, test.bedpe)
                    sample.summary.list[[si]] <-  r[['sample.summary']]
                    apobec.mut.list[[si]] <- r[['sample_muts_apobec']]
                    mut.list[[si]] <- r[['muts_df']]
                                      
                }
            }
    }
}

print(paste(lf, 'files loaded'))

save(sample.summary.list, apobec.mut.list, mut.list, file=paste0(output_folder, exp.name, '.RData'))



