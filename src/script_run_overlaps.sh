#!/bin/bash
#SBATCH -c 8                               # Request 16 cores
#SBATCH -t 0-2:05                         # Runtime in D-HH:MM format
#SBATCH -A park
#SBATCH -p short                              # Partition to run in
#SBATCH --mem=8G                         # Memory total in MB (for all cores)
#SBATCH -o log/overlaps_%A_%a.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e log/overlaps_%A_%a.err                 # File to which STDERR will be written, including job ID (%j)
#SBATCH --array=32-33                                     # 1-15 for channels, 1-21 for signatures


#module load gcc/9.2.0 python/3.10.11 R/4.3.1
#module load samtools
#export R_LIBS_USER="~/R-4.3.1-IRkernel/library"

module load gcc/14.2.0 python/3.13.1 R/4.4.2
module load samtools
export R_LIBS_USER="~/park_dglodzik/Renvs/R-4.4.2-IRkernel/library/"

config_fn=/home/dg204/projects/rsignatures/src/overlapConfigsSigs2.csv
#config_fn=/home/dg204/projects/rsignatures/src/overlapConfigs.csv
LNO=${SLURM_ARRAY_TASK_ID}
NTH_LINE=$(cat ${config_fn} | head -${LNO} | tail -1)
echo $NTH_LINE
# the line should be: tumor ID, tumor bam, normal ID, normal bam, (37|38)
SVCLASS=$(echo $NTH_LINE | awk -F"," '{ print $1}') 
echo "$SVCLASS"
SIZE_MIN=$(echo $NTH_LINE | awk -F"," '{ print $2}')  
echo "$SIZE_MIN"
SIZE_MAX=$(echo $NTH_LINE | awk -F"," '{ print $3}')  
echo "$SIZE_MAX"
MAX_SIG=$(echo $NTH_LINE | awk -F"," '{ print $4}')  
echo "$MAX_SIG"
SAMPLE_SUBSET=$(echo $NTH_LINE | awk -F"," '{ print $6}')  
echo "$SAMPLE_SUBSET"


#OUTPUT_FOLDER=/home/dg204/park_dglodzik/APOBEC_overlaps/
OUTPUT_FOLDER=/home/dg204/park_dglodzik/APOBEC_overlaps/sigs_all/

Rscript run_overlaps.R \
--svclass "$SVCLASS" \
--size_min "$SIZE_MIN" \
--size_max "$SIZE_MAX" \
--output "$OUTPUT_FOLDER" \
--max_sig "$MAX_SIG" \
--subset "$SAMPLE_SUBSET"

# look for notebooks/APOBEC/muts_dups_signatures_summary.ipynb