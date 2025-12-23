#!/bin/bash
#SBATCH -c 8                               # Request 16 cores
#SBATCH -t 0-2:05                         # Runtime in D-HH:MM format
#SBATCH -A park
#SBATCH -p short                              # Partition to run in
#SBATCH --mem=48G                         # Memory total in MB (for all cores)
#SBATCH -o log/sv_overlaps_%A_%a.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e log/sv_overlaps_%A_%a.err                 # File to which STDERR will be written, including job ID (%j)
#SBATCH --array=1-55                                     # 1-37, 1-15 for channels, 1-21 for signatures


module load gcc/14.2.0 python/3.13.1 R/4.4.2
module load samtools
export R_LIBS_USER="~/park_dglodzik/Renvs/R-4.4.2-IRkernel/library/"

config_fn=/home/dg204/projects/rsignatures/src/sv_overlap_config.csv
LNO=${SLURM_ARRAY_TASK_ID}
NTH_LINE=$(cat ${config_fn} | head -${LNO} | tail -1)
echo $NTH_LINE
# the line should be: tumor ID, tumor bam, normal ID, normal bam, (37|38)
EXP_NAME=$(echo $NTH_LINE | awk -F"," '{ print $1}') 
echo $EXP_NAME
FILTER_STR=$(echo $NTH_LINE | awk -F"," '{ print $2}') 
echo $FILTER_STR
SHUFFLE_STR=$(echo $NTH_LINE | awk -F"," '{ print $3}') 
echo $SHUFFLE_STR

Rscript run_sv_topography.R \
    "$EXP_NAME" \
    "$FILTER_STR" \
    "$SHUFFLE_STR"
