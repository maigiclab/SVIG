#!/bin/bash
#===============================================================================
# Purpose:
#   Run APOBEC overlap analysis for one parameter set defined in overlapConfigsSigs2.csv.
#   This script is designed to run as a Slurm job array, where each task reads
#   a specific line from the CSV config file and runs run_overlaps.R with those parameters.
#
# Inputs:
#   - overlapConfigsSigs2.csv : parameter table (one run per line)
#       Expected columns (comma-separated):
#         1) SVCLASS
#         2) SIZE_MIN
#         3) SIZE_MAX
#         4) MAX_SIG
#         5) (unused here)
#         6) SAMPLE_SUBSET
#
# Work script (R)
#   - run_overlaps.R : performs the overlap computation and writes results
#
# Outputs:
#   - ../data/processed/APOBEC/ : output directory for overlap results
#   - log/overlaps_<jobid>_<taskid>.out : STDOUT log
#   - log/overlaps_<jobid>_<taskid>.err : STDERR log
#
# Usage:
#   sbatch <this_script>.sh
#
# Notes:
#   - SLURM_ARRAY_TASK_ID is used as the CSV line number to read.
#   - After the array finishes, summarize results using:
#       notebooks/APOBEC/muts_dups_signatures_summary.ipynb
#===============================================================================

#SBATCH -c 8                               # Request 16 cores
#SBATCH -t 0-2:05                         # Runtime in D-HH:MM format
#SBATCH -A park
#SBATCH -p short                              # Partition to run in
#SBATCH --mem=8G                         # Memory total in MB (for all cores)
#SBATCH -o log/overlaps_%A_%a.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e log/overlaps_%A_%a.err                 # File to which STDERR will be written, including job ID (%j)
#SBATCH --array=32-33                                     # 1-15 for channels, 1-21 for signatures


module load gcc/14.2.0 python/3.13.1 R/4.4.2
module load samtools
export R_LIBS_USER="~/park_dglodzik/Renvs/R-4.4.2-IRkernel/library/"

config_fn=overlapConfigsSigs2.csv

OUTPUT_FOLDER="../data/processed/APOBEC/"
mkdir -p "$OUTPUT_FOLDER"

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

Rscript run_overlaps.R \
--svclass "$SVCLASS" \
--size_min "$SIZE_MIN" \
--size_max "$SIZE_MAX" \
--output "$OUTPUT_FOLDER" \
--max_sig "$MAX_SIG" \
--subset "$SAMPLE_SUBSET"

# summarize the results across the runs using this notebook: notebooks/APOBEC/muts_dups_signatures_summary.ipynb