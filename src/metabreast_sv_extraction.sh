#!/bin/bash

#SBATCH -c 18              
#SBATCH -t 20-23:50
#SBATCH -A park_contrib           
#SBATCH -p park            
#SBATCH --mem=32G            
#SBATCH -o logs/SV_SigExtraction_mvnmf_RFD_%A_%a.out 
#SBATCH -e logs/SV_SigExtraction_mvnmf__%A_%a.err 
#SBATCH --array=12-13     

module load gcc/14.2.0 python/3.13.1 R/4.4.2
source ~/park_dglodzik/envs/jupytervenv3.13/bin/activate

# Optional: Print some debug info
echo "Running on node: $(hostname)"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

# should be 100000
NO_ITERS=100000
SCRIPT_PATH="/home/dg204/projects/rsignatures/src/extraction/mb_extraction.py"

# Run different jobs depending on the array task ID
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_normalized" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 20 \
        --normalize-X

elif [ "${SLURM_ARRAY_TASK_ID}" -eq 2 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_rfd" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 20

elif [ "${SLURM_ARRAY_TASK_ID}" -eq 3 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "channels32_normalized" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 15 \
        --normalize-X

elif [ "${SLURM_ARRAY_TASK_ID}" -eq 4 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "channels32" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 15 
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 5 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_corr_normalized" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 15 
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 6 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_normalized" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 25 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/PCAWG/RFD/Selected_RFD.csv" \
        --cohort-name pancan_sel
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 7 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_rfd" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 25 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/PCAWG/RFD/Selected_RFD.csv" \
        --cohort-name pancan_sel
        elif [ "${SLURM_ARRAY_TASK_ID}" -eq 7 ]; then
python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_rfd" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 25 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/PCAWG/RFD/Selected_RFD.csv" \
        --cohort-name pancan_sel
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 8 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_corr_normalized" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 25 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/PCAWG/RFD/Selected_RFD.csv" \
        --cohort-name pancan_sel
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 9 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_replitime" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 25 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/PCAWG/repliTiming_mp/Selected_repliTiming_mp.csv" \
        --cohort-name pancan_sel
        # see the combinePCAWGcatalogues notebook
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 10 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_replitime_to_transpose" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 25 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/GEL/sel_tissue_cat_repliTiming_mp.csv" \
        --cohort-name gel_sel
        # see the combinePCAWGcatalogues notebook
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 11 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_rfd_to_transpose" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 25 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/GEL/sel_tissue_cat_RFD.csv" \
        --cohort-name gel_sel
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 12 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_replitime_to_transpose" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 30 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/GEL/all_tissue_cat_repliTiming_mp.csv" \
        --cohort-name gel_all
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 13 ]; then
    python "$SCRIPT_PATH" \
        --matrix-type "non_clustered_rfd_to_transpose" \
        --no-iters "$NO_ITERS" \
        --min-n-components 5 \
        --max-n-components 30 \
        --catalogue-matrix "/home/dg204/projects/rsignatures/data/processed/SVmatrices/GEL/all_tissue_cat_RFD.csv" \
        --cohort-name gel_all
else
    echo "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi