# SVIG: analysis of DNA replication at SV loci, and multi-class classifier of replication stress and DNA repair deficiencies 


Analysis of replication timing, origins and direction wrt. SV span
    src/script_run_sv_topography.sh
    src/run_sv_topography.R
    notebooks/2bp_summary.ipynb - generates a summary across experiments

Analysis of APOBEC strand
    src/script_run_overlaps.sh
    src/run_overlaps.R - PCAWG dataset, SV signature-focused experiments
    notebooks/APOBEC_summary.ipynb - generates a summary across experiments
[tested interactively; need to test to bash jobs + summary]


Plots and results regarding SV signatures in CCNE1 and CDK12 tumors
    notebooks/ccne1_cdk12_comparison.ipynb

Creation of SV catalogues for signature analysis
    src/utils/performOverlaps.R
    notebooks/cataloguePCAWG.ipynb
    notebooks/combinePCAWGcatalogues.ipynb (for PCAWG, the catalogues are combined here)

SV signature extraction
    src/metabreast_sv_extraction.sh
    src/mb_extraction.py

SV signature visualization and correction
    notebooks/Fig5_RFD_SigPlots.ipynb
    notebooks/RFD_SV_assignment.ipynb
    notebooks/FigSupp_RTime_SigPlots.ipynb

SVIG classifier training
    notebooks/classifier1.ipynb

Other stats used in the paper
    notebooks/paperStats.ipynb


