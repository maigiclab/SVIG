# SVIG
**SVIG**: Analysis of DNA replication at structural variant (SV) loci and a multi-class classifier of replication stress and DNA repair deficiencies.

This repository contains analysis pipelines, notebooks, and utilities used to study:
- Replication timing, origin activity, and fork directionality at SV loci
- APOBEC point mutation strand bias
- SV signature generation, extraction, and visualization
- A multi-class SV-based classifier of replication stress and DNA repair defects

---

## Repository structure

### 1. SV topography and replication features
Analysis of replication timing, origins, and fork directionality with respect to SV spans. PCAWG dataset; SV-signature-focused or genetic-alteration focused experiments.

**Scripts**
- `src/script_run_sv_topography.sh`
- `src/run_sv_topography.R`

**Notebooks**
- `notebooks/2bp_summary.ipynb`  
  Generates a summary across experiments.

---

### 2. APOBEC strand asymmetry analysis
Analysis of APOBEC mutagenesis strand bias at SV loci. PCAWG dataset; SV-signature-focused or genetic-alteration focused experiments.

**Scripts**
- `src/script_run_APOBEC_overlaps.sh`
- `src/run_APOBEC_overlaps.R`  
**Key function**
- `src/utils/performSVOverlaps.R`  

**Notebooks**
- `notebooks/APOBEC_summary.ipynb`  
  Generates a summary across experiments.

---

### 3. SV signatures in CCNE1 and CDK12 tumors
Plots and comparative analyses of SV signatures in CCNE1- and CDK12-amplified tumors.

**Notebooks**
- `notebooks/ccne1_cdk12_comparison.ipynb`

---

### 4. SV catalogue construction
Generation of SV catalogues for downstream signature analysis.

**Key function**
- `src/utils/performOverlaps.R`

**Notebooks**
- `notebooks/cataloguePCAWG.ipynb`
- `notebooks/combinePCAWGcatalogues.ipynb`  
  Combines PCAWG SV catalogues.

---

### 5. SV signature extraction
Extraction of SV signatures from curated catalogues.

**Scripts**
- `src/metabreast_sv_extraction.sh`
- `src/mb_extraction.py`

---

### 6. SV signature visualization and correction
Visualization and correction of SV signatures with respect to replication features.

**Notebooks**
- `notebooks/Fig5_RFD_SigPlots.ipynb`
- `notebooks/RFD_SV_assignment.ipynb`
- `notebooks/FigSupp_RTime_SigPlots.ipynb`

---

### 7. SVIG classifier training
Training of the SVIG multi-class classifier.

**Notebooks**
- `notebooks/classifier_svig.ipynb`

---

### 8. Additional statistics
Other analyses and statistics used in the manuscript.

**Notebooks**
- `notebooks/paperStats.ipynb`

---

## Notes
- Most analyses are notebook-driven and assume pre-computed inputs.
- Paths and job scripts may require adaptation for specific compute environments.
- This repository accompanies a research manuscript; code is provided for transparency and reproducibility.
