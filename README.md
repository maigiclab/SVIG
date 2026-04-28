# SVIG
**SVIG**: Analysis of DNA replication at structural variant (SV) loci and a multi-class classifier of replication stress and DNA repair deficiencies.

This repository contains analysis pipelines, notebooks, and utilities used to study:
- Replication timing, origin activity, and fork directionality at SV loci
- APOBEC point mutation strand bias
- SV signature generation, extraction, and visualization
- A multi-class SV-based classifier of replication stress and DNA repair defects

---

## Example runs

### 0. Download the data and set up the environment

Intermediate data files required to reproduce all analyses are deposited on Zenodo:

**Glodzik, D. (2026). Intermediate files for repository: maigiclab/SVIG (v2). Zenodo.**  
https://doi.org/10.5281/zenodo.19401853


Download and extract the archive into the repository root:
```bash
wget -O svig_data.tar.gz "https://zenodo.org/records/19601623/files/svig_data.tar.gz?download=1"
tar -xzvf svig_data.tar.gz
```
This will create a `data/` directory inside the repository.

A `Dockerfile` is provided with all R and Python dependencies listed.
a pre-built Docker image is available on Docker Hub:                                                                    
                                         
```bash
  docker pull dglodzik/svig
  docker run -v $(pwd):/app -m 2g -it dglodzik/svig /bin/bash        
```

Alternatively, **Python** dependencies are listed in `pyproject.toml`. To install:
```bash
pip install .
```
**R** package dependencies are listed in the Dockerfile. 



### 1. SV topography and replication features
```bash
cd src
R
```
```r
source('run_sv_topography.R')
```

**Expected output:** `../data/processed/RS1\ CDK12/` — including `.pdf` plots, `.csv` files and `.RData` objects that characterize overlap of SVs with replication features.

**SVs are loaded from:** `../data/interim/sample.rearrs.RData`

Other reference files in `../data/interim` are required.

### 2. APOBEC strand asymmetry analysis
```bash
cd src
```
```r
R
source('run_APOBEC_overlaps.R')
```

**Expected output:** `../data/processed/PCAWG_duplication_-1_300000_Ref.Sig.R1_CDK12.RData`

**Files required:**
- SVs loaded from: `../data/interim/sample.rearrs.RData`
- Point mutations, TimeR-processed: `../data/interim/TimeR_bb_all_sigs/`
- Other reference files in `../data/interim` are required.

## Repository structure

### A. SV topography and replication features
Analysis of replication timing, origins, and fork directionality with respect to SV spans. PCAWG dataset; SV-signature-focused or genetic-alteration focused experiments.

**Scripts**
- `src/script_run_sv_topography.sh`
- `src/run_sv_topography.R`

**Notebooks**
- `notebooks/2bp_summary.ipynb`  
  Generates a summary across experiments.

---

### B. APOBEC strand asymmetry analysis
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

### C. SV signatures in CCNE1 and CDK12 tumors
Plots and comparative analyses of SV signatures in CCNE1- and CDK12-amplified tumors.

**Notebooks**
- `notebooks/ccne1_cdk12_comparison.ipynb`

---

### D. SV catalogue construction
Generation of SV catalogues for downstream signature analysis.

**Key function**
- `src/utils/performOverlaps.R`

**Notebooks**
- `notebooks/cataloguePCAWG.ipynb`
- `notebooks/combinePCAWGcatalogues.ipynb`  
  Combines PCAWG SV catalogues.

---

### E. SV signature extraction
Extraction of SV signatures from curated catalogues.

**Scripts**
- `src/metabreast_sv_extraction.sh`
- `src/mb_extraction.py`

---

### F. SV signature visualization and correction
Visualization and correction of SV signatures with respect to replication features.

**Notebooks**
- `notebooks/Fig5_RFD_SigPlots.ipynb`
- `notebooks/RFD_SV_assignment.ipynb`
- `notebooks/FigSupp_RTime_SigPlots.ipynb`

---

### G. SVIG classifier training
Training of the SVIG multi-class classifier.

**Notebooks**
- `notebooks/classifier_svig.ipynb`

---

### H. Additional statistics
Other analyses and statistics used in the manuscript.

**Notebooks**
- `notebooks/paperStats.ipynb`

---

## Notes
- Most analyses are notebook-driven and assume pre-computed inputs.
- Paths and job scripts may require adaptation for specific compute environments.
- This repository accompanies a research manuscript; code is provided for transparency and reproducibility.
