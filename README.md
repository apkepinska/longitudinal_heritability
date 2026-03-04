# Longitudinal Heritability Pipeline

Longitudinal heritability pipeline for paediatric psychiatric phenotypes using twin-based ACE decomposition, SNP-based heritability estimation (GCTA-GREML), and linkage disequilibrium score regression (LDSC) across repeated assessments. Analyses are applied to [ABCD Study](https://abcdstudy.org/) data.

---

## Overview

This pipeline estimates the heritability of paediatric mental health phenotypes across development using three complementary approaches:

1. **Twin-based ACE decomposition** — decomposes phenotypic variance into additive genetic (A), shared environmental (C), and unique environmental (E) components using MZ and DZ twin pairs, run cross-sectionally at each assessment wave. Two implementations are used: conventional SEM via `mets::twinlm()`, and a Bayesian hierarchical linear model (HLM) via `brms` following Chen et al. (2025)
2. **GCTA-GREML** — estimates SNP-based heritability (h²_SNP) in unrelated individuals using genome-wide genetic relatedness matrices, applied to latent growth curve (LGC) intercepts representing each individual's stable symptom level across development
3. **LDSC** — estimates heritability from GWAS summary statistics (derived from LGC intercept phenotypes) using linkage disequilibrium score regression, and computes genetic correlations across phenotypes and with published adult GWAS

Analyses focus on **continuous mental health phenotypes** (e.g. CBCL subscales). Twin ACE models are run cross-sectionally at each assessment wave to track how heritability changes across development. GCTA-GREML and LDSC are applied to LGC intercepts extracted from the full longitudinal data in unrelated individuals, capturing the heritability of the stable component of each phenotype — following the approach of [Deng et al. (2024), *Molecular Psychiatry*](https://www.nature.com/articles/s41380-024-02704-4). The Bayesian HLM twin ACE implementation follows [Chen et al. (2025), *Frontiers in Genetics*](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2025.1522729/full).

---

## Directory Structure

```
longitudinal_heritability/
├── analysis_scripts/
│   ├── abcd_gwas_for_ldsc/
│   │   └── mental_health/          # GWAS pipeline (PLINK2-based) for LDSC input
│   ├── gcta_greml/
│   │   └── mental_health/          # GCTA-GREML heritability estimation
│   ├── latent_growth_curve/
│   │   └── mental_health/          # LGC models; extracts intercept/slope phenotypes
│   ├── ldsc/
│   │   └── abcd_gwas_ldsc/
│   │       └── mental_health/      # LDSC h² and genetic correlation analyses
│   └── twin_heritability/
│       └── mental_health/          # Twin ACE models
└── processing_scripts/
    ├── ldsc/                        # LDSC setup and reference panel downloads
    └── twin_heritability/           # Twin identification and phenotype QC
```

---

## Analysis Workflow

### Step 1: Data processing
- Identify twin pairs and cross-check zygosity (`processing_scripts/twin_heritability/`)
- Check phenotype availability and apply sample size filters
- Download and configure LDSC and LD reference panels (`processing_scripts/ldsc/`)

### Step 2: Latent growth curve modelling
- Fit LGC models to longitudinal mental health phenotypes
- Extract intercept and slope factor scores as phenotypes for downstream analyses

### Step 3: Twin heritability (ACE)
- Run ACE models for each phenotype using the twin subsample (~1,100 twins from ~550 families)
- Array job submission for batch processing across phenotypes

### Step 4: GCTA-GREML
- Build genetic relatedness matrix (GRM) from SNP data in unrelated individuals
- Run GREML for each phenotype; batch submission via array jobs

### Step 5: GWAS for LDSC
- Run GWAS on LGC intercept phenotypes using PLINK2
- Reformat summary statistics for LDSC input

### Step 6: LDSC
- Estimate h² from GWAS summary statistics
- Compute genetic correlations across phenotypes and with published adult GWAS

---

## Data

Analyses use data from the **Adolescent Brain Cognitive Development (ABCD) Study**, a longitudinal study of ~11,000 children recruited at ages 9–10 across 21 sites in the United States. Data access is managed through the [NIMH Data Archive](https://nda.nih.gov/).

> **Note:** No individual-level data are included in this repository. Scripts assume data are stored locally and paths are configured by the user (see Configuration below).

---

## Dependencies

**R packages:**
- `lavaan` — latent growth curve modelling
- `mets` — conventional twin ACE models (`twinlm()`)
- `brms` — Bayesian HLM twin ACE models (Chen et al. 2025 approach)
- `data.table`, `dplyr`, `tidyr` — data wrangling

**Command-line tools:**
- [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) (v1.94+)
- [PLINK2](https://www.cog-genomics.org/plink/2.0/)
- [LDSC](https://github.com/bulik/ldsc) (Python 2.7)

---

## Configuration

Each script contains a configuration section at the top where file paths and parameters should be set before running. No hardcoded paths are used outside of these sections.

---

## Citation

If you use this pipeline, please cite:

> Kępińska AP et al. (in preparation). Longitudinal heritability of paediatric psychiatric phenotypes in the ABCD Study.

Please also cite the key methodological papers this pipeline is based on:

> Deng WQ et al. (2024). Longitudinal characterization of impulsivity phenotypes boosts signal for genomic correlates and heritability. *Molecular Psychiatry*. https://doi.org/10.1038/s41380-024-02704-4

> Chen G et al. (2025). Heritability estimation through hierarchical linear modeling. *Frontiers in Genetics*. https://doi.org/10.3389/fgene.2025.1522729

Please also cite the ABCD Study, GCTA, PLINK2, and LDSC as appropriate.

---

## Contact

Adrianna Kępińska — [@apkepinska](https://github.com/apkepinska)
