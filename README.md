# Lumos Networks

**A transparent three-stage workflow for transcriptomic analysis, functional context, and prior-knowledge networks**

Lumos Networks is a Streamlit research application developed by EvoMol-Lab at
BioME/UFRN. The current suite contains three connected modules:

1. **DDEA — Diagonal Differential Expression Alley**
2. **APP — Arithmancy Pathway Profiler**
3. **PG — PrioriGraph**



## What each module computes

| Module | Implemented methods | Principal outputs |
|---|---|---|
| DDEA | PyDESeq2 Wald tests for raw integer RNA-seq counts; limma empirical-Bayes models for continuous expression; full-rank batch, covariate and paired-subject designs; random-effects meta-analysis | Complete DE tables, design/contrast provenance, diagnostic figures, forest plots, cross-dataset alignment, JSON manifest; self-contained interactive HTML report (offline) with searchable tables, CSV/TSV export and Save as PDF |
| APP | Explicit-background ORA through Speedrichr, separated by expression direction by default; STRING v12 networks; JASPAR/TRRUST target-set enrichment | KEGG/GO tables with query/universe sizes, STRING associations, candidate-regulator tables with provenance; self-contained interactive HTML report (offline) with searchable tables, CSV/TSV export and Save as PDF |
| PrioriGraph | Enriched-term → gene and candidate TF → target integration; HITS, PageRank, betweenness; personalized RWR, signed linear propagation and parameter-grid sensitivity | Interactive network, edge provenance, structural rankings, propagation and stability tables; self-contained interactive HTML report (offline) with searchable tables, CSV/TSV export and Save as PDF |

The README deliberately uses **candidate regulator**, **association**, and
**prior-knowledge network**. The implemented workflow does not establish
cohort-specific TF binding, regulatory causality, or pathway activation.

## Analytical workflow

```text
GEO matrix + user-defined groups
              │
              ▼
 DDEA: effect estimates + uncertainty
              │
              ▼
 APP: ORA + STRING + TF target-set enrichment
              │
              ▼
 PrioriGraph: inspectable topology + propagation
```

Streamlit session state connects modules during one browser session. CSV and
JSON downloads provide durable, reviewable hand-offs.

## Scientific behaviour

### DDEA

- Raw integer counts are validated before PyDESeq2. Fractional, negative,
  missing, or infinite values are rejected instead of silently truncated.
- Count features mapping to the same symbol are summed before PyDESeq2;
  microarray probes mapping to one symbol are summarized by per-sample median.
- Count filtering retains genes with at least 10 counts in at least the size of
  the smaller group.
- Continuous matrices retain finite, non-constant features. High-dynamic-range
  non-negative values can be log2-transformed when requested.
- The condition effect can be adjusted for batch, categorical/continuous
  covariates and paired subjects. Lumos rejects incomplete pairs and rank-
  deficient/confounded designs before model fitting.
- The publication route for microarrays and other continuous matrices uses
  Bioconductor limma with robust, trend-aware empirical-Bayes moderation.
  Welch is retained only as an explicitly labelled exploratory fallback.
- A DEG must pass both the selected FDR and absolute effect threshold.
- Single Dataset can retain several named contrasts from the same GSE. Results
  are aligned by gene and classified as consistent, shared, specific, or
  discordant; reused samples are not incorrectly treated as independent studies
  and no pooled meta-analytic effect is calculated.
- Multi-Dataset Meta-analysis is intended for independent datasets from the
  same tissue with comparable contrasts. Multi-Tissue Comparison keeps effects
  from different tissues separate and supports more than two datasets.
- The strict conserved multi-dataset signature requires every dataset to pass
  both thresholds with the same direction. It never substitutes minimum FDR or
  a mean effect for concordance.
- Random-effects meta-analysis reports confidence intervals, Q, I², tau²,
  study count, and contributing accessions.

### APP

- APP currently performs **ORA**, not ranked-list GSEA.
- Every ORA requires an explicit measured-gene universe. DDEA imports use all
  genes in the complete tested result table; users may instead paste a custom
  universe. Query genes outside the universe are excluded and counted.
- Default libraries are `KEGG_2021_Human` and
  `GO_Biological_Process_2023`; their names are retained in output tables.
- Up- and down-regulated lists are queried separately by default. Direction
  labels describe input lists; ORA alone does not infer activation/inhibition.
- JASPAR 2024 and TRRUST 2019 are queried independently and retained as an
  evidence-source column.
- STRING edges retain the returned score and are labelled associations.

### PrioriGraph

- Node rankings are calculated on the displayed TF → target subnetwork.
- TF-target exports preserve database, term, direction, overlap, and adjusted
  p-value when those fields are available.
- Personalized RWR conserves probability at dangling nodes and exposes the
  restart probability.
- Sensitivity analysis perturbs FDR, absolute-effect, RWR-restart and propagation-
  decay parameters and reports configuration coverage, top-k frequency, rank IQR
  and score variability for every node.

## Installation

### Docker

#### Obtain Docker image from Docker Hub


#### Build Docker image locally

### Local Python

Python 3.10 or 3.11 is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
streamlit run Lumos_Home.py
```

## Tests

Pure numerical functions live in `lumos_core.py` so they can be tested without
starting Streamlit.

```bash
python -m unittest discover -s tests -v
python -m py_compile Lumos_Home.py lumos_core.py lumos_report.py pages/*.py
```

The unit suite covers normalization, raw-count validation, experimental-design
rank checks, paired samples, GEO sample alignment, VST orientation, explicit ORA
universes, benchmark metrics, conserved signatures, RWR mass conservation and
network sensitivity. GitHub Actions additionally define:

- synthetic PyDESeq2-versus-R/DESeq2 and Lumos-versus-limma concordance gates;
- monthly/manual GSE117769 and GSE186063 count-matrix end-to-end analyses with archived evidence;
- a real Docker build, in-container unit test and live Streamlit health check.

## Project structure


## Reproducible reporting checklist


## Development team

- Dr. João Paulo M. S. Lima —
- MSc. Laís de Carvalho Gonçalves —
- Rodrigo Arruda Orvate — 

Laboratory of Molecular Evolution and Systems Biology (EvoMol-Lab), BioME,
Federal University of Rio Grande do Norte (UFRN), Natal, Brazil.

## Core references

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion
  for RNA-seq data with DESeq2. *Genome Biology*. 2014;15:550.
  <https://doi.org/10.1186/s13059-014-0550-8>
- Muzellec B, et al. PyDESeq2: a Python package for bulk RNA-seq differential
  expression analysis. *Bioinformatics*. 2023;39(9):btad547.
  <https://doi.org/10.1093/bioinformatics/btad547>
- Xie Z, et al. Gene Set Knowledge Discovery with Enrichr. *Current Protocols*.
  2021;1:e90. <https://doi.org/10.1002/cpz1.90>
- Szklarczyk D, et al. The STRING database in 2023: protein-protein association
  networks and functional enrichment analyses for any sequenced genome of
  interest. *Nucleic Acids Research*. 2023;51(D1):D638-D646.
  <https://doi.org/10.1093/nar/gkac1000>
- Fornes O, et al. JASPAR 2020: update of the open-access database of
  transcription factor binding profiles. *Nucleic Acids Research*.
  2020;48(D1):D87-D92. <https://doi.org/10.1093/nar/gkz1001>
- Han H, et al. TRRUST v2: an expanded reference database of human and mouse
  transcriptional regulatory interactions. *Nucleic Acids Research*.
  2018;46(D1):D380-D386. <https://doi.org/10.1093/nar/gkx1013>
- Hagberg AA, Schult DA, Swart PJ. Exploring network structure, dynamics, and
  function using NetworkX. *Proceedings of SciPy 2008*. 2008:11-15.

## License, release and DOI

Copyright © 2026 EvoMol-Lab/BioME/UFRN. Lumos Networks is distributed under
the OSI-approved Apache License 2.0; see `LICENSE` and `NOTICE`.

`CITATION.cff`, `.zenodo.json` and `RELEASE_CHECKLIST.md` prepare a tagged GitHub
release for Zenodo archival. 

Contact: jpmslima@gmail.com or lais.goncalves.034@ufrn.edu.br
