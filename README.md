# Lumos Networks

**A transparent three-stage workflow for transcriptomic analysis, functional context, and prior-knowledge networks**

Lumos Networks is a research application developed by EvoMol-Lab at BioME/UFRN. It runs on your own
computer as a local web interface (Flask, in the style of [Sauron](https://github.com/evomol-lab/Sauron)
and [The Room of Conformations](https://github.com/evomol-lab/roc)) and as a command line. The suite
contains three connected modules:

1. **DDEA — Diagonal Differential Expression Alley**
2. **APP — Arithmancy Pathway Profiler**
3. **PG — PrioriGraph**

Everything is explained in the **Documentation** tab of the interface ([DOCUMENTATION.md](DOCUMENTATION.md)).

## What's new in the Flask version (2.0 beta)

- **Local web interface** replacing the Streamlit pages: numbered steps, background jobs with
  progress, log and cancellation, light and dark themes, and the whole session (datasets, analyses,
  integrations, enrichment) kept in the address bar.
- **Groups chosen directly in the sample table**: tick rows (Shift+click for ranges), filter by any
  column, or create one group per value of a metadata column in one click. The expression column
  matched to every sample is shown.
- **Optional Ensembl → gene symbol conversion**: Ensembl/Entrez IDs via MyGene.info (Ensembl REST as
  fallback), probes via the GEO platform annotation; a switch shows every plot, table and network with
  gene symbols or with the original IDs.
- **Several datasets and analyses per session**, each with optional batch, covariates and paired
  design, and one **Integrate** tab for contrasts of one dataset, random-effects meta-analyses and
  multi-tissue comparisons.
- **Interactive HTML report** of any tab, plus a JSON manifest for every analysis, integration and
  enrichment.

**Groups chosen in the sample table**
![Lumos Networks - samples and groups](screenshots/samples.jpg)

**Differential expression, with gene symbols or Ensembl IDs**
![Lumos Networks - DDEA](screenshots/ddea.jpg)

**Integration of analyses**
![Lumos Networks - integration](screenshots/integrate.jpg)

**PrioriGraph**
![Lumos Networks - PrioriGraph](screenshots/priorigraph.jpg)

## What each module computes

| Module | Implemented methods | Principal outputs |
|---|---|---|
| DDEA | PyDESeq2 Wald tests for raw integer RNA-seq counts; limma empirical-Bayes models for continuous expression; full-rank batch, covariate and paired-subject designs; alignment of contrasts; random-effects meta-analysis; multi-tissue alignment | Complete DE tables, design/contrast provenance, diagnostic figures, forest plots, cross-dataset alignment, JSON manifest, interactive HTML report |
| APP | Explicit-background ORA through Speedrichr (local hypergeometric fallback), separated by expression direction by default; STRING v12 networks; TRRUST/JASPAR target-set enrichment | KEGG/GO tables with query/universe sizes and engine, STRING associations, candidate-regulator tables with provenance, JSON manifest, interactive HTML report |
| PrioriGraph | Enriched-term → gene and candidate TF → target integration; HITS, PageRank, betweenness; personalized RWR, signed linear propagation and parameter-grid sensitivity | Interactive network, edge provenance, structural rankings, propagation and stability tables, interactive HTML report |

The README deliberately uses **candidate regulator**, **association**, and **prior-knowledge
network**. The implemented workflow does not establish cohort-specific TF binding, regulatory
causality, or pathway activation.

## Analytical workflow

```text
GEO series or uploaded tables + groups chosen in the sample table
              │
              ▼
 DDEA: effect estimates + uncertainty (one or several analyses)
              │
              ├──▶ Integrate: contrasts · meta-analysis · multi-tissue
              ▼
 APP: ORA + STRING + TF target-set enrichment
              │
              ▼
 PrioriGraph: inspectable topology + propagation
```

Every step is a job stored on the server; CSV, JSON and HTML downloads provide durable, reviewable
hand-offs.

## Scientific behaviour

### DDEA

- Raw integer counts are validated before PyDESeq2. Fractional, negative, missing, or infinite values
  are rejected instead of silently truncated.
- The analysis unit is chosen per analysis: the features as provided (gene symbols are then labels
  that can be switched on and off), or gene symbols, in which case count features mapping to the same
  symbol are summed before PyDESeq2 and microarray probes are summarized by per-sample median.
- Count filtering retains genes with at least 10 counts in at least the size of the smaller group.
- Continuous matrices retain finite, non-constant features; non-log values are log2-transformed and
  microarrays quantile-normalized.
- The condition effect can be adjusted for batch, categorical/continuous covariates and paired
  subjects. Lumos rejects incomplete pairs and rank-deficient/confounded designs before model fitting.
- The publication route for microarrays and other continuous matrices uses Bioconductor limma with
  robust, trend-aware empirical-Bayes moderation. Welch is retained only as an explicitly labelled
  exploratory fallback (used automatically only when R/limma is not installed).
- A DEG must pass both the selected FDR and absolute effect threshold; thresholds can be changed on
  the results page without refitting.
- Several named analyses can be kept for one dataset. Integrated, their results are aligned by gene
  and classified as consistent, shared, specific, or discordant; reused samples are not treated as
  independent studies and no pooled effect is calculated.
- Meta-analysis is intended for independent datasets from the same tissue with comparable contrasts;
  it reports confidence intervals, Q, I², tau², study count, and contributing accessions.
  Multi-tissue comparison keeps effects from different tissues separate and supports more than two
  datasets.
- The strict conserved/consistent signature requires every analysis to pass both thresholds with the
  same direction. It never substitutes minimum FDR or a mean effect for concordance.

### APP

- APP currently performs **ORA**, not ranked-list GSEA.
- Every ORA requires an explicit measured-gene universe: by default all genes tested in the gene
  source; users may instead paste a custom universe. Query genes outside the universe are excluded and
  counted.
- Default libraries are `KEGG_2021_Human` and `GO_Biological_Process_2023`; any Enrichr library can
  be added, and library names and engines are retained in output tables.
- Up- and down-regulated lists are queried separately by default. Direction labels describe input
  lists; ORA alone does not infer activation/inhibition.
- `TRRUST_Transcription_Factors_2019` and `JASPAR_PWM_Human_2025` are queried independently and
  retained as an evidence-source column (ChEA 2022 and ENCODE/ChEA consensus are optional).
- Terms naming mouse, rat or murine experiments are removed by whole-word matching.
- STRING edges (functional or physical network) retain the returned score and are labelled
  associations.

### PrioriGraph

- Node rankings are calculated on the displayed TF → target subnetwork.
- TF-target exports preserve database, term, direction and adjusted p-value.
- Personalized RWR conserves probability at dangling nodes and exposes the restart probability.
- Sensitivity analysis perturbs FDR, absolute-effect, RWR-restart and propagation-decay parameters and
  reports configuration coverage, top-k frequency, rank IQR and score variability for every node.

## Installation

### Local Python

Python 3.11 or newer.

```bash
git clone https://github.com/evomol-lab/lumos-networks.git
cd lumos-networks
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python LumosGUI.py                 # then open http://127.0.0.1:5070
```

limma also needs R with the Bioconductor package `limma` (`BiocManager::install("limma")`, or
`r-bioc-limma` on Debian/Ubuntu). Without it, continuous data use the Welch t-test and the interface
says so.

### Docker / Podman (includes R and limma)

```bash
docker build -t lumos-networks .
docker run --rm -p 5070:5070 -v "$PWD/lumos_jobs":/app/lumos_jobs lumos-networks
```

### Apptainer (HPC)

```bash
apptainer build lumos.sif Apptainer.def
mkdir -p lumos_jobs
apptainer run --bind lumos_jobs:/opt/lumos/lumos_jobs lumos.sif
```

### Command line

```bash
python lumos.py dataset --geo GSE117769 --outdir run/data
python lumos.py de --dataset run/data --groups groups.json --reference Healthy --test RA \
    --covariate gender --outdir run/de
python lumos.py enrich --source run/de --genes-up up.txt --genes-down down.txt --outdir run/app
```

Jobs are stored in `lumos_jobs/` (`LUMOS_JOBS_DIR`); `LUMOS_HOST`, `LUMOS_PORT` and
`LUMOS_JOB_RETENTION_DAYS` change the address and the clean-up.

## Tests

Pure numerical functions live in `lumos_core.py`; the engine (`lumos.py`) is tested on synthetic data
without network access.

```bash
python -m pip install pytest
python -m pytest tests                      # all tests (limma tests run when R/limma is installed)
python -m unittest discover -s tests -v     # lumos_core unit tests
docker run --rm lumos-networks python -m pytest -q /app/tests
```

The unit suite covers normalization, raw-count validation, experimental-design rank checks, paired
samples, GEO sample alignment, VST orientation, explicit ORA universes, benchmark metrics, conserved
signatures, RWR mass conservation and network sensitivity; the engine tests cover table parsing,
sample matching, DESeq2/limma/Welch runs that recover planted genes, design formulas, the three
integrations, PrioriGraph construction and the HTML report.

## Project structure

```text
lumos-networks/
├── lumos.py             # engine and command line (dataset, de, integrate, enrich; PG, STRING, reports)
├── lumos_core.py        # pure, tested analytical helpers (design, alignment, ORA inputs, RWR, sensitivity)
├── lumos_report.py      # self-contained interactive HTML reports
├── LumosGUI.py          # Flask web interface (jobs, progress, downloads, on-demand networks)
├── r/limma_backend.R    # limma backend called through Rscript
├── templates/, static/  # page, style, script and images of the interface
├── DOCUMENTATION.md     # Documentation tab (methods, interpretation, outputs)
├── tests/               # unit and engine tests
├── Dockerfile, Apptainer.def, LumosGUI.spec   # container images and stand-alone executable
├── requirements.txt, pyproject.toml           # Python dependencies
├── CITATION.cff, zenodo.json, NOTICE, LICENSE
└── assets/, screenshots/                       # logos and README images
```

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
- Ritchie ME, et al. limma powers differential expression analyses for
  RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015;43(7):e47.
  <https://doi.org/10.1093/nar/gkv007>
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
- Wu C, MacLeod I, Su AI. BioGPS and MyGene.info: organizing online, gene-centric
  information. *Nucleic Acids Research*. 2013;41(D1):D561-D565.
  <https://doi.org/10.1093/nar/gks1114>
- Hagberg AA, Schult DA, Swart PJ. Exploring network structure, dynamics, and
  function using NetworkX. *Proceedings of SciPy 2008*. 2008:11-15.

## License, release and DOI

Copyright © 2026 EvoMol-Lab/BioME/UFRN. Lumos Networks is distributed under
the OSI-approved Apache License 2.0; see `LICENSE` and `NOTICE`.

`CITATION.cff`, `.zenodo.json` and `RELEASE_CHECKLIST.md` prepare a tagged GitHub
release for Zenodo archival. 

Contact: jpmslima@gmail.com or lais.goncalves.034@ufrn.edu.br
