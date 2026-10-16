# Lumos Networks: documentation

Lumos Networks takes transcriptome experiments from expression values to inspectable biological
hypotheses in three modules: **DDEA** estimates expression effects, **APP** adds functional and
association context, and **PG** (PrioriGraph) builds transparent prior-knowledge networks. This guide
describes every calculation and the biological interpretation it can (and cannot) support.

Lumos runs on your own computer as a local web application. Only accessions, gene lists and
identifiers are sent to the public services it uses (NCBI GEO, MyGene.info, Ensembl,
Enrichr/Speedrichr and STRING).

Contents

1. [Quick start](#quick-start)
2. [Choose the analysis](#choose-the-analysis)
3. [Step 1: datasets](#step-1-datasets)
4. [Step 2: samples and groups](#step-2-samples-and-groups)
5. [Step 3: differential expression (DDEA)](#step-3-differential-expression-ddea)
6. [Gene symbols for Ensembl, Entrez and probe IDs](#gene-symbols-for-ensembl-entrez-and-probe-ids)
7. [Integrating analyses](#integrating-analyses)
8. [Pathways and regulators (APP)](#pathways-and-regulators-app)
9. [PrioriGraph (PG)](#priorigraph-pg)
10. [Reports, manifests and files](#reports-manifests-and-files)
11. [Command line](#command-line)
12. [Installation and settings](#installation-and-settings)
13. [Reproducibility checklist and known limits](#reproducibility-checklist-and-known-limits)
14. [Changes in the Flask version](#changes-in-the-flask-version)
15. [References](#references)

---

## Quick start

```
pip install -r requirements.txt
python LumosGUI.py            # then open http://127.0.0.1:5070
```

1. **Dataset**: type a GEO series accession (e.g. `GSE117769`) or upload an expression table and,
   optionally, a sample table. Several datasets can be loaded in one session.
2. **Samples and groups**: tick samples in the table and add them to groups, or create one group per
   value of a column (e.g. *phenotype*) in one click. Choose the reference and the test group.
3. **Differential expression**: give the analysis a label, optionally add batch, covariates or a
   paired design, and run it. Repeat for other comparisons or models: every run is kept.
4. **Results**: volcano, MA and PCA plots, top genes, heatmap and the gene table. Thresholds,
   highlighted genes and **gene labels** (symbol or original ID) change at any time.
5. **Integrate** (optional): align contrasts of one dataset, meta-analyse independent datasets, or
   compare tissues gene by gene.
6. **APP**: choose the gene source (an analysis, a consistent signature, a profile or a meta-analysis),
   tick genes and run the enrichment; then build the STRING network.
7. **PG**: the PrioriGraph tab joins enriched terms, genes and candidate regulators, with rankings,
   standard network metrics, signal propagation and parameter sensitivity.

Every step runs as a separate job whose progress and log are shown and which can be cancelled.
The address bar keeps the whole session (datasets, analyses, integrations, enrichment), so a reload
or a bookmark reopens it. **HTML report** exports the tab being shown as a self-contained,
interactive file.

---

## Choose the analysis

| Question | What to run | Integration | Gene source for APP |
|---|---|---|---|
| One dataset, one comparison | one analysis | none | the analysis |
| One dataset, several comparisons or alternative models | several analyses of the same dataset | **Contrasts of one dataset**: aligned by gene; no pooled effect, because contrasts may reuse samples | consistent signature, or any contrast |
| Independent datasets of the **same tissue** and a comparable contrast | one analysis per dataset | **Meta-analysis**: random-effects combination of the study effects | meta-analytic effect |
| Datasets of **different tissues** | one analysis per dataset | **Multi-tissue comparison**: aligned by gene; no cross-tissue pooled effect | strict conserved signature, or any tissue |

Recommended order: inspect the sample metadata and the DDEA diagnostics first; run APP against the
complete tested-gene universe; then use PG to examine the evidence network and its sensitivity.
Downloads and manifests are the durable hand-off for a publication.

---

## Step 1: datasets

### From GEO

Lumos reads the **series matrix** of the accession for the sample descriptions (accession, title,
source and every `characteristics` field, such as *phenotype*, *treatment* or *age*). The expression
values are taken, in this order, from:

| Experiment | Expression source |
|---|---|
| Microarray | the table inside the series matrix (one value per probe and sample) |
| RNA-Seq | a supplementary matrix file (names with *count*, *raw* or *matrix* first), checked by how many of its columns match the samples |
| RNA-Seq | the per-sample files inside `*_RAW.tar`, joined into one matrix (only when no matrix file matched; the archive can be large) |
| either | the series-matrix table, as a last resort |

The column that holds each sample is matched by accession (`GSM…`) or title and shown in the sample
table; samples without expression values are greyed out. If no usable table exists, upload the count
files under **Count files for this series**. Series with several platforms have one series matrix per
platform: the metadata of all samples are listed and the expression table with most samples is used.
GEO layouts are heterogeneous: always compare the matched columns with the accession page.

### Uploaded tables

* **Expression**: one genes × samples matrix (CSV, TSV or whitespace-separated text, optionally
  gzipped) or one file per sample (HTSeq, featureCounts, STAR `ReadsPerGene`, salmon `quant.sf`,
  RSEM `genes.results`). Annotation columns are ignored; a gene-name column (`gene_name`, `symbol`…)
  is kept as gene symbols. HTSeq/STAR summary rows are removed and duplicated identifiers merged.
* **Sample table** (optional): one row per sample; the column that best matches the expression column
  names becomes the sample identifier. Other columns can define groups and design terms (batch,
  subject, age…).

### What Lumos detects

* **Identifiers**: Ensembl gene IDs (any species prefix; versions such as `.12` allowed), NCBI Gene
  (Entrez) IDs, microarray probe IDs or gene symbols.
* **Values**: raw counts (non-negative integers), log-scale values (any negative value or a maximum
  below 50) or other non-integer values (TPM, FPKM, intensities). This decides the default method.

---

## Step 2: samples and groups

Groups are defined **directly in the sample table**:

* Click a row to tick it; **Shift**+click ticks a range; the header box ticks every row visible after
  filtering. Filter with any word of any column.
* **Add ticked to** a group, to a **new group**, or **remove** them from their group. A sample belongs
  to at most one group. Group names can be edited in their chips.
* **One group per value of…** creates the groups in one step from a column with 2 to 12 values. A
  group whose name looks like a control (control, healthy, normal, WT, untreated, vehicle…) is
  proposed as the reference.
* Columns with a single value are hidden unless you ask for them.

The comparison is **Test vs Reference**: a positive log2 fold change means higher expression in the
test group. Each group needs at least two samples. Groups are remembered for each dataset.

---

## Step 3: differential expression (DDEA)

### Statistical routes

| Input | Engine | Behaviour |
|---|---|---|
| Raw integer counts | **PyDESeq2** negative-binomial Wald test | fractional, negative, missing or infinite counts are rejected, never rounded |
| Microarray / continuous expression | **limma** (R/Bioconductor) empirical Bayes, `trend = TRUE, robust = TRUE` | optional log2, quantile normalisation for microarrays |
| Continuous expression, exploratory | **Welch t-test** | labelled exploratory; ignores design terms; used automatically only when R/limma is not installed |

**DESeq2.** The read count $K_{ij}$ of gene $i$ in sample $j$ is modelled with a negative binomial
generalised linear model:

$$K_{ij} \sim NB(\mu_{ij}, \alpha_i), \qquad \mu_{ij} = s_j \, q_{ij}, \qquad \log_2 q_{ij} = x_j^\top \beta_i$$

* $s_j$ is the **size factor** of sample $j$ (median-of-ratios normalisation), which corrects for
  sequencing depth and composition bias.
* $\alpha_i$ is the gene-specific **dispersion** (variance beyond Poisson noise), shrunk towards a
  fitted mean–dispersion trend (empirical Bayes).
* $x_j$ holds the design terms (batch, covariates, subject) and the condition; the condition
  coefficient is the log2 fold change, tested with a **Wald test**. Cook's distances flag outliers and
  independent filtering removes genes with too little information; such genes keep no adjusted
  p-value and are omitted from the volcano plot.

**limma.** For log-scale values $y_{gj}$ a linear model $E(y_g) = X\beta_g$ is fitted per gene; the
gene variances are moderated towards a mean–variance trend (empirical Bayes, robust to outlier genes),
which stabilises inference with few samples. The condition coefficient is tested with a moderated t
statistic.

**Welch t-test.** Values are transformed (counts → log2(CPM + 1); other non-log values → log2(x + 1);
log values are kept) and, for microarrays, quantile-normalised. The log2 fold change is the difference
of the group means and the Welch test does not assume equal variances.

All p-values are adjusted with the **Benjamini–Hochberg** false discovery rate.

### Experimental design: batch, covariates and pairs

Under *Experimental design*, choose an additive **batch**, a **subject/pair** identifier (paired
design), **categorical covariates** and **continuous covariates**. The model formula is shown (for
example `~ batch + age + condition`) and the condition effect is estimated after these terms
(DESeq2 and limma). Lumos rejects incomplete pairs and rank-deficient or perfectly confounded designs
before fitting; interactions and splines are not offered.

### Analysis unit

* **Features as provided** (default): Ensembl IDs, probes… are tested as they are; gene symbols are
  labels that the plots can switch on and off.
* **Gene symbols**: counts of features with the same symbol are **summed** before DESeq2; other
  values are summarised by their **median** after normalisation (features without a symbol are
  dropped).

### Filtering, thresholds and several analyses

Counts need at least *Minimum count* (default 10) in as many samples as the smaller group; continuous
values need complete, non-constant rows. A gene is **up-** or **down-regulated** when FDR <
threshold and |log2FC| ≥ threshold; both can be changed on the results page without refitting.

Each run is an **analysis** with its own label. Run several comparisons of the same dataset, or the
same comparison with alternative models, and switch between them with *Analysis shown*.

### Plots

Volcano (click points to label them), MA, PCA of the 500 most variable genes (VST for DESeq2, the
transformed values otherwise), top up/down genes, heatmap of row z-scores (optional Ward clustering of
genes and samples) and library sizes.

---

## Gene symbols for Ensembl, Entrez and probe IDs

With **Convert … to gene symbols** on (default), every tested identifier gets a gene symbol once:

| Identifiers | Source |
|---|---|
| Ensembl gene IDs | MyGene.info (`ensembl.gene`), with the Ensembl REST API as fallback; versions are removed and the species follows the ID prefix (ENSG human, ENSMUSG mouse, ENSRNOG rat…) |
| Entrez IDs | MyGene.info (`entrezgene`), with the GEO taxon when known |
| Probe IDs (GEO microarrays) | the GEO annotation of the platform; a probe annotated to several genes takes the first symbol that is not a microRNA, LOC/LINC or antisense transcript |
| Uploaded table with a gene-name column | that column |

Symbols are cached per dataset. **Gene labels** switches every plot, table and network between
**gene symbols** and the **original IDs**. Enrichr, STRING and the integrations always use symbols.

---

## Integrating analyses

Open **Integrate**, tick two or more analyses and choose the integration.

**Contrasts of one dataset.** Results are joined by gene symbol and each gene is classified per
contrast (up, down, not significant) and across contrasts: *Consistent up/down (all contrasts)*,
*Shared*, *Specific …* and *Discordant*. Because contrasts of one dataset can reuse samples, they are
not treated as independent studies and no pooled effect is calculated.

**Meta-analysis (same tissue).** Each dataset is fitted independently. For gene $g$ with study
effects $\hat\beta_k$ and variances $v_k = SE_k^2$, the random-effects estimate is

$$\hat\beta_{RE} = \frac{\sum_k w_k^* \hat\beta_k}{\sum_k w_k^*}, \qquad w_k^* = \frac{1}{v_k + \hat\tau^2}$$

with the between-study variance $\hat\tau^2$ estimated by **Paule–Mandel** (default) or
**DerSimonian–Laird**. Lumos reports the 95% confidence interval, Cochran's $Q$, $I^2 = \max(0, (Q -
df)/Q)$, $\tau^2$, the number of studies and their accessions; p-values are BH-adjusted. When a gene
has several features in a study, the one with the smallest p-value represents it. Use datasets of the
same tissue with comparable contrasts, and inspect heterogeneity (forest plot, I²) before trusting a
pooled effect.

**Multi-tissue comparison.** Each tissue keeps its own effect; genes are aligned by symbol and
classified as *Conserved*, *Shared*, *Specific* or *Discordant*. A pooled cross-tissue effect is
intentionally not calculated.

**Strict consistent / conserved signature** (APP focus): every analysis must pass both thresholds
with the same direction. It never substitutes the minimum FDR or a mean effect for concordance; its
reported effect is the mean log2FC and its FDR the largest one.

---

## Pathways and regulators (APP)

### Over-representation with an explicit background

APP performs **over-representation analysis (ORA)**, not ranked-list GSEA. For a background of $N$
measured genes, a gene set with $M$ of them and a query of $n$ genes of which $k$ are in the set:

$$P(X \ge k) = \sum_{i=k}^{\min(n, M)} \frac{\binom{M}{i}\binom{N-M}{n-i}}{\binom{N}{n}}$$

* **Background**: by default every gene tested in the gene source (the universe of genes that could
  have been selected, not the whole genome); a custom list can be pasted. Query genes outside it are
  excluded and counted.
* **Engine**: Speedrichr (Enrichr with an explicit background) is used first; if it is unavailable,
  Lumos downloads the same Enrichr library and applies the one-sided hypergeometric test locally,
  restricting the universe to genes covered by the library, followed by Benjamini–Hochberg
  correction. The engine is recorded in every table.
* **Gene lists**: up- and down-regulated genes are queried **separately** by default; *Combined* uses
  one list. Direction labels describe the input list; ORA alone does not show that a pathway is
  activated or inhibited.
* **Libraries**: KEGG 2021 and GO Biological Process 2023 by default; any Enrichr library can be
  added. Terms naming mouse, rat or murine experiments are removed unless kept.
* **Gene source**: an analysis, the consistent/conserved signature or a single profile of an
  integration, or a meta-analysis.

### Candidate regulators

TRRUST 2019 (curated TF → target interactions) and JASPAR PWM 2025 (binding motifs) are tested by
default, each on its own and kept as evidence source; ChEA 2022 and ENCODE/ChEA consensus (ChIP-seq)
can be added. A reported factor is a **candidate regulator whose annotated target set is enriched**:
not proof of binding, activity, direction of regulation or causality.

### STRING associations

The selected genes with the largest effects are mapped to STRING v12 (species from the gene source)
and the **functional** or **physical** network is requested at a minimum combined score. Mapping
coverage, network type, score and query date are reported. STRING edges are prior-knowledge
associations, not interactions inferred from the cohort.

---

## PrioriGraph (PG)

### Knowledge-driven versus data-driven networks

Co-expression methods such as WGCNA build the adjacency from correlation across samples,
$a_{ij} = |\mathrm{cor}(x_i, x_j)|^\beta$, which is powerful for discovery but does not imply an
interaction and is sensitive to batch effects. PrioriGraph projects the selected genes on curated
knowledge, $V_{network} = V_{DEG} \cap V_{curated}$: an edge exists only where a database supports it.

### Construction

1. The **top enriched terms** (KEGG, GO…; FDR ≤ chosen threshold, merged across up/down lists) become
   hexagons.
2. Their selected genes that are targets of a **candidate regulator** (FDR ≤ chosen threshold) become
   circles, coloured by direction, linked to the term (dotted edges).
3. The regulators become squares with arrows to their targets; a regulator that is itself selected has
   a coloured border. Every edge keeps its database, term, list and FDR (**Edges with provenance**).

### Rankings and metrics

* **Regulators**: unique targets in the network (out-degree), targets among all selected genes, best
  FDR and evidence sources.
* **Targets**: number of regulators (in-degree) and terms, with the expression evidence.
* **Terms**: the regulator with most targets among each term's genes.
* **Metrics** on the directed regulator → target layer: in/out-degree, **HITS** hub and authority
  scores, **PageRank**, and **betweenness** on the undirected skeleton (bridges between processes).

### Signal propagation

* **Personalised random walk with restart (RWR)**: with column-normalised adjacency $W$ (regulator →
  target), seed vector $s \propto |log_2FC|$ and restart probability $r$,
  $p_{t+1} = (1 - r)\, W p_t + r\, s$, iterated to convergence; probability at nodes without
  outgoing edges is returned to the seeds, so the total stays 1.
* **Signed linear propagation**: $(I - \lambda W)^{-1} E$ with decay $\lambda$ and signed effects
  $E$. It is a transparent heuristic, intentionally not called SPIA, because signed regulatory
  coefficients and pathway-specific perturbation assumptions are absent.

### Sensitivity analysis

A grid of FDR and |log2FC| thresholds, restart probabilities and decays is evaluated; for every node
Lumos reports configuration coverage, top-k frequency, median rank, rank IQR and score variability.
This measures dependence on analytical choices; it is not a bootstrap and does not quantify
biological uncertainty.

---

## Reports, manifests and files

**HTML report** (results header) exports the tab being shown: DDEA, Integrate, APP or PG. The report
is a single offline file with interactive Plotly figures (drawn in a light theme), searchable,
sortable tables with CSV/TSV export, a table of contents and *Save as PDF*.

| Step | File | Content |
|---|---|---|
| dataset | `expression.tsv.gz`, `samples.csv` | expression table as read; sample table with the matched column |
| analysis | `de_results.csv` | every tested gene: ID, symbol, mean, log2FC, SE, statistic, p-value, FDR, status |
| analysis | `normalized_expression.tsv.gz` | values used for PCA and heatmap |
| analysis | `manifest.json` | dataset, samples and expression columns, model formula and design metadata, contrast, method, unit, filter, thresholds, gene-symbol source, package versions |
| integration | `integration.csv`, `manifest.json` | gene-level table; rule, analyses (with their manifests) and versions |
| APP | `enrichment_<library>.csv`, `regulators.csv` | every term with list, engine and genes; regulators with library, list and targets |
| APP | `evidence.csv`, `background.txt`, `manifest.json` | expression evidence, the explicit background, query genes, libraries and engines |

Every table on the page has its own download button.

---

## Command line

`lumos.py` runs every step without the web interface:

```
python lumos.py dataset --geo GSE117769 --outdir run/data
python lumos.py de --dataset run/data --groups groups.json --reference Healthy --test RA \
                   --covariate gender --label "RA vs Healthy (sex-adjusted)" --outdir run/de_ra
python lumos.py de --dataset run/data --groups groups.json --reference Healthy --test PsA --outdir run/de_psa
python lumos.py integrate --mode contrasts --analyses run/de_ra run/de_psa --outdir run/int
python lumos.py enrich --source run/int --focus consistent --genes-up up.txt --genes-down down.txt \
                   --outdir run/app
```

`groups.json` maps group names to sample IDs (the first column of `samples.csv`). Options include
`--engine auto|deseq2|limma|welch`, `--batch`, `--pair`, `--covariate`, `--continuous`,
`--unit features|symbols`, `--no-symbols`; `integrate --mode contrasts|meta|tissues --estimator
Paule-Mandel|DerSimonian-Laird`; `enrich --strategy separate|combined --background FILE
--libraries … --tf-libraries …`. Use `python lumos.py <command> -h` for all of them.

---

## Installation and settings

```
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python LumosGUI.py
```

limma needs R with the Bioconductor package (`BiocManager::install("limma")`, or `r-bioc-limma` on
Debian/Ubuntu). Without it, continuous data use the Welch t-test and the interface says so. The
container images include R and limma:

```
docker build -t lumos-networks .
docker run --rm -p 5070:5070 -v "$PWD/lumos_jobs":/app/lumos_jobs lumos-networks

apptainer build lumos.sif Apptainer.def
mkdir -p lumos_jobs && apptainer run --bind lumos_jobs:/opt/lumos/lumos_jobs lumos.sif
```

| Variable | Default | Meaning |
|---|---|---|
| `LUMOS_HOST` | `127.0.0.1` | address to listen on (`0.0.0.0` in containers) |
| `LUMOS_PORT` | `5070` | port |
| `LUMOS_JOBS_DIR` | `./lumos_jobs` | where jobs are stored |
| `LUMOS_CACHE_DIR` | `<jobs>/_cache` | downloaded Enrichr libraries (local ORA fallback) |
| `LUMOS_JOB_RETENTION_DAYS` | `0` | delete jobs older than this at start-up (0 keeps all) |
| `LUMOS_MAX_DOWNLOAD_MB` | `2048` | largest GEO supplementary file or archive downloaded |

---

## Reproducibility checklist and known limits

1. Archive the analysis table and `manifest.json` of every fitted contrast.
2. Report the GEO accession, matrix source, sample assignments, reference direction, model formula,
   filtering, FDR and absolute-effect thresholds.
3. For meta-analyses, report the estimator and inspect Q, I², τ² and the forest plots.
4. For multi-contrast or multi-tissue analyses, report the consistency rule and keep every
   profile-specific effect; do not describe the aligned mean as a pooled estimate.
5. For ORA, archive the query genes, background, library names, engine, date and results.
6. For STRING, report species, version, network type, score, mapping coverage and date.
7. For PG, archive the edges with provenance, metric definitions, propagation parameters and
   sensitivity results.
8. Validate central findings in an independent cohort or experimental system.

Known limits:

* Two-condition comparisons with additive design terms; interactions and splines are not offered.
* GEO metadata and supplementary layouts are heterogeneous: check the matched samples for every
  accession.
* External resources (GEO, Speedrichr/Enrichr, STRING, JASPAR, TRRUST, MyGene.info) change over
  time: keep the exported results and query dates. The default libraries are human; for other species
  add the corresponding Enrichr libraries.
* Knowledge-base coverage and filtering affect network rankings; sensitivity analysis does not replace
  independent validation.
* Lumos is hypothesis-generating research software and is not intended for clinical diagnosis.

---

## Changes in the Flask version

Lumos Networks 2.0 (beta) moves the 1.3 Streamlit suite to a local Flask interface in the style of
Sauron and The Room of Conformations, keeping its methods (`lumos_core.py`, `lumos_report.py`,
`r/limma_backend.R`):

* **Interface**: one page with numbered steps, background jobs with progress, log and cancellation,
  light/dark themes, Documentation tab, and the session in the address bar.
* **Groups chosen in the sample table**, with filtering, range selection and one-click grouping by a
  metadata column; the matched expression column is shown for every sample.
* **Gene symbols** for Ensembl/Entrez/probe IDs are optional, cached, and switchable in every plot.
* **Several datasets and analyses per session**, integrated from one tab (contrasts, meta-analysis,
  multi-tissue), with design terms available in every analysis.
* **Fixes**: the mouse/rat term filter matches whole words (a substring test also removed GO terms
  such as *proliferation* or *migration*); TF names are parsed correctly (`TP53 human` → TP53); the
  JASPAR library is `JASPAR_PWM_Human_2025` (`JASPAR_2024` does not exist in Enrichr); probes
  annotated to several genes no longer default to a microRNA symbol.

---

## References

* Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq
  data with DESeq2. *Genome Biology* 15:550.
* Muzellec B et al. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression
  analysis. *Bioinformatics* 39:btad547.
* Ritchie ME et al. (2015). limma powers differential expression analyses for RNA-sequencing and
  microarray studies. *Nucleic Acids Res* 43:e47.
* Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. *J R Stat Soc B* 57:289–300.
* DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. *Control Clin Trials* 7:177–188;
  Paule RC, Mandel J (1982). Consensus values and weighting factors. *J Res Natl Bur Stand* 87:377–385.
* Seabold S, Perktold J (2010). Statsmodels: econometric and statistical modeling with Python.
* Barrett T et al. (2013). NCBI GEO: archive for functional genomics data sets. *Nucleic Acids Res*
  41:D991–D995.
* Wu C, MacLeod I, Su AI (2013). BioGPS and MyGene.info. *Nucleic Acids Res* 41:D561–D565.
* Xie Z et al. (2021). Gene set knowledge discovery with Enrichr. *Current Protocols* 1:e90.
* Han H et al. (2018). TRRUST v2. *Nucleic Acids Res* 46:D380–D386; Rauluseviciute I et al. (2024).
  JASPAR 2024. *Nucleic Acids Res* 52:D174–D182.
* Szklarczyk D et al. (2023). The STRING database in 2023. *Nucleic Acids Res* 51:D638–D646.
* Köhler S et al. (2008). Walking the interactome for prioritization of candidate disease genes.
  *Am J Hum Genet* 82:949–958.
* Kleinberg JM (1999). Authoritative sources in a hyperlinked environment. *J ACM* 46:604–632;
  Page L et al. (1999). The PageRank citation ranking.
* Langfelder P, Horvath S (2008). WGCNA. *BMC Bioinformatics* 9:559.
* Hagberg AA, Schult DA, Swart PJ (2008). Exploring network structure, dynamics, and function using
  NetworkX. *Proc. SciPy 2008*.
