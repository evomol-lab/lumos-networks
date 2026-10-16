import os

import streamlit as st


VERSION = "1.3.0-beta"
st.set_page_config(page_title="Lumos Networks | Documentation", page_icon="📚", layout="wide")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)
LOGO_PATH = os.path.join(ROOT_DIR, "assets", "Lumos Networks.png")

st.markdown(
    """
    <style>
    [data-testid="stSidebarNav"] {display:none;}
    .stPageLink {background:rgba(151,166,195,.12);border-radius:18px;padding:8px;border:1px solid rgba(151,166,195,.3);}
    .doc-hero {padding:2.1rem 2.3rem;border-radius:24px;color:white;
      background:linear-gradient(120deg,#173b68 0%,#245f88 48%,#2f8f78 100%);
      box-shadow:0 12px 30px rgba(25,60,90,.16);margin-bottom:1.2rem;}
    .doc-hero h1 {margin:0 0 .35rem 0;font-size:2.35rem;}
    .doc-hero p {margin:0;opacity:.92;font-size:1.05rem;max-width:900px;}
    .doc-card {border:1px solid #e2e8f0;border-radius:18px;padding:1.15rem 1.25rem;
      min-height:150px;background:linear-gradient(180deg,#ffffff,#f8fafc);color:#1f2933;}
    .doc-card h3 {margin:.1rem 0 .45rem;color:#173b68;}
    .eyebrow {font-size:.78rem;letter-spacing:.08em;text-transform:uppercase;
      color:#4b718c;font-weight:700;}
    .good-note {border-left:5px solid #2f8f78;background:#eef9f5;padding:1rem 1.1rem;
      border-radius:8px;margin:.7rem 0;color:#1f2933;}
    .warn-note {border-left:5px solid #d99a23;background:#fff8e8;padding:1rem 1.1rem;
      border-radius:8px;margin:.7rem 0;color:#1f2933;}
    </style>
    """,
    unsafe_allow_html=True,
)

with st.sidebar:
    if os.path.exists(LOGO_PATH):
        st.image(LOGO_PATH, width=250)
    st.divider()
    st.markdown("### 🚀 Navigation")
    st.page_link("Lumos_Home.py", label="Home Page", icon="🏠")
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")
    st.divider()
    st.metric("Documentation release", f"v{VERSION}")
    st.caption("Methods, interpretation, provenance, and publication checklist.")

st.markdown(
    """
    <div class="doc-hero">
      <div class="eyebrow" style="color:#dcecf5">Lumos Networks · Technical guide</div>
      <h1>From expression data to inspectable biological hypotheses</h1>
      <p>DDEA estimates expression effects, APP adds functional and association context,
      and PG builds transparent prior-knowledge networks. This guide distinguishes every
      implemented calculation from the biological interpretation it can support.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

o1, o2, o3 = st.columns(3)
with o1:
    st.markdown(
        """<div class="doc-card"><div class="eyebrow">01 · DDEA</div><h3>Differential expression</h3>
        Effects, uncertainty, experimental designs, multiple contrasts, meta-analysis,
        and cross-tissue alignment.</div>""", unsafe_allow_html=True,
    )
with o2:
    st.markdown(
        """<div class="doc-card"><div class="eyebrow">02 · APP</div><h3>Functional context</h3>
        Explicit-background KEGG/GO ORA, STRING associations, and candidate-regulator
        target-set enrichment.</div>""", unsafe_allow_html=True,
    )
with o3:
    st.markdown(
        """<div class="doc-card"><div class="eyebrow">03 · PG</div><h3>Prior-knowledge networks</h3>
        Provenance-rich edges, topology, signal-propagation heuristics, and parameter
        sensitivity.</div>""", unsafe_allow_html=True,
    )

st.markdown("<br>", unsafe_allow_html=True)
tab_start, tab_ddea, tab_app, tab_pg, tab_repro = st.tabs([
    "Start here", "DDEA", "APP", "PG", "Reproducibility & limits",
])

with tab_start:
    st.header("Choose the analysis that matches the scientific question")
    st.dataframe(
        {
            "DDEA option": ["Single Dataset", "Multi-Dataset Meta-analysis", "Multi-Tissue Comparison"],
            "Use when": [
                "One GEO dataset; one or several contrasts among its sample groups",
                "Two or more independent datasets from the same tissue and a comparable contrast",
                "Two or more datasets representing different tissues",
            ],
            "Integration": [
                "Contrasts aligned by gene; no pooled effect when samples may be reused",
                "Study effects combined with a random-effects model",
                "Independent effects aligned by gene; no cross-tissue pooled effect",
            ],
            "Downstream APP focus": [
                "Consistent signature or any saved contrast",
                "Meta-analytic effect",
                "Strict conserved signature or any tissue profile",
            ],
        },
        width="stretch",
        hide_index=True,
    )
    st.markdown(
        """<div class="good-note"><b>Recommended order.</b> Inspect sample metadata and DDEA
        diagnostics first; then run APP with the complete tested-gene universe; finally use PG
        to examine the evidence network and its sensitivity. Downloads and manifests are the
        durable hand-off for a publication.</div>""",
        unsafe_allow_html=True,
    )
    st.subheader("What passes between modules")
    st.markdown(
        """
        - **DDEA → APP:** complete tested-gene results, effects, FDR values, analysis scope,
          and aligned profiles when applicable.
        - **APP → PG:** selected expression evidence, KEGG/GO results, STRING associations,
          and JASPAR/TRRUST candidate-regulator tables.
        - **Exports:** CSV tables and JSON manifests preserve thresholds, inputs, methods,
          versions, and mappings beyond the active Streamlit session.
        """
    )

with tab_ddea:
    st.header("DDEA — Diagonal Differential Expression Alley")
    st.subheader("1. Retrieve and verify")
    st.markdown(
        """
        DDEA retrieves GEO Series Matrix files and, when required, compatible supplementary
        tables. Select metadata columns that make sample identity and phenotype clear, then
        assign at least two unambiguous samples to each group. GEO layouts are heterogeneous:
        always compare the selected GSM accessions with the accession page and study design.
        """
    )
    st.subheader("2. Select the statistical route")
    st.dataframe(
        {
            "Input": ["Raw integer RNA-seq counts", "Microarray / continuous expression", "Exploratory continuous fallback"],
            "Engine": ["PyDESeq2 negative-binomial Wald test", "R/Bioconductor limma empirical Bayes", "Welch t-test"],
            "Key behavior": [
                "Rejects fractional/negative/non-finite counts; sums mapped count features",
                "Optional guarded log2, quantile normalization, mapped-probe median",
                "Explicitly labelled exploratory; not substituted silently for limma",
            ],
        },
        width="stretch",
        hide_index=True,
    )
    st.subheader("3. Model batch, covariates, and pairs")
    st.markdown(
        """
        Single Dataset exposes additive batch, categorical covariate, continuous covariate,
        and paired-subject terms. The displayed model formula and condition contrast are saved
        in the manifest. Lumos rejects incomplete pairs and rank-deficient or perfectly
        confounded designs before fitting. Arbitrary interactions and splines are not exposed.
        """
    )
    st.subheader("4. Save several analyses from one dataset")
    st.markdown(
        """
        Create sample groups, choose a reference and test group, give the analysis a unique
        label, and run it. Repeat for additional contrasts or alternative valid models. Use the
        saved-analysis selector to revisit each volcano plot, PCA, heatmap, table, and manifest.

        With two or more saved analyses, DDEA outer-joins results by gene and reports consistent,
        shared, specific, and discordant patterns. It also creates an effect heatmap and pairwise
        plots. Because contrasts from one dataset can reuse samples, **Lumos does not treat them
        as independent studies and does not pool them by meta-analysis**. APP can use either the
        strict consistent signature or any individual saved analysis; PG then receives that APP focus.
        """
    )
    st.subheader("5. Interpret the two multi-dataset options correctly")
    md1, md2 = st.columns(2)
    with md1:
        st.markdown("**Multi-Dataset Meta-analysis · same tissue**")
        st.markdown(
            """
            Each study is fitted independently. Gene-level log2 effects and standard errors are
            combined by a random-effects estimator (Paule–Mandel or DerSimonian–Laird). Outputs
            include confidence intervals, Q, I², tau², study count, and contributing accessions.
            Comparisons and biological tissue should be sufficiently comparable.
            """
        )
    with md2:
        st.markdown("**Multi-Tissue Comparison · different tissues**")
        st.markdown(
            """
            Add two or more tissue datasets. Effects remain tissue-specific and are aligned by
            symbol. Patterns distinguish conserved, shared, tissue-specific, and discordant genes.
            APP's strict conserved focus requires every tissue to pass both thresholds in the same
            direction. A pooled cross-tissue effect is intentionally not calculated.
            """
        )
    st.info(
        "A DEG must pass both the selected FDR and absolute log2-effect threshold. Changing "
        "thresholds changes classification, not the fitted effect estimate."
    )

with tab_app:
    st.header("APP — Arithmancy Pathway Profiler")
    st.subheader("KEGG and Gene Ontology")
    st.markdown(
        """
        APP performs **over-representation analysis (ORA)**, not ranked-list GSEA. It queries
        `KEGG_2021_Human` and `GO_Biological_Process_2023`, separating up- and down-regulated
        genes by default. Every query uses an explicit universe: normally all genes tested by
        DDEA, or a custom measured-gene universe supplied by the user.

        Speedrichr is attempted first because it supports a submitted background. If that service
        is unavailable, Lumos retrieves the same Enrichr gene-set library and performs the
        one-sided hypergeometric test locally, followed by Benjamini–Hochberg correction. The
        output records the engine, library, query size, effective universe, excluded genes, and UTC date.
        """
    )
    st.markdown(
        """<div class="warn-note"><b>Interpretation boundary.</b> “Up” and “Down” identify the input
        gene list. ORA does not, by itself, demonstrate activation or inhibition of a pathway.</div>""",
        unsafe_allow_html=True,
    )
    st.subheader("STRING associations")
    st.markdown(
        """
        APP maps submitted human symbols through STRING v12 before requesting the network.
        The export retains mapping coverage, network type, requested score threshold, returned
        association scores, and query provenance. Functional and physical networks are available.
        STRING edges are prior-knowledge associations; they are not interactions inferred from
        the expression cohort.
        """
    )
    st.subheader("Candidate regulators")
    st.markdown(
        """
        JASPAR 2024 and TRRUST 2019 target sets are tested independently using the same explicit-
        background ORA framework. A reported TF is a **candidate regulator whose annotated target
        set is enriched**. It is not proof of binding, activity, direction of regulation, or causality.
        """
    )

with tab_pg:
    st.header("PG — PrioriGraph")
    st.markdown(
        """
        PrioriGraph integrates enriched-term → member-gene edges and candidate-TF → annotated-
        target edges. Expression colors are applied only when the selected gene passes the active
        FDR/effect filter. Each export retains available database, term, direction, overlap, and
        adjusted-p-value provenance.
        """
    )
    st.subheader("Topology and propagation")
    st.markdown(
        """
        - Directed TF → target topology includes in/out-degree, HITS hubs/authorities, and PageRank.
        - Betweenness is computed on the undirected skeleton and labelled accordingly.
        - Personalized random walk with restart uses absolute expression evidence as seeds and
          redistributes dangling-node probability to the seed distribution.
        - Signed linear propagation computes `(I − λW)⁻¹E`. It is a transparent heuristic and is
          intentionally not called SPIA because signed regulatory coefficients and pathway-specific
          perturbation assumptions are absent.
        """
    )
    st.subheader("Sensitivity analysis")
    st.markdown(
        """
        The parameter grid perturbs FDR, absolute effect, restart probability, and propagation
        decay. It reports configuration coverage, top-k frequency, rank IQR, and score variability.
        This measures dependence on analytical choices; it is not a bootstrap and does not quantify
        biological uncertainty.
        """
    )

with tab_repro:
    st.header("Publication-oriented reproducibility")
    st.markdown(
        """
        1. Archive the DDEA complete table and JSON manifest for every fitted contrast.
        2. Report GEO accession, matrix source, GSM assignments, reference direction, model formula,
           filtering, FDR, and absolute-effect threshold.
        3. For meta-analysis, report the estimator and inspect Q, I², tau², and study influence.
        4. For multi-contrast or multi-tissue analyses, report the exact consistency rule and retain
           every profile-specific effect; do not describe the aligned mean as a pooled estimate.
        5. For ORA, archive query genes, effective background, library name, engine, date, and results.
        6. For STRING, report species, API version, network type, score, mapping coverage, and date.
        7. For PG, archive edges, evidence sources, topology definitions, propagation parameters,
           and sensitivity results.
        8. Validate central findings in an independent cohort or experimental system.
        """
    )
    st.subheader("Known limits")
    st.markdown(
        """
        - GEO metadata and supplementary layouts are heterogeneous; automated retrieval still
          requires human verification for every accession.
        - External GEO, Speedrichr/Enrichr, STRING, JASPAR, and TRRUST content can change. Preserve
          exported results and query dates.
        - Knowledge-base coverage and filtering affect network rankings. Sensitivity analysis
          does not replace independent validation.
        - Lumos is hypothesis-generating research software and is not intended for clinical diagnosis.
        """
    )
