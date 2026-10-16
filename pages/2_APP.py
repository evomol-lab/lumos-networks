import os
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import cloudscraper
import numpy as np
import pandas as pd
import plotly.express as px
import requests
import streamlit as st
from requests.adapters import HTTPAdapter
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
from streamlit_agraph import agraph, Node, Edge, Config
from urllib3.util.retry import Retry

from lumos_core import conserved_signature, prepare_ora_inputs
from lumos_report import build_html_report, table_spec, network_figure


# ============================================================
# PAGE / THEME
# ============================================================

st.set_page_config(
    page_title="Lumos Networks | Functional",
    page_icon="🧪",
    layout="wide"
)

st.markdown("""
<style>
[data-testid="stSidebarNav"] {display: none;}
.stPageLink {
    background-color: rgba(151, 166, 195, 0.12);
    border-radius: 20px;
    padding: 8px;
    border: 1px solid rgba(151, 166, 195, 0.3);
}
.section-header {
    color: #2E86C1;
    border-bottom: 2px solid #2E86C1;
    padding-bottom: 5px;
    font-weight: bold;
}
</style>
""", unsafe_allow_html=True)


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)

logo_path = os.path.join(PARENT_DIR, "assets", "APP.png")
lumos_logo = os.path.join(PARENT_DIR, "assets", "Lumos Networks.png")


# ============================================================
# SIDEBAR NAVIGATION
# ============================================================

with st.sidebar:
    if os.path.exists(lumos_logo):
        st.image(lumos_logo, width=250)

    st.divider()
    st.markdown("### 🚀 Navigation")
    st.page_link("Lumos_Home.py", label="Home Page", icon="🏠")
    
    st.markdown('<p style="color:#2E86C1; font-weight:bold; margin-bottom:0px; margin-top:10px;">📊 Analysis</p>', unsafe_allow_html=True)
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    
    st.markdown('<p style="color:#28B463; font-weight:bold; margin-bottom:0px; margin-top:10px;">🧬 Functional</p>', unsafe_allow_html=True)
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    
    st.markdown('<p style="color:#E67E22; font-weight:bold; margin-bottom:0px; margin-top:10px;">🕸️ Networks</p>', unsafe_allow_html=True)
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

    st.markdown('<p style="color:#E1AF12; font-weight:bold; margin-bottom:0px; margin-top:10px;">📚 Documentation</p>', unsafe_allow_html=True)
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")

    st.divider()

    if os.path.exists(logo_path):
        st.image(logo_path, width="stretch")

    st.info("You are in the functional module.")


# ============================================================
# HELPERS
# ============================================================

def _clean_genes(values):
    if values is None:
        return []

    if isinstance(values, pd.Series):
        values = values.tolist()

    clean = []
    for g in values:
        if pd.isna(g):
            continue
        g = str(g).strip().upper()
        if g and g != "NAN":
            clean.append(g)

    return sorted(set(clean))


def _detect_analysis_mode(df):
    if df is None or df.empty:
        return "none"

    if (
        "Analysis_Scope" in df.columns
        and df["Analysis_Scope"].astype(str).eq("single_dataset_multi_contrast").any()
    ):
        return "multi_contrast"
        
    effect_sets = {c.split("_", 1)[1] for c in df.columns if c.startswith("Log2FC_")}
    fdr_sets = {c.split("_", 1)[1] for c in df.columns if c.startswith("FDR_")}
    if len(effect_sets.intersection(fdr_sets)) >= 2:
        return "multi_tissue"

    if "Meta_Log2FC" in df.columns:
        return "meta"

    if "Log2FC" in df.columns:
        return "single"

    return "unknown"


def _standardize_result_dataframe(df, mode):
    df = df.copy()

    if "Symbol" not in df.columns:
        candidates = ["Gene", "gene", "symbol", "SYMBOL", "gene_symbol"]
        found = next((c for c in candidates if c in df.columns), None)
        if found:
            df["Symbol"] = df[found]

    if mode not in {"multi_tissue", "multi_contrast"}:
        if "FDR" not in df.columns:
            for candidate in ["Adjusted P-value", "padj", "FDR_adjusted", "adj.P.Val"]:
                if candidate in df.columns:
                    df["FDR"] = pd.to_numeric(df[candidate], errors="coerce")
                    break

        if mode == "meta":
            if "Meta_Log2FC" not in df.columns and "Log2FC" in df.columns:
                df["Meta_Log2FC"] = pd.to_numeric(df["Log2FC"], errors="coerce")
            df["Effect"] = pd.to_numeric(df["Meta_Log2FC"], errors="coerce")
        else:
            if "Log2FC" not in df.columns and "Meta_Log2FC" in df.columns:
                df["Log2FC"] = pd.to_numeric(df["Meta_Log2FC"], errors="coerce")
            df["Effect"] = pd.to_numeric(df["Log2FC"], errors="coerce")

        df["FDR"] = pd.to_numeric(df["FDR"], errors="coerce")

    if "Symbol" in df.columns:
        df["Symbol"] = (
            df["Symbol"]
            .astype(str)
            .str.strip()
            .str.upper()
        )

    if "Status" not in df.columns and mode not in {"multi_tissue", "multi_contrast"}:
        df["Status"] = "Not Significant"

    return df


def _selected_gene_symbols(df):
    if df.empty or "Symbol" not in df.columns:
        return []
    return _clean_genes(df["Symbol"].tolist())


# ============================================================
# ENRICHR
# ============================================================

def _retry_session(use_cloudscraper=False):
    retry = Retry(
        total=3, connect=3, read=3, status=3, backoff_factor=0.6,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}), raise_on_status=False,
    )
    session = (
        cloudscraper.create_scraper(
            browser={"browser": "chrome", "platform": "linux", "desktop": True}
        )
        if use_cloudscraper else requests.Session()
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=8, pool_maxsize=8)
    session.mount("https://", adapter)
    session.headers.update({
        "User-Agent": "LumosNetworks/1.3.0 (academic research; EvoMol-Lab/UFRN)"
    })
    return session


def _finalize_ora_table(df, clean_genes, clean_background, excluded, gene_sets, engine):
    if df.empty:
        return df
    df = df.copy()
    if "Genes_list" not in df.columns:
        df["Genes_list"] = df["Genes"].astype(str).str.split(";")
    df["Genes"] = df["Genes_list"].apply(lambda values: ";".join(map(str, values)))
    if "Overlap" not in df.columns:
        df["Overlap"] = df["Genes_list"].apply(
            lambda values: f"{len(values)}/{len(clean_genes)}"
        )
    df["Input_Gene_Count"] = len(clean_genes)
    df["Background_Gene_Count"] = len(clean_background)
    df["Excluded_Query_Gene_Count"] = len(excluded)
    df["Gene_Set_Library"] = gene_sets
    df["ORA_Engine"] = engine
    df = df[
        ~df["Term"].astype(str).str.contains(
            "Mouse|Mus musculus|Rat|Murine", case=False, na=False
        )
    ].copy()
    df["P-value"] = pd.to_numeric(df["P-value"], errors="coerce")
    df["Adjusted P-value"] = pd.to_numeric(df["Adjusted P-value"], errors="coerce")
    df["-log10(FDR)"] = -np.log10(
        df["Adjusted P-value"].fillna(1).clip(lower=1e-300)
    )
    df["Query_UTC"] = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    return df.sort_values(["Adjusted P-value", "P-value"], na_position="last")


def _run_speedrichr(clean_genes, clean_background, excluded, gene_sets):
    session = _retry_session(use_cloudscraper=True)
    response = session.post(
        "https://maayanlab.cloud/speedrichr/api/addList",
        files={
            "list": (None, "\n".join(clean_genes)),
            "description": (None, "Lumos_APP"),
        },
        timeout=(10, 35),
    )
    response.raise_for_status()
    user_list_id = response.json()["userListId"]

    response = session.post(
        "https://maayanlab.cloud/speedrichr/api/addbackground",
        data={"background": "\n".join(clean_background)}, timeout=(10, 35),
    )
    response.raise_for_status()
    background_id = response.json()["backgroundid"]

    response = session.post(
        "https://maayanlab.cloud/speedrichr/api/backgroundenrich",
        data={
            "userListId": user_list_id,
            "backgroundid": background_id,
            "backgroundType": gene_sets,
        },
        timeout=(10, 45),
    )
    response.raise_for_status()
    data = response.json()
    rows = data.get(gene_sets, [])
    if not rows:
        empty = pd.DataFrame()
        empty.attrs["diagnostic"] = "Speedrichr completed but returned no terms."
        return empty
    columns = [
        "Rank", "Term", "P-value", "Odds Ratio", "Combined Score",
        "Genes_list", "Adjusted P-value", "Old P-value", "Old Adjusted P-value",
    ]
    width = len(rows[0])
    if width > len(columns):
        columns.extend([f"Extra_{i}" for i in range(width - len(columns))])
    df = pd.DataFrame(rows, columns=columns[:width])
    return _finalize_ora_table(
        df, clean_genes, clean_background, excluded, gene_sets,
        "Speedrichr explicit-background API",
    )


@st.cache_data(show_spinner=False, ttl=86400)
def _fetch_enrichr_library(gene_sets):
    session = _retry_session()
    response = session.get(
        "https://maayanlab.cloud/Enrichr/geneSetLibrary",
        params={"mode": "text", "libraryName": gene_sets}, timeout=(10, 60),
    )
    response.raise_for_status()
    library = {}
    for line in response.text.splitlines():
        parts = [part.strip() for part in line.rstrip("\n").split("\t")]
        if len(parts) < 3:
            continue
        genes = _clean_genes(parts[2:])
        if genes:
            library[parts[0]] = genes
    if not library:
        raise ValueError(f"Enrichr returned an empty library for {gene_sets}.")
    return library


def _run_local_ora(clean_genes, clean_background, excluded, gene_sets):
    library = _fetch_enrichr_library(gene_sets)
    library_universe = set().union(*(set(genes) for genes in library.values()))
    background = set(clean_background).intersection(library_universe)
    query = set(clean_genes).intersection(background)
    if not query or len(background) < 2:
        empty = pd.DataFrame()
        empty.attrs["diagnostic"] = "No query genes were covered by the selected library/background."
        return empty

    population_size = len(background)
    query_size = len(query)
    tested = []
    for term, term_genes in library.items():
        term_background = background.intersection(term_genes)
        term_size = len(term_background)
        if term_size == 0:
            continue
        overlap = sorted(query.intersection(term_background))
        overlap_size = len(overlap)
        p_value = float(stats.hypergeom.sf(
            overlap_size - 1, population_size, term_size, query_size
        )) if overlap_size else 1.0
        odds_ratio = (
            (overlap_size + 0.5)
            * (population_size - term_size - query_size + overlap_size + 0.5)
            / ((query_size - overlap_size + 0.5) * (term_size - overlap_size + 0.5))
        )
        tested.append({
            "Term": term, "P-value": p_value, "Odds Ratio": odds_ratio,
            "Combined Score": np.nan, "Genes_list": overlap,
            "Overlap": f"{overlap_size}/{query_size}",
            "Effective_Background_Gene_Count": population_size,
        })
    if not tested:
        return pd.DataFrame()
    df = pd.DataFrame(tested)
    df["Adjusted P-value"] = fdrcorrection(df["P-value"].to_numpy(dtype=float))[1]
    df = df[df["Genes_list"].map(bool)].sort_values("Adjusted P-value").reset_index(drop=True)
    df.insert(0, "Rank", np.arange(1, len(df) + 1))
    return _finalize_ora_table(
        df, clean_genes, clean_background, excluded, gene_sets,
        "Local hypergeometric fallback using Enrichr gene-set library",
    )


@st.cache_data(show_spinner=False, ttl=86400)
def run_enrichr(gene_list, gene_sets, background_genes):
    clean_genes, clean_background, excluded = prepare_ora_inputs(
        gene_list, background_genes
    )
    if not clean_genes or not clean_background:
        empty = pd.DataFrame()
        empty.attrs["diagnostic"] = "The ORA query or measured background is empty."
        return empty
    speed_error = None
    try:
        result = _run_speedrichr(clean_genes, clean_background, excluded, gene_sets)
        if not result.empty:
            result.attrs["diagnostic"] = "Speedrichr query completed."
            return result
        if result.attrs.get("diagnostic"):
            return result
    except Exception as exc:
        speed_error = f"{type(exc).__name__}: {exc}"

    try:
        result = _run_local_ora(clean_genes, clean_background, excluded, gene_sets)
        result.attrs["diagnostic"] = (
            "Speedrichr was unavailable; explicit-background local ORA completed "
            f"from the Enrichr library. Speedrichr error: {speed_error}"
        )
        return result
    except Exception as exc:
        empty = pd.DataFrame()
        empty.attrs["diagnostic"] = (
            f"Speedrichr failed ({speed_error}); local fallback failed "
            f"({type(exc).__name__}: {exc})."
        )
        return empty


def run_directional_enrichr(genes_up, genes_down, database, background_genes):
    """Run ORA separately by expression direction and retain provenance."""
    frames = []
    for direction, genes in (("Up-regulated", genes_up), ("Down-regulated", genes_down)):
        if not _clean_genes(genes):
            continue
        result = run_enrichr(list(genes), database, background_genes)
        if isinstance(result, pd.DataFrame) and not result.empty:
            result = result.copy()
            result["Direction"] = direction
            frames.append(result)
    if frames:
        combined = pd.concat(frames, ignore_index=True)
        combined.attrs["diagnostic"] = " | ".join(
            frame.attrs.get("diagnostic", "completed") for frame in frames
        )
        return combined
    empty = pd.DataFrame()
    empty.attrs["diagnostic"] = "No directional ORA results were returned."
    return empty


def run_pathway_enrichment(genes_up, genes_down, combined_genes, strategy, background):
    libraries = ["KEGG_2021_Human", "GO_Biological_Process_2023"]

    def run_library(library):
        if strategy.startswith("Separate"):
            return run_directional_enrichr(genes_up, genes_down, library, background)
        return run_enrichr(combined_genes, library, background)

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(run_library, libraries))
    return dict(zip(libraries, results))

# ============================================================
# STRING
# ============================================================

@st.cache_data(show_spinner=False, ttl=86400)
def fetch_string_network(genes, confidence=400, max_nodes=150, network_type="functional"):
    genes = _clean_genes(genes)[:max_nodes]
    if len(genes) < 2:
        empty = pd.DataFrame()
        empty.attrs["diagnostic"] = "At least two unique genes are required for STRING."
        return empty

    try:
        session = _retry_session()
        base_url = "https://version-12-0.string-db.org/api/json"
        mapping_response = session.post(
            f"{base_url}/get_string_ids",
            data={
                "identifiers": "\r".join(genes), "species": 9606,
                "echo_query": 1, "caller_identity": "LumosNetworks",
            },
            timeout=(10, 45),
        )
        mapping_response.raise_for_status()
        mapping = pd.DataFrame(mapping_response.json())
        if mapping.empty or "stringId" not in mapping.columns:
            empty = pd.DataFrame()
            empty.attrs["diagnostic"] = "STRING could not map the submitted human identifiers."
            return empty
        mapped_ids = mapping["stringId"].dropna().astype(str).drop_duplicates().tolist()
        mapped_queries = set(mapping.get("queryItem", pd.Series(dtype=str)).astype(str).str.upper())
        unmapped = sorted(set(genes).difference(mapped_queries))
        if len(mapped_ids) < 2:
            empty = pd.DataFrame()
            empty.attrs["diagnostic"] = "Fewer than two identifiers mapped uniquely to STRING."
            return empty

        # STRING asks API clients to pause between calls; mapping is cached for reruns.
        time.sleep(1)
        params = {
            "identifiers": "\r".join(mapped_ids),
            "species": 9606,
            "required_score": confidence,
            "network_type": network_type,
            "caller_identity": "LumosNetworks",
        }
        response = session.post(
            f"{base_url}/network", data=params, timeout=(10, 60)
        )
        response.raise_for_status()
        result = pd.DataFrame(response.json())
        if not result.empty:
            result["STRING_Version"] = "12.0"
            result["STRING_Network_Type"] = network_type
            result["STRING_Required_Score"] = confidence
            result["STRING_Species"] = 9606
            result["STRING_Input_Gene_Count"] = len(genes)
            result["STRING_Mapped_Gene_Count"] = len(mapped_ids)
            result["STRING_Unmapped_Genes"] = ";".join(unmapped)
            result["Query_UTC"] = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
            result.attrs["diagnostic"] = (
                f"Mapped {len(mapped_ids)}/{len(genes)} genes to STRING 12.0."
            )
        else:
            result.attrs["diagnostic"] = (
                "STRING completed successfully but found no edges at the selected score threshold."
            )
        return result
    except Exception as exc:
        empty = pd.DataFrame()
        empty.attrs["diagnostic"] = f"STRING query failed: {type(exc).__name__}: {exc}"
        return empty


# ============================================================
# FUNCTIONAL MAPPING
# ============================================================

def build_functional_mapping(enrichment_frames, network_genes):
    mapping_data = []

    if not enrichment_frames:
        return pd.DataFrame()

    current_network_genes = set(_clean_genes(network_genes))

    for df_enrich in enrichment_frames:
        if not isinstance(df_enrich, pd.DataFrame) or df_enrich.empty:
            continue

        for _, row in df_enrich.iterrows():
            via_name = row.get("Term", "")
            genes_in_via = [
                g for g in str(row.get("Genes", "")).split(";")
                if g.upper() in current_network_genes
            ]

            if genes_in_via:
                mapping_data.append({
                    "Process / Biological Pathway": via_name,
                    "Related Genes in Network": ", ".join(
                        sorted(set(genes_in_via))
                    ),
                    "Count": len(set(genes_in_via)),
                    "Significance (FDR)": (
                        f"{float(row.get('Adjusted P-value', np.nan)):.2e}"
                        if pd.notna(row.get("Adjusted P-value", np.nan))
                        else "NA"
                    )
                })

    return pd.DataFrame(mapping_data)


# ============================================================
# LOAD DATA
# ============================================================

st.title("Arithmancy Pathway Profiler 🕸️🍲")

with st.sidebar:
    st.header("1. Data Entry")

    data_source = st.radio(
        "Data Source:",
        ["Import from DDEA", "Upload CSVs"]
    )

    df_full = pd.DataFrame()
    imported_from_meta = False
    meta_studies = []
    mt_gses = []
    contrast_ids = []
    contrast_name_map = {}

    if data_source == "Import from DDEA":
        
        # ----------------------------------------------------
        # MULTI-TISSUE HAS PRIORITY OVER META
        # ----------------------------------------------------
        if (
            "mt_results" in st.session_state
            and isinstance(st.session_state["mt_results"], pd.DataFrame)
            and not st.session_state["mt_results"].empty
        ):
            df_full = st.session_state["mt_results"].copy()
            mt_gses = st.session_state.get("mt_valid_gses", [])
            st.success("Multi-Tissue alignment results successfully imported from DDEA.")

        # ----------------------------------------------------
        # META-ANALYSIS 
        # ----------------------------------------------------
        elif (
            "meta_results" in st.session_state
            and isinstance(st.session_state["meta_results"], pd.DataFrame)
            and not st.session_state["meta_results"].empty
        ):
            df_full = st.session_state["meta_results"].copy()
            imported_from_meta = True
            meta_studies = st.session_state.get("meta_gse_ids", [])
            st.success("Meta-analysis results successfully imported from DDEA.")

        # ----------------------------------------------------
        # MULTIPLE CONTRASTS FROM ONE DATASET
        # ----------------------------------------------------
        elif (
            "single_contrast_results" in st.session_state
            and isinstance(st.session_state["single_contrast_results"], pd.DataFrame)
            and not st.session_state["single_contrast_results"].empty
        ):
            df_full = st.session_state["single_contrast_results"].copy()
            contrast_ids = st.session_state.get("single_contrast_ids", [])
            contrast_name_map = st.session_state.get("single_contrast_name_map", {})
            st.success("Multiple analyses from one dataset were successfully imported from DDEA.")

        # ----------------------------------------------------
        # FALLBACK: SINGLE DATASET
        # ----------------------------------------------------
        elif (
            ("res" in st.session_state and isinstance(st.session_state["res"], pd.DataFrame))
            or ("df_diff" in st.session_state and isinstance(st.session_state["df_diff"], pd.DataFrame))
        ):
            source = st.session_state.get("res", st.session_state.get("df_diff", pd.DataFrame()))
            df_full = source.copy()
            st.success("Complete single-dataset result table imported from DDEA (including the tested universe).")

        else:
            st.warning("No results found. Run the analysis in the DDEA module first.")

    else:
        uploaded_files = st.file_uploader(
            "Upload CSV result files",
            type=["csv"],
            accept_multiple_files=True
        )

        if uploaded_files:
            frames = []
            for uploaded in uploaded_files:
                temp = pd.read_csv(uploaded)
                temp["Source_File"] = uploaded.name
                frames.append(temp)

            if frames:
                df_full = pd.concat(frames, ignore_index=True, sort=False)

    # --------------------------------------------------------
    # MODE DETECTION
    # --------------------------------------------------------

    analysis_mode = _detect_analysis_mode(df_full)

    if analysis_mode == "meta":
        imported_from_meta = True
    elif analysis_mode == "multi_tissue" and not mt_gses:
        mt_gses = sorted({c.split("_", 1)[1] for c in df_full.columns if c.startswith("Log2FC_")})
    elif analysis_mode == "multi_contrast" and not contrast_ids:
        contrast_ids = sorted({
            c.split("_", 1)[1] for c in df_full.columns if c.startswith("Log2FC_")
        })
        contrast_name_map = {contrast_id: contrast_id for contrast_id in contrast_ids}

    if analysis_mode == "unknown" and not df_full.empty:
        st.error("The imported file does not contain a recognized effect-size column.")
        df_full = pd.DataFrame()
        analysis_mode = "none"

    # --------------------------------------------------------
    # SETTINGS
    # --------------------------------------------------------

    st.divider()
    st.header("2. Filters and Settings")

    fdr_cut = st.number_input(
        "FDR Cutoff:",
        min_value=0.0001,
        max_value=1.0,
        value=0.05,
        format="%.3f"
    )

    fc_cut = st.number_input(
        "Abs Log2FC Cutoff:",
        min_value=0.0,
        max_value=20.0,
        value=1.0,
        format="%.2f"
    )

    fdr_viz = st.slider("FDR Cut for Charts:", 0.01, 1.0, 0.05)

    st.divider()
    st.header("📐 Graph Dimensions")

    g_width = st.slider("Net width:", 600, 1800, 1100, 50)
    g_height = st.slider("Net height:", 400, 1200, 800, 50)


# ============================================================
# MAIN ANALYSIS
# ============================================================

if not df_full.empty:

    df_full = _standardize_result_dataframe(df_full, analysis_mode)
    required_columns = (
        ["Symbol"]
        if analysis_mode in {"multi_tissue", "multi_contrast"}
        else ["Symbol", "FDR", "Effect"]
    )
    missing_columns = [column for column in required_columns if column not in df_full.columns]
    if missing_columns:
        st.error("Missing required result columns: " + ", ".join(missing_columns))
        st.stop()

    analysis_universe = _selected_gene_symbols(df_full)

    # --------------------------------------------------------
    # MODE BANNER AND DYNAMIC MULTI-TISSUE ROUTER
    # --------------------------------------------------------

    if analysis_mode in {"multi_tissue", "multi_contrast"}:
        profile_ids = mt_gses if analysis_mode == "multi_tissue" else contrast_ids
        profile_names = (
            {profile_id: profile_id for profile_id in profile_ids}
            if analysis_mode == "multi_tissue"
            else contrast_name_map
        )
        integration_label = "Multi-Tissue Integration" if analysis_mode == "multi_tissue" else "Multi-Contrast Integration"
        consistent_label = (
            "Conserved genes (same direction in every dataset)"
            if analysis_mode == "multi_tissue"
            else "Consistent genes (same direction in every analysis)"
        )
        with st.sidebar:
            st.divider()
            st.header(f"2b. {integration_label}")
            if len(profile_ids) >= 2:
                focus = st.radio(
                    "Functional Focus:", 
                    [consistent_label] + [f"Profile: {profile_names.get(item, item)}" for item in profile_ids]
                )
            else:
                st.error("Integrated mode requires at least two analysis-specific effect columns.")
                st.stop()
        
        if focus == consistent_label:
            scope_word = "dataset" if analysis_mode == "multi_tissue" else "analysis"
            st.success(
                f"🧬 **Strict consistent-signature mode:** genes must pass both thresholds "
                f"and have the same direction in every {scope_word}."
            )
            try:
                eligible_columns = [
                    column
                    for profile_id in profile_ids
                    for column in (f"Log2FC_{profile_id}", f"FDR_{profile_id}")
                ]
                eligible_for_consistency = df_full.dropna(subset=eligible_columns)
                analysis_universe = _selected_gene_symbols(eligible_for_consistency)
                df_full, conserved_by_direction = conserved_signature(
                    df_full, profile_ids, fdr_cut, fc_cut
                )
            except ValueError as exc:
                st.error(str(exc))
                st.stop()
            if df_full.empty:
                st.warning("No genes satisfy the strict conserved-signature definition.")
                
        else:
            tg = next(
                item for item in profile_ids
                if focus == f"Profile: {profile_names.get(item, item)}"
            )
            profile_scope = "tissue dataset" if analysis_mode == "multi_tissue" else "saved analysis"
            st.info(f"📊 **Specific profile:** using {profile_scope} {profile_names.get(tg, tg)}.")
            df_full["Effect"] = pd.to_numeric(df_full[f"Log2FC_{tg}"], errors="coerce")
            df_full["FDR"] = pd.to_numeric(df_full[f"FDR_{tg}"], errors="coerce")
            df_full["Status"] = df_full[f"Status_{tg}"]
            analysis_universe = _selected_gene_symbols(
                df_full.dropna(subset=["Effect", "FDR"])
            )

        st.session_state["analysis_mode"] = analysis_mode
        state_table_key = "df_mt" if analysis_mode == "multi_tissue" else "df_multi_contrast"
        st.session_state[state_table_key] = df_full.copy()

    elif analysis_mode == "meta":
        st.success("🧬 **Meta-analysis mode:** APP is using the combined effect size (`Meta_Log2FC`) generated by DDEA.")
        if meta_studies: st.caption("Datasets included: " + ", ".join(map(str, meta_studies)))
        st.session_state["df_meta"] = df_full.copy()
        st.session_state["analysis_mode"] = "meta"
    else:
        st.info("📊 **Single-dataset mode:** APP is using the `Log2FC` generated by DDEA.")
        st.session_state["analysis_mode"] = "single"

    # --------------------------------------------------------
    # DATA SELECTION & SIGNIFICANCE FILTER
    # --------------------------------------------------------

    with st.sidebar:
        st.divider()
        st.header("3. Data Selection")

        all_symbols = df_full["Symbol"].dropna().astype(str).str.strip().str.upper().unique().tolist()
        selected_genes = st.multiselect("Select genes for functional analysis:", sorted(all_symbols), default=sorted(all_symbols))
        df_full = df_full[df_full["Symbol"].isin(selected_genes)].copy()

    df_sig = df_full[
        (df_full["FDR"] < fdr_cut) & (df_full["Effect"].abs() >= fc_cut)
    ].copy()

    raw_symbols = _selected_gene_symbols(df_sig)
    unique_symbols = sorted(set(raw_symbols))

    genes_up = set(df_sig.loc[df_sig["Effect"] > 0, "Symbol"].astype(str).str.upper())
    genes_down = set(df_sig.loc[df_sig["Effect"] < 0, "Symbol"].astype(str).str.upper())

    with st.sidebar:
        enrichment_strategy = st.radio(
            "ORA gene-list strategy:",
            ["Separate up/down (recommended)", "Combined DEGs"],
            help="Separate analysis preserves effect direction. ORA itself does not infer pathway activation.",
        )
        st.markdown("#### Explicit ORA background")
        background_source = st.radio(
            "Universe source",
            ["All genes tested in the result table", "Paste a custom measured-gene universe"],
            help="The universe must represent genes that could have been selected, not all genes in the genome.",
        )
        if background_source.startswith("Paste"):
            background_text = st.text_area(
                "Measured genes (one per line, comma or semicolon separated)", height=120
            )
            ora_background = _clean_genes([x for x in re.split(r"[\s,;]+", background_text) if x])
        else:
            ora_background = analysis_universe
        if not ora_background:
            st.error("An explicit ORA background is required.")
            st.stop()
        st.caption(f"Universe size: {len(ora_background):,} unique genes")

    # --------------------------------------------------------
    # SUMMARY METRICS
    # --------------------------------------------------------

    if analysis_mode == "meta":
        median_i2 = pd.to_numeric(df_full["I2_percent"], errors="coerce").median() if "I2_percent" in df_full.columns else np.nan
        metrics = [
            ("Datasets", len(meta_studies) if meta_studies else (int(df_full["N_Studies"].max()) if "N_Studies" in df_full.columns and not df_full.empty else "—")),
            ("Genes Evaluated", len(unique_symbols)),
            ("Meta-DEGs", len(df_sig)),
            ("Median I²", f"{median_i2:.1f}%" if pd.notna(median_i2) else "—"),
        ]
        cols = st.columns(len(metrics))
        for col, (label, val) in zip(cols, metrics):
            col.metric(label, val)

    elif analysis_mode in {"multi_tissue", "multi_contrast"}:
        first_label = "Datasets Aligned" if analysis_mode == "multi_tissue" else "Analyses Aligned"
        first_val = len(mt_gses) if analysis_mode == "multi_tissue" else len(contrast_ids)
        metrics = [
            (first_label, first_val),
            ("Union/Focus Genes Evaluated", len(unique_symbols)),
            ("DEGs in Focus", len(df_sig)),
        ]
        cols = st.columns(len(metrics))
        for col, (label, val) in zip(cols, metrics):
            col.metric(label, val)

    else:
        metrics = [
            ("Total DEGs Evaluated", len(raw_symbols)),
            ("Unique DEGs", len(unique_symbols)),
            ("Target Scope", len(df_sig)),
        ]
        cols = st.columns(len(metrics))
        for col, (label, val) in zip(cols, metrics):
            col.metric(label, val)

    # --------------------------------------------------------
    # TABS
    # --------------------------------------------------------

    report = {}
    tab1, tab2, tab3, tab4 = st.tabs([
        "📊 Pathway (KEGG/GO)",
        "🕸️ STRING Interactive & Synergy",
        "👑 Candidate Regulators",
        "🔬 Cross-Evaluation Evidence"
    ])

    # ========================================================
    # ENRICHMENT
    # ========================================================

    with st.spinner("Running functional enrichment..."):
        pathway_results = run_pathway_enrichment(
            genes_up, genes_down, unique_symbols,
            enrichment_strategy, ora_background,
        )
        k_res = pathway_results["KEGG_2021_Human"]
        g_res = pathway_results["GO_Biological_Process_2023"]

    st.session_state["k_res"] = k_res
    st.session_state["g_res"] = g_res

    # ========================================================
    # TAB 1 — PATHWAYS
    # ========================================================

    with tab1:

        title_suffix = ""
        if analysis_mode == "meta": title_suffix = " — Meta-analysis genes"
        if analysis_mode == "multi_tissue": title_suffix = f" — {focus}"
        if analysis_mode == "multi_contrast": title_suffix = f" — {focus}"
        st.subheader("Functional Enrichment" + title_suffix)
        st.caption(
            "Explicit-background over-representation analysis (ORA). Speedrichr is used first; "
            "if unavailable, Lumos downloads the same Enrichr gene-set library and applies a "
            "one-sided hypergeometric test with Benjamini–Hochberg correction locally. "
            "Adjusted p-values are library/query-specific; "
            f"an explicit measured-gene background of {len(ora_background):,} genes is used. "
            "Direction labels describe the input list and do not prove pathway activation."
        )

        c1, c2 = st.columns(2)

        if isinstance(k_res, pd.DataFrame) and not k_res.empty:
            k_plot = k_res[pd.to_numeric(k_res["Adjusted P-value"], errors="coerce") <= fdr_cut].head(15)
            with c1:
                if k_plot.empty:
                    st.info("No KEGG term passed the selected FDR threshold.")
                else:
                    fig_kegg = px.bar(
                        k_plot,
                        x="-log10(FDR)",
                        y="Term",
                        orientation="h",
                        title="KEGG",
                        color="Adjusted P-value",
                        color_continuous_scale="Reds_r",
                        hover_data=[c for c in ["Direction", "Overlap"] if c in k_plot.columns],
                    )
                    report["kegg"] = fig_kegg
                    st.plotly_chart(fig_kegg, width="stretch")
        else:
            c1.info("No KEGG terms were returned.")
            diagnostic = getattr(k_res, "attrs", {}).get("diagnostic")
            if diagnostic: c1.caption(diagnostic)

        if isinstance(g_res, pd.DataFrame) and not g_res.empty:
            g_plot = g_res[pd.to_numeric(g_res["Adjusted P-value"], errors="coerce") <= fdr_cut].head(15).copy()
            g_plot["Count"] = g_plot["Overlap"].apply(lambda x: int(str(x).split("/")[0]) if "/" in str(x) else 10)

            with c2:
                if g_plot.empty:
                    st.info("No GO Biological Process term passed the selected FDR threshold.")
                else:
                    fig_go = px.scatter(
                        g_plot,
                        x="-log10(FDR)",
                        y="Term",
                        size="Count",
                        title="GO Biological Process",
                        color="Adjusted P-value",
                        color_continuous_scale="Blues_r",
                        hover_data=[c for c in ["Direction", "Overlap"] if c in g_plot.columns],
                    )
                    report["go"] = fig_go
                    st.plotly_chart(fig_go, width="stretch")
        else:
            c2.info("No GO Biological Process terms were returned.")
            diagnostic = getattr(g_res, "attrs", {}).get("diagnostic")
            if diagnostic: c2.caption(diagnostic)

        if isinstance(k_res, pd.DataFrame):
            st.download_button("📥 Download KEGG enrichment", k_res.to_csv(index=False).encode("utf-8"), "Lumos_KEGG_enrichment.csv", mime="text/csv")
        if isinstance(g_res, pd.DataFrame):
            st.download_button("📥 Download GO enrichment", g_res.to_csv(index=False).encode("utf-8"), "Lumos_GO_enrichment.csv", mime="text/csv")
        ora_engines = sorted({
            str(engine)
            for result in (k_res, g_res)
            if isinstance(result, pd.DataFrame) and "ORA_Engine" in result.columns
            for engine in result["ORA_Engine"].dropna().unique()
        })
        app_manifest = {
            "analysis_mode": analysis_mode,
            "fdr_threshold": fdr_cut,
            "absolute_log2fc_threshold": fc_cut,
            "enrichment_strategy": enrichment_strategy,
            "input_genes": unique_symbols,
            "up_regulated_genes": sorted(genes_up),
            "down_regulated_genes": sorted(genes_down),
            "gene_set_libraries": ["KEGG_2021_Human", "GO_Biological_Process_2023"],
            "ora_engines": ora_engines,
            "kegg_diagnostic": getattr(k_res, "attrs", {}).get("diagnostic"),
            "go_diagnostic": getattr(g_res, "attrs", {}).get("diagnostic"),
            "background_source": background_source,
            "background_gene_count": len(ora_background),
            "background_genes": ora_background,
        }
        st.download_button(
            "📥 Download APP query manifest",
            json.dumps(app_manifest, indent=2),
            "Lumos_APP_query_manifest.json",
            mime="application/json",
        )

    # ========================================================
    # TAB 2 — STRING
    # ========================================================

    with tab2:

        title_suffix = ""
        if analysis_mode == "meta": title_suffix = " — Meta-analysis"
        if analysis_mode == "multi_tissue": title_suffix = " — Aligned Tissues"
        if analysis_mode == "multi_contrast": title_suffix = " — Aligned Analyses"
        st.subheader("STRING Association Network" + title_suffix)
        st.caption("Edges are STRING associations under the selected network type and score threshold; they are not inferred from this expression dataset.")

        c_p, c_m = st.columns([1, 4])

        conf = c_p.selectbox("STRING Confidence:", [150, 400, 700, 900], index=1)
        string_network_type = c_p.selectbox("Network type:", ["functional", "physical"], index=0)
        max_n = c_p.number_input("Maximum Number of Genes:", min_value=2, max_value=500, value=150)

        if c_p.button("🌐 Create Master Network"):
            top_genes_list = (
                df_sig.sort_values(by="Effect", key=lambda s: s.abs(), ascending=False)
                .head(int(max_n))["Symbol"].tolist()
            )

            string_df = fetch_string_network(top_genes_list, conf, int(max_n), string_network_type)

            if not string_df.empty:
                st.session_state["string_df"] = string_df
                diagnostic = string_df.attrs.get("diagnostic")
                if diagnostic:
                    c_m.caption(diagnostic)
                nodes, edges, nodes_added = [], [], set()
                top_vias_list = []

                if isinstance(k_res, pd.DataFrame): top_vias_list.extend(k_res.head(12)["Term"].tolist())
                if isinstance(g_res, pd.DataFrame): top_vias_list.extend(g_res.head(13)["Term"].tolist())

                for _, row in string_df.iterrows():
                    g1 = str(row["preferredName_A"]).upper()
                    g2 = str(row["preferredName_B"]).upper()
                    shared_paths = []

                    if isinstance(k_res, pd.DataFrame):
                        k_m = k_res[k_res["Term"].isin(top_vias_list) & k_res["Genes"].apply(lambda x: g1 in str(x).upper() and g2 in str(x).upper())]
                        shared_paths.extend(k_m["Term"].tolist())
                    if isinstance(g_res, pd.DataFrame):
                        g_m = g_res[g_res["Term"].isin(top_vias_list) & g_res["Genes"].apply(lambda x: g1 in str(x).upper() and g2 in str(x).upper())]
                        shared_paths.extend(g_m["Term"].tolist())

                    score = pd.to_numeric(row.get("score"), errors="coerce")
                    tooltip_text = "Association"
                    if pd.notna(score): tooltip_text += f" | STRING score: {score:.3f}"
                    if shared_paths: tooltip_text += " | Shared terms: " + " | ".join(sorted(set(shared_paths)))

                    for g in [g1, g2]:
                        if g not in nodes_added:
                            if g in genes_up: node_color = "#FF4B4B"
                            elif g in genes_down: node_color = "#1C83E1"
                            else: node_color = "#D3D3D3"

                            gene_info = df_full[df_full["Symbol"] == g]
                            title = g

                            if not gene_info.empty:
                                eff = gene_info.iloc[0]["Effect"]
                                fdr = gene_info.iloc[0]["FDR"]
                                if pd.notna(eff): title += f" | Effect: {eff:.2f}"
                                if pd.notna(fdr): title += f" | FDR: {fdr:.2e}"

                                if analysis_mode == "meta" and "I2_percent" in gene_info.columns:
                                    i2 = gene_info.iloc[0]["I2_percent"]
                                    if pd.notna(i2): title += f" | I²: {i2:.1f}%"
                                elif analysis_mode == "multi_tissue":
                                    for gse in mt_gses:
                                        eff_t = gene_info.iloc[0].get(f"Log2FC_{gse}", np.nan)
                                        if pd.notna(eff_t): title += f" | {gse} FC: {eff_t:.2f}"
                                elif analysis_mode == "multi_contrast":
                                    for contrast_id in contrast_ids:
                                        eff_c = gene_info.iloc[0].get(f"Log2FC_{contrast_id}", np.nan)
                                        if pd.notna(eff_c):
                                            title += f" | {contrast_name_map.get(contrast_id, contrast_id)} FC: {eff_c:.2f}"

                            nodes.append(Node(id=g, label=g, size=15, color=node_color, title=title))
                            nodes_added.add(g)

                    edges.append(Edge(source=g1, target=g2, title=tooltip_text))

                config = Config(width=g_width, height=g_height, directed=False, physics=False, hierarchical=False, stabilization=True)
                with c_m: agraph(nodes=nodes, edges=edges, config=config)
                st.session_state['app_network'] = (nodes, edges)
                st.divider()
                st.subheader("📋 Functional Mapping: Network Processes vs. Genes")
                mapping_data = build_functional_mapping([k_res, g_res], list(nodes_added))

                if not mapping_data.empty: st.dataframe(mapping_data, width="stretch", hide_index=True)
                else: st.info("No functional mapping was detected for the current network.")
                st.download_button(
                    "📥 Download STRING associations",
                    string_df.to_csv(index=False).encode("utf-8"),
                    "Lumos_STRING_associations.csv", mime="text/csv",
                )
            else:
                st.info("No STRING interactions detected.")
                diagnostic = string_df.attrs.get("diagnostic")
                if diagnostic:
                    st.caption(diagnostic)

    # ========================================================
    # TAB 3 — MASTER REGULATORS
    # ========================================================

    with tab3:
        st.subheader("Candidate Transcriptional Regulators (JASPAR/TRRUST)")
        st.caption("Candidates are identified by target-set over-representation. This is prior-knowledge enrichment, not causal inference or motif scanning of the uploaded samples.")
        st.session_state["tf_regulators"] = pd.DataFrame()
        tf_frames = []
        for db in ["JASPAR_2024", "TRRUST_Transcription_Factors_2019"]:
            if enrichment_strategy.startswith("Separate"):
                current = run_directional_enrichr(genes_up, genes_down, db, ora_background)
            else:
                current = run_enrichr(unique_symbols, db, ora_background)
            if isinstance(current, pd.DataFrame) and not current.empty:
                current = current.copy()
                current["Regulatory_Database"] = db
                tf_frames.append(current)

        res_tf = pd.concat(tf_frames, ignore_index=True) if tf_frames else pd.DataFrame()

        if isinstance(res_tf, pd.DataFrame) and not res_tf.empty:
            res_tf = res_tf.copy()
            res_tf["TF_Symbol"] = res_tf["Term"].apply(
                lambda x: str(x).split(" (")[0].strip().upper()
            )
            df_final = res_tf[res_tf["Adjusted P-value"] <= fdr_viz].copy()

            if df_final.empty:
                st.info("No candidate regulator passed the selected visualization FDR.")

            fig_tf = px.bar(df_final.head(20), x="-log10(FDR)", y="TF_Symbol", orientation="h", color="Adjusted P-value", color_continuous_scale="Viridis_r", title="Top candidate regulators", hover_data=[c for c in ["Regulatory_Database", "Direction", "Overlap"] if c in df_final.columns])
            st.plotly_chart(fig_tf, width="stretch")

            export_columns = ["TF_Symbol", "Term", "Genes", "Adjusted P-value", "Overlap", "Regulatory_Database", "Direction"]
            export_df = df_final[[c for c in export_columns if c in df_final.columns]].copy()
            st.session_state["tf_regulators"] = export_df
            report["regulators"] = (fig_tf, export_df)
            st.dataframe(export_df, width="stretch", hide_index=True)
            st.download_button("📥 Download CSV for PG", export_df.to_csv(index=False).encode("utf-8"), "lumos_regulation_table.csv", mime="text/csv")
        else:
            st.info("No candidate-regulator enrichment was returned for the current gene lists.")

    # ========================================================
    # TAB 4 — CROSS-EVALUATION EVIDENCE
    # ========================================================

    with tab4:
        if analysis_mode in {"multi_tissue", "multi_contrast"}:
            profile_ids = mt_gses if analysis_mode == "multi_tissue" else contrast_ids
            scope_title = "Cross-Tissue Evidence" if analysis_mode == "multi_tissue" else "Evidence Across Analyses"
            scope_caption = (
                "These results reflect the alignment between tissues imported from the DDEA."
                if analysis_mode == "multi_tissue"
                else "These are aligned contrasts from one dataset; effects remain separate and are not meta-analyzed."
            )
            export_name = "Lumos_APP_MultiTissue.csv" if analysis_mode == "multi_tissue" else "Lumos_APP_MultiContrast.csv"
            st.subheader(f"🔬 {scope_title}")
            st.caption(scope_caption)
            mt_cols = [
                c for c in df_full.columns
                if c.startswith("Log2FC_") or c.startswith("FDR_")
                or c.startswith("Status_") or c == "Symbol"
                or c in {"Cross_Tissue_Pattern", "Cross_Contrast_Pattern", "N_Significant_Datasets", "N_Significant_Contrasts"}
            ]
            report["evidence"] = ([], df_full[mt_cols], "Symbol")
            st.dataframe(df_full[mt_cols].sort_values("Symbol"), width="stretch", hide_index=True)
            st.download_button(
                f"📥 Download APP {scope_title} Table",
                df_full[mt_cols].to_csv(index=False).encode("utf-8"),
                export_name, mime="text/csv",
            )

        elif analysis_mode == "meta":
            st.subheader("🔬 Evidence Across Datasets")
            st.caption("These results are taken directly from the DDEA meta-analysis output.")
            
            if meta_studies:
                st.markdown("### Datasets included")
                study_cols = st.columns(min(len(meta_studies), 4))
                for i, study in enumerate(meta_studies): study_cols[i % len(study_cols)].metric(f"Dataset {i + 1}", str(study))

            if "I2_percent" in df_full.columns:
                st.markdown("### Heterogeneity")
                i2_values = pd.to_numeric(df_full["I2_percent"], errors="coerce").dropna()
                if not i2_values.empty:
                    hc1, hc2, hc3 = st.columns(3)
                    hc1.metric("Median I²", f"{i2_values.median():.1f}%")
                    hc2.metric("Mean I²", f"{i2_values.mean():.1f}%")
                    hc3.metric("Genes with I² > 50%", f"{(i2_values > 50).sum()}")

            st.markdown("### Meta-analysis results")
            preferred_cols = ["Symbol", "Meta_Log2FC", "Meta_SE", "CI95_Low", "CI95_High", "PValue", "FDR", "I2_percent", "Tau2", "N_Studies", "Studies", "Status"]
            display_cols = [c for c in preferred_cols if c in df_full.columns]
            st.dataframe(df_full[display_cols].sort_values("FDR", na_position="last"), width="stretch", hide_index=True)
            st.download_button("📥 Download APP Meta-analysis Table", df_full.to_csv(index=False).encode("utf-8"), "Lumos_APP_Meta_Analysis.csv", mime="text/csv")

            meta_plot = df_full.copy()
            meta_plot["-log10(FDR)"] = -np.log10(meta_plot["FDR"].fillna(1).clip(lower=1e-300))
            meta_plot["Status"] = "Not Significant"
            meta_plot.loc[(meta_plot["FDR"] < fdr_cut) & (meta_plot["Meta_Log2FC"] >= fc_cut), "Status"] = "Up-regulated"
            meta_plot.loc[(meta_plot["FDR"] < fdr_cut) & (meta_plot["Meta_Log2FC"] <= -fc_cut), "Status"] = "Down-regulated"

            st.markdown("### Meta-analysis effect landscape")
            meta_volcano = px.scatter(
                meta_plot, x="Meta_Log2FC", y="-log10(FDR)", color="Status", hover_name="Symbol", hover_data=[c for c in ["N_Studies", "I2_percent", "Studies"] if c in meta_plot.columns], title="Meta-analysis Volcano Plot"
            )
            st.plotly_chart(meta_volcano, width="stretch")
            report["evidence"] = ([meta_volcano], df_full, "FDR")
            top_meta = df_full.sort_values("FDR", na_position="last").head(20).copy()
            if not top_meta.empty:
                st.markdown("### Top meta-analysis genes")
                st.dataframe(top_meta[display_cols], width="stretch", hide_index=True)

        else:
            st.info("This tab becomes active when integrated results from DDEA are imported.")

    # --------------------------------------------------------
    # HTML REPORT
    # --------------------------------------------------------

    if report:
        st.sidebar.header("Export Report")
        if st.sidebar.button("🌐 Generate HTML Report"):
            fmts = {
                "Adjusted P-value": "exp2",
                "P-value": "exp2",
                "PValue": "exp2",
                "FDR": "exp2",
                "Significance (FDR)": "exp2",
                "Combined Score": "fixed3",
                "Odds Ratio": "fixed3",
                "Meta_Log2FC": "fixed3",
                "Meta_SE": "fixed3",
                "CI95_Low": "fixed3",
                "CI95_High": "fixed3",
                "I2_percent": "fixed3",
                "Tau2": "fixed3",
                "Count": "int",
                "N_Studies": "int",
            }
            fmt = lambda c: fmts.get(c) or ("fixed3" if c.startswith("Log2FC_") else "exp2" if c.startswith("FDR_") else "text")

            def spec(df, sort, cols=None, **kw):
                if not isinstance(df, pd.DataFrame) or df.empty:
                    return None
                return table_spec(df, [(c, c, fmt(c)) for c in (cols or df.columns) if c in df.columns], sort, **kw)

            sdf = st.session_state.get("string_df")
            has_net = isinstance(sdf, pd.DataFrame) and {"preferredName_A", "preferredName_B"} <= set(sdf.columns)
            mapping = build_functional_mapping([k_res, g_res], pd.concat([sdf["preferredName_A"], sdf["preferredName_B"]]).tolist()) if has_net else pd.DataFrame()
            if not mapping.empty:
                mapping["Significance (FDR)"] = pd.to_numeric(mapping["Significance (FDR)"], errors="coerce")
            tf_fig, tf_df = report.get("regulators", (None, None))
            ev_figs, ev_df, ev_sort = report.get("evidence", ([], None, None))
            status_kw = {"filter_col": "Status", "badge_col": "Status", "badge_map": {"Up-regulated": "up", "Down-regulated": "down"}}
            enr_cols = ["Term", "Direction", "Overlap", "P-value", "Adjusted P-value", "Combined Score", "Genes", "ORA_Engine"]
            dir_kw = {"filter_col": "Direction"}
            sections = [
                {"id": "kegg", "title": "KEGG Enrichment", "desc": "Over-representation analysis with explicit background (KEGG_2021_Human; Speedrichr or local fallback).", "figs": [report.get("kegg")], "table": spec(k_res, "Adjusted P-value", enr_cols, **dir_kw)},
                {"id": "go", "title": "GO Enrichment", "desc": "Over-representation analysis with explicit background (GO_Biological_Process_2023; Speedrichr or local fallback).", "figs": [report.get("go")], "table": spec(g_res, "Adjusted P-value", enr_cols, **dir_kw)},
                {"id": "network", "title": "Interaction Network", "desc": f"Last STRING network built in this session ({len(sdf) if has_net else 0} interactions): the same nodes and edges as the page graph, drawn with a force-directed layout; hover a node for its details. Functional mapping of the network genes below.",
                 "figs": [network_figure(*st.session_state['app_network']) if 'app_network' in st.session_state else None], "table": spec(mapping, "Significance (FDR)")},
                {"id": "regulators", "title": "Master Regulators", "desc": f"Candidate regulators by target-set over-representation in JASPAR_2024 and TRRUST_Transcription_Factors_2019 (adjusted p ≤ {fdr_viz}). Prior-knowledge enrichment, not causal inference.", "figs": [tf_fig], "table": spec(tf_df, "Adjusted P-value", filter_col="Regulatory_Database")},
                {"id": "evidence", "title": "Cross-Evaluation Evidence", "desc": "Per-dataset evidence imported from DDEA (multi-tissue, multi-contrast or meta-analysis).", "figs": ev_figs,
                 "table": spec(ev_df, ev_sort, **(status_kw if analysis_mode == "meta" else {}))},
            ]
            mode_names = {
                "single": "Single dataset",
                "meta": "Meta-analysis",
                "multi_tissue": "Multi-tissue",
                "multi_contrast": "Multi-contrast",
            }
            info = [
                ("Data source", data_source),
                ("Analysis mode", mode_names.get(analysis_mode, analysis_mode)),
            ]
            if analysis_mode in {"multi_tissue", "multi_contrast"}:
                info.append(("Functional focus", focus))

            if analysis_mode == "meta" and meta_studies:
                info.append(("Datasets", ", ".join(map(str, meta_studies))))
            elif analysis_mode == "multi_tissue" and mt_gses:
                info.append(("Datasets", ", ".join(map(str, mt_gses))))
            elif analysis_mode == "multi_contrast" and contrast_ids:
                info.append(("Datasets", ", ".join(map(str, contrast_ids))))

            info += [
                ("Input genes", len(unique_symbols)),
                ("Significant genes", f"{len(df_sig)} ({len(genes_up)} up, {len(genes_down)} down)"),
                ("Enrichment strategy", enrichment_strategy),
                ("ORA background", f"{background_source} ({len(ora_background):,} genes)"),
                ("Databases", "Enrichr libraries KEGG_2021_Human, GO_Biological_Process_2023, JASPAR_2024, TRRUST_Transcription_Factors_2019 (Speedrichr or local ORA); STRING v12 (Homo sapiens)"),
                ("Thresholds", f"FDR < {fdr_cut}, |Log2FC| ≥ {fc_cut}; regulators adjusted p ≤ {fdr_viz}"),
            ]
            footer = ("Arithmancy Pathway Profiler (APP) &mdash; EvoMol-Lab, BioME, UFRN<br>"
                      "Enrichment: Enrichr (Kuleshov et al., 2016). Protein interactions: STRING v12 (Szklarczyk et al., 2023).")
            with st.spinner("Generating HTML report..."):
                html_str = build_html_report(
                    "APP Enrichment Report",
                    logo_path,
                    info,
                    [(l, v, "") for l, v in metrics],
                    sections,
                    footer,
                    "APP",
                    stats_desc="Analysis mode, gene counts, databases and thresholds used in this enrichment.",
                )
                st.sidebar.download_button(
                    "📥 Download HTML Report",
                    data=html_str.encode("utf-8"),
                    file_name=f"APP_Report_{datetime.now().strftime('%Y%m%d')}.html",
                    mime="text/html",
                    on_click="ignore",
                )

else:
    st.info("Waiting for data entry. Run a Single Dataset, Multi-Dataset Meta-analysis, or Multi-Tissue in DDEA first.")

# ============================================================
# CONNECTION NOTES FOR PG
# ============================================================
st.session_state["k_res"] = st.session_state.get("k_res", pd.DataFrame())
st.session_state["g_res"] = st.session_state.get("g_res", pd.DataFrame())
