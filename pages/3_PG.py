import os
from datetime import datetime
import streamlit as st
import pandas as pd
import numpy as np
from streamlit_agraph import agraph, Node, Edge, Config

from lumos_core import network_sensitivity_analysis, random_walk_with_restart
from lumos_report import build_html_report, table_spec, network_figure

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    nx = None
    NETWORKX_AVAILABLE = False


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="Lumos Networks | Networks",
    page_icon="🧬",
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

LOGO_PATH = os.path.join(PARENT_DIR, "assets", "Lumos Networks.png")
PG_LOGO_PATH = os.path.join(PARENT_DIR, "assets", "PG.png")


# ============================================================
# SIDEBAR NAVIGATION
# ============================================================

with st.sidebar:
    if os.path.exists(LOGO_PATH):
        st.image(LOGO_PATH, width=250)

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
    if os.path.exists(PG_LOGO_PATH):
        st.image(PG_LOGO_PATH, width="stretch")

    st.info("You are in the networking module.")


# ============================================================
# HELPERS
# ============================================================

def clean_gene(value):
    if pd.isna(value): return ""
    return str(value).strip().upper()

def get_evidence_table():
    candidates = ["df_multi_contrast", "df_mt", "df_meta", "meta_results", "df_diff"]
    for key in candidates:
        obj = st.session_state.get(key)
        if isinstance(obj, pd.DataFrame) and not obj.empty:
            df = obj.copy()
            if "Symbol" in df.columns:
                df["Symbol"] = df["Symbol"].astype(str).str.strip().str.upper()
            return df, key
    return pd.DataFrame(), None

def get_effect_column(df):
    for col in ["Effect", "Meta_Log2FC", "Log2FC"]:
        if col in df.columns: return col
    return None

def get_fdr_column(df):
    for col in ["FDR", "Adjusted P-value", "padj", "adj.P.Val"]:
        if col in df.columns: return col
    return None

def prepare_regulatory_edges():
    raw = st.session_state.get("tf_regulators")
    if not isinstance(raw, pd.DataFrame) or raw.empty: return pd.DataFrame()
    df = raw.copy()
    tf_col = next((c for c in ["TF_Symbol", "TF_Name", "TF"] if c in df.columns), None)
    target_col = next((c for c in ["Genes", "Target_Gene", "Target"] if c in df.columns), None)
    if not tf_col or not target_col: return pd.DataFrame()
    df = df.assign(Target=df[target_col].astype(str).str.split(";")).explode("Target")
    df["TF"] = df[tf_col].astype(str).map(clean_gene)
    df["Target"] = df["Target"].astype(str).map(clean_gene)
    df = df[(df["TF"] != "") & (df["Target"] != "") & (df["TF"] != "NAN") & (df["Target"] != "NAN")]
    provenance = [
        c for c in ["Regulatory_Database", "Direction", "Adjusted P-value", "Term", "Overlap"]
        if c in df.columns
    ]
    return df[["TF", "Target"] + provenance].drop_duplicates().reset_index(drop=True)

def prepare_enrichment():
    k_res = st.session_state.get("k_res")
    g_res = st.session_state.get("g_res")
    frames = []
    if isinstance(k_res, pd.DataFrame) and not k_res.empty: frames.append(k_res.copy())
    if isinstance(g_res, pd.DataFrame) and not g_res.empty: frames.append(g_res.copy())
    if not frames: return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)

def enrich_gene_set(df_enrichment, term, genes):
    if df_enrichment.empty: return set()
    rows = df_enrichment[df_enrichment["Term"].astype(str).str.upper() == str(term).upper()]
    if rows.empty: return set()
    pathway_genes = set()
    for value in rows["Genes"]:
        pathway_genes.update(clean_gene(g) for g in str(value).split(";") if clean_gene(g))
    return pathway_genes.intersection(set(clean_gene(g) for g in genes))

@st.cache_data(show_spinner=False, max_entries=32)
def build_effect_map(df_evidence, mode):
    if df_evidence.empty or "Symbol" not in df_evidence.columns: return {}
    effect_col = get_effect_column(df_evidence)
    fdr_col = get_fdr_column(df_evidence)
    effect_map = {}
    for _, row in df_evidence.iterrows():
        gene = clean_gene(row.get("Symbol"))
        if not gene: continue
        effect = pd.to_numeric(row.get(effect_col), errors="coerce") if effect_col else np.nan
        fdr = pd.to_numeric(row.get(fdr_col), errors="coerce") if fdr_col else np.nan
        i2 = np.nan
        for col in ["I2_percent", "I2", "I_squared"]:
            if col in df_evidence.columns:
                i2 = pd.to_numeric(row.get(col), errors="coerce")
                break
        info = {"effect": effect, "fdr": fdr, "i2": i2}
        if mode in ["multi_tissue", "multi_contrast"]:
            for c in df_evidence.columns:
                if c.startswith("Log2FC_") or c.startswith("FDR_") or c.startswith("Status_"):
                    info[c] = row.get(c)
        effect_map[gene] = info
    return effect_map

def node_gene_color(gene, effect_map, fdr_threshold=0.05):
    info = effect_map.get(gene, {})
    effect = info.get("effect", np.nan)
    fdr = info.get("fdr", np.nan)
    if pd.notna(effect) and (pd.isna(fdr) or fdr <= fdr_threshold):
        if effect > 0: return "#FF4B4B"
        if effect < 0: return "#1C83E1"
    return "#D3D3D3"

def node_gene_title(gene, effect_map, mode):
    info = effect_map.get(gene, {})
    effect = info.get("effect", np.nan)
    fdr = info.get("fdr", np.nan)
    i2 = info.get("i2", np.nan)
    parts = [gene]
    if mode in ["multi_tissue", "multi_contrast"]:
        if pd.notna(effect): parts.append(f"Focus Effect: {effect:.2f}")
        if pd.notna(fdr): parts.append(f"Focus FDR: {fdr:.2e}")
        for k, v in info.items():
            if k.startswith("Log2FC_") and pd.notna(v) and isinstance(v, (int, float)):
                gse = k.replace("Log2FC_", "")
                parts.append(f"{gse} FC: {v:.2f}")
    else:
        if pd.notna(effect): parts.append(f"Effect: {effect:.2f}")
        if pd.notna(fdr): parts.append(f"FDR: {fdr:.2e}")
        if pd.notna(i2): parts.append(f"I²: {i2:.1f}%")
    return " | ".join(parts)


cached_network_sensitivity = st.cache_data(show_spinner=False, max_entries=16)(network_sensitivity_analysis)


# --- SIGNAL PROPAGATION ---
@st.cache_data(show_spinner=False, max_entries=32)
def calculate_signal_propagation(
    df_reg, all_network_genes, df_evidence, mode, restart_probability=0.7, propagation_decay=0.5
):
    if df_reg.empty or not all_network_genes:
        return pd.DataFrame()
    effect_map = build_effect_map(df_evidence, mode)
        
    nodes = sorted(all_network_genes.union(set(df_reg['TF'])))
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    n = len(nodes)
    
    # Adjacency matrix (Directed: TF -> Target)
    A = np.zeros((n, n))
    for _, row in df_reg.iterrows():
        tf, target = row['TF'], row['Target']
        if tf in node_to_idx and target in node_to_idx:
            A[node_to_idx[target], node_to_idx[tf]] = 1.0 
            
    # RWR
    col_sums = A.sum(axis=0)
    W = np.zeros_like(A)
    non_dangling = col_sums > 0
    W[:, non_dangling] = A[:, non_dangling] / col_sums[non_dangling]
    
    seeds = np.zeros(n)
    for i, node in enumerate(nodes):
        eff = effect_map.get(node, {}).get('effect', 0)
        if pd.notna(eff):
            seeds[i] = abs(eff)
    
    if seeds.sum() > 0:
        seeds = seeds / seeds.sum()
    else:
        seeds = np.ones(n) / n

    p, iterations, converged = random_walk_with_restart(
        A, seeds, restart_probability=restart_probability
    )

    # Linear, direction-aware propagation. This is intentionally not labelled
    # SPIA because the network has no signed activation/inhibition coefficients.
    E = np.zeros(n)
    for i, node in enumerate(nodes):
        eff = effect_map.get(node, {}).get('effect', 0)
        E[i] = eff if pd.notna(eff) else 0

    I = np.eye(n)
    try:
        PF = np.linalg.solve(I - propagation_decay * W, E)
    except np.linalg.LinAlgError:
        PF = E 
        
    return pd.DataFrame({
        'Gene': nodes,
        'Original_Log2FC': E,
        'Personalized_RWR_Score': p,
        'Linear_Network_Propagation': PF,
        'RWR_Converged': converged,
        'RWR_Iterations': iterations,
        'Restart_Probability': restart_probability,
        'Propagation_Decay': propagation_decay,
    }).sort_values('Personalized_RWR_Score', ascending=False)


@st.cache_data(show_spinner=False, max_entries=32)
def calculate_standard_topology(df_reg, all_network_genes):
    """Calculate complementary graph metrics for the displayed TF-target layer."""
    if not NETWORKX_AVAILABLE or df_reg.empty or not all_network_genes:
        return pd.DataFrame()

    relevant_edges = (
        df_reg[df_reg["Target"].isin(all_network_genes)][["TF", "Target"]]
        .dropna()
        .drop_duplicates()
    )
    if relevant_edges.empty:
        return pd.DataFrame()

    graph = nx.DiGraph()
    graph.add_edges_from(relevant_edges.itertuples(index=False, name=None))

    pagerank = nx.pagerank(graph, alpha=0.85)
    betweenness = nx.betweenness_centrality(graph.to_undirected(), normalized=True)

    try:
        hub_scores, authority_scores = nx.hits(
            graph,
            max_iter=1000,
            tol=1.0e-08,
            normalized=True,
        )
    except (nx.PowerIterationFailedConvergence, nx.NetworkXException):
        hub_scores = {node: np.nan for node in graph.nodes}
        authority_scores = {node: np.nan for node in graph.nodes}

    tf_nodes = set(relevant_edges["TF"])
    target_nodes = set(relevant_edges["Target"])
    rows = []
    for node in graph.nodes:
        if node in tf_nodes and node in target_nodes:
            node_type = "TF and Target"
        elif node in tf_nodes:
            node_type = "Transcription Factor"
        else:
            node_type = "Target Gene"

        rows.append({
            "Node": node,
            "Node Type": node_type,
            "Out-Degree Count": int(graph.out_degree(node)),
            "In-Degree Count": int(graph.in_degree(node)),
            "HITS Hub Score": float(hub_scores.get(node, np.nan)),
            "HITS Authority Score": float(authority_scores.get(node, np.nan)),
            "PageRank Score": float(pagerank.get(node, np.nan)),
            "Betweenness Centrality (Undirected)": float(betweenness.get(node, np.nan)),
        })

    return pd.DataFrame(rows).sort_values(
        ["HITS Hub Score", "PageRank Score", "Out-Degree Count"],
        ascending=[False, False, False],
        na_position="last",
    )


# ============================================================
# MAIN INTERFACE
# ============================================================

st.title("PrioriGraph 👑")
st.markdown("### Integration of Regulatory, Functional, and Topologic Flow Networks")

k_res = st.session_state.get("k_res")
g_res = st.session_state.get("g_res")

df_evidence, source_key = get_evidence_table()
df_reg = prepare_regulatory_edges()
df_enrichment = prepare_enrichment()

analysis_mode = st.session_state.get("analysis_mode", "single")
if analysis_mode == "single":
    if source_key == "df_multi_contrast": analysis_mode = "multi_contrast"
    elif source_key == "df_mt": analysis_mode = "multi_tissue"
    elif source_key == "df_meta": analysis_mode = "meta"


# ============================================================
# SIDEBAR CONTROLS
# ============================================================

with st.sidebar:
    st.header("🔍 Search and Filters")
    search_query = st.text_input("Locate Gene or Factor (Highlighted in Green):", "").strip().upper()
    st.divider()
    st.header("📊 Analysis Mode")

    if analysis_mode == "meta": st.success("🧬 Meta-analysis network")
    elif analysis_mode == "multi_tissue": st.success("🧬 Multi-Tissue Alignment network")
    elif analysis_mode == "multi_contrast": st.success("🧬 Multi-Contrast Alignment network")
    else: st.info("📊 Single-dataset network")

    st.divider()
    tf_source = st.radio("Regulatory Data Source:", ["Import from APP", "Upload Manual Tables"])
    num_vias = st.slider("Number of enriched terms:", 2, 15, 6)

    if analysis_mode in ["meta", "multi_tissue", "multi_contrast"] and not df_evidence.empty:
        st.divider()
        st.header(f"🧬 {analysis_mode.replace('_',' ').title()} Filters")
        meta_fdr = st.number_input("Maximum FDR:", min_value=0.0001, max_value=1.0, value=0.05, format="%.3f")
        meta_fc = st.number_input("Minimum |Effect| (Log2FC):", min_value=0.0, max_value=20.0, value=1.0, format="%.2f")
        min_studies = st.number_input("Min supporting datasets:", min_value=1, max_value=100, value=1, step=1) if analysis_mode == "meta" else 1
    else:
        meta_fdr = 0.05
        meta_fc = 1.0
        min_studies = 1

# ============================================================
# REGULATORY DATA
# ============================================================

if tf_source == "Import from APP":
    if not df_reg.empty:
        st.sidebar.success(f"Loaded {len(df_reg):,} TF-target interactions.")
    else:
        st.sidebar.error("Regulatory data is missing. Open the Candidate Regulators tab in APP first.")
else:
    uploaded_files = st.sidebar.file_uploader("Upload Tables JASPAR/TRRUST (CSVs)", type=["csv"], accept_multiple_files=True)
    if uploaded_files:
        all_data = []
        for file in uploaded_files:
            temp = pd.read_csv(file)
            c_tf = next((c for c in ["TF_Symbol", "TF_Name", "TF"] if c in temp.columns), None)
            c_tg = next((c for c in ["Genes", "Target_Gene", "Target"] if c in temp.columns), None)
            if c_tf and c_tg:
                temp = temp.assign(Target=temp[c_tg].astype(str).str.split(";")).explode("Target")
                keep = [c_tf, "Target"] + [c for c in ["Regulatory_Database", "Direction", "Adjusted P-value", "Term", "Overlap"] if c in temp.columns]
                all_data.append(temp[keep].rename(columns={c_tf: "TF"}))
        if all_data:
            df_reg = pd.concat(all_data, ignore_index=True)
            df_reg["TF"] = df_reg["TF"].astype(str).map(clean_gene)
            df_reg["Target"] = df_reg["Target"].astype(str).map(clean_gene)
            df_reg = df_reg[(df_reg["TF"] != "") & (df_reg["Target"] != "")].drop_duplicates()


# ============================================================
# LEGEND & ANALYSIS SUMMARY
# ============================================================

st.markdown("""
<div style="background-color:rgba(151, 166, 195, 0.08); padding:15px; border-radius:10px; border:1px solid rgba(151, 166, 195, 0.3); margin-bottom:20px;">
<h4 style="margin-top:0;">📍 Network Legend</h4>
<span style="color:#1C83E1; font-weight:bold;">● Blue Circle:</span> Target Gene with negative effect |
<span style="color:#FF4B4B; font-weight:bold;">● Red Circle:</span> Target Gene with positive effect |
<span style="color:#FF4B4B; font-weight:bold;">■ Red Square:</span> Candidate transcriptional regulator |
<span style="color:#FFD700; font-weight:bold;">⬢ Yellow Hexagon:</span> Pathways/Processes |
<span style="color:#28B463; font-weight:bold;">★ Green Highlight:</span> Searched Item
</div>
""", unsafe_allow_html=True)


# ============================================================
# NETWORK ELIGIBILITY & RENDERING
# ============================================================

if df_reg.empty:
    st.info("💡 Please process the Candidate Regulators tab in APP before building the PrioriGraph network.")
elif df_enrichment.empty:
    st.info("💡 No KEGG/GO enrichment data are currently available from APP.")
else:
    enrichment_fdr = pd.to_numeric(df_enrichment["Adjusted P-value"], errors="coerce")
    top_vias = (
        df_enrichment.loc[enrichment_fdr <= meta_fdr]
        .dropna(subset=["Term"])
        .sort_values("Adjusted P-value", na_position="last")
        .head(num_vias)
        .copy()
    )
    if top_vias.empty:
        st.warning("No enriched term passes the active FDR threshold. Adjust the APP/PG filters or inspect the complete APP tables.")
        st.stop()
    effect_map = build_effect_map(df_evidence, analysis_mode)
    focus_gene_set = set()
    use_focus = analysis_mode in ["meta", "multi_tissue", "multi_contrast"] and "Symbol" in df_evidence.columns

    if use_focus:
        meta_filtered = df_evidence.copy()
        effect_col = get_effect_column(meta_filtered)
        fdr_col = get_fdr_column(meta_filtered)
        if effect_col: meta_filtered[effect_col] = pd.to_numeric(meta_filtered[effect_col], errors="coerce")
        if fdr_col: meta_filtered[fdr_col] = pd.to_numeric(meta_filtered[fdr_col], errors="coerce")
        
        mask = pd.Series(True, index=meta_filtered.index)
        if fdr_col: mask &= (meta_filtered[fdr_col] <= meta_fdr)
        if effect_col: mask &= (meta_filtered[effect_col].abs() >= meta_fc)
        if analysis_mode == "meta" and "N_Studies" in meta_filtered.columns and min_studies > 1:
            n_studies_numeric = pd.to_numeric(meta_filtered["N_Studies"], errors="coerce")
            mask &= (n_studies_numeric >= min_studies)
        focus_gene_set = set(meta_filtered.loc[mask, "Symbol"].map(clean_gene))

    nodes = []
    edges = []
    added = set()
    all_network_genes = set()
    represented_term_names = set()

    for _, row in top_vias.iterrows():
        via_name = str(row["Term"]).strip().upper()
        via_genes = set(clean_gene(g) for g in str(row.get("Genes", "")).split(";") if clean_gene(g))
        target_genes = via_genes.intersection(set(df_reg["Target"]))
        
        if use_focus:
            target_genes = target_genes.intersection(focus_gene_set)
            
        if not target_genes: continue
        represented_term_names.add(via_name)
        is_via_search = (search_query != "" and search_query in via_name)

        if via_name not in added:
            nodes.append(Node(id=via_name, label=via_name[:35], size=50 if is_via_search else 45, color="#28B463" if is_via_search else "#FFD700", shape="hexagon"))
            added.add(via_name)

        for gene in sorted(target_genes):
            gene = clean_gene(gene)
            if not gene: continue
            all_network_genes.add(gene)
            is_gene_search = (search_query != "" and search_query == gene)
            if gene not in added:
                gene_color = node_gene_color(gene, effect_map, meta_fdr)
                gene_title = node_gene_title(gene, effect_map, analysis_mode)
                nodes.append(Node(id=gene, label=gene, size=45 if is_gene_search else 20, color="#28B463" if is_gene_search else gene_color, shape="dot", title=gene_title))
                added.add(gene)
            edges.append(Edge(source=via_name, target=gene, color="#FFD700", width=1, dashed=True))
            
            rel_tfs = df_reg[df_reg["Target"] == gene]["TF"].unique()
            for tf in rel_tfs:
                tf = clean_gene(tf)
                if not tf: continue
                is_tf_search = (search_query != "" and search_query == tf)
                if tf not in added:
                    nodes.append(Node(id=tf, label=tf, size=55 if is_tf_search else 30, color="#28B463" if is_tf_search else "#FF4B4B", shape="square", title=f"{tf} | Candidate regulator from target-set enrichment"))
                    added.add(tf)
                edges.append(Edge(source=tf, target=gene, directed=True, color="#999"))

    if nodes and edges:
        config = Config(width=1100, height=800, directed=True, physics=False, hierarchical=False)
        st.subheader(f"🕸️ Integrated Prior-Knowledge Network ({len(represented_term_names)} represented terms)")
        st.caption("TF → target edges originate from JASPAR/TRRUST target-set enrichment; pathway → gene edges indicate term membership. The graph is hypothesis-generating and does not establish causality.")
        agraph(nodes=nodes, edges=edges, config=config)

        displayed_edges = df_reg[df_reg["Target"].isin(all_network_genes)].drop_duplicates().copy()
        st.download_button(
            "📥 Download displayed TF-target edges with provenance",
            displayed_edges.to_csv(index=False).encode("utf-8"),
            "Lumos_PrioriGraph_edges.csv",
            mime="text/csv",
        )

        st.divider()
        standard_topology_df = pd.DataFrame()
        sensitivity_summary = pd.DataFrame()

        t1, t2, t3, t4, t5, t6, t7 = st.tabs([
            "👑 Candidate Regulators",
            "🎯 Highly Regulated Targets",
            "🧬 TF-Via Mapping",
            "🔬 Analysis Evidence",
            "📐 Standard Metrics",
            "🌊 Topology & Propagation",
            "🧭 Sensitivity",
        ])

        with t1:
            st.subheader("Candidate Regulators by Observed Out-Degree")
            m_rank = df_reg[df_reg["Target"].isin(all_network_genes)].groupby("TF")["Target"].nunique().sort_values(ascending=False).reset_index()
            m_rank.columns = ["Transcription Factor", "Unique Targets in Network"]
            st.dataframe(m_rank, width="stretch", hide_index=True)

        with t2:
            st.subheader("Target Genes by Observed In-Degree")
            h_rank = df_reg[df_reg["Target"].isin(all_network_genes)].groupby("Target")["TF"].nunique().sort_values(ascending=False).reset_index()
            h_rank.columns = ["Target Gene", "Total Number of Regulators"]
            st.dataframe(h_rank, width="stretch", hide_index=True)

        with t3:
            st.subheader("Hierarchical Analysis: Key Regulators by Pathway")
            mapping_rows = []
            for _, row in top_vias.iterrows():
                via = row["Term"]
                v_genes = set(clean_gene(g) for g in str(row.get("Genes", "")).split(";") if clean_gene(g)).intersection(all_network_genes)
                reg_in_via = df_reg[df_reg["Target"].isin(v_genes)]
                tf_counts = reg_in_via["TF"].value_counts()
                if not tf_counts.empty:
                    mapping_rows.append({
                        "Biological Pathway / Process": via,
                        "Top Candidate Regulator": tf_counts.index[0],
                        "TF Target Genes in This Pathway": int(tf_counts.iloc[0]),
                        "Genes in Pathway Represented in Network": reg_in_via["Target"].nunique()
                    })
            mapping_df = pd.DataFrame(mapping_rows)
            st.dataframe(mapping_df, width="stretch", hide_index=True)

        with t4:
            st.subheader(f"🔬 Gene-level Evidence")
            evidence_rows = []
            for gene in sorted(all_network_genes):
                info = effect_map.get(gene, {})
                evidence_rows.append({"Gene": gene, "Focus Effect": info.get("effect", np.nan), "Focus FDR": info.get("fdr", np.nan)})
            evidence_df = pd.DataFrame(evidence_rows).sort_values("Focus FDR", na_position="last")
            st.dataframe(evidence_df, width="stretch", hide_index=True)
            
        with t5:
            st.subheader("📐 Standard Directed-Network Metrics")
            st.caption(
                "HITS and PageRank are calculated on the directed TF → target layer. "
                "Betweenness is calculated on its undirected skeleton to identify bridge nodes."
            )
            if not NETWORKX_AVAILABLE:
                st.error("NetworkX is required for the standard topology metrics. Add networkx to the application requirements.")
            else:
                standard_topology_df = calculate_standard_topology(df_reg, all_network_genes)
                if not standard_topology_df.empty:
                    st.dataframe(standard_topology_df, width="stretch", hide_index=True)
                    st.download_button(
                        "📥 Download Standard Metrics",
                        standard_topology_df.to_csv(index=False).encode("utf-8"),
                        "Lumos_Standard_Topology_Metrics.csv",
                        mime="text/csv",
                    )
                else:
                    st.info("Insufficient TF-target connectivity for standard topology metrics.")

        with t6:
            st.subheader("🌊 Signal Propagation & Topological Impact")
            st.caption("Personalized RWR propagates absolute expression evidence downstream. Linear propagation retains the sign of the expression effect. Neither score is a causal effect estimate.")
            pc1, pc2 = st.columns(2)
            restart_probability = pc1.slider("RWR restart probability", 0.05, 0.95, 0.70, 0.05)
            propagation_decay = pc2.slider("Linear propagation decay", 0.05, 0.95, 0.50, 0.05)
            
            topological_df = calculate_signal_propagation(
                df_reg, all_network_genes, df_evidence, analysis_mode,
                restart_probability=restart_probability,
                propagation_decay=propagation_decay,
            )
            if not topological_df.empty:
                st.dataframe(topological_df, width="stretch", hide_index=True)
                st.download_button("📥 Download Topology Data", topological_df.to_csv(index=False).encode("utf-8"), "Lumos_Topological_Impact.csv", mime="text/csv")
            else:
                st.info("Insufficient network connectivity for propagation analysis.")

        with t7:
            st.subheader("🧭 Network Parameter Sensitivity")
            st.caption(
                "Perturbs gene-selection thresholds, RWR restart and linear-propagation decay. "
                "Stable nodes recur across configurations; this is not a biological bootstrap or causal validation."
            )
            sc1, sc2 = st.columns(2)
            sensitivity_fdr = sc1.multiselect(
                "FDR thresholds", [0.01, 0.025, 0.05, 0.10], default=[0.01, 0.05, 0.10]
            )
            sensitivity_effect = sc2.multiselect(
                "|Log2FC| thresholds", [0.0, 0.5, 1.0, 1.5, 2.0], default=[0.5, 1.0, 1.5]
            )
            sc3, sc4, sc5 = st.columns(3)
            sensitivity_restart = sc3.multiselect(
                "RWR restart", [0.3, 0.5, 0.7, 0.9], default=[0.5, 0.7, 0.9]
            )
            sensitivity_decay = sc4.multiselect(
                "Propagation decay", [0.25, 0.5, 0.75], default=[0.25, 0.5, 0.75]
            )
            sensitivity_top_k = sc5.number_input("Top-k stability", 5, 100, 20, 5)

            evidence_for_sensitivity = pd.DataFrame([
                {
                    "Symbol": str(gene),
                    "Effect": values.get("effect", np.nan),
                    "FDR": values.get("fdr", np.nan),
                }
                for gene, values in effect_map.items()
            ])
            if all([sensitivity_fdr, sensitivity_effect, sensitivity_restart, sensitivity_decay]):
                try:
                    sensitivity_rankings, sensitivity_summary = cached_network_sensitivity(
                        df_reg[["TF", "Target"]], evidence_for_sensitivity,
                        fdr_thresholds=sensitivity_fdr,
                        effect_thresholds=sensitivity_effect,
                        restart_probabilities=sensitivity_restart,
                        propagation_decays=sensitivity_decay,
                        top_k=int(sensitivity_top_k),
                    )
                except ValueError as exc:
                    st.error(str(exc))
                    sensitivity_rankings, sensitivity_summary = pd.DataFrame(), pd.DataFrame()
                if not sensitivity_summary.empty:
                    st.metric("Configurations evaluated", sensitivity_rankings["Configuration"].nunique())
                    st.dataframe(sensitivity_summary, width="stretch", hide_index=True)
                    sdl1, sdl2 = st.columns(2)
                    sdl1.download_button(
                        "📥 Download stability summary",
                        sensitivity_summary.to_csv(index=False).encode("utf-8"),
                        "Lumos_Network_Sensitivity_Summary.csv", mime="text/csv",
                    )
                    sdl2.download_button(
                        "📥 Download all configurations",
                        sensitivity_rankings.to_csv(index=False).encode("utf-8"),
                        "Lumos_Network_Sensitivity_Rankings.csv", mime="text/csv",
                    )
                else:
                    st.info("No TF-target subnetwork survived the selected sensitivity grid.")
            else:
                st.warning("Select at least one value for every sensitivity dimension.")

        # ========================================================
        # HTML REPORT
        # ========================================================

        st.sidebar.header("Export Report")
        if st.sidebar.button("🌐 Generate HTML Report"):
            fmts = {
                "Focus Effect": "fixed3", "Focus FDR": "exp2", "Pathway Adjusted P-value": "exp2",
                "Original_Log2FC": "fixed3", "Linear_Network_Propagation": "fixed3",
                "Personalized_RWR_Score": "exp2", "Restart_Probability": "fixed3",
                "Propagation_Decay": "fixed3", "RWR_Iterations": "int",
                "Betweenness Centrality (Undirected)": "fixed3", "HITS Hub Score": "exp2",
                "HITS Authority Score": "exp2", "PageRank Score": "exp2",
                "Unique Targets in Network": "int", "Total Number of Regulators": "int",
                "TF Target Genes in This Pathway": "int", "Genes in Pathway Represented in Network": "int",
                "Out-Degree Count": "int", "In-Degree Count": "int",
                "Configurations_Present": "int", "TopK_Frequency": "fixed3", "Configuration_Coverage": "fixed3",
                "Median_RWR_Rank": "fixed3", "IQR_RWR_Rank": "fixed3", "Mean_RWR_Score": "exp2", "SD_RWR_Score": "exp2",
                "Median_Linear_Propagation": "fixed3"
            }

            def spec(df, sort, asc=False, **kw):
                if not isinstance(df, pd.DataFrame) or df.empty:
                    return None
                return table_spec(df, [(c, c, fmts.get(c, "text")) for c in df.columns], sort, asc, **kw)

            pathway_padj = top_vias.drop_duplicates("Term").set_index("Term")["Adjusted P-value"]
            if not mapping_df.empty and "Biological Pathway / Process" in mapping_df.columns:
                mapping_df = mapping_df.assign(**{"Pathway Adjusted P-value": mapping_df["Biological Pathway / Process"].map(pathway_padj)})

            ev_kw = {}
            if {"Symbol", "Status"} <= set(df_evidence.columns) and df_evidence["Status"].isin(["Up-regulated", "Down-regulated"]).any():
                evidence_df = evidence_df.assign(Status=evidence_df["Gene"].map(df_evidence.drop_duplicates("Symbol", keep="last").set_index("Symbol")["Status"]))
                ev_kw = {"filter_col": "Status", "badge_col": "Status", "badge_map": {"Up-regulated": "up", "Down-regulated": "down"}}

            info = [
                ("Analysis mode", {"single": "Single dataset", "meta": "Meta-analysis", "multi_tissue": "Multi-tissue", "multi_contrast": "Multi-contrast"}.get(analysis_mode, analysis_mode)),
                ("Regulatory data source", tf_source),
                ("Represented terms", f"{len(represented_term_names)} in network (top {num_vias} enriched terms screened)"),
            ]
            if use_focus:
                info += [
                    ("Filters", f"FDR ≤ {meta_fdr}, |Effect| ≥ {meta_fc}" + (f", supporting datasets ≥ {min_studies}" if analysis_mode == "meta" else "")),
                    ("Genes passing filters", len(focus_gene_set)),
                ]
            if search_query:
                info.append(("Highlighted search", search_query))

            n_edges = len({(e.source, e.to) for e in edges})
            cards = [
                ("Nodes", len(nodes), ""),
                ("Edges", n_edges, ""),
                ("Represented terms", len(represented_term_names), ""),
                ("Candidate regulators", len(m_rank), ""),
                ("Target genes", len(all_network_genes), ""),
            ]
            sections = [
                {
                    "id": "network", "title": "Integrated Network",
                    "figs": [network_figure(nodes, edges) if NETWORKX_AVAILABLE else None],
                    "desc": "Yellow hexagons: pathways/processes; red squares: candidate transcriptional regulators; red/blue circles: target genes with positive/negative effect passing the FDR threshold (grey: not significant or no effect data); green: searched item. Grey arrows: TF → target; dotted yellow edges: term membership."
                },
                {
                    "id": "candidate-regulators", "title": "Candidate Regulators",
                    "desc": "Transcription factors ranked by unique target genes in the network.",
                    "table": spec(m_rank, "Unique Targets in Network", False)
                },
                {
                    "id": "highly-regulated-targets", "title": "Highly Regulated Targets",
                    "desc": "Target genes ranked by total number of regulating transcription factors.",
                    "table": spec(h_rank, "Total Number of Regulators", False)
                },
                {
                    "id": "tf-via-mapping", "title": "TF-Via Mapping",
                    "desc": "Key regulators mapped to biological pathways and processes.",
                    "table": spec(mapping_df, "Pathway Adjusted P-value" if "Pathway Adjusted P-value" in mapping_df.columns else "Biological Pathway / Process", True)
                },
                {
                    "id": "analysis-evidence", "title": "Analysis Evidence",
                    "desc": "Gene-level effect size and FDR from the focus analysis.",
                    "table": spec(evidence_df, "Focus FDR", True, **ev_kw)
                },
                {
                    "id": "standard-metrics", "title": "Standard Network Metrics",
                    "desc": "Complementary standard network metrics computed on the directed TF-target graph.",
                    "table": spec(standard_topology_df if NETWORKX_AVAILABLE else None, "HITS Hub Score", False, filter_col="Node Type")
                },
                {
                    "id": "topology-propagation", "title": "Topology & Propagation",
                    "desc": f"Personalized Random Walk with Restart (RWR, restart probability = {restart_probability}) and linear signal propagation (decay = {propagation_decay}).",
                    "table": spec(topological_df, "Personalized_RWR_Score", False)
                },
                {
                    "id": "sensitivity", "title": "Parameter Sensitivity",
                    "desc": "Node stability across FDR/effect thresholds, RWR restart and propagation decay configurations (top-k frequency, rank spread). Not a biological bootstrap or causal validation.",
                    "table": spec(sensitivity_summary, "TopK_Frequency", False)
                },
            ]
            footer = (
                "PrioriGraph (PG) &mdash; EvoMol-Lab, BioME, UFRN<br>"
                "Network metrics: NetworkX (Hagberg et al., 2008). Regulatory interactions: JASPAR (Rauluseviciute et al., 2024) / TRRUST (Han et al., 2018) via Enrichr. "
                "Propagation: RWR (Köhler et al., 2008), SPIA (Tarca et al., 2009)."
            )
            with st.spinner("Generating HTML report..."):
                html_str = build_html_report(
                    "PrioriGraph Network Report", PG_LOGO_PATH, info, cards, sections, footer, "PG",
                    stats_desc="Analysis mode, filters and size of the integrated regulatory network."
                )
                st.sidebar.download_button(
                    "📥 Download HTML Report", data=html_str.encode("utf-8"),
                    file_name=f"PG_Report_{datetime.now().strftime('%Y%m%d')}.html",
                    mime="text/html", on_click="ignore"
                )

st.session_state["pg_regulatory_edges"] = df_reg
st.session_state["pg_analysis_mode"] = analysis_mode
if analysis_mode in ["meta", "multi_tissue", "multi_contrast"]:
    st.session_state["pg_evidence"] = df_evidence.copy()
