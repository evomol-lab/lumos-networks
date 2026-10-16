import streamlit as st
import os

# 1. Configuração da página
st.set_page_config(
    page_title="Lumos Networks | Home",
    page_icon="🧬",
    layout="wide"
)

# --- CSS PARA ESCONDER O MENU PADRÃO E ESTILIZAR OS BALÕES ---
st.markdown("""
    <style>
    [data-testid="stSidebarNav"] {display: none;} /* Esconde o menu cinza padrão */
    .stPageLink {
        background-color: rgba(151, 166, 195, 0.12);
        border-radius: 20px;
        padding: 8px;
        border: 1px solid rgba(151, 166, 195, 0.3);
    }
    </style>
    """, unsafe_allow_html=True)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# --- CABEÇALHO ---
caminho_do_logo = os.path.join(BASE_DIR, "assets", "Lumos Networks.png") 
if os.path.exists(caminho_do_logo):
    col_espaço_1, col_logo, col_espaço_2 = st.columns([1, 1, 1])
    with col_logo:
        st.image(caminho_do_logo, width=500) 

titulo_magico = os.path.join(BASE_DIR, "titulo_lumos_hp.png")
if os.path.exists(titulo_magico):
    col_t1, col_t2, col_t3 = st.columns([1, 4, 1])
    with col_t2:
        st.image(titulo_magico, width="stretch")
else:
    st.markdown("<h1 style='text-align: center; font-size: 50px; color: #E1AF12;'>LUMOS NETWORKS</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center; color: #5D6D7E;'>Systems Biology Suite</h3>", unsafe_allow_html=True)

st.markdown('<p class="subtitle" style="text-align: center; font-size: 20px; color: #5D6D7E; margin-bottom: 30px;">Integrated Analysis Systems: Transcriptomics and Networks</p>', unsafe_allow_html=True)
st.divider()

   # --- SIDEBAR ---
with st.sidebar:
    if os.path.exists(caminho_do_logo):
        st.image(caminho_do_logo, width=300)
    
    st.divider()
    st.markdown("### 🚀 Modules")
    
    st.page_link("Lumos_Home.py", label="Home", icon="🏠")
    
    st.markdown('<p style="color:#2E86C1; font-weight:bold; margin-bottom:0px; margin-top:10px;">📊 Analysis</p>', unsafe_allow_html=True)
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    
    st.markdown('<p style="color:#28B463; font-weight:bold; margin-bottom:0px; margin-top:10px;">🧬 Functional</p>', unsafe_allow_html=True)
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    
    st.markdown('<p style="color:#E67E22; font-weight:bold; margin-bottom:0px; margin-top:10px;">🕸️ Networks</p>', unsafe_allow_html=True)
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

    st.markdown('<p style="color:#E1AF12; font-weight:bold; margin-bottom:0px; margin-top:10px;">📚 Documentation</p>', unsafe_allow_html=True)
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")

    st.divider()
    
    st.info("Select a module above to begin your analysis.")  
    
    st.markdown("### 🧠 Quick Guide to Concepts")
    
    with st.expander("🔎 What is FDR (P-adj)?"):
        st.caption("False Positive Rate. Corrects the p-value for multiple tests.")
    with st.expander("🕸️ What is a Gene Hub?"):
        st.caption("Defined by the high degree of connectivity in a network.")
    with st.expander("👑 What is a candidate regulator?"):
        st.caption("A transcription factor whose known target set is enriched among the selected genes; experimental validation is still required.")
    with st.expander("🧪 What is ORA?"):
        st.caption("Over-representation analysis tests whether a gene set contains more pathway members than expected under a reference universe.")
    
    st.divider()
    st.metric(label="Release", value="v1.3.0-beta")
    st.markdown("### Institutional Support")
    st.write("🔬 **EvoMol-Lab** e 🎓 **BioME - UFRN**")
    st.markdown("### Development Team")
    st.write("Dr. João Paulo M. S. Lima; MSc. Laís de Carvalho Gonçalves; Rodrigo Arruda Orvate")

st.divider()

# --- CONTEÚDO PRINCIPAL ---
st.markdown("""
    <style>
    .section-header { color: #2E86C1; border-bottom: 2px solid #2E86C1; padding-bottom: 5px; font-weight: bold; }
    .module-text { text-align: justify; font-size: 14px; }
    </style>
    """, unsafe_allow_html=True)

col1, col2 = st.columns([2, 1])

with col1:
    st.markdown('<h2 class="section-header">About the project</h2>', unsafe_allow_html=True)
    st.write("""
    **Lumos Networks** is a three-module research application that connects differential-expression
    analysis, functional over-representation analysis, and transparent prior-knowledge networks.
    Results are hypothesis-generating and retain the statistical and database provenance needed for review.
    """) 
    st.info("👈 **Navigation:** Use the menu on the left to switch between the analysis modules.")
    st.subheader("Suite Modules")
    
    # --- MÓDULO 1: DDEA ---
    with st.expander("📊 1. DDEA (Diagonal Differential Expression Alley)", expanded=True):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            <div class="module-text">
            <b>Identification of Differentially Expressed Genes (DEGs).</b><br>
            Raw integer RNA-seq counts are modeled with <b>PyDESeq2</b>; continuous-expression
            inputs use limma empirical-Bayes models in the publication image, with Welch retained only as an explicitly exploratory fallback. The module provides:
            <ul>
                <li>Control of the false positive rate using the Benjamini-Hochberg correction (FDR).</li>
                <li>Diagnostic visualization using Volcano Plots and MA Plots.</li>
                <li>Complete result tables and a machine-readable analysis manifest.</li>
            </ul>
            </div>
            """, unsafe_allow_html=True)
        with c2:
            caminho_ddea = os.path.join(BASE_DIR, "assets", "DDEA.png")
            if os.path.exists(caminho_ddea):
                st.image(caminho_ddea, width=80)
            else:
                st.write("📈")

    # --- MÓDULO 2: APP ---
    with st.expander("🧬 2. APP (Arithmancy Pathway Profiler)", expanded=False):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            <div class="module-text">
            <b>Biological Contextualization and Ontology.</b><br>
            Converts gene lists into biological meaning through multidimensional integration:
            <ul>
                <li><b>Enrichment:</b> Enrichr-based ORA against versioned KEGG and Gene Ontology libraries, optionally separated by expression direction.</li>
                <li><b>Associations:</b> Functional or physical STRING networks with retained confidence scores.</li>
                <li><b>Regulators:</b> Candidate TF target-set enrichment using JASPAR and TRRUST libraries.</li>
            </ul>
            </div>
            """, unsafe_allow_html=True)
        with c2:
            caminho_app = os.path.join(BASE_DIR, "assets", "APP.png")
            if os.path.exists(caminho_app):
                st.image(caminho_app, width=80)
            else:
                st.write("🧪")

    # --- MÓDULO 3: PG ---
    with st.expander("🕸️ 3. PG (PrioriGraph)", expanded=False):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            <div class="module-text">
            <b>Prior-knowledge network integration.</b><br>
            Connects enriched terms, genes, and candidate TFs without claiming causal inference:
            <ul>
                <li><b>Topology:</b> In/out-degree, HITS, PageRank, and bridge centrality.</li>
                <li><b>Provenance:</b> Export of displayed TF-target edges with source database and enrichment evidence.</li>
                <li><b>Propagation:</b> Personalized random walk and a clearly labelled linear propagation heuristic.</li>
            </ul>
            </div>
            """, unsafe_allow_html=True)
        with c2:
            caminho_pg = os.path.join(BASE_DIR, "assets", "PG.png")
            if os.path.exists(caminho_pg):
                st.image(caminho_pg, width=80)
            else:
                st.write("🕸️")

with col2:
    st.markdown('<h2 class="section-header">Technical Summary</h2>', unsafe_allow_html=True)
    st.success("""
    **Integrated Pipeline:**
    1. **DDEA:** estimates differential-expression effects and uncertainty.
    2. **APP:** contextualizes selected genes using versioned knowledge bases.
    3. **PG:** integrates evidence into an inspectable, exportable network.
    """)

# Footer final
st.divider()
footer_html = """<div style="text-align: center; color: #5D6D7E; padding: 10px;">
    <p><b>Lumos Networks © 2026</b> | Developed by members of the <b>EvoMol-Lab</b> (BioME - UFRN)</p></div>"""
st.markdown(footer_html, unsafe_allow_html=True)
