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
        background-color: #f0f2f6;
        border-radius: 20px;
        padding: 8px;
        border: 1px solid #e0e4eb;
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
        st.image(titulo_magico, use_container_width=True)
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
    # AJUSTE OS NOMES ABAIXO PARA CASAR COM OS ARQUIVOS QUE VOCÊ VIU NA LISTA:
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    
    st.markdown('<p style="color:#28B463; font-weight:bold; margin-bottom:0px; margin-top:10px;">🧬 Functional</p>', unsafe_allow_html=True)
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    
    st.markdown('<p style="color:#E67E22; font-weight:bold; margin-bottom:0px; margin-top:10px;">🕸️ Networks</p>', unsafe_allow_html=True)
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

    st.markdown('<p style="color:#E1AF12; font-weight:bold; margin-bottom:0px; margin-top:10px;">📚 Documentation</p>', unsafe_allow_html=True)
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")

    st.divider()
    
    # O restante do seu código (Guia de Conceitos, Metric, etc) 
    # deve seguir este mesmo alinhamento (indentação)
    st.info("Select a module above to begin your analysis.")  
    
    st.markdown("### 🧠 Quick Guide to Concepts")
    
    with st.expander("🔎 What is FDR (P-adj)?"):
        st.caption("False Positive Rate. Corrects the p-value for multiple tests.")
    with st.expander("🕸️ What is a Gene Hub?"):
        st.caption("Defined by the high degree of connectivity in a network.")
    with st.expander("👑 What is a Master Regulator?"):
        st.caption("A transcription factor that orchestrates large-scale genetic programs.")
    with st.expander("🧪 What is GSEA?"):
        st.caption("Gene Set Enrichment Analysis. Analyzes gene sets without fixed cutoffs.")
    
    st.divider()
    st.metric(label="Release", value="v1.0.0-Beta")
    st.markdown("### Institutional Support")
    st.write("🔬 **EvoMol-Lab** e 🎓 **BioME - UFRN**")
    st.markdown("### Development Team")
    st.write("Dr. João Paulo M. S. Lima; MSc. Laís de Carvalho Gonçalves; Rodrigo Arruda Orvate")

st.divider()

# --- CONTEÚDO PRINCIPAL (MANTIDO IGUAL) ---
st.markdown("""
    <style>
    .section-header { color: #2E86C1; border-bottom: 2px solid #2E86C1; padding-bottom: 5px; font-weight: bold; }
    .module-text { text-align: justify; font-size: 14px; }
    </style>
    """, unsafe_allow_html=True)

col1, col2 = st.columns([2, 1])

with col1:
    st.markdown('<h2 class="section-header">Sobre o Projeto</h2>', unsafe_allow_html=True)
    st.write("""
    **Lumos Networks** is a bioinformatics suite designed to streamline the transition from raw sequencing 
    data to systemic biological interpretation. Developed at **EvoMol-Lab**, the system integrates statistical rigor, 
    functional contextualization, and complex network modeling.
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
            This module uses the <b>PyDESeq2</b> framework to perform library-size normalization 
            and genomic dispersion estimation. Designed for RNA-Seq and microarray experiments, it allows you to:
            <ul>
                <li>Control of the false positive rate using the Benjamini-Hochberg correction (FDR).</li>
                <li>Diagnostic visualization using Volcano Plots and MA Plots.</li>
                <li>Identification of candidate biomarkers for downstream analyses.</li>
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
                <li><b>Enrichment:</b> Over-representation analysis (ORA) and GSEA in KEGG and Gene Ontology (GO) databases.</li>
                <li><b>Interatoma:</b> Integration with the <b>STRING</b> API for constructing protein-protein interaction (PPI) networks.</li>
                <li><b>Master Regulators:</b> Prediction of transcription factors using JASPAR and TRRUST to identify the hierarchical control of the system.</li>
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
            <b>Modeling of Systems Biology Networks.</b><br>
            Performs topological modeling and robustness analysis of biological systems:
            <ul>
                <li><b>Centrality:</b> Calculation of degree and betweenness metrics to identify <b>Hub Genes</b>.</li>
                <li><b>Topology:</b> Construction of directed networks (TF-Target) and detection of functional communities/modules.</li>
                <li><b>Prioritization:</b> Ranking of therapeutic targets or regulatory switches based on the network architecture.</li>
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
    1. **DDEA:** Filters the signal from statistical noise.
    2. **APP:** Find out ‘what’ genes do and ‘who’ controls them.
    3. **PG:** It reveals *how* the system is organized.
    """)
    # st.warning("⚠️ Certifique-se de que os nomes dos arquivos de imagem no diretório `src` coincidem com os chamados no código.")


# Footer final
st.divider()
footer_html = """<div style="text-align: center; color: #5D6D7E; padding: 10px;">
    <p><b>Lumos Networks © 2026</b> | Developed by members of the <b>EvoMol-Lab</b> (BioME - UFRN)</p></div>"""
st.markdown(footer_html, unsafe_allow_html=True)
