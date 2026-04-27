import streamlit as st
import os

# 1. Configuração da página (ajustei o título para cada módulo)
st.set_page_config(page_title="Lumos Networks | Analysis", page_icon="🧬", layout="wide")

# 2. CSS para manter o padrão visual (Igual à Home)
st.markdown("""
    <style>
    [data-testid="stSidebarNav"] {display: none;} /* Esconde o menu original */
    .stPageLink {
        background-color: #f0f2f6;
        border-radius: 20px;
        padding: 8px;
        border: 1px solid #e0e4eb;
    }
    .section-header { color: #2E86C1; border-bottom: 2px solid #2E86C1; padding-bottom: 5px; font-weight: bold; }
    </style>
    """, unsafe_allow_html=True)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# Como os módulos estão dentro de 'pages', subimos um nível para achar a logo
LOGO_PATH = os.path.join(os.path.dirname(BASE_DIR), "assets", "Lumos Networks.png")

# --- SIDEBAR PADRONIZADA ---
with st.sidebar:
    if os.path.exists(LOGO_PATH):
        st.image(LOGO_PATH, width=250)
    
    st.divider()
    st.markdown("### 🚀 Navigation")
    
    # IMPORTANTE: Caminhos saindo da pasta 'pages'
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
    st.info("You are in the Documentation module.")
    
st.set_page_config(page_title="Lumos | Documentation", layout="wide", page_icon="📚")

st.title("📚 Technical Documentation and Fundamentals")
st.markdown("---")

tab1, tab2, tab3 = st.tabs(["🧬 DDEA (Statistics)", "🧪 APP (Functional)", "🕸️ PG (Systems)"])

with tab1:
    st.header("Diagonal Differential Expression Alley (DDEA)")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(r"""
        ### 🔬 Statistical Foundations
        
        The **DDEA** module implements differential expression analysis based on generalized linear models (**GLM**).
        The choice of the **Negative Binomial (NB)** distribution is essential for addressing the discrete nature of sequencing data 
        and the characteristic of **overspread**, where the variance is greater than the mean.
        #### 1. Stochastic Modeling
        The read count $K_{ij}$ for gene $i$ in sample $j$ is modeled as:
        $$K_{ij} \sim NB(\mu_{ij}, \alpha_i)$$
        Where:
        * **$\mu_{ij}$**: represents the expected average of the counts.
        * **$\alpha_i$**: is the gene-specific dispersion parameter, which captures extra-Poisson biological variability.
        #### 2. Normalization and Robustness
        To enable comparisons between samples with different sequencing depths, we used the **Median-of-Ratios** method. 
        This method is superior to simple normalization per million (RPM/TPM) because:
        
        * **Corrects compositional biases:** Ensures that genes highly expressed in one condition do not cause false negatives in the others.
        * **Stability:** It uses a virtual reference sample based on the geometric mean, making the subsequent Wald test much more reliable.
        “”")

    with col2:
        st.subheader("🛠️ Workflow")
        st.info("""
        1. **Pre-processing:** Filtering genes with low copy numbers (Low Counts).
        2. **Estimation:** Calculation of size factors (Size Factors).
        3. **Testing:** Wald's test for statistical significance.
        4. **Correction:** P-value adjustment using the Benjamini-Hochberg method (FDR).
        """)

with tab2:
    st.header("Arithmancy Pathway Profiler (APP)")
    
    st.markdown("""
    This module performs **functional annotation** of the genes identified in the DDEA. It does not merely list terms; it calculates the statistical probability that a pathway is active.
    """)
    
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("### 📊 Methods")
        st.write("- **ORA (Over-Representation Analysis):** It uses the hypergeometric test to evaluate the intersection between your list and the database.")
        st.latex(r"P(X > k) = \sum_{i=k+1}^{n} \frac{\binom{M}{i} \binom{N-M}{n-i}}{\binom{N}{n}}")
        st.write("- **GSEA:** Analyze the complete gene ranking without threshold cuts.")

    with c2:
        st.markdown("### 🗃️ Databases")
        st.markdown("""
        **1. Kyoto Encyclopedia of Genes and Genomes (KEGG):**
        Focus on metabolic pathways and classical biochemical signaling cascades.
        
        **2. Gene Ontology (GO):**
        Hierarchical classification of genes and their products. Lumos Networks prioritizes Biological Process (BP) analysis, enabling the identification of specific metabolic pathways and signaling cascades activated in the biological system under study.
      
        **3. STRING (PPI Networks):**
        Mapping of protein-protein interactions. This allows you to see whether the products of differentially expressed genes interact physically or functionally.
        
        **4. JASPAR & TRRUST:**
        Enrichment analysis of **transcription factor (TF)** binding sites, identifying potential master regulators of the system.
        """)

with tab3:
    st.header("PrioriGraph (PG)")
    st.markdown("PG models the system topology by integrating omics data with curated biological knowledge.")
    
    st.subheader("🕸️ Network Architecture and Hierarchy")
    
    c_h1, c_h2 = st.columns(2)
    with c_h1:
        st.markdown("**1. Gene Hubs**")
        st.caption("Based on degree centrality. Identifies genes with high local connectivity that are essential for the stability of the immediate module.")
    
    with c_h2:
        st.markdown("**2. Master Regulators**")
        st.caption("Based on hierarchical influence. Genes that control transcriptional cascades and modulate the system’s overall response.")

    st.divider()
    
    st.subheader("📊 Centrality Metrics")
    cc1, cc2, cc3 = st.columns(3)
    with cc1:
        st.markdown("**Betweenness**")
        st.caption("Identifies “links” between different biological processes.")
    with cc2:
        st.markdown("**Clustering Coefficient**")
        st.caption("Measures the functional cohesion of the identified modules.")
    with cc3:
        st.markdown("**Eigenvector Centrality**")
        st.caption("Evaluates the influence of a gene based on the importance of its neighbors.")
        
    st.warning("⚠️ **Research Note:** PrioriGraph prioritizes interactions supported by experimental evidence (Curated Databases), thereby reducing the noise associated with purely computational predictions.")

st.markdown("---")
st.caption("Lumos Networks Documentation - EvoMol-Lab 2026")
