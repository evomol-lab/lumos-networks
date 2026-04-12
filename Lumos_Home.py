import streamlit as st
import os

# Configuração da página
st.set_page_config(
    page_title="Lumos Networks | Home",
    page_icon="🧬",
    layout="wide"
)

# --- 1. LOGO NA BARRA LATERAL ---
with st.sidebar:
    logo_lumos = "assets/logos/Lumos Networks.png"
    if os.path.exists(logo_lumos):
        st.image(logo_lumos, use_container_width=True)
    st.divider()
    st.info("Selecione um módulo acima para iniciar sua análise.")

# Estilização CSS personalizada
st.markdown("""
    <style>
    .main-title { font-size: 45px; font-weight: bold; color: #2E86C1; text-align: center; margin-top: -20px; }
    .subtitle { font-size: 20px; text-align: center; color: #5D6D7E; margin-bottom: 30px; }
    .section-header { color: #2E86C1; border-bottom: 2px solid #2E86C1; padding-bottom: 5px; }
    </style>
    """, unsafe_allow_html=True)

# Cabeçalho Principal
st.markdown('<p class="main-title">LUMOS NETWORKS: Systems Biology Suite </p>', unsafe_allow_html=True)
st.markdown('<p class="subtitle">Sistemas de Análise Integrada: Transcriptômica e Redes</p>', unsafe_allow_html=True)

st.divider()

col1, col2 = st.columns([2, 1])

with col1:
    st.markdown('<h2 class="section-header">Sobre o Projeto</h2>', unsafe_allow_html=True)
    st.write("""
    O **Lumos Networks** é uma suíte bioinformática projetada para otimizar a transição entre dados brutos de sequenciamento 
    e a interpretação biológica sistêmica. Desenvolvido no EvoMol-Lab, o sistema integra análise estatística, 
    contextualização funcional e modelagem de redes.
    """)
    
    st.subheader("Módulos da Suíte")
    
    # --- MÓDULO 1: DDEA ---
    with st.expander("1. DDEA (Diagonal Differential Expression Alley)", expanded=True):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            Análise de expressão gênica diferencial utilizando **PyDESeq2** e modelos lineares. 
            Focado em rigor estatístico e normalização robusta para Microarray e RNA-Seq.
            """)
        with c2:
            if os.path.exists("assets/logos/DDEA.png"):
                st.image("assets/logos/DDEA.png", width=80)
        
    # --- MÓDULO 2: APP ---
    with st.expander("2. APP (Arithmancy Pathway Profiler)", expanded=False):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            Integração multidimensional com bases KEGG, GO, STRING e circuitos regulatórios (TFs).
            """)
        with c2:
            if os.path.exists("assets/logos/Arithmancy Pathway Profiler (APP) .png"):
                st.image("assets/logos/Arithmancy Pathway Profiler (APP) .png", width=80)

    # --- MÓDULO 3: PRIORIGRAPH (O QUE FALTAVA) ---
    with st.expander("3. PG (PrioriGraph)", expanded=False):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            **Modelagem de Redes de Biologia de Sistemas.** Constrói grafos de interação baseados em conhecimento prévio, integrando DEGs e Fatores de 
            Transcrição para identificar hubs regulatórios e comunidades moleculares.
            """)
        with c2:
            if os.path.exists("assets/logos/PrioriGraph(PG).jpeg"):
                st.image("assets/logos/PrioriGraph(PG).jpeg", width=80)

with col2:
    st.info("👈 **Navegação:** Use o menu à esquerda para alternar entre os módulos.")
    
    if os.path.exists("assets/logos/Lumos Networks.png"):
        st.image("assets/logos/Lumos Networks.png", use_container_width=True)
    
    st.metric(label="Release", value="v1.0.0-Beta")
    
    st.markdown("### Apoio Institucional")
    st.write("🔬 **EvoMol-Lab**")
    st.write("🎓 **BioME - UFRN**")

st.divider()

# Rodapé centralizado e sem nomes individuais (conforme solicitado)
footer_html = """
<div style="text-align: center; color: #5D6D7E; padding: 10px;">
    <p><b>Lumos Networks © 2026</b> | Desenvolvido por integrantes do <b>EvoMol-Lab</b> (BioME - UFRN)</p>
</div>
"""
st.markdown(footer_html, unsafe_allow_html=True)
