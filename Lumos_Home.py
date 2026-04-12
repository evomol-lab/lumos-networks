import streamlit as st

# Configuração da página
st.set_page_config(
    page_title="Lumos Networks | Home",
    page_icon="🧬",
    layout="wide"
)

# Estilização CSS personalizada
st.markdown("""
    <style>
    .main-title {
        font-size: 45px;
        font-weight: bold;
        color: #2E86C1;
        text-align: center;
        margin-top: -20px;
    }
    .subtitle {
        font-size: 20px;
        text-align: center;
        color: #5D6D7E;
        margin-bottom: 30px;
    }
    .section-header {
        color: #2E86C1;
        border-bottom: 2px solid #2E86C1;
        padding-bottom: 5px;
    }
    </style>
    """, unsafe_allow_html=True)

# Cabeçalho Principal
st.markdown('<p class="main-title">LUMOS NETWORKS: Systems Biology Suite </p>', unsafe_allow_html=True)
st.markdown('<p class="subtitle">Sistemas de Análise Integrada: Transcriptômica e Redes</p>', unsafe_allow_html=True)

st.divider()

# Colunas para organizar o conteúdo principal
col1, col2 = st.columns([2, 1])

with col1:
    st.markdown('<h2 class="section-header">Sobre o Projeto</h2>', unsafe_allow_html=True)
    st.write("""
    O **Lumos Networks** é uma suíte bioinformática projetada para otimizar a transição entre dados brutos de sequenciamento 
    e a interpretação biológica sistêmica. Desenvolvido no contexto de pesquisas de alto desempenho, o Lumos foca em 
    reprodutibilidade, precisão estatística e visualização dinâmica de redes complexas.
    """)
    
    st.subheader("Módulos da Suíte")
    
    with st.expander("1. DDEA (Differential Expression)", expanded=True):
        st.markdown("""
        Realiza a análise de expressão gênica diferencial utilizando o framework **PyDESeq2**, 
        garantindo rigor estatístico para experimentos com réplicas biológicas.
        """)
        
    with st.expander("2. APP (Enrichment Analysis)", expanded=False):
        st.markdown("""
        Módulo voltado para análise de enriquecimento funcional. Identifica vias metabólicas e processos 
        biológicos super-representados em listas de genes de interesse.
        """)
        
    with st.expander("3. PrioriGraph (Network Biology)", expanded=False):
        st.markdown("""
        Ferramenta de construção de redes baseada em conhecimento prévio e integração de dados, 
        permitindo a visualização de interações moleculares e identificação de hubs biológicos.
        """)

with col2:
    # Card de Informações
    st.info("👈 **Navegação:** Use o menu à esquerda para alternar entre os módulos.")
    
    st.metric(label="Release", value="v1.0.0-Beta")
    
    st.markdown("### Desenvolvimento")
    st.write("🔬 **Laboratório:** EvoMol-Lab (UFRN)")
    st.write("💻 **Tecnologias:** Python, Streamlit, PyDESeq2")
    
    st.markdown("### Apoio")
    st.write("🎓 **BioME** - Programa de Pós-Graduação em Bioinformática")
    # Se tiver os logos no GitHub, descomente abaixo:
    # st.image("assets/evomol_logo.png", width=150)

st.divider()

# Rodapé formatado
col_f1, col_f2 = st.columns([3, 1])
with col_f1:
    st.caption("Lumos Networks © 2026 | Desenvolvido por Integrantes do EvoMol-Lab no Centro de Biociências - UFRN")
with col_f2:
    st.caption("📍 Natal, RN - Brasil")
