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

st.markdown('<p class="subtitle" style="text-align: center; font-size: 20px; color: #5D6D7E; margin-bottom: 30px;">Sistemas de Análise Integrada: Transcriptômica e Redes</p>', unsafe_allow_html=True)
st.divider()

   # --- SIDEBAR ---
with st.sidebar:
    if os.path.exists(caminho_do_logo):
        st.image(caminho_do_logo, width=300)
    
    st.divider()
    st.markdown("### 🚀 Módulos")
    
    st.page_link("Lumos_Home.py", label="Home", icon="🏠")
    
    st.markdown('<p style="color:#2E86C1; font-weight:bold; margin-bottom:0px; margin-top:10px;">📊 Análise</p>', unsafe_allow_html=True)
    # AJUSTE OS NOMES ABAIXO PARA CASAR COM OS ARQUIVOS QUE VOCÊ VIU NA LISTA:
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    
    st.markdown('<p style="color:#28B463; font-weight:bold; margin-bottom:0px; margin-top:10px;">🧬 Funcional</p>', unsafe_allow_html=True)
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    
    st.markdown('<p style="color:#E67E22; font-weight:bold; margin-bottom:0px; margin-top:10px;">🕸️ Redes</p>', unsafe_allow_html=True)
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

    st.markdown('<p style="color:#E1AF12; font-weight:bold; margin-bottom:0px; margin-top:10px;">📚 Documentação</p>', unsafe_allow_html=True)
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")

    st.divider()
    
    # O restante do seu código (Guia de Conceitos, Metric, etc) 
    # deve seguir este mesmo alinhamento (indentação)
    st.info("Selecione um módulo acima para iniciar sua análise.")  
    
    st.markdown("### 🧠 Guia Rápido de Conceitos")
    
    with st.expander("🔎 O que é FDR (P-adj)?"):
        st.caption("Falsa Taxa de Descoberta. Corrige o valor de p para múltiplos testes.")
    with st.expander("🕸️ O que é um Gene Hub?"):
        st.caption("Definido pela alta conectividade (grau) em uma rede.")
    with st.expander("👑 O que é um Regulador Mestre?"):
        st.caption("Fator de Transcrição que orquestra grandes programas genéticos.")
    with st.expander("🧪 O que é GSEA?"):
        st.caption("Gene Set Enrichment Analysis. Analisa conjuntos de genes sem cortes fixos.")
    
    st.divider()
    st.metric(label="Release", value="v1.0.0-Beta")
    st.markdown("### Apoio Institucional")
    st.write("🔬 **EvoMol-Lab** e 🎓 **BioME - UFRN**")
    st.markdown("### Equipe de Desenvolvimento")
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
    O **Lumos Networks** é uma suíte bioinformática projetada para otimizar a transição entre dados brutos de sequenciamento 
    e a interpretação biológica sistêmica. Desenvolvido no **EvoMol-Lab**, o sistema integra rigor estatístico, 
    contextualização funcional e modelagem de redes complexas.
    """) 
    st.info("👈 **Navegação:** Use o menu à esquerda para alternar entre os módulos de análise.")
    st.subheader("Módulos da Suíte")
    
    # --- MÓDULO 1: DDEA ---
    with st.expander("📊 1. DDEA (Diagonal Differential Expression Alley)", expanded=True):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            <div class="module-text">
            <b>Identificação de Genes Diferencialmente Expressos (DEGs).</b><br>
            Este módulo utiliza o framework <b>PyDESeq2</b> para realizar a normalização por tamanho de biblioteca 
            e estimativa de dispersão genômica. Focado em experimentos de RNA-Seq e Microarray, permite:
            <ul>
                <li>Controle de taxa de falsos positivos via correção de Benjamini-Hochberg (FDR).</li>
                <li>Visualização diagnóstica através de Volcano Plots e MA Plots.</li>
                <li>Identificação de biomarcadores candidatos para análises a jusante.</li>
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
            <b>Contextualização Biológica e Ontologia.</b><br>
            Converte listas de genes em significado biológico através de integração multidimensional:
            <ul>
                <li><b>Enriquecimento:</b> Análise de sobre-representação (ORA) e GSEA em bancos KEGG e Gene Ontology (GO).</li>
                <li><b>Interatoma:</b> Integração com a API do <b>STRING</b> para construção de redes de interação proteína-proteína (PPI).</li>
                <li><b>Reguladores Mestres:</b> Predição de Fatores de Transcrição utilizando JASPAR e TRRUST para identificar o controle hierárquico do sistema.</li>
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
            <b>Modelagem de Redes de Biologia de Sistemas.</b><br>
            Realiza a modelagem topológica e análise de robustez de sistemas biológicos:
            <ul>
                <li><b>Centralidade:</b> Cálculo de métricas de grau (Degree) e intermediação (Betweenness) para identificar <b>Genes Hub</b>.</li>
                <li><b>Topologia:</b> Construção de redes direcionadas (TF-Target) e detecção de comunidades/módulos funcionais.</li>
                <li><b>Priorização:</b> Ranqueamento de alvos terapêuticos ou chaves regulatórias baseadas na arquitetura da rede.</li>
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
    st.markdown('<h2 class="section-header">Resumo Técnico</h2>', unsafe_allow_html=True)
    st.success("""
    **Pipeline Integrado:**
    1. **DDEA:** Filtra o sinal do ruído estatístico.
    2. **APP:** Descobre 'o quê' os genes fazem e 'quem' os manda.
    3. **PG:** Revela 'como' o sistema está organizado.
    """)
    # st.warning("⚠️ Certifique-se de que os nomes dos arquivos de imagem no diretório `src` coincidem com os chamados no código.")


# Footer final
st.divider()
footer_html = """<div style="text-align: center; color: #5D6D7E; padding: 10px;">
    <p><b>Lumos Networks © 2026</b> | Desenvolvido por integrantes do <b>EvoMol-Lab</b> (BioME - UFRN)</p></div>"""
st.markdown(footer_html, unsafe_allow_html=True)
