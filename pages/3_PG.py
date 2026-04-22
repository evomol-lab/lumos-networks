import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from streamlit_agraph import agraph, Node, Edge, Config

st.set_page_config(layout="wide", page_title="PrioriGraph")

import os
import streamlit as st


# 1. Configuração da página (ajuste o título para cada módulo)
st.set_page_config(page_title="Lumos Networks | Análise", page_icon="🧬", layout="wide")

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
    st.markdown("### 🚀 Navegação")
    
    # IMPORTANTE: Caminhos saindo da pasta 'pages'
    st.page_link("Lumos_Home.py", label="Página Inicial", icon="🏠")
    
    st.markdown('<p style="color:#2E86C1; font-weight:bold; margin-bottom:0px; margin-top:10px;">📊 Análise</p>', unsafe_allow_html=True)
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    
    st.markdown('<p style="color:#28B463; font-weight:bold; margin-bottom:0px; margin-top:10px;">🧬 Funcional</p>', unsafe_allow_html=True)
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    
    st.markdown('<p style="color:#E67E22; font-weight:bold; margin-bottom:0px; margin-top:10px;">🕸️ Redes</p>', unsafe_allow_html=True)
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

    st.markdown('<p style="color:#E1AF12; font-weight:bold; margin-bottom:0px; margin-top:10px;">📚 Documentação</p>', unsafe_allow_html=True)
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")

    st.divider()
    st.info("Você está no módulo de redes.")


# 1. Localização atual: /code/src/pages/seu_script.py
FILE_DIR = os.path.dirname(os.path.abspath(__file__))

# 2. Subir UM nível para chegar na pasta 'src'
PARENT_DIR = os.path.dirname(FILE_DIR)

# 3. Apontar para o arquivo que está solto na 'src'
logo_path = os.path.join(PARENT_DIR, "assets", "PG.png")

with st.sidebar:
    if os.path.exists(logo_path):
        st.image(logo_path, use_container_width=True)

# ============================================================
# INTERFACE PRINCIPAL
# ============================================================

st.title("PrioriGraph 👑")
st.markdown("### Integração de Redes Regulatórias (TF -> Genes)")

with st.sidebar:
    st.header("📂 1. Carga de Dados")
    uploaded_files = st.file_uploader("Upload Tabelas JASPAR/TRRUST (CSVs)", type=['csv'], accept_multiple_files=True)
    
    st.divider()
    st.header("🔍 2. Localizar na Rede")
    search_query = st.text_input("Digite o nome do TF ou Gene:", "").strip().upper()
    
    st.divider()
    st.header("📐 Ajustes de Visualização")
    g_width = st.slider("Largura da Área:", 600, 1800, 1100, 50)
    g_height = st.slider("Altura da Área:", 400, 1200, 800, 50)
    
    st.divider()
    st.header("🎨 Legenda da Rede")
    st.markdown(f"""
    <div style="display: flex; align-items: center; margin-bottom: 10px;">
        <div style="width: 20px; height: 20px; background-color: #FF4B4B; border-radius: 50%; margin-right: 10px; border: 1px solid #333;"></div>
        <b>Fator de Transcrição (TF)</b>
    </div>
    <div style="display: flex; align-items: center; margin-bottom: 10px;">
        <div style="width: 20px; height: 20px; background-color: #1C83E1; border-radius: 50%; margin-right: 10px; border: 1px solid #333;"></div>
        <b>Gene Alvo (Target)</b>
    </div>
    <div style="display: flex; align-items: center; margin-bottom: 10px;">
        <div style="width: 20px; height: 20px; background-color: #FFD700; border-radius: 50%; margin-right: 10px; border: 2px solid #000;"></div>
        <b>Destaque: {search_query if search_query else 'Nenhum'}</b>
    </div>
    """, unsafe_allow_html=True)

if uploaded_files:
    all_interactions = []
    # DataFrame para o gráfico de barras (precisa do FDR/Adj P-value)
    plot_data_list = []

    for f in uploaded_files:
        df_temp = pd.read_csv(f)
        col_tf = next((c for c in ['TF_Symbol', 'TF_Name'] if c in df_temp.columns), None)
        col_genes = next((c for c in ['Genes', 'Target_Gene'] if c in df_temp.columns), None)
        
        if col_tf and col_genes:
            # Dados para o gráfico (Top TFs por significância)
            if 'Adjusted P-value' in df_temp.columns:
                plot_data_list.append(df_temp[[col_tf, 'Adjusted P-value']].copy())

            # Dados para a rede
            df_temp[col_genes] = df_temp[col_genes].astype(str)
            df_temp = df_temp.assign(Gene=df_temp[col_genes].str.split(';')).explode('Gene')
            df_temp = df_temp[[col_tf, 'Gene']].rename(columns={col_tf: 'TF', 'Gene': 'Target'})
            all_interactions.append(df_temp)
    
    if not all_interactions:
        st.error("Nenhuma coluna compatível encontrada.")
        st.stop()
        
    # === 1. CONSOLIDAÇÃO DOS DADOS ===
    df_final = pd.concat(all_interactions).drop_duplicates()
    df_final['TF'] = df_final['TF'].astype(str).str.strip().str.upper()
    df_final['Target'] = df_final['Target'].astype(str).str.strip().str.upper()
    tfs_list = set(df_final['TF'].unique())

    # === 2. MOTOR DE INTEGRAÇÃO (A "MAGIA" ÚNICA) ===
    # Recuperamos os resultados salvos no st.session_state pelos módulos anteriores
    k_res = st.session_state.get('k_res')       # Vias KEGG (do APP)
    g_res = st.session_state.get('g_res')       # Processos GO (do APP)
    s_df = st.session_state.get('string_df')    # Rede STRING (do APP)

    def get_integrated_weight(tf_name, interactions_df):
        # Genes alvos que este TF regula nesta rede
        my_targets = set(interactions_df[interactions_df['TF'] == tf_name]['Target'])
        if not my_targets: return 1
        
        score = 0
        # Bônus KEGG: Se os alvos do TF estão em vias metabólicas
        if k_res is not None:
            for _, row in k_res.iterrows():
                if my_targets.intersection(set(str(row['Genes']).split(';'))):
                    score += 2
        
        # Bônus GO: Se os alvos do TF participam de processos biológicos
        if g_res is not None:
            for _, row in g_res.iterrows():
                if my_targets.intersection(set(str(row['Genes']).split(';'))):
                    score += 1
        
        # Bônus STRING: Se o TF tem interação física comprovada
        if s_df is not None:
            string_genes = set(s_df['preferredName_A']).union(set(s_df['preferredName_B']))
            if tf_name in string_genes:
                score += 3
                
        return max(1, score)

    # Calculamos o peso de impacto para cada TF
    tf_impact_map = {tf: get_integrated_weight(tf, df_final) for tf in tfs_list}

    # === 3. CONSTRUÇÃO DA REDE VISUAL (AGRAPH) ===
    nodes, edges, nodes_added = [], [], set()

    # Primeiro, criamos as conexões (Arestas)
    for _, row in df_final.iterrows():
        tf, target = row['TF'], row['Target']
        if tf == 'NAN' or target == 'NAN': continue

        # Criar Nó do TF (Diamante)
        if tf not in nodes_added:
            impacto = tf_impact_map.get(tf, 1)
            # Dourado se impacto > 10 (Mestre Regulador), senão Vermelho
            color = "#FFD700" if impacto > 10 else "#FF4B4B" 
            size = 25 + (impacto * 2) # Tamanho proporcional à importância biológica
            nodes.append(Node(id=tf, label=f"{tf} (Impacto: {impacto})", 
                              size=size, color=color, shape="diamond"))
            nodes_added.add(tf)

        # Criar Nó do Gene Alvo (Círculo)
        if target not in nodes_added:
            # Destaque se for o gene pesquisado na sidebar
            is_search = search_query and search_query == target
            color = "#FFD700" if is_search else "#1C83E1"
            size = 40 if is_search else 15
            nodes.append(Node(id=target, label=target, size=size, color=color, shape="dot"))
            nodes_added.add(target)

        # Criar a linha de conexão
        edges.append(Edge(source=tf, target=target, directed=True, color="#999"))

    # === 4. EXIBIÇÃO DA REDE ===
    st.subheader("🕸️ Rede Regulatória com Carga Funcional (KEGG + GO + STRING)")
    config = Config(width=g_width, height=g_height, directed=True, physics=False, hierarchical=False)
    agraph(nodes=nodes, edges=edges, config=config)
    
    # 3. Rankings Decrescentes
    st.divider()
    col_rank1, col_rank2 = st.columns(2)
    with col_rank1:
        st.subheader("🏆 Master Regulators (Out-Degree)")
        master_rank = df_final['TF'].value_counts().reset_index()
        master_rank.columns = ['Fator de Transcrição', 'Nº de Alvos']
        st.dataframe(master_rank, use_container_width=True, hide_index=True)
    with col_rank2:
        st.subheader("🎯 Genetic Hubs (In-Degree)")
        hub_rank = df_final['Target'].value_counts().reset_index()
        hub_rank.columns = ['Gene Alvo', 'Nº de Reguladores']
        st.dataframe(hub_rank, use_container_width=True, hide_index=True)

    st.divider()
    st.subheader("📄 Tabela Consolidada de Interações")
    st.dataframe(df_final.sort_values('TF'), use_container_width=True, hide_index=True)

else:
    st.info("Aguardando upload dos arquivos CSV.")
