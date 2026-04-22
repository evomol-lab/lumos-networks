import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from streamlit_agraph import agraph, Node, Edge, Config
import os

# 1. Configuração ÚNICA da página (Removendo duplicidade)
st.set_page_config(page_title="Lumos Networks | PrioriGraph", page_icon="🕸️", layout="wide")

# 2. CSS Padronizado
st.markdown("""
    <style>
    [data-testid="stSidebarNav"] {display: none;}
    .stPageLink { background-color: #f0f2f6; border-radius: 20px; padding: 8px; border: 1px solid #e0e4eb; }
    </style>
    """, unsafe_allow_html=True)

# --- 3. INICIALIZAÇÃO DE VARIÁVEIS DE ESTADO (Resolve o NameError) ---
if 'search_term' not in st.session_state:
    st.session_state['search_term'] = ""

# --- SIDEBAR: CONFIGURAÇÕES ---
with st.sidebar:
    st.header("🔍 Filtros da Rede")
    # Usamos o session_state diretamente no widget
    st.session_state['search_term'] = st.text_input("Localizar Gene/TF:", value=st.session_state['search_term']).strip().upper()
    
    st.divider()
    uploaded_files = st.file_uploader("Upload CSVs (JASPAR/TRRUST)", type=['csv'], accept_multiple_files=True)
    
    st.divider()
    g_width = st.slider("Largura do Grafo:", 600, 1800, 1100)
    g_height = st.slider("Altura do Grafo:", 400, 1200, 800)
    
    st.divider()
    st.markdown("### 🚀 Navegação")
    st.page_link("Lumos_Home.py", label="Início", icon="🏠")
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

# --- RECUPERAÇÃO DE DADOS DO APP ---
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')

st.title("PrioriGraph 👑")

if uploaded_files:
    # Processamento de Dados
    all_interactions = []
    for f in uploaded_files:
        df_temp = pd.read_csv(f)
        col_tf = next((c for c in ['TF_Symbol', 'TF_Name', 'TF'] if c in df_temp.columns), None)
        col_genes = next((c for c in ['Genes', 'Target_Gene', 'Target'] if c in df_temp.columns), None)
        if col_tf and col_genes:
            df_temp = df_temp.assign(Gene=df_temp[col_genes].astype(str).str.split(';')).explode('Gene')
            df_temp = df_temp[[col_tf, 'Gene']].rename(columns={col_tf: 'TF', 'Gene': 'Target'})
            all_interactions.append(df_temp)
    
    df_final = pd.concat(all_interactions).drop_duplicates()
    df_final['TF'] = df_final['TF'].astype(str).str.upper()
    df_final['Target'] = df_final['Target'].astype(str).str.upper()
    tfs_list = set(df_final['TF'].unique())

    # --- CÁLCULO DE IMPACTO FUNCIONAL ---
    def calculate_impact(tf_name):
        targets = set(df_final[df_final['TF'] == tf_name]['Target'])
        score = 0
        if k_res is not None:
            for _, r in k_res.iterrows():
                if targets.intersection(set(str(r['Genes']).split(';'))): score += 2
        if g_res is not None:
            for _, r in g_res.iterrows():
                if targets.intersection(set(str(r['Genes']).split(';'))): score += 1
        return max(1, score)

    impact_map = {tf: calculate_impact(tf) for tf in tfs_list}

    # --- CONSTRUÇÃO DO GRAFO (Sem movimento) ---
    nodes, edges, nodes_added = [], [], set()
    s_query = st.session_state['search_term']

    for _, row in df_final.iterrows():
        u, v = row['TF'], row['Target']
        
        # Nó TF (Diamante)
        if u not in nodes_added:
            imp = impact_map.get(u, 1)
            color = "#FFD700" if imp > 8 or u == s_query else "#FF4B4B"
            nodes.append(Node(id=u, label=f"{u} (Imp:{imp})", size=25+(imp*2), color=color, shape="diamond"))
            nodes_added.add(u)

        # Nó Target (Círculo)
        if v not in nodes_added:
            is_match = (v == s_query)
            nodes.append(Node(id=v, label=v, size=45 if is_match else 15, color="#FFD700" if is_match else "#1C83E1"))
            nodes_added.add(v)
        
        edges.append(Edge(source=u, target=v, directed=True, color="#999"))

    # CONFIGURAÇÃO: physics=False trava os nós no lugar
    config = Config(width=g_width, height=g_height, directed=True, physics=False, hierarchical=False)
    agraph(nodes=nodes, edges=edges, config=config)

    # --- RANKINGS ORIGINAIS ---
    st.divider()
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("🏆 Master Regulators")
        st.dataframe(df_final['TF'].value_counts().reset_index().rename(columns={'count':'Alvos'}), use_container_width=True)
    with c2:
        st.subheader("🎯 Genetic Hubs")
        st.dataframe(df_final['Target'].value_counts().reset_index().rename(columns={'count':'Reguladores'}), use_container_width=True)

else:
    st.info("Aguardando upload dos ficheiros CSV para reconstruir a rede.")
