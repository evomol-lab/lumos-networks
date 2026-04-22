import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from streamlit_agraph import agraph, Node, Edge, Config
import os

# 1. Configuração ÚNICA e Inicial
st.set_page_config(page_title="Lumos Networks | PrioriGraph", page_icon="🕸️", layout="wide")

# 2. Inicialização de Variáveis de Busca (Evita NameError)
if 'search_term' not in st.session_state:
    st.session_state['search_term'] = ""

# 3. CSS para manter o padrão visual
st.markdown("""
    <style>
    [data-testid="stSidebarNav"] {display: none;}
    .stPageLink { background-color: #f0f2f6; border-radius: 20px; padding: 8px; border: 1px solid #e0e4eb; }
    </style>
    """, unsafe_allow_html=True)

# --- SIDEBAR: CONFIGURAÇÕES ---
with st.sidebar:
    st.header("🔍 Filtros e Busca")
    st.session_state['search_term'] = st.text_input("Localizar Gene ou TF:", value=st.session_state['search_term']).strip().upper()
    
    st.divider()
    st.header("📂 Entrada de Dados")
    uploaded_files = st.file_uploader("Upload CSVs (JASPAR/TRRUST)", type=['csv'], accept_multiple_files=True)
    
    st.divider()
    st.header("📐 Área do Grafo")
    g_width = st.slider("Largura:", 600, 1800, 1100)
    g_height = st.slider("Altura:", 400, 1200, 800)
    
    st.divider()
    st.markdown("### 🚀 Navegação")
    st.page_link("Lumos_Home.py", label="Início", icon="🏠")
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

# --- RECUPERAÇÃO DE DADOS DO MÓDULO APP ---
k_res = st.session_state.get('k_res')   # KEGG
g_res = st.session_state.get('g_res')   # GO
s_df = st.session_state.get('string_df') # STRING

st.title("PrioriGraph 👑")

# Painel de Diagnóstico de Integração
with st.expander("🔍 Status da Integração Biológica", expanded=False):
    c1, c2, c3 = st.columns(3)
    c1.metric("KEGG Pathways", len(k_res) if k_res is not None else "0")
    c2.metric("GO Processes", len(g_res) if g_res is not None else "0")
    c3.metric("STRING Phys.", "Ativo" if s_df is not None else "Inativo")

# --- MOTOR DE IMPACTO (A LOGICA ÚNICA) ---
def calculate_integrated_impact(tf_name, df_int):
    targets = set(df_int[df_int['TF'] == tf_name]['Target'])
    if not targets: return 1
    
    score = 0
    # Bônus KEGG (Vias)
    if k_res is not None:
        for _, row in k_res.iterrows():
            if targets.intersection(set(str(row['Genes']).split(';'))): score += 3
    # Bônus GO (Processos)
    if g_res is not None:
        for _, row in g_res.iterrows():
            if targets.intersection(set(str(row['Genes']).split(';'))): score += 1
    # Bônus STRING (Físico)
    if s_df is not None:
        string_genes = set(s_df['preferredName_A']).union(set(s_df['preferredName_B']))
        if tf_name in string_genes: score += 5
        
    return max(1, score)

# --- PROCESSAMENTO ---
if uploaded_files:
    all_data = []
    for f in uploaded_files:
        temp = pd.read_csv(f)
        c_tf = next((c for c in ['TF_Symbol', 'TF_Name', 'TF'] if c in temp.columns), None)
        c_tg = next((c for c in ['Genes', 'Target_Gene', 'Target'] if c in temp.columns), None)
        
        if c_tf and c_tg:
            temp = temp.assign(Target=temp[c_tg].astype(str).str.split(';')).explode('Target')
            temp = temp[[c_tf, 'Target']].rename(columns={c_tf: 'TF'})
            all_data.append(temp)

    df_final = pd.concat(all_data).drop_duplicates()
    df_final['TF'] = df_final['TF'].str.upper().str.strip()
    df_final['Target'] = df_final['Target'].str.upper().str.strip()
    tfs_list = set(df_final['TF'].unique())

    # Mapa de Pesos
    impact_map = {tf: calculate_integrated_impact(tf, df_final) for tf in tfs_list}

    # Construção da Rede
    nodes, edges, added = [], [], set()
    s_query = st.session_state['search_term']

    for _, row in df_final.iterrows():
        u, v = row['TF'], row['Target']
        
        if u not in added:
            imp = impact_map.get(u, 1)
            # Dourado se impacto for alto (Mestre Regulador)
            is_gold = imp > 10 or u == s_query
            color = "#FFD700" if is_gold else "#FF4B4B"
            nodes.append(Node(id=u, label=f"{u} (Imp: {imp})", size=25+(imp*1.5), color=color, shape="diamond"))
            added.add(u)
            
        if v not in added:
            is_search = (v == s_query)
            nodes.append(Node(id=v, label=v, size=45 if is_search else 15, color="#FFD700" if is_search else "#1C83E1"))
            added.add(v)
            
        edges.append(Edge(source=u, target=v, directed=True, color="#999"))

    # CONFIGURAÇÃO: physics=False para travar o movimento adoidado
    config = Config(width=g_width, height=g_height, directed=True, physics=False, hierarchical=False)
    agraph(nodes=nodes, edges=edges, config=config)

    # --- RANKINGS ORIGINAIS RESTAURADOS ---
    st.divider()
    col_r1, col_r2 = st.columns(2)
    with col_r1:
        st.subheader("🏆 Master Regulators (Out-Degree)")
        st.dataframe(df_final['TF'].value_counts().reset_index().rename(columns={'count':'Nº de Alvos'}), use_container_width=True)
    with col_r2:
        st.subheader("🎯 Genetic Hubs (In-Degree)")
        st.dataframe(df_final['Target'].value_counts().reset_index().rename(columns={'count':'Reguladores'}), use_container_width=True)

else:
    st.info("Aguardando upload dos arquivos CSV regulatórios.")
