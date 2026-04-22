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
# INTERFACE PRINCIPAL LÓGICA DE INTEGRAÇÃO 
# ============================================================
# --- RECUPERAÇÃO DE DADOS DO APP ---
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')
s_df = st.session_state.get('string_df')

st.title("PrioriGraph 👑")
if k_res is not None or g_res is not None:
    st.success("🔗 Integração Ativa: Impacto calculado via KEGG/GO/STRING.")
else:
    st.info("💡 Dica: Rode o módulo APP para ativar o cálculo de Impacto Funcional.")

# --- CARGA DE ARQUIVOS ---
uploaded_files = st.file_uploader("Upload Tabelas JASPAR/TRRUST (CSVs)", type=['csv'], accept_multiple_files=True)

if uploaded_files:
    all_interactions = []
    for f in uploaded_files:
        df_temp = pd.read_csv(f)
        col_tf = next((c for c in ['TF_Symbol', 'TF_Name'] if c in df_temp.columns), None)
        col_genes = next((c for c in ['Genes', 'Target_Gene'] if c in df_temp.columns), None)
        if col_tf and col_genes:
            df_temp[col_genes] = df_temp[col_genes].astype(str)
            df_temp = df_temp.assign(Gene=df_temp[col_genes].str.split(';')).explode('Gene')
            df_temp = df_temp[[col_tf, 'Gene']].rename(columns={col_tf: 'TF', 'Gene': 'Target'})
            all_interactions.append(df_temp)
    
    df_final = pd.concat(all_interactions).drop_duplicates()
    df_final['TF'] = df_final['TF'].astype(str).str.upper()
    df_final['Target'] = df_final['Target'].astype(str).str.upper()
    tfs_list = set(df_final['TF'].unique())

    # --- MOTOR DE IMPACTO ---
    def get_impact_score(tf_name, df_int):
        my_targets = set(df_int[df_int['TF'] == tf_name]['Target'])
        score = 0
        if k_res is not None:
            for _, r in k_res.iterrows():
                if my_targets.intersection(set(str(r['Genes']).split(';'))): score += 2
        if g_res is not None:
            for _, r in g_res.iterrows():
                if my_targets.intersection(set(str(r['Genes']).split(';'))): score += 1
        return max(1, score)

    impact_map = {tf: get_impact_score(tf, df_final) for tf in tfs_list}

    # --- CONSTRUÇÃO DOS NÓS E ARESTAS ---
    nodes, edges, nodes_added = [], [], set()
    for _, row in df_final.iterrows():
        tf, target = row['TF'], row['Target']
        
        # Nó do TF
        if tf not in nodes_added:
            imp = impact_map.get(tf, 1)
            color = "#FFD700" if imp > 8 else "#FF4B4B"
            nodes.append(Node(id=tf, label=f"{tf} (Imp: {imp})", size=25+(imp*2), color=color, shape="diamond"))
            nodes_added.add(tf)

        # Nó do Target
        if target not in nodes_added:
            is_search = search_query == target
            nodes.append(Node(id=target, label=target, size=40 if is_search else 15, color="#FFD700" if is_search else "#1C83E1"))
            nodes_added.add(target)
        
        edges.append(Edge(source=tf, target=target, directed=True, color="#999"))

    # Renderização Estática (Sem movimento louco)
    config = Config(width=g_width, height=g_height, directed=True, physics=False, hierarchical=False)
    agraph(nodes=nodes, edges=edges, config=config)

    # --- RESTORE DAS FUNÇÕES LEGAIS (Rankings) ---
    st.divider()
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("🏆 Master Regulators (Out-Degree)")
        m_rank = df_final['TF'].value_counts().reset_index()
        m_rank.columns = ['TF', 'Genes Regulados']
        st.dataframe(m_rank, use_container_width=True)
    with c2:
        st.subheader("🎯 Genetic Hubs (In-Degree)")
        h_rank = df_final['Target'].value_counts().reset_index()
        h_rank.columns = ['Gene', 'Reguladores']
        st.dataframe(h_rank, use_container_width=True)

else:
    st.info("Aguardando upload de CSVs para gerar a rede.")
