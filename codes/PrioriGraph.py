import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from streamlit_agraph import agraph, Node, Edge, Config

st.set_page_config(layout="wide", page_title="PrioriGraph")

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
        
    df_final = pd.concat(all_interactions).drop_duplicates()
    df_final['TF'] = df_final['TF'].astype(str).str.strip().str.upper()
    df_final['Target'] = df_final['Target'].astype(str).str.strip().str.upper()

    # 1. Gráfico de Interação (Master Regulators por FDR)
    if plot_data_list:
        st.subheader("📊 Significância Regulatória (Top TFs)")
        df_plot = pd.concat(plot_data_list).drop_duplicates()
        df_plot['-log10(FDR)'] = -np.log10(df_plot['Adjusted P-value'] + 1e-10)
        df_plot = df_plot.sort_values('-log10(FDR)', ascending=False).head(15)
        fig = px.bar(df_plot, x='-log10(FDR)', y=df_plot.columns[0], orientation='h', 
                     color='Adjusted P-value', color_continuous_scale='Viridis_r')
        st.plotly_chart(fig, use_container_width=True)

    # 2. Rede Agraph
    nodes, edges, nodes_added = [], [], set()
    tfs_list = set(df_final['TF'].unique())

    for _, row in df_final.iterrows():
        tf, target = row['TF'], row['Target']
        if tf == 'NAN' or target == 'NAN' or not tf or not target: continue

        for n in [tf, target]:
            if n not in nodes_added:
                color = "#FF4B4B" if n in tfs_list else "#1C83E1"
                size, stroke_w, stroke_c = 20, 1, "#333"
                if search_query and search_query == n:
                    color, size, stroke_w, stroke_c = "#FFD700", 45, 3, "#000"
                nodes.append(Node(id=n, label=n, size=size, color=color, strokeWidth=stroke_w, strokeColor=stroke_c))
                nodes_added.add(n)
        edges.append(Edge(source=tf, target=target, directed=True, color="#999"))

    st.subheader("🕸️ Rede Regulatória Integrada")
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
    
