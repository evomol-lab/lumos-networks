import streamlit as st
import pandas as pd
import polars as pl  # ADICIONADO
import numpy as np
import plotly.express as px
import os
import gc  # ADICIONADO para gestão de memória
from streamlit_agraph import agraph, Node, Edge, Config

st.set_page_config(layout="wide", page_title="PrioriGraph")

# 1. Encontrar o diretório base do projeto
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 2. Construir o caminho para o logo específico
logo_path = os.path.join(BASE_DIR, "assets", "logos", "PG.png")

with st.sidebar:
    if os.path.exists(logo_path):
        st.image(logo_path, use_container_width=True)
    else:
        st.error(f"Erro: Logo não encontrado em {logo_path}")

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
    all_interactions_pl = []
    plot_data_list = []

    for f in uploaded_files:
        # ATUALIZAÇÃO: Leitura ultrarrápida com Polars
        df_pl = pl.read_csv(f, infer_schema_length=1000)
        
        col_tf = next((c for c in ['TF_Symbol', 'TF_Name'] if c in df_pl.columns), None)
        col_genes = next((c for c in ['Genes', 'Target_Gene'] if c in df_pl.columns), None)
        
        if col_tf and col_genes:
            # Dados para o gráfico (Top TFs)
            if 'Adjusted P-value' in df_pl.columns:
                plot_data_list.append(df_pl.select([col_tf, 'Adjusted P-value']).to_pandas())

            # ATUALIZAÇÃO: Explosão da lista de genes usando Polars (muito mais eficiente em RAM)
            # Transforma "GENE1;GENE2" em linhas separadas
            df_inter = (
                df_pl.select([col_tf, col_genes])
                .with_columns(pl.col(col_genes).str.split(";"))
                .explode(col_genes)
                .rename({col_tf: "TF", col_genes: "Target"})
            )
            all_interactions_pl.append(df_inter)
    
    if not all_interactions_pl:
        st.error("Nenhuma coluna compatível encontrada (TF_Symbol/Genes).")
        st.stop()
        
    # Consolidação dos dados com Polars
    df_final_pl = pl.concat(all_interactions_pl).unique()
    
    # Limpeza de strings e remoção de nulos
    df_final_pl = df_final_pl.filter(
        (pl.col("TF").is_not_null()) & (pl.col("Target").is_not_null())
    ).with_columns([
        pl.col("TF").cast(pl.Utf8).str.strip_chars().str.to_uppercase(),
        pl.col("Target").cast(pl.Utf8).str.strip_chars().str.to_uppercase()
    ])

    # Convertemos para Pandas apenas o necessário para a visualização
    df_final = df_final_pl.to_pandas()
    del all_interactions_pl, df_final_pl; gc.collect()

    # 1. Gráfico de Interação
    if plot_data_list:
        st.subheader("📊 Significância Regulatória (Top TFs)")
        df_plot = pd.concat(plot_data_list).drop_duplicates()
        df_plot['-log10(FDR)'] = -np.log10(df_plot['Adjusted P-value'].astype(float) + 1e-10)
        df_plot = df_plot.sort_values('-log10(FDR)', ascending=False).head(15)
        fig = px.bar(df_plot, x='-log10(FDR)', y=df_plot.columns[0], orientation='h', 
                     color='Adjusted P-value', color_continuous_scale='Viridis_r')
        st.plotly_chart(fig, use_container_width=True)

    # 2. Rede Agraph (Otimizada)
    nodes, edges, nodes_added = [], [], set()
    tfs_list = set(df_final['TF'].unique())

    # Limitador de segurança para não travar o navegador com redes gigantes
    if len(df_final) > 2000:
        st.warning(f"Rede muito grande ({len(df_final)} conexões). Mostrando apenas os top 50 TFs para evitar travamento.")
        top_tfs = df_final['TF'].value_counts().head(50).index
        df_viz = df_final[df_final['TF'].isin(top_tfs)]
    else:
        df_viz = df_final

    for _, row in df_viz.iterrows():
        tf, target = str(row['TF']), str(row['Target'])
        if tf == 'NAN' or target == 'NAN': continue

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
    config = Config(width=g_width, height=g_height, directed=True, physics=False, hierarchical=False, stabilization=True)
    agraph(nodes=nodes, edges=edges, config=config)

    # 3. Rankings (Cálculo rápido com Pandas já que o volume aqui é menor)
    st.divider()
    col_rank1, col_rank2 = st.columns(2)
    with col_rank1:
        st.subheader("👑 Master Regulators (Out-Degree)")
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
