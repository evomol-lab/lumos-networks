import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import gseapy as gp
import requests
from streamlit_agraph import agraph, Node, Edge, Config

st.set_page_config(layout="wide", page_title="Arithmancy Pathway Profiler")

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
    st.info("Você está no módulo funcional.")
    

# 1. Localização atual: /code/src/pages/seu_script.py
FILE_DIR = os.path.dirname(os.path.abspath(__file__))

# 2. Subir UM nível para chegar na pasta 'src'
PARENT_DIR = os.path.dirname(FILE_DIR)

# 3. Apontar para o arquivo que está solto na 'src'
logo_path = os.path.join(PARENT_DIR, "assets", "APP.png")

with st.sidebar:
    if os.path.exists(logo_path):
        st.image(logo_path, use_container_width=True)

# ============================================================
# FUNÇÕES DE PROCESSAMENTO
# ============================================================

@st.cache_data(show_spinner=False)
def run_enrichr(gene_list, gene_sets):
    try:
        clean_genes = list(set([str(g).strip().upper() for g in gene_list if str(g).strip().lower() != 'nan' and g != ""]))
        if not clean_genes: return None
        enr = gp.enrichr(gene_list=clean_genes, gene_sets=gene_sets, organism='human', outdir=None)
        df = enr.results
        if df is None or df.empty: return None
        df = df[~df['Term'].str.contains('Mouse|Mus musculus|Rat|Murine', case=False, na=False)]
        df['-log10(FDR)'] = -np.log10(df['Adjusted P-value'] + 1e-10)
        return df.sort_values('Adjusted P-value')
    except: return None

@st.cache_data(show_spinner=False)
def fetch_string_network(genes, confidence=400, max_nodes=150):
    try:
        url = "https://version-12-0.string-db.org/api/json/network"
        params = {"identifiers": "\r".join(genes[:max_nodes]), "species": 9606, "required_score": confidence}
        resp = requests.post(url, data=params, timeout=30)
        return pd.DataFrame(resp.json()) if resp.status_code == 200 else pd.DataFrame()
    except: return pd.DataFrame()

# ============================================================
# INTERFACE PRINCIPAL
# ============================================================

st.title("Arithmancy Pathway Profiler 🕸️🍲")

with st.sidebar:
    st.header("1. Entrada de Dados")
    uploaded_files = st.file_uploader("Upload CSVs (App 1)", type=['csv'], accept_multiple_files=True)
    st.divider()
    st.header("2. Filtros e Ajustes")
    fdr_cut = st.number_input("FDR Cutoff:", value=0.05, format="%.3f")
    fc_cut = st.number_input("Abs Log2FC Cutoff:", value=1.0, format="%.2f")
    fdr_viz = st.slider("Corte de FDR para Gráficos:", 0.01, 1.0, 1.0)
    
    st.divider()
    st.header("📐 Dimensões do Grafo")
    g_width = st.slider("Largura da Rede:", 600, 1800, 1100, 50)
    g_height = st.slider("Altura da Rede:", 400, 1200, 800, 50)

if uploaded_files:
    all_dfs = [pd.read_csv(f) for f in uploaded_files]
    df_full = pd.concat(all_dfs)
    raw_symbols = df_full['Symbol'].dropna().astype(str).str.strip().str.upper().tolist()
    unique_symbols = sorted(list(set(raw_symbols)))
    
    df_sig = df_full[(df_full['FDR'] < fdr_cut) & (df_full['Log2FC'].abs() >= fc_cut)].copy()
    genes_up = set(df_sig[df_sig['Log2FC'] > 0]['Symbol'].str.upper())
    genes_down = set(df_sig[df_sig['Log2FC'] < 0]['Symbol'].str.upper())

    m1, m2 = st.columns(2)
    m1.metric("Soma Total de Genes", len(raw_symbols))
    m2.metric("Genes Únicos", len(unique_symbols))

    tab1, tab2, tab3 = st.tabs(["📊 Vias (KEGG/GO)", "🕸️ STRING Interativo & Sinergia", "👑 Master Regulators"])

    k_res = run_enrichr(unique_symbols, 'KEGG_2021_Human')
    g_res = run_enrichr(unique_symbols, 'GO_Biological_Process_2023')

    with tab1:
        c1, c2 = st.columns(2)
        if isinstance(k_res, pd.DataFrame):
            with c1: st.plotly_chart(px.bar(k_res.head(15), x='-log10(FDR)', y='Term', orientation='h', title="KEGG", color='Adjusted P-value', color_continuous_scale='Reds_r'), use_container_width=True)
        if isinstance(g_res, pd.DataFrame):
            with c2: 
                g_res['Count'] = g_res['Overlap'].apply(lambda x: int(x.split('/')[0]))
                st.plotly_chart(px.scatter(g_res.head(15), x='-log10(FDR)', y='Term', size='Count', title="GO BP", color='Adjusted P-value', color_continuous_scale='Blues_r'), use_container_width=True)

    with tab2:
        st.subheader("Rede de Interações Físicas (Análise de Tooltips)")
        st.caption("Passe o mouse sobre as LINHAS para descobrir os processos compartilhados.")
        c_p, c_m = st.columns([1, 4])
        conf = c_p.selectbox("Confiança do STRING:", [150, 400, 700, 900], index=1)
        max_n = c_p.number_input("Máximo de Genes na Rede:", value=150)
        
        if c_p.button("🌐 Gerar Rede Master"):
            top_genes_list = df_sig.sort_values(by='Log2FC', key=abs, ascending=False).head(max_n)['Symbol'].tolist()
            string_df = fetch_string_network(top_genes_list, conf, max_n)
            
            if not string_df.empty:
                nodes, edges, nodes_added = [], [], set()
                
                # Selecionar Top 25 Vias/Processos
                top_vias_list = []
                if isinstance(k_res, pd.DataFrame): top_vias_list.extend(k_res.head(12)['Term'].tolist())
                if isinstance(g_res, pd.DataFrame): top_vias_list.extend(g_res.head(13)['Term'].tolist())

                # Montar Nós e Arestas
                for _, row in string_df.iterrows():
                    g1, g2 = row['preferredName_A'], row['preferredName_B']
                    
                    # LOGICA DE TOOLTIP PARA LINHAS:
                    shared_paths = []
                    if isinstance(k_res, pd.DataFrame):
                        k_m = k_res[k_res['Term'].isin(top_vias_list) & k_res['Genes'].apply(lambda x: g1 in str(x) and g2 in str(x))]
                        shared_paths.extend(k_m['Term'].tolist())
                    if isinstance(g_res, pd.DataFrame):
                        g_m = g_res[g_res['Term'].isin(top_vias_list) & g_res['Genes'].apply(lambda x: g1 in str(x) and g2 in str(x))]
                        shared_paths.extend(g_m['Term'].tolist())
                    
                    tooltip_text = "Interação: " + (" | ".join(list(set(shared_paths))) if shared_paths else "Física")
                    
                    for g in [g1, g2]:
                        if g not in nodes_added:
                            color = '#FF4B4B' if g in genes_up else ('#1C83E1' if g in genes_down else '#D3D3D3')
                            nodes.append(Node(id=g, label=g, size=15, color=color))
                            nodes_added.add(g)
                    
                    # Adiciona a linha com o título que aparece ao passar o mouse
                    edges.append(Edge(source=g1, target=g2, title=tooltip_text))

                config = Config(width=g_width, height=g_height, directed=False, physics=False, hierarchical=False, stabilization=True) 
                with c_m: agraph(nodes=nodes, edges=edges, config=config)
                
                # Tabela de Mapeamento (Conforme aprovado anteriormente)
                st.divider()
                st.subheader("📋 Mapeamento Funcional: Processos vs Genes da Rede")
                mapping_data = []
                current_network_genes = list(nodes_added)
                full_enrichment = pd.concat([k_res.head(12), g_res.head(13)]) if (k_res is not None and g_res is not None) else pd.DataFrame()
                
                if not full_enrichment.empty:
                    for _, row in full_enrichment.iterrows():
                        via_name = row['Term']
                        genes_in_via = [g for g in str(row['Genes']).split(';') if g in current_network_genes]
                        if genes_in_via:
                            mapping_data.append({
                                "Processo / Via Biológica (Top 25)": via_name,
                                "Genes Relacionados na Rede": ", ".join(genes_in_via),
                                "Contagem": len(genes_in_via),
                                "Significância (FDR)": f"{row['Adjusted P-value']:.2e}"
                            })
                st.table(pd.DataFrame(mapping_data))
            else:
                st.info("Nenhuma interação detectada.")

    with tab3:
        st.subheader("Master Regulators (JASPAR/TRRUST)")
        res_tf = None
        for db in ['JASPAR_2022', 'JASPAR_2024', 'TRRUST_Transcription_Factors_2019']:
            if res_tf is None: res_tf = run_enrichr(unique_symbols, db)
        
        if isinstance(res_tf, pd.DataFrame):
            res_tf['TF_Symbol'] = res_tf['Term'].apply(lambda x: str(x).split(' (')[0].split('_')[0].upper())
            df_final = res_tf[res_tf['Adjusted P-value'] <= fdr_viz].copy()
            st.plotly_chart(px.bar(df_final.head(20), x='-log10(FDR)', y='TF_Symbol', orientation='h', color='Adjusted P-value', color_continuous_scale='Viridis_r'), use_container_width=True)
            export_df = df_final[['TF_Symbol', 'Genes', 'Adjusted P-value', 'Overlap']].copy()
            st.dataframe(export_df, use_container_width=True)
            st.download_button("📥 Baixar CSV para App 3", export_df.to_csv(index=False).encode('utf-8'), "tabela_regulacao_lumos.csv")
        else:
            st.error("Falha ao carregar reguladores.")
else:
    st.info("Aguardando CSVs.")
