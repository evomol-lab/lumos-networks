import streamlit as st
import pandas as pd
import numpy as np
from streamlit_agraph import agraph, Node, Edge, Config

st.set_page_config(layout="wide", page_title="PrioriGraph")

import os
import streamlit as st

# 1. Configuração da página 
st.set_page_config(page_title="Lumos Networks | Networks", page_icon="🧬", layout="wide")

# 2. CSS para manter o padrão visual
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
    st.markdown("### 🚀 Navigation")
    
    # IMPORTANTE: Caminhos saindo da pasta 'pages'
    st.page_link("Lumos_Home.py", label="Home Page", icon="🏠")
    
    st.markdown('<p style="color:#2E86C1; font-weight:bold; margin-bottom:0px; margin-top:10px;">📊 Analysis</p>', unsafe_allow_html=True)
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    
    st.markdown('<p style="color:#28B463; font-weight:bold; margin-bottom:0px; margin-top:10px;">🧬 Functional</p>', unsafe_allow_html=True)
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    
    st.markdown('<p style="color:#E67E22; font-weight:bold; margin-bottom:0px; margin-top:10px;">🕸️ Networks</p>', unsafe_allow_html=True)
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

    st.markdown('<p style="color:#E1AF12; font-weight:bold; margin-bottom:0px; margin-top:10px;">📚 Documentation</p>', unsafe_allow_html=True)
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")

    st.divider()
    st.info("You are in the networking module.")


# 1. Localização atual: /code/src/pages/seu_script.py
FILE_DIR = os.path.dirname(os.path.abspath(__file__))

# 2. Subir UM nível
PARENT_DIR = os.path.dirname(FILE_DIR)

# 3. Apontar para o arquivo 
logo_path = os.path.join(PARENT_DIR, "assets", "PG.png")

with st.sidebar:
    if os.path.exists(logo_path):
        st.image(logo_path, use_container_width=True)

# ============================================================
# INTERFACE PRINCIPAL
# ============================================================

st.title("PrioriGraph 👑")
st.markdown("### Integration of Regulatory Networks (TF -> Genes)")

# Recuperação de Dados do APP (Sincronização entre módulos)
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')

# --- SIDEBAR: CONTROLES ---
with st.sidebar:
    st.header("🔍 Search and Filters")
    # A busca agora foca no destaque VERDE
    search_query = st.text_input("Locate Gene or Factor (Highlighted in Green):", "").strip().upper()
    
    st.divider()
    uploaded_files = st.file_uploader("Upload Tables JASPAR/TRRUST (CSVs)", type=['csv'], accept_multiple_files=True)
    num_vias = st.slider("Number of Clusters (Paths):", 2, 15, 6)

# --- LEGENDA VISUAL FIXA ---
st.markdown("""
    <div style="background-color: #f9f9f9; padding: 15px; border-radius: 10px; border: 1px solid #ddd; margin-bottom: 20px;">
        <h4 style="margin-top: 0;">📍 Network Legend</h4>
        <span style="color: #1C83E1; font-weight: bold;">● Blue Circle:</span> Target Gene | 
        <span style="color: #FF4B4B; font-weight: bold;">■ Red Square:</span> Transcription Factor | 
        <span style="color: #FFD700; font-weight: bold;">⬢ Yellow Hexagon:</span> Pathways/Processes | 
        <span style="color: #28B463; font-weight: bold;">★ Green Highlight:</span> Searched Item
    </div>
    """, unsafe_allow_html=True)

if uploaded_files and (k_res is not None or g_res is not None):
    # --- 1. PROCESSAMENTO DE DADOS ---
    all_data = []
    for f in uploaded_files:
        temp = pd.read_csv(f)
        c_tf = next((c for c in ['TF_Symbol', 'TF_Name', 'TF'] if c in temp.columns), None)
        c_tg = next((c for c in ['Genes', 'Target_Gene', 'Target'] if c in temp.columns), None)
        if c_tf and c_tg:
            # Explode genes caso estejam separados por ';'
            temp = temp.assign(Target=temp[c_tg].astype(str).str.split(';')).explode('Target')
            all_data.append(temp[[c_tf, 'Target']].rename(columns={c_tf: 'TF'}))
    
    df_reg = pd.concat(all_data).drop_duplicates()
    df_reg['TF'] = df_reg['TF'].str.upper().str.strip()
    df_reg['Target'] = df_reg['Target'].str.upper().str.strip()

    # Seleção das Top Vias/Processos para os Clusters
    top_vias = pd.concat([k_res.head(num_vias), g_res.head(num_vias)]).sort_values('Adjusted P-value')

# --- CONSTRUÇÃO DA REDE ---
    nodes, edges, added = [], [], set()

    for _, row in top_vias.iterrows():
        via_name = row['Term'].strip().upper() # Padroniza para busca
        via_genes = set([g.strip().upper() for g in str(row['Genes']).split(';')])
        target_genes_in_via = via_genes.intersection(set(df_reg['Target']))

        if target_genes_in_via:
            # 1. Nó de Via: Hexágono Amarelo (Destaque Verde se pesquisado)
            if via_name not in added:
                # Verificação de busca blindada
                is_via_search = (search_query != "" and search_query in via_name)
                nodes.append(Node(id=via_name, 
                                  label=via_name[:35], 
                                  size=50 if is_via_search else 45, 
                                  color="#28B463" if is_via_search else "#FFD700", 
                                  shape="hexagon"))
                added.add(via_name)

            for gene in target_genes_in_via:
                gene = gene.strip().upper()
                # 2. Nó de Gene: Bola Azul (Destaque Verde se pesquisado)
                if gene not in added:
                    is_g_search = (search_query != "" and search_query == gene)
                    nodes.append(Node(id=gene, 
                                      label=gene, 
                                      size=45 if is_g_search else 20, 
                                      color="#28B463" if is_g_search else "#1C83E1", 
                                      shape="dot"))
                    added.add(gene)
                
                edges.append(Edge(source=via_name, target=gene, color="#FFD700", width=1, dashed=True))

                # 3. Nó de TF: Quadrado Vermelho (Destaque Verde se pesquisado)
                rel_tfs = df_reg[df_reg['Target'] == gene]['TF'].unique()
                for tf in rel_tfs:
                    tf = tf.strip().upper()
                    if tf not in added:
                        is_tf_search = (search_query != "" and search_query == tf)
                        nodes.append(Node(id=tf, 
                                          label=tf, 
                                          size=55 if is_tf_search else 30, 
                                          color="#28B463" if is_tf_search else "#FF4B4B", 
                                          shape="square"))
                        added.add(tf)
                    edges.append(Edge(source=tf, target=gene, directed=True, color="#999"))

    # Configuração da Área do Grafo (Física Desligada)
    config = Config(width=1100, height=800, directed=True, physics=False, hierarchical=False)
    
    st.subheader(f"🕸️ Integrated Cluster Network ({num_vias} Processes)")
    agraph(nodes=nodes, edges=edges, config=config)

    # --- 3. QUADROS DE ANÁLISE E RANKINGS ---
    st.divider()
    t1, t2, t3 = st.tabs(["👑 Master Regulators", "🎯 Genetic Hubs", "🧬 TF-Via Mapping"])

    with t1:
        st.subheader("Transcription Factors by Out-Degree")
        m_rank = df_reg['TF'].value_counts().reset_index()
        m_rank.columns = ['Transcription Factor', 'Total Targets in the Network']
        st.dataframe(m_rank, use_container_width=True, hide_index=True)

    with t2:
        st.subheader("Target Genes by In-Degree")
        h_rank = df_reg['Target'].value_counts().reset_index()
        h_rank.columns = ['Target Gene', 'Total Number of Regulators']
        st.dataframe(h_rank, use_container_width=True, hide_index=True)

    with t3:
        st.subheader("Hierarchical Analysis: Key Regulators by Pathway.")
        mapping_rows = []
        for _, row in top_vias.iterrows():
            via = row['Term']
            v_genes = set(str(row['Genes']).split(';'))
            
            # Filtra as interações que pertencem a esta via específica
            reg_in_via = df_reg[df_reg['Target'].isin(v_genes)]
            tf_counts = reg_in_via['TF'].value_counts()
            
            if not tf_counts.empty:
                main_tf = tf_counts.index[0]
                n_targets_in_via = tf_counts.values[0]
                total_genes_via_rede = len(reg_in_via['Target'].unique())
                
                mapping_rows.append({
                    "Biological Pathway / Process": via,
                    "Master Regulator": main_tf,
                    "TF Target Genes in This Pathway": n_targets_in_via,
                    "Total number of genes in the pathway compared to the entire network.": total_genes_via_rede,
                    "Significance (FDR)": f"{row['Adjusted P-value']:.2e}"
                })
        
        st.table(pd.DataFrame(mapping_rows))

else:
    st.info("💡 Please upload the regulatory CSV files. Make sure you have processed the data in the APP module before uploading the clusters.")
