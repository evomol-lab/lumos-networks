import streamlit as st
import pandas as pd
import numpy as np
from streamlit_agraph import agraph, Node, Edge, Config

# 1. Configuração ÚNICA
st.set_page_config(page_title="Lumos Networks | PrioriGraph", page_icon="🕸️", layout="wide")

# 2. Recuperação de Dados do APP
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')

# --- SIDEBAR: CONTROLES E BUSCA ---
with st.sidebar:
    st.header("🔍 Central de Comando")
    search_query = st.text_input("Localizar Gene ou Fator (TF):", "").strip().upper()
    
    st.divider()
    uploaded_files = st.file_uploader("Upload JASPAR/TRRUST (CSVs)", type=['csv'], accept_multiple_files=True)
    num_vias = st.slider("Clusters de Vias (Top N):", 2, 12, 5)
    
    st.divider()
    st.markdown("### 🚀 Navegação")
    st.page_link("Lumos_Home.py", label="Início", icon="🏠")
    st.page_link("pages/2_APP.py", label="APP (Funcional)", icon="🧪")
    st.page_link("pages/3_PG.py", label="PG (Redes)", icon="🕸️")

st.title("PrioriGraph 👑")

if uploaded_files and (k_res is not None or g_res is not None):
    # --- 1. PROCESSAMENTO DE DADOS ---
    all_data = []
    for f in uploaded_files:
        temp = pd.read_csv(f)
        c_tf = next((c for c in ['TF_Symbol', 'TF'] if c in temp.columns), None)
        c_tg = next((c for c in ['Genes', 'Target'] if c in temp.columns), None)
        if c_tf and c_tg:
            temp = temp.assign(Target=temp[c_tg].astype(str).str.split(';')).explode('Target')
            all_data.append(temp[[c_tf, 'Target']].rename(columns={c_tf: 'TF'}))
    
    df_reg = pd.concat(all_data).drop_duplicates()
    df_reg['TF'] = df_reg['TF'].str.upper().str.strip()
    df_reg['Target'] = df_reg['Target'].str.upper().str.strip()

    # Seleção das Top Vias
    top_vias = pd.concat([k_res.head(num_vias), g_res.head(num_vias)]).sort_values('Adjusted P-value')

    # --- 2. CONSTRUÇÃO DA REDE (CLUSTERIZADA) ---
    nodes, edges, added = [], [], set()

    for _, row in top_vias.iterrows():
        via_name = row['Term']
        via_genes = set(str(row['Genes']).split(';'))
        target_genes_in_via = via_genes.intersection(set(df_reg['Target']))

        if target_genes_in_via:
            # Nó da Via (Hexágono)
            if via_name not in added:
                nodes.append(Node(id=via_name, label=via_name[:35], size=40, color="#FFD700", shape="hexagon"))
                added.add(via_name)

            for gene in target_genes_in_via:
                # Nó do Gene (Círculo)
                if gene not in added:
                    is_search = (gene == search_query)
                    nodes.append(Node(id=gene, label=gene, size=40 if is_search else 18, 
                                      color="#FFD700" if is_search else "#1C83E1", shape="dot"))
                    added.add(gene)
                
                edges.append(Edge(source=via_name, target=gene, color="#FFD700", width=1, dashed=True))

                # Conectar TFs que regulam este gene
                relevant_tfs = df_reg[df_reg['Target'] == gene]['TF'].unique()
                for tf in relevant_tfs:
                    if tf not in added:
                        is_tf_search = (tf == search_query)
                        nodes.append(Node(id=tf, label=tf, size=50 if is_tf_search else 30, 
                                          color="#FFD700" if is_tf_search else "#FF4B4B", shape="diamond"))
                        added.add(tf)
                    edges.append(Edge(source=tf, target=gene, directed=True, color="#999"))

    # Configuração Estática (Sem movimento)
    config = Config(width=1100, height=800, directed=True, physics=False, hierarchical=False)
    
    st.subheader(f"🕸️ Rede de Clusters: {num_vias} Principais Processos")
    agraph(nodes=nodes, edges=edges, config=config)

    # --- 3. TABELAS DE ANÁLISE (O QUE VOCÊ PEDIU) ---
    st.divider()
    tab_m1, tab_m2, tab_m3 = st.tabs(["👑 Master Regulators", "🎯 Genetic Hubs", "🧬 Top TFs por Via"])

    with tab_m1:
        st.subheader("Fatores de Transcrição com Maior Alcance (Out-Degree)")
        m_rank = df_reg['TF'].value_counts().reset_index()
        m_rank.columns = ['Fator de Transcrição', 'Nº de Alvos na Rede']
        st.dataframe(m_rank, use_container_width=True, hide_index=True)

    with tab_m2:
        st.subheader("Genes Alvos com Maior Pressão Regulatória (In-Degree)")
        h_rank = df_reg['Target'].value_counts().reset_index()
        h_rank.columns = ['Gene Alvo', 'Nº de Reguladores']
        st.dataframe(h_rank, use_container_width=True, hide_index=True)

    with tab_m3:
        st.subheader("Mapeamento: TFs que mais regulam cada Via")
        tf_via_data = []
        for _, row in top_vias.iterrows():
            via = row['Term']
            via_genes = set(str(row['Genes']).split(';'))
            # Filtra TFs que regulam genes desta via
            tfs_in_via = df_reg[df_reg['Target'].isin(via_genes)]['TF'].value_counts()
            if not tfs_in_via.empty:
                top_tf_name = tfs_in_via.index[0]
                count = tfs_in_via.values[0]
                tf_via_data.append({
                    "Via Biológica": via,
                    "Principal Regulador": top_tf_name,
                    "Genes Alvos na Via": count,
                    "Significância da Via (FDR)": f"{row['Adjusted P-value']:.2e}"
                })
        st.table(pd.DataFrame(tf_via_data))

else:
    st.warning("⚠️ Aguardando dados: Certifique-se de carregar os CSVs e ter rodado o APP primeiro para gerar os clusters.")
