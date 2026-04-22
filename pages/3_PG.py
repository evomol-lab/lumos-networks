import streamlit as st
import pandas as pd
import numpy as np
from streamlit_agraph import agraph, Node, Edge, Config

# Configuração ÚNICA
st.set_page_config(layout="wide", page_title="PrioriGraph  👑")
# Configuração da página (ajuste o título para cada módulo)
st.set_page_config(page_title="Lumos Networks | REDES", page_icon="🕸️", layout="wide")

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) # Pasta 'pages'
PARENT_DIR = os.path.dirname(BASE_DIR) # Raiz do projeto
LOGO_PATH = os.path.join(PARENT_DIR, "assets", "Lumos Networks.png")

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
    st.info("Você está no módulo de REDES.")


# 1. Localização atual: /code/src/pages/seu_script.py
FILE_DIR = os.path.dirname(os.path.abspath(__file__))

# 2. Subir UM nível para chegar na pasta 'src'
PARENT_DIR = os.path.dirname(FILE_DIR)

# 3. Apontar para o arquivo que está solto na 'src'
logo_path = os.path.join(PARENT_DIR, "assets", "PG.png")

with st.sidebar:
    if os.path.exists(logo_path):
        st.image(logo_path, use_container_width=True)


# 1. Configuração da página (ajuste o título para cada módulo)
st.set_page_config(page_title="Lumos Networks | Análise", page_icon="🧬", layout="wide")
# 2. Recuperação de Dados do APP (Sincronização entre módulos)
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')

# --- SIDEBAR: CONTROLES ---
with st.sidebar:
    st.header("🔍 Busca e Filtros")
    search_query = st.text_input("Localizar Gene ou Fator (TF):", "").strip().upper()
    
    st.divider()
    uploaded_files = st.file_uploader("Upload Tabelas JASPAR/TRRUST (CSVs)", type=['csv'], accept_multiple_files=True)
    num_vias = st.slider("Quantidade de Clusters (Vias):", 2, 15, 6)


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

    # --- 2. CONSTRUÇÃO DA REDE ---
    nodes, edges, added = [], [], set()

    for _, row in top_vias.iterrows():
        via_name = row['Term']
        via_genes = set(str(row['Genes']).split(';'))
        target_genes_in_via = via_genes.intersection(set(df_reg['Target']))

        if target_genes_in_via:
            # Nó do Cluster (Via)
            if via_name not in added:
                nodes.append(Node(id=via_name, label=via_name[:35], size=40, color="#FFD700", shape="hexagon"))
                added.add(via_name)

            for gene in target_genes_in_via:
                # Nó do Gene Alvo
                if gene not in added:
                    is_s = (gene == search_query)
                    nodes.append(Node(id=gene, label=gene, size=45 if is_s else 18, 
                                      color="#FFD700" if is_s else "#1C83E1", shape="dot"))
                    added.add(gene)
                
                edges.append(Edge(source=via_name, target=gene, color="#FFD700", width=1, dashed=True))

                # Conectar TFs reguladores
                rel_tfs = df_reg[df_reg['Target'] == gene]['TF'].unique()
                for tf in rel_tfs:
                    if tf not in added:
                        is_tf_s = (tf == search_query)
                        nodes.append(Node(id=tf, label=tf, size=55 if is_tf_s else 30, 
                                          color="#FFD700" if is_tf_s else "#FF4B4B", shape="diamond"))
                        added.add(tf)
                    edges.append(Edge(source=tf, target=gene, directed=True, color="#999"))

    # Configuração da Área do Grafo (Física Desligada)
    config = Config(width=1100, height=800, directed=True, physics=False, hierarchical=False)
    
    st.subheader(f"🕸️ Rede de Clusters Integrada ({num_vias} Processos)")
    agraph(nodes=nodes, edges=edges, config=config)

    # --- 3. QUADROS DE ANÁLISE E RANKINGS ---
    st.divider()
    t1, t2, t3 = st.tabs(["👑 Master Regulators", "🎯 Genetic Hubs", "🧬 Mapeamento TF-Via"])

    with t1:
        st.subheader("Fatores de Transcrição por Grau de Saída (Out-Degree)")
        m_rank = df_reg['TF'].value_counts().reset_index()
        m_rank.columns = ['Fator de Transcrição', 'Total de Alvos na Rede']
        st.dataframe(m_rank, use_container_width=True, hide_index=True)

    with t2:
        st.subheader("Genes Alvos por Grau de Entrada (In-Degree)")
        h_rank = df_reg['Target'].value_counts().reset_index()
        h_rank.columns = ['Gene Alvo', 'Total de Reguladores']
        st.dataframe(h_rank, use_container_width=True, hide_index=True)

    with t3:
        st.subheader("Análise Hierárquica: Reguladores Principais por Via")
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
                    "Via Biológica / Processo": via,
                    "Regulador Principal": main_tf,
                    "Genes Alvos do TF nesta Via": n_targets_in_via,
                    "Total de Genes da Via na Rede": total_genes_via_rede,
                    "Significância (FDR)": f"{row['Adjusted P-value']:.2e}"
                })
        
        st.table(pd.DataFrame(mapping_rows))

else:
    st.info("💡 Por favor, carregue os arquivos CSV regulatórios. Certifique-se de ter processado os dados no módulo APP para carregar os clusters.")
