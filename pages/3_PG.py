import streamlit as st
import pandas as pd
import numpy as np
from streamlit_agraph import agraph, Node, Edge, Config

# 1. Configuração ÚNICA
st.set_page_config(page_title="Lumos Networks | PrioriGraph", page_icon="🕸️", layout="wide")

# 2. Inicialização de Estado
if 'search_term' not in st.session_state: st.session_state['search_term'] = ""

# 3. Recuperação de Dados do APP
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')
s_df = st.session_state.get('string_df')

# --- SIDEBAR ---
with st.sidebar:
    st.header("🧬 Controle da Rede")
    st.session_state['search_term'] = st.text_input("Buscar Gene/TF:", value=st.session_state['search_term']).upper()
    
    st.divider()
    uploaded_files = st.file_uploader("Upload CSVs Regulatórios", type=['csv'], accept_multiple_files=True)
    
    # NOVO: Filtro para limpar a estranheza visual
    top_n = st.slider("Mostrar Top X Reguladores (por impacto):", 5, 50, 15)
    
    st.divider()
    st.markdown("### 🚀 Navegação")
    st.page_link("Lumos_Home.py", label="Início", icon="🏠")
    st.page_link("pages/2_APP.py", label="APP (Funcional)", icon="🧪")
    st.page_link("pages/3_PG.py", label="PG (Redes)", icon="🕸️")

st.title("PrioriGraph 👑")

if uploaded_files:
    # --- PROCESSAMENTO ---
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

    # --- MOTOR DE IMPACTO ---
    def get_score(tf_name):
        targets = set(df_final[df_final['TF'] == tf_name]['Target'])
        score = 1
        if k_res is not None:
            for _, r in k_res.iterrows():
                if targets.intersection(set(str(r['Genes']).split(';'))): score += 5
        if g_res is not None:
            for _, r in g_res.iterrows():
                if targets.intersection(set(str(r['Genes']).split(';'))): score += 2
        return score

    tfs_list = df_final['TF'].unique()
    impact_map = {tf: get_score(tf) for tf in tfs_list}
    
    # Filtra os Top N TFs para a rede ficar limpa
    top_tfs = sorted(impact_map, key=impact_map.get, reverse=True)[:top_n]
    df_filtered = df_final[df_final['TF'].isin(top_tfs)]

    # --- CONSTRUÇÃO DA REDE ---
    nodes, edges, added = [], [], set()
    s_query = st.session_state['search_term']

    for _, row in df_filtered.iterrows():
        u, v = row['TF'], row['Target']
        
        if u not in added:
            imp = impact_map.get(u, 1)
            is_master = imp > 15 or u == s_query
            nodes.append(Node(id=u, label=f"{u} (Impacto: {imp})", 
                              size=30 + (imp * 1.5), 
                              color="#FFD700" if is_master else "#FF4B4B", 
                              shape="diamond", shadow=True))
            added.add(u)
            
        if v not in added:
            is_target_match = (v == s_query)
            nodes.append(Node(id=v, label=v, size=20 if is_target_match else 12, 
                              color="#FFD700" if is_target_match else "#1C83E1", 
                              shape="dot"))
            added.add(v)
            
        edges.append(Edge(source=u, target=v, directed=True, color="#CCCCCC", arrowStrikethrough=False))

    # Configuração Estática e Organizada
    config = Config(width=1100, height=800, directed=True, physics=False, 
                    hierarchical=False, borderWidth=2)
    
    st.info(f"Exibindo os {top_n} principais reguladores integrados com KEGG/GO.")
    agraph(nodes=nodes, edges=edges, config=config)

    # --- RANKINGS ---
    st.divider()
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("🏆 Reguladores por Impacto Biológico")
        rank_df = pd.DataFrame(impact_map.items(), columns=['TF', 'Carga Funcional']).sort_values('Carga Funcional', ascending=False)
        st.dataframe(rank_df, use_container_width=True, hide_index=True)
    with c2:
        st.subheader("🎯 Genes Alvos mais Regulados")
        st.dataframe(df_filtered['Target'].value_counts().reset_index(), use_container_width=True)

else:
    st.info("Aguardando os CSVs para gerar a PrioriGraph.")
