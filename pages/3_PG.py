import streamlit as st
import pandas as pd
from streamlit_agraph import agraph, Node, Edge, Config

# 1. Recuperação de Dados do APP (Sincronização)
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')

st.title("PrioriGraph 👑 | Clusters Biológicos")

with st.sidebar:
    st.header("📂 Entrada e Filtros")
    uploaded_files = st.file_uploader("Upload JASPAR/TRRUST", type=['csv'], accept_multiple_files=True)
    num_vias = st.slider("Quantidade de Clusters (Vias):", 2, 10, 4)
    st.divider()
    st.info("A física foi desligada para evitar movimentação excessiva.")

if uploaded_files and (k_res is not None or g_res is not None):
    # Processamento de Dados dos Reguladores
    all_data = []
    for f in uploaded_files:
        temp = pd.read_csv(f)
        c_tf = next((c for c in ['TF_Symbol', 'TF'] if c in temp.columns), None)
        c_tg = next((c for c in ['Genes', 'Target'] if c in temp.columns), None)
        if c_tf and c_tg:
            temp = temp.assign(Target=temp[c_tg].astype(str).str.split(';')).explode('Target')
            all_data.append(temp[[c_tf, 'Target']].rename(columns={c_tf: 'TF'}))
    
    df_reg = pd.concat(all_data).drop_duplicates()
    
    # Seleção das Top Vias (Cluster Centers)
    top_vias = pd.concat([k_res.head(num_vias), g_res.head(num_vias)]).sort_values('Adjusted P-value')

    nodes, edges, added = [], [], set()

    # --- LÓGICA DE CONSTRUÇÃO DOS CLUSTERS ---
    for i, row in top_vias.iterrows():
        via_name = row['Term']
        via_genes = set(str(row['Genes']).split(';'))
        
        # 1. Nó do Cluster (O Processo Biológico)
        if via_name not in added:
            nodes.append(Node(id=via_name, 
                              label=via_name[:40], 
                              size=45, 
                              color="#FFD700", 
                              shape="hexagon"))
            added.add(via_name)

        # 2. Genes que pertencem a este cluster
        # Filtramos apenas genes que aparecem na rede de regulação
        target_genes_in_via = via_genes.intersection(set(df_reg['Target']))
        
        for gene in target_genes_in_via:
            if gene not in added:
                nodes.append(Node(id=gene, label=gene, size=18, color="#1C83E1", shape="dot"))
                added.add(gene)
            
            # Ligação Funcional (Via -> Gene)
            edges.append(Edge(source=via_name, target=gene, color="#FFD700", width=1.5, dashed=True))

            # 3. Encontrar TFs que regulam este gene específico
            tfs_regulators = df_reg[df_reg['Target'] == gene]['TF'].unique()
            for tf in tfs_regulators:
                if tf not in added:
                    nodes.append(Node(id=tf, label=tf, size=30, color="#FF4B4B", shape="diamond"))
                    added.add(tf)
                
                # Ligação Regulatória (TF -> Gene)
                edges.append(Edge(source=tf, target=gene, directed=True, color="#999"))

    # CONFIGURAÇÃO CRÍTICA: physics=False para parar o movimento
    config = Config(
        width=1100, 
        height=800, 
        directed=True, 
        physics=False,  # <--- ISOLAMENTO DO MOVIMENTO
        hierarchical=False,
        stabilization=True
    )
    
    st.subheader(f"Visualização de {num_vias} Clusters de Vias e Processos")
    agraph(nodes=nodes, edges=edges, config=config)

    # --- TABELA DE APOIO ---
    st.divider()
    st.subheader("📋 Resumo dos Clusters")
    summary = []
    for i, row in top_vias.iterrows():
        v_genes = set(str(row['Genes']).split(';')).intersection(set(df_reg['Target']))
        summary.append({"Processo": row['Term'], "Genes na Rede": len(v_genes), "FDR": row['Adjusted P-value']})
    st.table(pd.DataFrame(summary))

else:
    st.warning("⚠️ Dados insuficientes. Certifique-se de carregar os CSVs e ter rodado o APP primeiro.")
