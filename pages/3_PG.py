import streamlit as st
import pandas as pd
from streamlit_agraph import agraph, Node, Edge, Config

# 1. Recuperação de Dados
k_res = st.session_state.get('k_res')
g_res = st.session_state.get('g_res')

st.title("PrioriGraph 👑 | Cluster Edition")

with st.sidebar:
    st.header("📂 Entrada")
    uploaded_files = st.file_uploader("Upload JASPAR/TRRUST", type=['csv'], accept_multiple_files=True)
    num_vias = st.slider("Quantidade de Vias para Agrupar:", 3, 15, 5)

if uploaded_files and (k_res is not None or g_res is not None):
    # Processamento dos Reguladores
    all_data = []
    for f in uploaded_files:
        temp = pd.read_csv(f)
        c_tf = next((c for c in ['TF_Symbol', 'TF'] if c in temp.columns), None)
        c_tg = next((c for c in ['Genes', 'Target'] if c in temp.columns), None)
        if c_tf and c_tg:
            temp = temp.assign(Target=temp[c_tg].astype(str).str.split(';')).explode('Target')
            all_data.append(temp[[c_tf, 'Target']].rename(columns={c_tf: 'TF'}))
    
    df_reg = pd.concat(all_data).drop_duplicates()
    
    # Seleção das Top Vias para criar os Clusters
    top_vias = pd.concat([k_res.head(num_vias), g_res.head(num_vias)]).sort_values('Adjusted P-value')

    nodes, edges, added = [], [], set()

    # --- CRIAR CLUSTERS (Nós de Via/Processo) ---
    for i, row in top_vias.iterrows():
        via_id = row['Term']
        via_genes = set(str(row['Genes']).split(';'))
        
        # Nó da Via (O Sol do Cluster)
        if via_id not in added:
            nodes.append(Node(id=via_id, label=via_id[:30]+"...", size=40, 
                              color="#FFD700", shape="hexagon", shadow=True))
            added.add(via_id)

        # Conectar Genes da Via ao Centro do Cluster
        # Filtramos para mostrar apenas genes que também estão na sua rede de TFs
        genes_na_rede = via_genes.intersection(set(df_reg['Target']))
        
        for gene in genes_na_rede:
            if gene not in added:
                nodes.append(Node(id=gene, label=gene, size=15, color="#1C83E1", shape="dot"))
                added.add(gene)
            # Aresta: Via -> Gene (Mostra a pertinência funcional)
            edges.append(Edge(source=via_id, target=gene, color="#FFD700", width=2, label="pertence"))

    # --- CONECTAR TFs AOS CLUSTERS ---
    for _, row in df_reg.iterrows():
        tf, target = row['TF'].upper(), row['Target'].upper()
        
        if target in added: # Só desenha TFs que regulam genes das vias selecionadas
            if tf not in added:
                nodes.append(Node(id=tf, label=tf, size=30, color="#FF4B4B", shape="diamond"))
                added.add(tf)
            
            # Aresta: TF -> Gene (Ação regulatória)
            edges.append(Edge(source=tf, target=target, directed=True, color="#999"))

    # Configuração de Layout
    # Usamos hierarchical=False para permitir que os clusters se formem naturalmente
    config = Config(width=1100, height=800, directed=True, physics=True, 
                    stabilization=True, hierarchical=False)
    
    st.info(f"Rede clusterizada em torno das {num_vias} principais vias detectadas no APP.")
    agraph(nodes=nodes, edges=edges, config=config)

else:
    st.warning("⚠️ Para ver os clusters, você precisa: 1. Rodar o APP e 2. Subir os CSVs aqui.")
