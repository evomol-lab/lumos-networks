import streamlit as st
import pandas as pd
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(layout="wide", page_title="PrioriGraph Master")

# ============================================================
# FUNÇÕES DE PROCESSAMENTO
# ============================================================

def load_and_combine_degs(files):
    all_dfs = []
    for f in files:
        df = pd.read_csv(f)
        df.columns = [c.strip() for c in df.columns]
        if 'Symbol' in df.columns:
            all_dfs.append(df)
    if not all_dfs: return pd.DataFrame()
    return pd.concat(all_dfs).sort_values('FDR').drop_duplicates('Symbol')

def build_network(files_tfs, df_degs):
    G = nx.DiGraph()
    # Criamos um set de busca rápida para garantir que só conectamos aos nossos DEGs
    valid_genes = set(df_degs['Symbol'].astype(str).str.strip().str.upper().unique())
    
    for f in files_tfs:
        df_i = pd.read_csv(f)
        df_i.columns = [c.strip() for c in df_i.columns]
        
        # O App 2/Enrichr geralmente usa 'Term' para TF e 'Genes' para a lista
        tf_col = next((c for c in ['Transcription Factor', 'Term', 'Transcription Factors'] if c in df_i.columns), None)
        gene_col = next((c for c in ['Genes', 'Target', 'Gene'] if c in df_i.columns), None)
        
        if tf_col and gene_col:
            for _, row in df_i.iterrows():
                # 1. Limpeza do TF (exatamente como no R)
                tf_raw = str(row[tf_col]).split(' (')[0].split('_')[0].strip().upper()
                
                # 2. TRATAMENTO CRÍTICO: Quebrar a célula de genes do App 2
                # O Enrichr usa ';' ou ',' dependendo da versão
                raw_text = str(row[gene_col])
                if ';' in raw_text:
                    targets = raw_text.split(';')
                else:
                    targets = raw_text.split(',')
                
                for g in targets:
                    g_clean = g.strip().upper()
                    
                    # 3. Filtro de Interseção (Só conta se o gene alvo estiver nos seus DEGs)
                    if g_clean in valid_genes:
                        G.add_edge(tf_raw, g_clean)
                        
    G.remove_edges_from(nx.selfloop_edges(G))
    return G

# ============================================================
# INTERFACE
# ============================================================

st.title("PrioriGraph Master: Bioinformática de Sistemas 🧬")

with st.sidebar:
    st.header("📂 Gerenciador de Arquivos")
    up_degs = st.file_uploader("Upload DEGs CSV(s)", type=['csv'], accept_multiple_files=True)
    up_tfs = st.file_uploader("Upload TFs/JASPAR CSV(s)", type=['csv'], accept_multiple_files=True)
    
    st.divider()
    st.header("⚙️ Ajustes")
    layout_type = st.selectbox("Layout da Rede:", ["Spring", "Kamada-Kawai", "Circular"])
    top_n = st.slider("Quantidade de Hubs/MRs nos gráficos:", 5, 50, 10)

if up_degs and up_tfs:
    df_degs = load_and_combine_degs(up_degs)
    G = build_network(up_tfs, df_degs)

    if len(G.nodes) == 0:
        st.error("Erro: Nenhuma conexão encontrada. Verifique se os nomes dos genes no App 1 e App 2 coincidem.")
        st.stop()

    # Cálculo de Métricas
    metrics = pd.DataFrame({
        'Node': list(G.nodes()),
        'Out_Degree': [G.out_degree(n) for n in G.nodes()],
        'In_Degree': [G.in_degree(n) for n in G.nodes()],
        'Total_Degree': [G.degree(n) for n in G.nodes()]
    })
    
    # Define Tipos para cores (igual ao seu código R)
    metrics['Type'] = 'Gene'
    tfs_list = [n for n, d in G.out_degree() if d > 0]
    metrics.loc[metrics['Node'].isin(tfs_list), 'Type'] = 'TF'
    
    mr_list = metrics[metrics['Type'] == 'TF'].nlargest(top_n, 'Out_Degree')['Node'].tolist()
    hub_list = metrics.nlargest(top_n, 'Total_Degree')['Node'].tolist()

    # --- TABELAS E GRÁFICOS (ESTILO RSTUDIO) ---
    tab1, tab2 = st.tabs(["🕸️ Network Visualization", "📊 Hubs & Master Regulators"])

    with tab1:
        st.subheader("Interactive Gene Regulation Network")
        
        # Escolha do Layout
        if layout_type == "Spring": pos = nx.spring_layout(G, k=0.2)
        elif layout_type == "Circular": pos = nx.circular_layout(G)
        else: pos = nx.kamada_kawai_layout(G)

        # 1. Desenhar Arestas (Linhas)
        edge_traces = []
        for edge in G.edges():
            x0, y0 = pos[edge[0]]; x1, y1 = pos[edge[1]]
            edge_traces.append(go.Scatter(x=[x0, x1, None], y=[y0, y1, None], 
                                          line=dict(width=0.4, color='#bdbdbd'), mode='lines', hoverinfo='none', showlegend=False))

        # 2. Lógica de Cores e Categorias para a Legenda
        categories = {
            "Master Regulator": {"color": "#984ea3", "size": 28},
            "Hub Gene": {"color": "#FFD700", "size": 22},
            "Transcription Factor": {"color": "#377eb8", "size": 15},
            "Target Gene": {"color": "#e41a1c", "size": 10}
        }

        # Criar traços invisíveis apenas para a LEGENDA
        legend_traces = []
        for cat, style in categories.items():
            legend_traces.append(go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(size=style['size'], color=style['color']),
                legendgroup=cat, showlegend=True, name=cat
            ))

        # 3. Desenhar os Nós Reais
        node_x, node_y, node_color, node_size, node_text = [], [], [], [], []
        for node in G.nodes():
            x, y = pos[node]; node_x.append(x); node_y.append(y)
            
            if node in mr_list:
                node_color.append(categories["Master Regulator"]["color"])
                node_size.append(categories["Master Regulator"]["size"])
            elif node in hub_list:
                node_color.append(categories["Hub Gene"]["color"])
                node_size.append(categories["Hub Gene"]["size"])
            elif node in tfs_list:
                node_color.append(categories["Transcription Factor"]["color"])
                node_size.append(categories["Transcription Factor"]["size"])
            else:
                node_color.append(categories["Target Gene"]["color"])
                node_size.append(categories["Target Gene"]["size"])
            node_text.append(node)

        node_trace = go.Scatter(
            x=node_x, y=node_y, mode='markers+text',
            text=[n if n in mr_list or n in hub_list else '' for n in G.nodes()],
            textposition="top center",
            marker=dict(color=node_color, size=node_size, line_width=1.5, line_color='white'),
            hoverinfo='text', hovertext=node_text, showlegend=False
        )

        # Combinar tudo no gráfico
        fig_net = go.Figure(data=edge_traces + legend_traces + [node_trace])
        
        fig_net.update_layout(
            legend=dict(title="Legenda da Rede:", orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            template='simple_white', height=750,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
        )
        
        st.plotly_chart(fig_net, use_container_width=True)

    with tab2:
        col1, col2 = st.columns(2)
        
        # 1. Gráfico Master Regulators (Out-degree)
        with col1:
            st.markdown("### 🏆 Master Regulators")
            # Filtrar TFs e ordenar
            df_mr = metrics[metrics['Type'] == 'TF'].nlargest(top_n, 'Out_Degree').sort_values('Out_Degree')
            
            fig_mr = px.bar(
                df_mr, 
                x='Out_Degree', 
                y='Node', 
                orientation='h',
                title=f"Top {top_n} MRs (by Out-Degree)",
                color_discrete_sequence=['#377eb8'], # Azul para TFs
                labels={'Node': 'Transcription Factor', 'Out_Degree': 'Grau de Saída'}
            )
            fig_mr.update_layout(template="simple_white", yaxis={'categoryorder':'total ascending'})
            st.plotly_chart(fig_mr, use_container_width=True)

        # 2. Gráfico Hub Genes (Total-degree)
        with col2:
            st.markdown("### 🎯 Hub Nodes")
            # Pega os principais hubs idependente do tipo
            df_hub = metrics.nlargest(top_n, 'Total_Degree').sort_values('Total_Degree')
            
            # Mapeamento de cores para os hubs
            color_map_hubs = {"TF": "#377eb8", "Gene": "#e41a1c"}
            
            fig_hub = px.bar(
                df_hub, 
                x='Total_Degree', 
                y='Node', 
                orientation='h',
                title=f"Top {top_n} Hubs (by Total Degree)",
                color='Type', # Define a cor baseada no tipo (Gene ou TF)
                color_discrete_map=color_map_hubs,
                labels={'Node': 'Node Name', 'Total_Degree': 'Grau Total', 'Type': 'Categoria'}
            )
            fig_hub.update_layout(template="simple_white", yaxis={'categoryorder':'total ascending'})
            st.plotly_chart(fig_hub, use_container_width=True)
            
        st.divider()
        st.subheader("📋 Centrality Data Table")
        # Exibe a tabela completa com a coluna de tipo para conferência
        st.dataframe(metrics.sort_values('Total_Degree', ascending=False), use_container_width=True)

else:
    st.info("Aguardando upload dos arquivos CSV para gerar a rede e os gráficos de centralidade.")