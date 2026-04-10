import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import gseapy as gp
import requests
import networkx as nx

st.set_page_config(layout="wide", page_title="Cauldron Networks Leaky")

# ============================================================
# FUNÇÕES DE ENRIQUECIMENTO (GSEAPY) E REDES (STRING)
# ============================================================

@st.cache_data(show_spinner=False)
def run_enrichr(gene_list, gene_sets):
    try:
        # Acessa a API do Enrichr silenciosamente
        enr = gp.enrichr(gene_list=gene_list, 
                         gene_sets=gene_sets, 
                         organism='human', 
                         outdir=None)
        df = enr.results
        if df.empty: return None
        # Calcula o score de plotagem (-log10 do Adjusted P-value)
        df['-log10(FDR)'] = -np.log10(df['Adjusted P-value'] + 1e-10)
        return df.sort_values('Adjusted P-value')
    except Exception as e:
        return str(e)

@st.cache_data(show_spinner=False)
def fetch_string_network(genes, confidence=400, max_nodes=50):
    try:
        string_api_url = "https://version-12-0.string-db.org/api/json/network"
        params = {
            "identifiers": "\r".join(genes[:max_nodes]), # Limita para não bugar o grafo
            "species": 9606, # Homo sapiens
            "required_score": confidence,
            "caller_identity": "DDEA_Streamlit_App"
        }
        response = requests.post(string_api_url, data=params, timeout=30)
        response.raise_for_status()
        return pd.DataFrame(response.json())
    except Exception as e:
        return str(e)

def plot_network(df_string, genes_up, genes_down):
    if df_string is None or df_string.empty: return None
    
    G = nx.Graph()
    for _, row in df_string.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])
        
    pos = nx.spring_layout(G, k=0.5, iterations=50)
    
    edge_x, edge_y = [], []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')
    
    node_x, node_y, node_colors, node_text = [], [], [], []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(node)
        
        # Pinta o nó baseado se é UP ou DOWN no dataset original
        if node in genes_up: node_colors.append('red')
        elif node in genes_down: node_colors.append('blue')
        else: node_colors.append('grey') # Proteínas adicionadas pela rede (se houver)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=node_text,
        textposition="top center",
        marker=dict(size=20, color=node_colors, line_width=2))

    fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                showlegend=False,
                hovermode='closest',
                margin=dict(b=0, l=0, r=0, t=0),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                plot_bgcolor='white'
             ))
    return fig

# ============================================================
# APP PRINCIPAL
# ============================================================

st.title("Cauldron Networks Leaky 🕸️🍲 ")

st.sidebar.header("1. Upload de DEGs")
st.sidebar.markdown("Faça upload de um ou mais CSVs gerados no DDEA. O app irá combinar os genes únicos de todos os arquivos.")

# Alteração: accept_multiple_files=True
uploaded_files = st.sidebar.file_uploader("Upload DEGs (CSV)", type=['csv'], accept_multiple_files=True)

if uploaded_files:
    dfs = []
    for uploaded_file in uploaded_files:
        df_temp = pd.read_csv(uploaded_file)
        
        # Validação Básica por arquivo
        required_cols = {'Symbol', 'Log2FC', 'FDR'}
        if not required_cols.issubset(df_temp.columns):
            st.error(f"Erro no arquivo {uploaded_file.name}: Deve conter as colunas: {required_cols}")
            continue
        dfs.append(df_temp)

    if not dfs:
        st.stop()

    df_raw = pd.concat(dfs).sort_values('FDR').drop_duplicates('Symbol', keep='first')

    st.sidebar.divider()
    st.sidebar.header("2. Filtros de Análise")
    fdr_cut = st.sidebar.number_input("FDR Cutoff:", value=0.05, format="%.3f")
    fc_cut = st.sidebar.number_input("Abs Log2FC Cutoff:", value=1.0, format="%.2f")
    
    # Filtragem
    df_sig = df_raw[(df_raw['FDR'] < fdr_cut) & (df_raw['Log2FC'].abs() >= fc_cut)].copy()
    df_sig['Direction'] = np.where(df_sig['Log2FC'] > 0, 'UP', 'DOWN')
    
    genes_up = df_sig[df_sig['Direction'] == 'UP']['Symbol'].dropna().tolist()
    genes_down = df_sig[df_sig['Direction'] == 'DOWN']['Symbol'].dropna().tolist()
    all_genes = df_sig['Symbol'].dropna().tolist()

    st.success(f"Carregado com sucesso! Genes Filtrados: **{len(all_genes)}** (UP: {len(genes_up)} | DOWN: {len(genes_down)})")

    if len(all_genes) < 3:
        st.warning("Poucos genes significativos para realizar análise de vias e redes. Ajuste os filtros na barra lateral.")
        st.stop()

    # ============================================================
    # ABAS DE ANÁLISE
    # ============================================================
    tab1, tab2, tab3 = st.tabs(["🧬 Vias (KEGG / GO)", "🕸️ Rede Proteína-Proteína (STRING)", "👑 Master Regulators"])

    # ----------------------------------------------------
    # TAB 1: KEGG e GO
    # ----------------------------------------------------
    with tab1:
        st.markdown("### Enriquecimento Funcional (KEGG 2021 Human)")
        with st.spinner("Consultando servidor Enrichr (KEGG)..."):
            kegg_res = run_enrichr(all_genes, 'KEGG_2021_Human')
            
            if isinstance(kegg_res, str):
                st.error(f"Falha de conexão com a API: {kegg_res}")
            elif kegg_res is None or kegg_res[kegg_res['Adjusted P-value'] < 0.05].empty:
                st.info("Nenhuma via do KEGG apresentou significância estatística (FDR < 0.05).")
            else:
                top_kegg = kegg_res[kegg_res['Adjusted P-value'] < 0.05].head(15)
                fig_kegg = px.bar(top_kegg, x='-log10(FDR)', y='Term', orientation='h', 
                                  color='Adjusted P-value', color_continuous_scale='Reds_r',
                                  title="Top 15 Vias do KEGG Enriquecidas")
                fig_kegg.update_layout(yaxis={'categoryorder':'total ascending'}, template="simple_white")
                st.plotly_chart(fig_kegg, use_container_width=True)

        st.divider()
        st.markdown("### Processos Biológicos (GO Biological Process 2023)")
        with st.spinner("Consultando servidor Enrichr (GO)..."):
            go_res = run_enrichr(all_genes, 'GO_Biological_Process_2023')
            
            if isinstance(go_res, str):
                st.error(f"Falha de conexão com a API: {go_res}")
            elif go_res is None or go_res[go_res['Adjusted P-value'] < 0.05].empty:
                st.info("Nenhum termo do GO apresentou significância estatística (FDR < 0.05).")
            else:
                top_go = go_res[go_res['Adjusted P-value'] < 0.05].head(15).copy()
                
                # Limpa os IDs complicados do GO para o gráfico ficar legível
                top_go['Term'] = top_go['Term'].apply(lambda x: x.split(' (GO:')[0])
                
                # CORREÇÃO: Converte a string "12/200" para o número inteiro 12
                top_go['Overlap_Count'] = top_go['Overlap'].apply(lambda x: int(x.split('/')[0]))
                
                # Usa 'Overlap_Count' no parâmetro 'size'
                fig_go = px.scatter(top_go, x='-log10(FDR)', y='Term', size='Overlap_Count', 
                                    color='Adjusted P-value', color_continuous_scale='Blues_r',
                                    title="Dotplot: Processos Biológicos GO")
                fig_go.update_layout(yaxis={'categoryorder':'total ascending'}, template="simple_white")
                st.plotly_chart(fig_go, use_container_width=True)

    # ----------------------------------------------------
    # TAB 2: STRING PPI NETWORK
    # ----------------------------------------------------
    with tab2:
        st.markdown("### Rede de Interação Física (STRING DB)")
        st.markdown("*Nós vermelhos = Up-regulados | Nós azuis = Down-regulados*")
        
        c1, c2 = st.columns([1, 3])
        confidence = c1.selectbox("Grau de Confiança (STRING):", [("Baixa (0.150)", 150), ("Média (0.400)", 400), ("Alta (0.700)", 700), ("Altíssima (0.900)", 900)], index=1)[1]
        max_n = c1.number_input("Máximo de Genes na Rede:", min_value=10, max_value=200, value=50, step=10)
        
        if c1.button("🌐 Gerar Rede", use_container_width=True):
            with st.spinner(f"Construindo rede com os {min(len(all_genes), max_n)} genes mais significativos..."):
                # Pega os genes mais drásticos pelo Log2FC para a rede não virar uma "bola de pelo" preta inútil
                top_drastic_genes = df_sig.sort_values(by='Log2FC', key=abs, ascending=False).head(int(max_n))['Symbol'].tolist()
                
                string_df = fetch_string_network(top_drastic_genes, confidence, int(max_n))
                
                if isinstance(string_df, str):
                    c2.error(f"Erro na API do STRING: {string_df}")
                elif string_df.empty:
                    c2.info("Nenhuma interação física conhecida encontrada para estes genes neste nível de confiança.")
                else:
                    fig_net = plot_network(string_df, genes_up, genes_down)
                    c2.plotly_chart(fig_net, use_container_width=True)
                    with st.expander("Ver Tabela de Interações (Arestas)"):
                        st.dataframe(string_df[['preferredName_A', 'preferredName_B', 'score']])

    # ----------------------------------------------------
    # TAB 3: MASTER REGULATORS (Fatores de Transcrição)
    # ----------------------------------------------------
    with tab3:
        st.markdown("### Fatores de Transcrição (ChEA 2022)")
        st.markdown("Identifica quais fatores de transcrição estão governando a lista de DEGs identificada no seu projeto.")
        
        with st.spinner("Consultando banco de dados ChEA..."):
            chea_res = run_enrichr(all_genes, 'ChEA_2022')
            
            if isinstance(chea_res, str):
                st.error(f"Falha de conexão com a API: {chea_res}")
            elif chea_res is None or chea_res[chea_res['Adjusted P-value'] < 0.05].empty:
                st.info("Nenhum Fator de Transcrição significativo encontrado.")
            else:
                top_chea = chea_res[chea_res['Adjusted P-value'] < 0.05].head(15)
                fig_chea = px.bar(top_chea, x='-log10(FDR)', y='Term', orientation='h', 
                                  color='Adjusted P-value', color_continuous_scale='Purples_r',
                                  title="Top Master Regulators (ChEA)")
                fig_chea.update_layout(yaxis={'categoryorder':'total ascending'}, template="simple_white")
                st.plotly_chart(fig_chea, use_container_width=True)
                
                st.markdown("**Detalhes de Interseção dos TFs:**")
                st.dataframe(top_chea[['Term', 'Adjusted P-value', 'Overlap', 'Genes']], use_container_width=True)

else:
    st.info("👈 Faça o upload da sua tabela de DEGs na barra lateral para iniciar a análise sistêmica.")
