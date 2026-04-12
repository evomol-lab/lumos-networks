import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import gseapy as gp
import requests
import networkx as nx

st.set_page_config(layout="wide", page_title="Arithmancy Pathway Profiler")

# ============================================================
# FUNÇÕES DE PROCESSAMENTO
# ============================================================

@st.cache_data(show_spinner=False)
def run_enrichr(gene_list, gene_sets):
    try:
        # Limpeza para padrão Humano: remove nans e coloca em caixa alta
        clean_genes = list(set([str(g).strip().upper() for g in gene_list if str(g).strip().lower() != 'nan']))
        
        enr = gp.enrichr(gene_list=clean_genes, 
                         gene_sets=gene_sets, 
                         organism='human', 
                         outdir=None)
        df = enr.results
        if df is None or df.empty: return None
        
        # FILTRO STRICT HUMAN: Remove termos de outras espécies que possam vazar das bibliotecas
        df = df[~df['Term'].str.contains('Mouse|Mus musculus|Rat|Murine', case=False, na=False)]
        
        df['-log10(FDR)'] = -np.log10(df['Adjusted P-value'] + 1e-10)
        return df.sort_values('Adjusted P-value')
    except Exception as e:
        return str(e)

@st.cache_data(show_spinner=False)
def fetch_string_network(genes, confidence=400, max_nodes=50):
    try:
        url = "https://version-12-0.string-db.org/api/json/network"
        params = {
            "identifiers": "\r".join(genes[:max_nodes]),
            "species": 9606, 
            "required_score": confidence,
            "caller_identity": "PrioriGraph_App"
        }
        resp = requests.post(url, data=params, timeout=30)
        return pd.DataFrame(resp.json()) if resp.status_code == 200 else pd.DataFrame()
    except:
        return pd.DataFrame()

def plot_network(df_string, genes_up, genes_down):
    if df_string is None or not isinstance(df_string, pd.DataFrame) or df_string.empty: return None
    G = nx.Graph()
    for _, row in df_string.iterrows():
        G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])
    
    pos = nx.spring_layout(G, k=0.5)
    edge_x, edge_y = [], []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]; x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None]); edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')
    
    node_x, node_y, node_colors, node_text = [], [], [], []
    for node in G.nodes():
        x, y = pos[node]; node_x.append(x); node_y.append(y); node_text.append(node)
        if node in genes_up: node_colors.append('red')
        elif node in genes_down: node_colors.append('blue')
        else: node_colors.append('grey')

    node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text,
                            textposition="top center", marker=dict(size=15, color=node_colors, line_width=1.5))

    return go.Figure(data=[edge_trace, node_trace], layout=go.Layout(showlegend=False, template='simple_white',
                                                                  xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                                                  yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

# ============================================================
# INTERFACE PRINCIPAL
# ============================================================

st.title("Arithmancy Pathway Profiler 🕸️🍲")

with st.sidebar:
    st.header("1. Upload de DEGs")
    uploaded_files = st.file_uploader("Upload DEGs (CSV)", type=['csv'], accept_multiple_files=True)
    
    st.divider()
    st.header("2. Filtros de Análise")
    fdr_cut = st.number_input("FDR Cutoff:", value=0.05, format="%.3f")
    fc_cut = st.number_input("Abs Log2FC Cutoff:", value=1.0, format="%.2f")
    fdr_viz = st.slider("Corte de FDR para Gráficos:", 0.01, 1.0, 1.0)

if uploaded_files:
    # Consolidação de Dados de múltiplos arquivos
    all_dfs = [pd.read_csv(f) for f in uploaded_files]
    for df in all_dfs: df.columns = [c.strip() for c in df.columns]
    
    df_full = pd.concat(all_dfs)
    
    # SOMA TOTAL DE GENES (Variável que mantém as repetições para a estatística)
    all_symbols = df_full['Symbol'].dropna().astype(str).str.strip().str.upper().tolist()
    
    # Definições para a REDE (Filtragem baseada no FDR e Fold Change)
    df_sig = df_full[(df_full['FDR'] < fdr_cut) & (df_full['Log2FC'].abs() >= fc_cut)].copy()
    genes_up = df_sig[df_sig['Log2FC'] > 0]['Symbol'].str.upper().tolist()
    genes_down = df_sig[df_sig['Log2FC'] < 0]['Symbol'].str.upper().tolist()
    
    col_stat1, col_stat2 = st.columns(2)
    col_stat1.metric("Total de Genes Citados (Soma)", len(all_symbols))
    col_stat2.metric("Genes Únicos", len(set(all_symbols)))

    tab1, tab2, tab3 = st.tabs(["📊 Vias (KEGG / GO)", "🕸️ Rede STRING", "👑 Master Regulators"])

    # --- TAB 1: KEGG e GO ---
    with tab1:
        st.subheader("Enriquecimento Funcional (Strictly Human)")
        c1, c2 = st.columns(2)
        with c1:
            k_res = run_enrichr(all_symbols, 'KEGG_2021_Human')
            if isinstance(k_res, pd.DataFrame):
                st.plotly_chart(px.bar(k_res.head(15), x='-log10(FDR)', y='Term', orientation='h', color='Adjusted P-value', color_continuous_scale='Reds_r'), use_container_width=True)
        with c2:
            g_res = run_enrichr(all_symbols, 'GO_Biological_Process_2023')
            if isinstance(g_res, pd.DataFrame):
                g_res['Count'] = g_res['Overlap'].apply(lambda x: int(x.split('/')[0]))
                st.plotly_chart(px.scatter(g_res.head(15), x='-log10(FDR)', y='Term', size='Count', color='Adjusted P-value', color_continuous_scale='Blues_r'), use_container_width=True)

    # --- TAB 2: REDE STRING ---
    with tab2:
        st.markdown("### Rede de Interação Física (STRING DB)")
        st.markdown("*Nós vermelhos = Up-regulados | Nós azuis = Down-regulados*")
        
        c1, c2 = st.columns([1, 3])
        confidence = c1.selectbox("Grau de Confiança (STRING):", [("Baixa (0.150)", 150), ("Média (0.400)", 400), ("Alta (0.700)", 700), ("Altíssima (0.900)", 900)], index=1)[1]
        max_n = c1.number_input("Máximo de Genes na Rede:", min_value=10, max_value=200, value=50, step=10)
        
        if c1.button("🌐 Gerar Rede", use_container_width=True):
            # AJUSTE AQUI: Trocado all_genes por all_symbols
            with st.spinner(f"Construindo rede com os {min(len(all_symbols), max_n)} genes citados..."):
                top_drastic_genes = df_sig.sort_values(by='Log2FC', key=abs, ascending=False).head(int(max_n))['Symbol'].tolist()
                
                string_df = fetch_string_network(top_drastic_genes, confidence, int(max_n))
                
                if isinstance(string_df, str):
                    c2.error(f"Erro na API do STRING: {string_df}")
                elif string_df.empty:
                    c2.info("Nenhuma interação conhecida encontrada com estes parâmetros.")
                else:
                    fig_net = plot_network(string_df, genes_up, genes_down)
                    c2.plotly_chart(fig_net, use_container_width=True)
                    with st.expander("Ver Tabela de Interações"):
                        st.dataframe(string_df[['preferredName_A', 'preferredName_B', 'score']])

    # --- TAB 3: JASPAR ---
    with tab3:
        st.subheader("Reguladores Mestres (JASPAR Human)")
        
        with st.spinner("Consultando bases de dados de regulação..."):
            # Tentativa Sequencial (JASPAR 2022 -> 2024 -> TRRUST)
            res_tf = run_enrichr(all_symbols, 'JASPAR_2022')
            if res_tf is None or isinstance(res_tf, str):
                res_tf = run_enrichr(all_symbols, 'JASPAR_2024')
            if res_tf is None or isinstance(res_tf, str):
                res_tf = run_enrichr(all_symbols, 'TRRUST_Transcription_Factors_2019')

            if isinstance(res_tf, pd.DataFrame):
                # Limpeza e Formatação
                res_tf['TF'] = res_tf['Term'].apply(lambda x: str(x).split(' (')[0].split('_')[0].upper())
                
                # Aplicação do Filtro do Slider
                df_p = res_tf[res_tf['Adjusted P-value'] <= fdr_viz].head(30)
                
                if not df_p.empty:
                    fig_tf = px.bar(df_p, x='-log10(FDR)', y='TF', orientation='h', 
                                    color='Adjusted P-value', color_continuous_scale='Viridis_r', height=700)
                    fig_tf.update_layout(template='simple_white', yaxis={'categoryorder':'total ascending'})
                    st.plotly_chart(fig_tf, use_container_width=True)
                    
                    st.markdown("### Tabela de Alvos (TF -> Genes da sua lista)")
                    st.dataframe(res_tf[['TF', 'Adjusted P-value', 'Overlap', 'Genes']], use_container_width=True)
                    st.download_button("📥 Baixar Resultados JASPAR", res_tf.to_csv(index=False), "jaspar_results.csv")
                else:
                    st.warning(f"Nenhum fator passou pelo filtro de FDR <= {fdr_viz}. Tente aumentar o slider na barra lateral.")
            else:
                st.error("Falha ao obter dados. Verifique sua conexão com a internet ou se os Gene Symbols são válidos.")

else:
    st.info("Aguardando upload dos arquivos CSV.")