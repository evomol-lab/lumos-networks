import streamlit as st

st.set_page_config(page_title="Lumos | Documentação", layout="wide", page_icon="📚")

st.title("📚 Documentação Técnica e Fundamentos")
st.markdown("---")

tab1, tab2, tab3 = st.tabs(["🧬 DDEA (Estatística)", "🧪 APP (Funcional)", "🕸️ PG (Sistemas)"])

with tab1:
    st.header("Diagonal Differential Expression Alley (DDEA)")
    
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("🔬 Fundamentos")
        st.markdown("""
        O DDEA é baseado no modelo linear generalizado (GLM) da distribuição **Binomial Negativa**, ideal para dados de contagem de RNA-Seq que apresentam sobre-dispersão:
        
        $$K_{ij} \sim NB(\mu_{ij}, \alpha_i)$$
        
        Onde $\mu_{ij}$ é a média e $\alpha_i$ é o parâmetro de dispersão. A normalização é feita por **Median-of-Ratios**, garantindo que a variação na profundidade de sequenciamento não gere falsos positivos.
        """)

    with col2:
        st.subheader("🛠️ Workflow")
        st.info("""
        1. **Pre-processing:** Filtragem de genes com baixa contagem (Low Counts).
        2. **Estimation:** Cálculo dos fatores de tamanho (Size Factors).
        3. **Testing:** Teste de Wald para significância estatística.
        4. **Correction:** Ajuste de P-value via Benjamini-Hochberg (FDR).
        """)

with tab2:
    st.header("Arithmancy Pathway Profiler (APP)")
    
    st.markdown("""
    Este módulo realiza a **Contextualização Funcional** dos genes identificados no DDEA. Ele não apenas lista termos, ele calcula a probabilidade estatística de uma via estar ativa.
    """)
    
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("### 📊 Métodos")
        st.write("- **ORA (Over-Representation Analysis):** Utiliza o teste hipergeométrico para avaliar a intersecção entre sua lista e o banco de dados.")
        st.latex(r"P(X > k) = \sum_{i=k+1}^{n} \frac{\binom{M}{i} \binom{N-M}{n-i}}{\binom{N}{n}}")
        st.write("- **GSEA:** Analisa o ranqueamento completo dos genes sem cortes de threshold.")

    with c2:
        st.markdown("### 🗃️ Bases de Dados")
        st.markdown("""
        - **KEGG/GO:** Vias metabólicas e processos biológicos.
        - **STRING:** Interações proteína-proteína (PPI) com score de evidência.
        - **JASPAR/TRRUST:** Fatores de Transcrição e alvos regulatórios validados.
        """)

with tab3:
    st.header("PrioriGraph (PG)")
    
    st.markdown("""
    O PG integra os resultados ômicos com o conhecimento biológico prévio para modelar a topologia do sistema.
    """)
    
    st.subheader("🕸️ Métricas de Rede (Graph Theory)")
    
    cc1, cc2, cc3 = st.columns(3)
    with cc1:
        st.markdown("**Degree Centrality**")
        st.caption("Identifica os 'Hubs' (genes com maior número de conexões).")
    with cc2:
        st.markdown("**Betweenness**")
        st.caption("Identifica genes que funcionam como 'pontes' entre diferentes módulos biológicos.")
    with cc3:
        st.markdown("**Clustering Coefficient**")
        st.caption("Mede o quão conectada está a vizinhança de um gene.")

    st.warning("⚠️ **Nota de Pesquisa:** O PrioriGraph prioriza interações com suporte experimental (Curated Databases), reduzindo o ruído de predições computacionais puras.")

st.markdown("---")
st.caption("Lumos Networks Documentation - EvoMol-Lab 2026")
