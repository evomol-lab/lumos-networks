import streamlit as st

st.set_page_config(page_title="Lumos | Documentação", layout="wide", page_icon="📚")

st.title("📚 Documentação Técnica e Fundamentos")
st.markdown("---")

tab1, tab2, tab3 = st.tabs(["🧬 DDEA (Estatística)", "🧪 APP (Funcional)", "🕸️ PG (Sistemas)"])

with tab1:
    st.header("Diagonal Differential Expression Alley (DDEA)")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(r"""
        ### 🔬 Fundamentos Estatísticos
        
        O módulo **DDEA** implementa a análise de expressão diferencial fundamentada em modelos lineares generalizados (**GLM**). 
        A escolha da distribuição **Binomial Negativa (NB)** é essencial para lidar com a natureza discreta dos dados de 
        sequenciamento e a característica de **sobre-dispersão**, onde a variância é superior à média.

        #### 1. Modelagem Estocástica
        A contagem de leituras $K_{ij}$ para o gene $i$ na amostra $j$ é modelada como:

        $$K_{ij} \sim NB(\mu_{ij}, \alpha_i)$$

        Onde:
        * **$\mu_{ij}$**: representa a média esperada das contagens.
        * **$\alpha_i$**: é o parâmetro de dispersão específico do gene, que captura a variabilidade biológica extra-Poisson.

        #### 2. Normalização e Robustez
        Para permitir a comparação entre amostras com diferentes profundidades de sequenciamento, utilizamos o método **Median-of-Ratios**. 
        Este processo é superior à simples normalização por milhão (RPM/TPM) pois:
        
        * **Corrige vieses de composição:** Garante que genes altamente expressos em uma condição não criem falsos negativos nos demais.
        * **Estabilidade:** Utiliza uma amostra de referência virtual baseada na média geométrica, tornando o teste de Wald subsequente muito mais confiável.
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
