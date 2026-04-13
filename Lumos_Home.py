import streamlit as st
import os

# Configuração da página
st.set_page_config(
    page_title="Lumos Networks | Home",
    page_icon="🧬",
    layout="wide"
)

# --- NOVO CABEÇALHO MÁGICO ---

# 1. Centralizar o Logo Principal (Opcional, mas fica lindo)
logo_principal = "assets/logos/Lumos Networks.png" # Se tiver um logo grande, use aqui
if os.path.exists(logo_principal):
    # Usamos colunas para centralizar a imagem
    col_l1, col_l2, col_l3 = st.columns([1, 2, 1])
    with col_l2:
        st.image(logo_principal, use_container_width=True)

# 2. O Título em Letras Clássicas (A Mágica acontece aqui)
# Em vez de st.title(), usamos st.image() com o texto já formatado.
# Isso garante que a fonte 'Harry Potter' apareça para TODOS os usuários.

titulo_magico = "assets/logo/titulo_lumos_hp.png" # Imagem que vamos criar

if os.path.exists(titulo_magico):
    # Centraliza o título gerado
    st.image(titulo_magico, use_container_width=True)
else:
    # Caso a imagem não exista (fallback), exibe o texto normal bem grande
    st.markdown("<h1 style='text-align: center; font-size: 50px; color: #E1AF12;'>LUMOS NETWORKS</h1>", unsafe_allow_html=True)
    st.markdown("<h3 style='text-align: center; color: #5D6D7E;'>Systems Biology Suite</h3>", unsafe_allow_html=True)

# Subtítulo (opcional, já que está no título da imagem)
# st.markdown('<p class="main-title">LUMOS NETWORKS: Systems Biology Suite </p>', unsafe_allow_html=True)
st.markdown('<p class="subtitle">Sistemas de Análise Integrada: Transcriptômica e Redes</p>', unsafe_allow_html=True)

st.divider()


# --- 1. LOGO NA BARRA LATERAL ---
with st.sidebar:
    logo_lumos = "assets/logos/Lumos Networks.png"
    if os.path.exists(logo_lumos):
        st.image(logo_lumos, use_container_width=True)
    st.divider()
    st.info("Selecione um módulo acima para iniciar sua análise.")

# Estilização CSS personalizada
st.markdown("""
    <style>
    .main-title { font-size: 45px; font-weight: bold; color: #2E86C1; text-align: center; margin-top: -20px; }
    .subtitle { font-size: 20px; text-align: center; color: #5D6D7E; margin-bottom: 30px; }
    .section-header { color: #2E86C1; border-bottom: 2px solid #2E86C1; padding-bottom: 5px; }
    </style>
    """, unsafe_allow_html=True)

st.divider()

col1, col2 = st.columns([2, 1])

with col1:
    st.markdown('<h2 class="section-header">Sobre o Projeto</h2>', unsafe_allow_html=True)
    st.write("""
    O **Lumos Networks** é uma suíte bioinformática projetada para otimizar a transição entre dados brutos de sequenciamento 
    e a interpretação biológica sistêmica. Desenvolvido no EvoMol-Lab, o sistema integra análise estatística, 
    contextualização funcional e modelagem de redes.
    """)
    
    st.subheader("Módulos da Suíte")
    
    # --- MÓDULO 1: DDEA ---
    with st.expander("1. DDEA (Diagonal Differential Expression Alley)", expanded=True):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            Análise de expressão gênica diferencial utilizando **PyDESeq2** e modelos lineares. 
            Focado em rigor estatístico e normalização robusta para Microarray e RNA-Seq.
            """)
        with c2:
            if os.path.exists("assets/logos/DDEA.png"):
                st.image("assets/logos/DDEA.png", width=80)
        
    # --- MÓDULO 2: APP ---
    with st.expander("2. APP (Arithmancy Pathway Profiler)", expanded=False):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            Realiza a contextualização funcional do transcriptoma através da integração multidimensional com bases de 
            dados globais como KEGG, GO e STRING. O módulo utiliza algoritmos de sobre-representação (ORA) e GSEA para 
            converter listas de genes em cascatas bioquímicas e complexos de interação proteína-proteína (PPI). 
            Seu diferencial estratégico reside no mapeamento de circuitos regulatórios via JASPAR/TRRUST, permitindo a 
            identificação dos fatores de transcrição que atuam como reguladores mestres do sistema biológico estudado.
            """)
        with c2:
            if os.path.exists("assets/logos/Arithmancy Pathway Profiler (APP) .png"):
                st.image("assets/logos/Arithmancy Pathway Profiler (APP) .png", width=80)

    # --- MÓDULO 3: PRIORIGRAPH (O QUE FALTAVA) ---
    with st.expander("3. PG (PrioriGraph)", expanded=False):
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown("""
            **Modelagem de Redes de Biologia de Sistemas.** Realiza a modelagem topológica de sistemas 
            biológicos através da construção de redes direcionadas. Ao integrar fatores de transcrição aos 
            seus genes alvos, o módulo aplica métricas de centralidade para priorizar hubs regulatórios 
            e detectar módulos funcionais conservados.
            """)
        with c2:
            if os.path.exists("assets/logos/PrioriGraph(PG).jpeg"):
                st.image("assets/logos/PrioriGraph(PG).jpeg", width=80)

with col2:
    # Este card de informação fica no corpo principal (Home)
    st.info("👈 **Navegação:** Use o menu à esquerda para alternar entre os módulos.")
    
    # Métrica de Release para destaque na Home
    st.metric(label="Release", value="v1.0.0-Beta")
    
    st.markdown("### Apoio Institucional")
    st.write("🔬 **EvoMol-Lab** e 🎓 **BioME - UFRN**")
    st.markdown("### Development Team")
    st.write("Dr. João Paulo M. S. Lima; MSc. Laís de Carvalho Gonçalves; Rodrigo Arruda Orvate")

st.divider()

# --- ORGANIZAÇÃO DA SIDEBAR (BARRA LATERAL) ---
with st.sidebar:
    st.markdown("### 🧠 Guia Rápido de Conceitos")
    st.divider()
    
    # Dica dinâmica
    st.info("💡 **Lumos Tip:** Comece sua jornada pelo módulo **DDEA** para identificar os genes com expressão diferencial estatisticamente significativa.")

    # Glossário Compacto - Útil para consulta rápida durante o uso
    with st.expander("🔎 O que é FDR (P-adj)?"):
        st.caption("""
        **Falsa Taxa de Descoberta.** Corrige o valor de p para múltiplos testes. 
        Em bioinformática, evita que identifiquemos genes como diferenciais por puro acaso.
        """)
        
    with st.expander("🕸️ O que é um Gene Hub?"):
        st.caption("""
        Definido pela **alta conectividade** (grau) em uma rede. 
        Funciona como um nó central de interações, sendo essencial para a 
        estabilidade estrutural da rede (proteínas de scaffold, por exemplo).
        """)

    with st.expander("👑 O que é um Regulador Mestre?"):
        st.caption("""
        Definido pela **função hierárquica**. Geralmente é um Fator de Transcrição 
        que orquestra grandes programas genéticos, controlando diretamente 
        ou indiretamente diversos genes a jusante (*downstream*).
        """)
        
    with st.expander("🧪 O que é GSEA?"):
        st.caption("""
        **Gene Set Enrichment Analysis.** Analisa se um conjunto de genes (vias) tende a estar no topo 
        ou no fundo de uma lista ranqueada, sem precisar de um corte fixo de P-value.
        """)

    st.divider()

# Rodapé centralizado e sem nomes individuais (conforme solicitado)
footer_html = """
<div style="text-align: center; color: #5D6D7E; padding: 10px;">
    <p><b>Lumos Networks © 2026</b> | Desenvolvido por integrantes do <b>EvoMol-Lab</b> (BioME - UFRN)</p>
</div>
"""
st.markdown(footer_html, unsafe_allow_html=True)
