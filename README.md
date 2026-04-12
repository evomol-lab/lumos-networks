

## Lumos Networks 🧬
An Integrated Suite for Transcriptomics and Biological Network Analysis

Desenvolvido por Integrantes do EvoMol-Lab no Centro de Biociências - UFRN"

_______________________________________________________________________________________________________________________________________________________________________________________________________
🌟 **1. Introduction**

Lumos Networks is a modular Python web application designed to bridge the gap between raw transcriptomic data and systems biology insights. Developed at EvoMol-Lab (UFRN), the suite provides a streamlined workflow for Differential Expression Analysis, Functional Enrichment, and Knowledge-based Network construction.

________________________________________________________________________________________________________________________________________________________________________________________________________
🛠 **2. The Lumos Suite (Modules)**

Lumos is organized into three specialized modules, accessible via the sidebar in the multipage application:

<img width="480" height="590" alt="DDEA" src="https://github.com/user-attachments/assets/a7afe547-b1e9-4dd1-8407-aa5a47aef8ed" />

📊 *1. DDEA (Differential Differential Expression Analysis)*

Powered by PyDESeq2, this module handles the statistical heavy lifting of RNA-Seq.

Input: Raw count matrices and metadata.

Features: Normalization, dispersion estimation, and Wald tests for differential expression.

Output: Volcano plots, MA plots, and interactive result tables.

<img width="2760" height="1504" alt="Arithmancy Pathway Profiler (APP) " src="https://github.com/user-attachments/assets/33fd5ac2-f63a-416f-b214-9bbe4e45ef45" />

🧬 *2. APP (Analysis of Pathways & Processes)*

A functional enrichment module using GSEApy.

Features: Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA).

Databases: Integrated support for KEGG, Gene Ontology (GO), and Reactome.



🕸 *3. PrioriGraph (Network Biology)*

Our flagship tool for building networks based on prior biological knowledge.

Features: Integration of DEGs into interaction networks.

Metrics: Identification of hub genes through degree and centrality analysis.

Visualization: High-performance interactive graph rendering.

__________________________________________________________________________________________________________________________________________________________________________________________________________
⚙ **3. Installation & Local Execution**

Bash

Clone the repository
git clone https://github.com/evomol-lab/lumos-networks.git
cd lumos-networks

Install dependencies
pip install -r requirements.txt

Run the application
streamlit run Lumos_Home.py

_________________________________________________________________________________________________________________________________________________________________________________________________________
📄 **4. Core Requirements**

Lumos relies on the following state-of-the-art libraries:

Statistics: pydeseq2, scipy, statsmodels.

Bioinformatics: gseapy, biopython.

Networks: networkx, streamlit-agraph.

Visualization: plotly, matplotlib, seaborn.

__________________________________________________________________________________________________________________________________________________________________________________________________________
📁 **5. Project Structure**
The repository follows the Streamlit Multipage pattern:

Plaintext

lumos-networks/
├── Lumos_Home.py          # Landing page & Global Config
├── requirements.txt       # Version-controlled dependencies
├── pages/                 # Module directory
│   ├── 1_DDEA.py          # Differential Expression logic
│   ├── 2_APP.py           # Enrichment Analysis logic
│   └── 3_PrioriGraph.py   # Network construction logic
└── assets/                # Logos and documentation images

__________________________________________________________________________________________________________________________________________________________________________________________________________
🤝 **6. Credits & Support**
The Lumos Networks suite is an ongoing collaborative effort developed at the EvoMol-Lab (Laboratory of Molecular Evolution and Bioinformatics), part of the Bioinformatics Multidisciplinary Environment (BioME) at the Federal University of Rio Grande do Norte (UFRN), Brazil.

👥 Development Team
Dr. João Paulo M. S. Lima – Principal Investigator (PI)
MSc. Laís de Carvalho Gonçalves – PhD Student & Lead Developer
Djorkaeff Oliveira Fontinele – Master’s Student & Developer
Rodrigo Arruda Orvate – Master’s Student & Developer

🏛 Institutions & Partners
UFRN: Universidade Federal do Rio Grande do Norte.
BioME: Multi-user Bioinformatics Center UFRN.
EvoMol-Lab: Laboratory of Molecular Evolution and Systems Biology.

💰 Financial Support
This project is supported by the following Brazilian research agencies:
CAPES (Coordination for the Improvement of Higher Education Personnel)
UFRN (Institutional Support)

_______________________________________________________________________________________________________________________________________________________________________________________________________
⚖ **7. Disclaimer**
This software is provided for research purposes. The developer team utilized generative AI for UI/UX optimization, PDF reporting architecture, and documentation refinement to ensure the highest code quality and user experience.

_______________________________________________________________________________________________________________________________________________________________________________________________________
### Contato
Dúvidas ou sugestões? [bioinfo.imd.ufrn.br](http://bioinfo.imd.ufrn.br).
_______________________________________________________________________________________________________________________________________________________________________________________________________
Lumos Networks © 2026 | Natal, RN - Brazil





















___________________________________________________________________________________________________________________________________________________________________________________________________________
## Diagonal Differential Expression Alley (DDEA) 🧬

O **DDEA** é uma plataforma analítica de alta performance desenvolvida para simplificar e elevar o rigor científico na análise de expressão diferencial. A ferramenta automatiza o processamento de dados de **Microarray** e **RNA-Seq** provenientes do **NCBI Gene Expression Omnibus (GEO)**, garantindo resultados prontos para publicação.

Desenvolvido pelo **EvoMol-Lab** (BioME, UFRN, Brasil).

---

## 💎 Pilares de Validade Científica

O DDEA implementa um pipeline de bioinformática de nível de publicação, fundamentado em cinco pilares rigorosos:

### 1. Modelagem Estatística Específica por Tecnologia
O código segrega o processamento para respeitar a natureza físico-química de cada dado:
* **Para RNA-Seq (PyDESeq2):** Utiliza um modelo baseado na **distribuição Binomial Negativa** para lidar com a sobredispersão de dados de contagem bruta. É o *gold standard* exigido por revistas de alto impacto (Nature/Cell).
* **Para Microarray (Modelos Lineares):** Aplica **Normalização por Quantis** para eliminar ruídos técnicos do chip, seguida de regressão de mínimos quadrados ordinários (OLS).

### 2. Controle Rigoroso de Falsos Positivos (FDR)
O código não confia no P-valor bruto. Ele aplica a **Correção de Benjamini-Hochberg** para gerar a **Taxa de Falsa Descoberta (FDR - False Discovery Rate)**. Isso garante que, se definido um limite de 5%, no máximo 5% dos genes declarados diferencialmente expressos (DEGs) serão falsos positivos, garantindo a integridade genômica da análise.

### 3. Redução de Dimensionalidade (PCA) para Controle de Qualidade
A **Análise de Componentes Principais (PCA)** incluída funciona como uma prova de sanidade. Ela valida se os grupos (Controle vs. Teste) se separam biologicamente antes da análise individual de genes. Falhas na separação do PCA indicam possíveis problemas experimentais ou efeitos de lote (*batch effects*).

### 4. Mapeamento Ontológico de Identificadores
O código traduz identificadores técnicos (Probes Affymetrix, IDs Ensembl) para **Gene Symbols (HGNC)** oficiais através de integração com arquivos de plataforma (GPL) e APIs genômicas (**MyGene.info**). Isso torna o dado biologicamente interpretável e pronto para análises de enriquecimento (GO/KEGG).

### 5. Resgate Profundo de Dados Brutos
O mecanismo de extração varre servidores FTP do NCBI para localizar arquivos suplementares compactados quando as matrizes principais são inconsistentes. Ele garante que o motor estatístico receba a matéria-prima correta, mitigando falhas de curadoria do próprio GEO.

---

## 🚀 Fluxo de Trabalho

* **Upstream:** Puxa dados brutos, limpa, normaliza, aplica a matemática pesada (DESeq2/Limma) e identifica os Genes Diferencialmente Expressos (DEGs).
* **Downstream:** Exporta listas prontas para descoberta funcional (Vias, Redes e Fatores de Transcrição).

---

## 🛠️ Tecnologias Utilizadas

* **Streamlit:** Interface de usuário.
* **Pandas & NumPy:** Processamento de matrizes.
* **Plotly:** Gráficos interativos (Volcano, Heatmap, PCA).
* **PyDESeq2:** Modelagem estatística avançada para RNA-Seq.
* **Biopython & MyGene:** Anotação e interface com NCBI.
* **FPDF:** Geração de relatórios técnicos em PDF.

---

## 📋 Como Utilizar

1.  **Input:** Insira o GSE ID e selecione a tecnologia (Microarray ou RNA-Seq).
2.  **Grupos:** Defina os nomes dos grupos e selecione as amostras correspondentes.
3.  **Parâmetros:** Ajuste os limiares de FDR e Log2 Fold Change.
4.  **Análise:** Clique em `Run Analysis` para processar o pipeline estatístico.
5.  **Relatório:** Baixe o relatório PDF completo com métricas e gráficos de nível de publicação.

---

## 🧬 Citação e Créditos

> **DDEA Analytics - EvoMol-Lab (BioME, UFRN, Natal-RN, Brazil).**
> Desenvolvido no Instituto do Cérebro (ICe/UFRN) para garantir reprodutibilidade e aderência às melhores práticas de estatística genômica.

---

### Contato
Dúvidas ou sugestões? [bioinfo.imd.ufrn.br](http://bioinfo.imd.ufrn.br).


---
r
