

<img width="2048" height="1401" alt="Lumos Networks" src="https://github.com/user-attachments/assets/126e70c5-159f-47cd-b2fe-a7533de62d89" />




## Lumos Networks 🧬
An Integrated Suite for Transcriptomics and Biological Network Analysis

Developed by members of the EvoMol-Lab at the BioME - UFRN

_______________________________________________________________
🌟 **1. Introduction**

Lumos Networks is a modular Python web application designed to bridge the gap between raw transcriptomic data and systems biology insights. Developed at EvoMol-Lab (UFRN), the suite provides a streamlined workflow for Differential Expression Analysis, Functional Enrichment, and Knowledge-based Network construction.

_______________________________________________________________
🛠 **2. The Lumos Suite (Modules)**

Lumos is organized into three specialized modules, accessible via the sidebar in the multipage application:


<img width="480" height="590" alt="DDEA" src="https://github.com/user-attachments/assets/a7afe547-b1e9-4dd1-8407-aa5a47aef8ed" />



📊 **a. DDEA (Diagonal Differential Expression Analysis)**

Powered by PyDESeq2, this module handles the statistical heavy lifting of RNA-Seq. Differential expression analysis utilizes the generalized linear model (GLM) of the Negative Binomial family:

<img width="200" height="78" alt="calculo" src="https://github.com/user-attachments/assets/92c1af98-8873-4805-b116-bbe6d7d19b54" />


**Shrinkage Estimation:** We implemented empirical Bayesian dispersion estimation to stabilize fold change in genes with low counts, reducing technical noise and false positives.

**Multiple Testing Correction:** The software automatically applies the Benjamini-Hochberg (FDR) method to control the rate of false discoveries in large-scale experiments.

**Input:** Raw count matrices and metadata.

**Features:** Normalization, dispersion estimation, and Wald tests for differential expression.

**Output:** Volcano plots, MA plots, and interactive result tables.


<img width="2760" height="1504" alt="Arithmancy Pathway Profiler (APP) " src="https://github.com/user-attachments/assets/33fd5ac2-f63a-416f-b214-9bbe4e45ef45" />



🧬 **b. APP (Arithmancy Pathway Profiler)**

A functional enrichment module using GSEApy.

**Features:** Gene Set Enrichment Analysis (GSEA).

**Databases:** Integrated support for KEGG, Gene Ontology (GO) String and JASPAR/TRRUST.


![PrioriGraph(PG)](https://github.com/user-attachments/assets/48efb59a-0c1d-4a4e-8c6b-081a9d54dd1f)



🕸 **c. PG (PrioriGraph)**

Our flagship tool for building networks based on prior biological knowledge.

**Features:** Integration of DEGs and transcription factors in interaction networks.

**Metrics:** Identification of core genes and transcription factors through in-degree/out-degree and centrality analysis..

**Visualization:** High-performance interactive graph rendering.

_______________________________________________
⚙ **3. Installation & Local Execution**

```bash
# Bash
________________________________________________
# 1. Clone the repository
git clone [https://github.com/evomol-lab/lumos-networks.git](https://github.com/evomol-lab/lumos-networks.git)

# 2. Enter the directory
cd lumos-networks

# 3. Create a virtual environment (Optional but Recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# 4. Install dependencies
pip install -r requirements.txt

# 5. Run the application
streamlit run Lumos_Home.py
```
___________________________________________________
📄 **4. Core Requirements**

Lumos relies on the following state-of-the-art libraries:

**Statistics:** pydeseq2, scipy, statsmodels.

**Bioinformatics:** gseapy, biopython.

**Networks:** networkx, streamlit-agraph.

**Visualization:** plotly, matplotlib, seaborn.

______________________________________________________
📁 **5. Project Structure**

The repository follows the Streamlit Multipage pattern:

```text
text
______________________________________________________
lumos-networks/
├── .devcontainer/         # Configurations for an isolated development environment (Docker)
├── assets/                # Logos, images and visual documentation
├── fonts/                 # Custom font files for the Streamlit UI
├── packages               # Core logic, statistical functions, and reusable classes
├── pages/                 # Streamlit multipage application pages
├── LICENSE                # File containing the text of the Apache License 2.0
├── Lumos_Home.py          # Main file (App entry point)
├── README.md              # Main documentation
└── requirements.txt       # List of dependencies and versions
```

________________________________________________________
🤝 **6. Credits & Support**

The Lumos Networks suite is an ongoing collaborative effort developed at the EvoMol-Lab (Laboratory of Molecular Evolution and Bioinformatics), part of the Bioinformatics Multidisciplinary Environment (BioME) at the Federal University of Rio Grande do Norte (UFRN), Brazil.

👥 **Development Team**

**Dr. João Paulo M. S. Lima** – Principal Investigator (PI) - https://github.com/jpmslima

**MSc. Laís de Carvalho Gonçalves** – PhD Student & Lead Developer - https://github.com/laisdcg

**Rodrigo Arruda Orvate** – Master’s Student & Developer - https://github.com/RodrigoOrvate

**Djorkaeff Oliveira Fontinele** – Master’s Student & Developer - https://github.com/djkfof


🏛 **Institutions & Partners**

**UFRN:** Universidade Federal do Rio Grande do Norte.

**BioME:** Multi-user Bioinformatics Center UFRN.

**EvoMol-Lab:** Laboratory of Molecular Evolution and Systems Biology.


💰 **Financial Support**

This project is supported by the following Brazilian research agencies:

**CAPES** (Coordination for the Improvement of Higher Education Personnel)

**UFRN** (Institutional Support)

_____________________________________________________________________________________
⚖ **7. Disclaimer**

This software is provided for research purposes. The developer team utilized generative AI for UI/UX optimization, PDF reporting architecture, and documentation refinement to ensure the highest code quality and user experience.

_____________________________________________________________________________________
### Contato
Dúvidas ou sugestões? jpmslima@gmail.com
_____________________________________________________________________________________




