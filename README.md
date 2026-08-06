<img width="540" height="450" alt="Lumos Networks" src="https://github.com/evomol-lab/lumos-networks/blob/111ccd4b2d08c3af9d81d88181dbbd0ba81cd94c/assets/Lumos%20Nexus.png" />


<img width="540" height="450" alt="Lumos Networks" src="https://github.com/evomol-lab/lumos-networks/blob/25d9efc53274244b2a08f49b8673dac719bff44f/assets/LumosPy.png" />


<img width="540" height="450" alt="Lumos Networks" src="https://github.com/user-attachments/assets/126e70c5-159f-47cd-b2fe-a7533de62d89" />


## Lumos Networks 🧬
An Integrated Suite for Transcriptomics and Biological Network Analysis

Developed by members of the EvoMol-Lab at the BioME - UFRN

_______________________________________________________________
🌟 **1. Introduction**

Lumos Networks is a modular Python web application designed to bridge the gap between raw transcriptomic data and systems biology insights. Developed at EvoMol-Lab (UFRN), the suite provides a streamlined workflow for Differential Expression Analysis, Functional Enrichment, and Knowledge-based Network construction.

_______________________________________________________________
🛠 **2. The Lumos Suite (Modules)**

Lumos is organized into three specialized modules, accessible via the sidebar in the multipage application:


<img width="250" height="350" alt="DDEA" src="https://github.com/user-attachments/assets/a7afe547-b1e9-4dd1-8407-aa5a47aef8ed" />



📊 **a. DDEA (Diagonal Differential Expression Analysis)**

Powered by PyDESeq2, this module handles the statistical heavy lifting of RNA-Seq. Differential expression analysis utilizes the generalized linear model (GLM) of the Negative Binomial family:

<img width="200" height="78" alt="calculo" src="https://github.com/user-attachments/assets/92c1af98-8873-4805-b116-bbe6d7d19b54" />

* **Shrinkage Estimation:** We implemented empirical Bayesian dispersion estimation to stabilize fold change in genes with low counts, reducing technical noise and false positives.

* **Multiple Testing Correction:** The software automatically applies the Benjamini-Hochberg (FDR) method to control the rate of false discoveries in large-scale experiments.
  
* **Input:** Raw count matrices and metadata.

* **Features:** Normalization, dispersion estimation, and Wald tests for differential expression.
  
* **Output:** Volcano plots, MA plots, interactive result tables, and automated PDF reporting via fpdf2.


<img width="320" height="200" alt="Arithmancy Pathway Profiler (APP) " src="https://github.com/user-attachments/assets/33fd5ac2-f63a-416f-b214-9bbe4e45ef45" />



🧬 **b. APP (Arithmancy Pathway Profiler)**

A functional enrichment module using GSEApy. Unlike simple list tools, the APP module uses the weighted Gene Set Enrichment Analysis (GSEA) algorithm:

- The Enrichment Score (ES) calculation reflects the degree to which a set of genes $S$ is represented at the extremes (top or bottom) of a ranked list $L$.

- This allows for the identification of activated or repressed biological pathways even when individual genes do not reach isolated statistical significance.

**Features:** Gene Set Enrichment Analysis (GSEA).

**Databases:** Integrated support for KEGG, Gene Ontology (GO) String and JASPAR/TRRUST.


<img src="https://github.com/user-attachments/assets/48efb59a-0c1d-4a4e-8c6b-081a9d54dd1f" width="200" alt="PrioriGraph Module Preview">
</details>



🕸 **c. PG (PrioriGraph)**

Our flagship tool for building networks based on prior biological knowledge. The distinguishing feature of PrioriGraph is the transition from co-expression to regulatory causality.

- **Eigenvector Centrality:** We identify Master Regulators not only by the number of connections (degree), but by the quality of those connections (being connected to other important nodes).

- **Integration by Prior Knowledge:** Unlike de novo networks, PrioriGraph overlays its RNA-Seq data with validated interactions (curated in STRING/JASPAR), mitigating the problem of spurious correlations.

**Features:** Integration of DEGs and transcription factors in interaction networks.

**Metrics:** Identification of core genes and transcription factors through in-degree/out-degree and centrality analysis..

**Visualization:** High-performance interactive graph rendering.

_______________________________________________
🛠️ **3. Professional Engineering Standards**
To ensure reproducibility (a fundamental pillar of Open Science), the project adopts:

### 🐳 Reproducibility with Docker & .devcontainer
The repository is fully containerized. The use of Docker and `.devcontainer` ensures that the execution environment (Python versions, C++ compilers for PyDESeq2, and system drivers) is identical for any researcher. This eliminates the common error of "conflicting dependencies" on Linux/Ubuntu systems and allows for platform-agnostic deployment.

### 🧪 Unit Testing & Validation
The repository includes a test suite (`/tests`) that validates:
* The integrity of count normalization.
* The convergence of network layout algorithms.
* The consistency of metadata during module parsing.

______________________________________________________
📊 **4. Use Case: From Raw Counts to Discovery**

**Normalization:** The raw data is transformed to a logarithmic stable scale (VST or rlog).

**Contrast Matrix:** The user defines the experimental vs. control condition.

**Cross-Validation:** The DEGs identified in DDEA are automatically injected into PrioriGraph for system visualization.

______________________________________________________
⚙ **5. Installation & Execution**

### Docker Deployment
To guarantee environment reproducibility, build and run the application via Docker:

```bash
# 1. Clone the repository
git clone [https://github.com/evomol-lab/lumos-networks.git](https://github.com/evomol-lab/lumos-networks.git)
cd lumos-networks

# 2. Build the Docker image
docker build -t lumos-networks .

# 3. Run the container
docker run -p 8501:8501 lumos-networks
```
___________________________________________________
📄 **6. Core Requirements**

Lumos relies on the following state-of-the-art libraries:

**Statistics:** pydeseq2, scipy, statsmodels.

**Bioinformatics:** gseapy, biopython.

**Networks:** networkx, streamlit-agraph.

**Visualization:** plotly, matplotlib, seaborn.

______________________________________________________
📁 **7. Project Structure**

The repository follows the Streamlit Multipage pattern:

```text
text
______________________________________________________
lumos-networks/
├── Dockerfile             # Containerization and deployment configuration
├── assets/                # Logos, images, and visual documentation
├── fonts/                 # Custom font files for the Streamlit UI
├── packages               # Core logic, statistical functions, and reusable classes
├── pages/                 # Streamlit multipage application pages (DDEA, APP, PrioriGraph)
├── LICENSE                # Apache License 2.0 text
├── Lumos_Home.py          # Main file (App entry point)
├── README.md              # Main documentation
├── packages.txt           # System-level OS packages required for deployment
└── requirements.txt       # Python dependency version control
```

________________________________________________________
🤝 **8. Credits & Support**

The Lumos Networks suite is an ongoing collaborative effort developed at the EvoMol-Lab (Laboratory of Molecular Evolution and Bioinformatics), part of the Bioinformatics Multidisciplinary Environment (BioME) at the Federal University of Rio Grande do Norte (UFRN), Brazil.

👥 **Development Team**

**Dr. João Paulo M. S. Lima** – Principal Investigator (PI) - https://github.com/jpmslima

**MSc. Laís de Carvalho Gonçalves** – PhD Student & Lead Developer - https://github.com/laisdcg

**Rodrigo Arruda Orvate** – Master’s Student & Developer - https://github.com/RodrigoOrvate

🏛 **Institutions & Partners**

**UFRN:** Universidade Federal do Rio Grande do Norte.

**BioME:** Multi-user Bioinformatics Center UFRN.

**EvoMol-Lab:** Laboratory of Molecular Evolution and Systems Biology.


💰 **Financial Support**

This project is supported by the following Brazilian research agencies:

**CAPES** (Coordination for the Improvement of Higher Education Personnel)

**UFRN** (Institutional Support)

____________________________________________________________________________________
**9.References & Citations**

To support the ecosystem, please credit the open-source libraries integrated into this tool when publishing your research:

Streamlit Team. (2023). Streamlit: The fastest way to build and share data apps. https://streamlit.io

Muzellec, B., et al. (2023). PyDESeq2: a python implementation of DESeq2. https://github.com/owkin/PyDESeq2

Fang, Z., et al. (2023). GSEApy: Gene Set Enrichment Analysis in Python. https://github.com/zqfang/GSEApy

Hagberg, A., et al. (2008). Exploring Network Structure, Dynamics, and Function using NetworkX. https://networkx.org

Cock, P. J., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics.

Harris, C. R., et al. (2020). Array programming with NumPy. Nature; McKinney, W. (2010). Data Structures for Statistical Computing in Python.

Pedregosa, F., et al. (2011). Scikit-learn: Machine Learning in Python. Journal of Machine Learning Research.

Virtanen, P., et al. (2020). SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods.

Seabold, S., & Perktold, J. (2010). Statsmodels: Econometric and statistical modeling with python.

Plotly Technologies Inc. (2015). Collaborative data science. https://plot.ly

Thommes, M. (2023). Streamlit-agraph: Interactive Graph Visualizations in Streamlit.

Reingold, J., et al. (2023). fpdf2: Simple PDF generation for Python. https://github.com/py-pdf/fpdf2


_____________________________________________________________________________________
⚖ **10. Research Software Disclaimer & Development Notice**

This software is provided exclusively for research purposes. To ensure peak code quality and an optimized user experience, the development team integrated generative AI tools throughout the engineering process.

AI assistance was specifically utilized for:

**System Architecture & UI/UX**: Optimization of user interfaces and the structural design of the PDF reporting architecture.

**Code Engineering**: Technical revision, logic optimization, and performance refinement.

**Documentation & Language**: Elaborating complex topic structures and conducting comprehensive English language reviews for clarity and precision.

**Disclaimer**: This software is intended for academic and scientific use. While generative AI was employed to enhance technical efficiency and documentation standards, the final output remains under the oversight of the development team to ensure research integrity.

_____________________________________________________________________________________
### Contact
Questions or suggestions? Reach out to jpmslima@gmail.com

________________________________________________________________________________

**Lumos Networks** © 2026 | Natal, RN - Brazil


