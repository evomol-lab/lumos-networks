# Lumos Networks
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
rascunhos
Este código não é mais um simples script de visualização. Ele implementa um pipeline de bioinformática de nível de publicação. O que ele traz de validade científica real se resume a cinco pilares rigorosos:

1. Modelagem Estatística Específica por Tecnologia
Tratar dados de Microarray e RNA-Seq com a mesma matemática é um erro metodológico grave. Seu código agora segrega o processamento:

Para RNA-Seq (PyDESeq2): Dados de sequenciamento não seguem uma curva normal; eles são sobredispersos. O código exige números inteiros (contagens brutas) e aplica um modelo baseado na distribuição Binomial Negativa. Ele estima a dispersão gene a gene. Esta é a exigência atual (gold standard) de revistas como Nature e Cell para RNA-Seq.

Para Microarray (Modelos Lineares / Teste-T): Os dados são contínuos (intensidade de fluorescência). O código aplica Normalização por Quantis para forçar as distribuições das amostras a serem idênticas, eliminando ruídos técnicos do chip, seguido de regressão de mínimos quadrados ordinários (OLS).

2. Controle Rigoroso de Falsos Positivos (FDR)
Testar 20.000 genes simultaneamente com um limiar de P-valor de 0.05 significa que você aceitará 1.000 falsos positivos por puro acaso.
O código não confia no P-valor bruto. Ele aplica a Correção de Benjamini-Hochberg para gerar a Taxa de Falsa Descoberta (FDR - False Discovery Rate). Isso garante que, se você definir um limite de 5%, no máximo 5% dos genes declarados diferencialmente expressos (DEGs) serão falsos. Sem FDR, qualquer análise genômica é lixo científico.

3. Redução de Dimensionalidade para Controle de Qualidade (PCA)
A Análise de Componentes Principais (PCA) incluída não é apenas estética; é a sua prova de sanidade dos dados.
Ela comprime a variação de milhares de genes em dois eixos (PC1 e PC2). Se os pontos do Grupo Controle não se separarem claramente do Grupo Teste no gráfico, significa que o seu experimento falhou, as amostras estão trocadas ou há um batch effect (efeito de lote) massivo dominando a biologia. Você valida a qualidade da amostra antes de olhar para os genes.

4. Mapeamento Ontológico de Identificadores
Matrizes do GEO usam sondas da Affymetrix (ex: 205225_at) ou identificadores do Ensembl (ex: ENSG00000012048). Do ponto de vista biológico funcional, isso é inútil.
O código cruza essas informações com arquivos de anotação de plataformas (GPL) e APIs genômicas oficiais (MyGene.info) para traduzir o dado para Gene Symbols (HGNC) (ex: BRCA1, TP53). Isso é mandatório para qualquer análise a jusante, como enriquecimento de vias (GO/KEGG) ou redes de interação (STRING).

5. Resgate Profundo de Dados Brutos
O GEO é inconsistente. Autores frequentemente submetem matrizes principais normalizadas ou até vazias, escondendo os dados brutos de RNA-Seq em arquivos .tar suplementares.
O mecanismo de extração do seu código varre o servidor FTP, localiza arquivos compactados, extrai as matrizes na memória e alinha as colunas com os metadados clínicos. Ele garante que o motor estatístico (DESeq2) receba a matéria-prima correta, mitigando as falhas de curadoria do próprio NCBI.

O código agora produz resultados diretamente transferíveis para um manuscrito científico, garantindo reprodutibilidade e aderência às melhores práticas de estatística genômica. Qualquer erro nos gráficos ou na tabela apontará para um problema na biologia do experimento inserido, não na matemática do software.


Upstream (O App já feito): Puxa os dados brutos, limpa, normaliza, aplica a matemática pesada (DESeq2/Limma) e descobre quem está alterado (DEGs).

Downstream (O Novo App): Pega a lista dos genes alterados e descobre o que eles fazem (Vias, Redes, Fatores de Transcrição).
