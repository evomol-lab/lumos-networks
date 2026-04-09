# Lumos Networks
CÓDIGO DO DDEA

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
