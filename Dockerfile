# Usa uma imagem oficial do Python, versão slim para economizar espaço
FROM python:3.10-slim

# Define o diretório de trabalho dentro do contêiner
WORKDIR /app

# Instala dependências do sistema necessárias para compilar pacotes científicos em C++
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copia os arquivos de dependência primeiro para otimizar o cache do Docker
COPY requirements.txt .

# Instala as bibliotecas do Python
RUN pip install --no-cache-dir -r requirements.txt

# Copia todo o restante do código-fonte para o contêiner
COPY . .

# Expõe a porta padrão que o Streamlit usa
EXPOSE 8501

# Define a variável de ambiente para o Streamlit rodar em modo headless (sem abrir navegador no servidor)
ENV STREAMLIT_SERVER_HEADLESS=true

# Comando obrigatório para iniciar a aplicação
CMD ["streamlit", "run", "Lumos_Home.py"]
