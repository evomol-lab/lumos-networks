import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial import distance
import requests
import gzip
import io
import re
import gc
import os
import tarfile
import json
import platform
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from Bio import Entrez
from statsmodels.stats.multitest import fdrcorrection

# Meta-analysis support (statsmodels >= 0.14)
try:
    from statsmodels.stats.meta_analysis import combine_effects
    HAS_META = True
except ImportError:
    HAS_META = False
from sklearn.decomposition import PCA

from lumos_report import build_html_report, table_spec
from lumos_core import (
    __version__ as LUMOS_VERSION,
    align_named_contrasts,
    align_supplementary_columns,
    annotate_cross_dataset_patterns,
    assemble_geo_tar_tables,
    build_design_metadata,
    classify_differential_expression,
    collapse_symbol_results,
    make_unique_sample_labels,
    package_versions,
    orient_gene_sample_matrix,
    quantile_normalize,
    resolve_sample_columns,
    validate_count_matrix,
)

# --- BLOCO DE IMPORTAÇÃO UNIFICADO E PROTEGIDO ---
try:
    # Para pydeseq2 >= 0.4.0, o DeseqDataSet fica em .dds
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    HAS_DESEQ2 = True
except ImportError:
    HAS_DESEQ2 = False

# 1. Configuração da página
st.set_page_config(page_title="Lumos Networks | Analysis", page_icon="🧬", layout="wide")

# 2. Configurações de API e Headers
HEADERS = {'User-Agent': 'DDEA/4.0 (Streamlit App; Academic Version)'}
Entrez.email = os.getenv("NCBI_EMAIL", "jpmslima@gmail.com")

# 3. Gerenciamento de Caminhos e Logo
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT_DIR = os.path.dirname(BASE_DIR)

LOGO_PATH = os.path.join(PARENT_DIR, "assets", "Lumos Networks.png")
DDEA_LOGO_PATH = os.path.join(PARENT_DIR, "assets", "DDEA.png")

# 4. CSS Estilizado
st.markdown("""
    <style>
    [data-testid="stSidebarNav"] {display: none;}
    .stPageLink {
        background-color: rgba(151, 166, 195, 0.12);
        border-radius: 20px;
        padding: 8px;
        border: 1px solid rgba(151, 166, 195, 0.3);
    }
    </style>
    """, unsafe_allow_html=True)

# --- SIDEBAR PADRONIZADA ---
with st.sidebar:
    if os.path.exists(LOGO_PATH):
        st.image(LOGO_PATH, width=250)
    
    st.divider()
    st.markdown("### 🚀 Navigation")
    st.page_link("Lumos_Home.py", label="Home Page", icon="🏠")
    
    st.markdown('<p style="color:#2E86C1; font-weight:bold; margin-bottom:0px; margin-top:10px;">📊 Analysis</p>', unsafe_allow_html=True)
    st.page_link("pages/1_DDEA.py", label="DDEA", icon="📈")
    
    st.markdown('<p style="color:#28B463; font-weight:bold; margin-bottom:0px; margin-top:10px;">🧬 Functional</p>', unsafe_allow_html=True)
    st.page_link("pages/2_APP.py", label="APP", icon="🧪")
    
    st.markdown('<p style="color:#E67E22; font-weight:bold; margin-bottom:0px; margin-top:10px;">🕸️ Networks</p>', unsafe_allow_html=True)
    st.page_link("pages/3_PG.py", label="PG", icon="🕸️")

    st.markdown('<p style="color:#E1AF12; font-weight:bold; margin-bottom:0px; margin-top:10px;">📚 Documentation</p>', unsafe_allow_html=True)
    st.page_link("pages/Documentation.py", label="Documentation", icon="📚")

    st.divider()
    if os.path.exists(DDEA_LOGO_PATH):
        st.image(DDEA_LOGO_PATH, width="stretch")
    st.info("You are in the analysis module.")


# ============================================================
# HTML REPORT GENERATION
# ============================================================

def generate_html_report(ctx: dict) -> str:
    df = ctx.get('df_diff')
    df = df if df is not None else pd.DataFrame(columns=['Symbol', 'Log2FC', 'FDR', 'Status'])
    figs = ctx.get('figures') or {}
    n_up, n_down = int((df['Status'] == 'Up-regulated').sum()), int((df['Status'] == 'Down-regulated').sum())
    tn, rn = ctx.get('tn', ''), ctx.get('rn', '')
    info = [("GSE ID", ctx.get('gse_id', 'N/A')), ("Experiment type", ctx.get('mode', '')), ("Comparison", f"{tn} vs {rn}"),
            ("Samples per group", f"{tn}: {int(ctx.get('n_test', 0))} · {rn}: {int(ctx.get('n_ref', 0))}"),
            ("Analysis method", ", ".join(ctx.get('methods') or []) or "N/A"), ("Design formula", ctx.get('design_formula') or "~condition"),
            ("Thresholds", f"FDR < {ctx.get('p_thr')}, |Log2FC| ≥ {ctx.get('fc_thr')}"), ("Lumos version", ctx.get('version', ''))]
    cards = [('Total DEGs', len(df), ''), ('Up-regulated', n_up, 'up'), ('Down-regulated', n_down, 'down')]
    deg = table_spec(df, [('Symbol', 'Symbol', 'text'), ('Log2FC', 'Log2FC', 'fixed3'), ('FDR', 'FDR', 'exp2'), ('Status', 'Status', 'text')], 'FDR',
                     filter_col='Status', badge_col='Status', badge_map={'Up-regulated': 'up', 'Down-regulated': 'down'}) if len(df) else None
    sections = [{'id': 'volcano', 'title': 'Volcano Plot', 'desc': 'Global distribution of fold change versus adjusted significance.', 'figs': [figs.get('fig_v')]},
                {'id': 'ma', 'title': 'MA Plot', 'desc': 'Effect size versus mean abundance.', 'figs': [figs.get('fig_ma')]},
                {'id': 'pca', 'title': 'PCA', 'desc': 'Sample clustering on the first two principal components.', 'figs': [figs.get('fig_pca')]},
                {'id': 'top-degs', 'title': 'Top Differentially Expressed Genes', 'desc': 'Strongest up- and down-regulated genes.', 'figs': [figs.get('fig_up'), figs.get('fig_down')]},
                {'id': 'deg-table', 'title': 'DEG Table', 'desc': 'All differentially expressed genes: search, filter, sort and export.', 'table': deg, 'table_id': 'DEGs'},
                {'id': 'heatmap', 'title': 'Expression Heatmap', 'desc': f"Top {int(ctx.get('max_plot') or 0)} DEGs, row z-scores.", 'figs': [figs.get('fig_h')]}]
    cites = {'PyDESeq2': 'PyDESeq2 (Muzellec et al., 2023)', 'limma': 'limma (Ritchie et al., 2015)'}
    used = [c for k, c in cites.items() if any(k.lower() in m.lower() for m in ctx.get('methods') or [])]
    cite = f"<br>Differential expression computed with {' and '.join(used)}." if used else ''
    return build_html_report('DDEA Analysis Report', DDEA_LOGO_PATH, info, cards, sections,
                             f'Diagonal Differential Expression Alley (DDEA) &mdash; EvoMol-Lab, BioME, UFRN{cite}',
                             f"DDEA_{ctx.get('gse_id') or 'report'}", doc_title=f"DDEA Report - {ctx.get('gse_id', '')}")


# ============================================================
# NORMALIZAÇÃO E UTILITÁRIOS
# ============================================================

def detect_index_type(index_values):
    samples = [str(v).strip() for v in index_values[:50] if str(v).strip()]
    if not samples: return 'unknown'
    ensembl = sum(1 for s in samples if re.match(r'^ENS[A-Z]*G\d{11}', s))
    entrez = sum(1 for s in samples if re.match(r'^\d+$', s))
    probe = sum(1 for s in samples if re.match(r'^\d+_[a-z]', s) or re.match(r'^[A-Z]{1,3}\d{6,}', s))
    symbol = sum(1 for s in samples if re.match(r'^[A-Za-z][A-Za-z0-9\-\.]{1,15}$', s) and not re.match(r'^\d+$', s) and '_' not in s)
    scores = {'ensembl': ensembl, 'entrez': entrez, 'probe': probe, 'symbol': symbol}
    best = max(scores, key=scores.get)
    return best if scores[best] > 0 else 'unknown'

def _strip_ensembl_version(ensembl_id: str) -> str:
    return ensembl_id.split('.')[0] if '.' in ensembl_id else ensembl_id
    
# ============================================================
# MAPEAMENTO 
# ============================================================

@st.cache_data(show_spinner=False)
def get_gene_mapping_microarray(gse_id):
    try:
        search_handle = Entrez.esearch(db="gds", term=f"{gse_id}[ACCN]")
        uid = Entrez.read(search_handle)["IdList"][0]
        record = Entrez.read(Entrez.esummary(db="gds", id=uid))
        gpl_id = record[0]['GPL']
        gpl_prefix = f"GPL{gpl_id[:-3]}nnn" if len(gpl_id) > 3 else "GPLnnn"
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_prefix}/GPL{gpl_id}/soft/GPL{gpl_id}_family.soft.gz"
        response = requests.get(url, stream=True, headers=HEADERS)
        with gzip.open(response.raw, 'rt', encoding='utf-8', errors='ignore') as f:
            table_lines, in_table = [], False
            for line in f:
                if line.startswith('!platform_table_begin'): in_table = True; continue
                if line.startswith('!platform_table_end'): break
                if in_table: table_lines.append(line)
                if len(table_lines) > 250000: break
            map_df = pd.read_csv(io.StringIO("".join(table_lines)), sep='\t', low_memory=False)
            target_cols = ['Gene.Symbol', 'Gene Symbol', 'GENE_SYMBOL', 'Symbol', 'SYMBOL']
            symbol_col = next((c for c in map_df.columns if any(k.upper() == c.upper() for k in target_cols)), None)
            if symbol_col:
                map_df['ID'] = map_df['ID'].astype(str).str.strip().str.replace('"', '')
                map_df[symbol_col] = map_df[symbol_col].astype(str).apply(lambda x: x.split(' /// ')[0])
                return map_df[['ID', symbol_col]].rename(columns={'ID': 'Probe_ID', symbol_col: 'Symbol'}), map_df.head(100)
        return None, None
    except: return None, None

@st.cache_data(show_spinner=False)
def get_gene_mapping_rnaseq(index_ids: tuple, id_type: str):
    if id_type == 'symbol': return pd.DataFrame({'Probe_ID': list(index_ids), 'Symbol': list(index_ids)}), "Já contém Gene Symbols."
    if id_type not in ('entrez', 'ensembl'): return None, f"Tipo '{id_type}' não suportado."
    
    ids_raw = [str(i).strip() for i in index_ids if str(i).strip()]
    scope = "entrezgene" if id_type == 'entrez' else "ensembl.gene"
    original_to_clean = {i: i if id_type == 'entrez' else _strip_ensembl_version(i) for i in ids_raw}
    clean_to_original = {v: k for k, v in original_to_clean.items()}
    ids_clean = list(clean_to_original.keys())

    results = []
    for i in range(0, len(ids_clean), 1000):
        try:
            resp = requests.post("https://mygene.info/v3/query", data={"q": ",".join(ids_clean[i:i + 1000]), "scopes": scope, "fields": "symbol", "species": "human"}, timeout=30)
            if resp.status_code == 200:
                for item in resp.json():
                    if not item.get('notfound'):
                        q_id = str(item.get('query', ''))
                        results.append({"Probe_ID": clean_to_original.get(q_id, q_id), "Symbol": item.get('symbol', None)})
        except: continue

    if not results: return None, "Sem resultados MyGene."
    mapping_df = pd.DataFrame(results).drop_duplicates('Probe_ID')
    mapping_df["Symbol"] = mapping_df["Symbol"].fillna(mapping_df["Probe_ID"])
    return mapping_df, f"{(mapping_df['Symbol'] != mapping_df['Probe_ID']).sum()}/{len(mapping_df)} mapeados."

# ============================================================
# PARSE E EXTRAÇÃO GEO
# ============================================================

def _parse_matrix_bytes(raw_bytes, filename=""):
    try:
        content = gzip.decompress(raw_bytes) if filename.endswith('.gz') else raw_bytes
        for sep in ['\t', ',']:
            try:
                df = pd.read_csv(io.BytesIO(content), sep=sep, header=None, low_memory=False)
                
                # Se a primeira linha contiver texto, usar como header
                if any(df.iloc[0].astype(str).str.contains(r'[A-Za-z]', regex=True)):
                    df.columns = df.iloc[0]
                    df = df[1:].reset_index(drop=True)
                
                # A primeira coluna quase sempre é a coluna de genes (Probe/Symbol/Ensembl)
                gene_col = df.columns[0]
                
                # Procura por uma coluna de contagem clara
                count_col = next((c for c in df.columns if isinstance(c, str) and ('count' in c.lower() or 'read' in c.lower())), None)
                
                # Se não encontrar uma coluna de "count" explícita, mas o dataframe for bidimensional (Genes x Amostras)
                if not count_col and len(df.columns) > 1:
                    df_clean = df.copy()
                    df_clean[gene_col] = df_clean[gene_col].astype(str).str.strip().str.replace('"', '')
                    df_clean = df_clean.dropna(subset=[gene_col]).set_index(gene_col)
                    
                    for col in df_clean.columns:
                        df_clean[col] = pd.to_numeric(df_clean[col], errors='coerce')
                    
                    df_clean = df_clean.dropna(axis=1, how='all')
                    
                    if not df_clean.empty:
                        df_clean = df_clean.loc[:, ~df_clean.columns.duplicated()]
                        df_clean = df_clean.groupby(df_clean.index).sum()
                        return df_clean
                
                # Fallback para o layout antigo
                if count_col:
                    df_clean = df[[gene_col, count_col]].copy()
                    df_clean.columns = ['gene_id', 'raw_counts']
                    df_clean['gene_id'] = df_clean['gene_id'].astype(str).str.strip().str.replace('"', '')
                    df_clean['raw_counts'] = pd.to_numeric(df_clean['raw_counts'], errors='coerce')
                    df_clean = df_clean.dropna().set_index('gene_id')
                    
                    df_clean = df_clean[~df_clean.index.duplicated(keep='first')]
                    return df_clean['raw_counts']

            except Exception:
                continue
    except Exception:
        return None

def _try_series_matrix(gse_id):
    num = gse_id.replace("GSE", "")
    prefix = f"GSE{num[:-3]}nnn" if len(num) > 3 else "GSEnnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    try:
        r = requests.get(url, timeout=60, headers=HEADERS)
        if r.status_code != 200: return None, None, None, [], 'unknown'
        with gzip.open(io.BytesIO(r.content), 'rt') as f:
            titles, gsms, gsm_order, stype = [], [], [], ""
            meta_dict = {}
            df = pd.DataFrame()
            for line in f:
                line = line.strip()
                if line.startswith('!Series_type'): 
                    stype = line.lower()
                elif line.startswith('!Sample_title'): 
                    meta_dict['Title'] = [t.strip().replace('"', '') for t in line.split('\t')[1:]]
                elif line.startswith('!Sample_geo_accession'):
                    gsms = [t.strip().replace('"', '') for t in line.split('\t')[1:]]
                    meta_dict['Accession'] = gsms
                    gsm_order = gsms[:]
                elif line.startswith('!Sample_source_name_ch1'):
                    meta_dict['source_name_ch1'] = [t.strip().replace('"', '') for t in line.split('\t')[1:]]
                elif line.startswith('!Sample_characteristics_ch1'):
                    vals = [v.strip().replace('"', '') for v in line.split('\t')[1:]]
                    if vals:
                        key = vals[0].split(':')[0] if ':' in vals[0] else f"Char_{len(meta_dict)}"
                        meta_dict[key] = [v.split(': ')[1] if ': ' in v else v for v in vals]
                elif line.startswith('ID_REF') or line.startswith('"ID_REF"'):
                    df = pd.read_csv(f, sep='\t', header=None, low_memory=True)
                    break

            det_type = 'RNASeq' if 'sequencing' in stype else 'Microarray' if 'array' in stype else 'unknown'
            meta_df = pd.DataFrame(meta_dict)
            if "Accession" in meta_df.columns:
                meta_df["Accession"] = meta_df["Accession"].astype(str).str.extract(
                    r"(GSM\d+)", expand=False
                )
                meta_df = meta_df.dropna(subset=["Accession"]).drop_duplicates(
                    "Accession", keep="first"
                ).reset_index(drop=True)
                gsm_order = meta_df["Accession"].tolist()
                gsms = gsm_order.copy()

            if df.empty or len(df.columns) < 2: return None, meta_df, gsms, gsm_order, det_type
            df = df.set_index(0)
            df.index = df.index.astype(str).str.strip().str.replace('"', '')
            if len(df.columns) != len(gsm_order):
                return None, meta_df, gsms, gsm_order, det_type
            if len(set(gsm_order)) != len(gsm_order):
                return None, meta_df, gsms, gsm_order, det_type
            df.columns = gsm_order
            df = df.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
            return df, meta_df, gsms, gsm_order, det_type
    except: return None, None, None, [], 'unknown'

def _sync_suppl_columns_with_gsms(df_suppl, gsm_order):
    return align_supplementary_columns(df_suppl, gsm_order)


def get_geo_full_data(gse_id, mode, log_cb=None):
    gse_id = str(gse_id).strip().upper()
    if not re.fullmatch(r'GSE\d+', gse_id):
        return None, None, None, [], None, "Invalid GEO Series accession. Expected GSE followed by digits.", 'unknown'
    if log_cb: log_cb("Buscando metadados...")
    
    try:
        res_matrix = _try_series_matrix(gse_id)
        if res_matrix and len(res_matrix) == 5:
            df_matrix, meta_df, gsms, gsm_order, detected_type = res_matrix
        else:
            return None, None, None, [], None, "Formato inválido do Series Matrix.", 'unknown'
    except Exception as e:
        return None, None, None, [], None, f"Erro ao ler Series Matrix: {str(e)}", 'unknown'
        
    if meta_df is None: 
        return None, None, None, [], None, "Falha GEO: Metadados não encontrados.", 'unknown'
    
    if df_matrix is not None and df_matrix.shape[1] >= 2 and (mode != "RNASeq" or df_matrix.max().max() > 50):
        return df_matrix, meta_df, gsms, gsm_order, "Series Matrix", None, detected_type

    num = gse_id.replace("GSE", "")
    prefix = f"GSE{num[:-3]}nnn" if len(num) > 3 else "GSEnnn"
    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/suppl/"
    
    try:
        r = requests.get(base_url, timeout=20, headers=HEADERS)
        if r.status_code == 200:
            urls = [base_url + f[0] for f in re.findall(r'href="([^"]+\.(txt|tsv|csv|tar)(\.gz)?)"', r.text, re.IGNORECASE)]
            tar_url = next((u for u in urls if '_RAW.tar' in u or '.tar' in u.lower()), None)

            # Prefer a study-level count matrix.  It is smaller, already
            # consolidated and avoids downloading a potentially large RAW TAR.
            direct_errors = []
            non_tar_urls = [u for u in urls if '.tar' not in u.lower()]
            for u in sorted(
                non_tar_urls,
                key=lambda x: (
                    'count' in x.lower(),
                    'matrix' in x.lower(),
                    'raw' in x.lower(),
                ),
                reverse=True,
            ):
                try:
                    response = requests.get(u, timeout=300, headers=HEADERS)
                    response.raise_for_status()
                    df_suppl = _parse_matrix_bytes(response.content, u)
                    if isinstance(df_suppl, pd.DataFrame) and df_suppl.shape[1] >= 2:
                        df_synced, msg = _sync_suppl_columns_with_gsms(df_suppl, gsm_order)
                        return df_synced, meta_df, gsms, gsm_order, f"Suppl ({msg})", None, detected_type
                except Exception as direct_exc:
                    direct_errors.append(f"{os.path.basename(u)}: {direct_exc}")
            
            if tar_url:
                local_tar_path = f"{gse_id}_RAW_temp.tar"
                progress_bar = st.progress(0, text="🧬 Conectando ao GEO para baixar o arquivo TAR...")
                
                with requests.get(tar_url, stream=True, headers=HEADERS, timeout=300) as stream_resp:
                    stream_resp.raise_for_status()
                    total_size = int(stream_resp.headers.get('content-length', 0))
                    bytes_downloaded = 0
                    
                    with open(local_tar_path, 'wb') as f_out:
                        for chunk in stream_resp.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                f_out.write(chunk)
                                bytes_downloaded += len(chunk)
                                if total_size > 0:
                                    percent = min(1.0, bytes_downloaded / total_size)
                                    progress_bar.progress(percent, text=f"📥 Baixando arquivo RAW.tar ({bytes_downloaded // (1024*1024)}MB baixados)...")

                progress_bar.progress(1.0, text="📦 Descompactando e consolidando a matriz de contagens...")
                with tarfile.open(local_tar_path, 'r') as tar:
                    members = [m for m in tar.getmembers() if m.isfile() and (m.name.endswith('.gz') or m.name.endswith('.txt'))]

                    def parsed_members():
                        for i, member in enumerate(members):
                            f_extracted = tar.extractfile(member)
                            if f_extracted:
                                filename = member.name.split('/')[-1]
                                parsed = _parse_matrix_bytes(f_extracted.read(), filename)
                                if parsed is not None:
                                    yield filename, parsed
                            if i % 5 == 0 or i == len(members) - 1:
                                progress_bar.progress(1.0, text=f"🧬 Processando amostras: {i+1} de {len(members)} concluídas...")

                    combined_df, tar_method = assemble_geo_tar_tables(parsed_members(), gsm_order)

                if os.path.exists(local_tar_path):
                    os.remove(local_tar_path)
                progress_bar.empty()

                if not combined_df.empty:
                    combined_df = combined_df.fillna(0).astype(int)
                    df_synced, msg = _sync_suppl_columns_with_gsms(combined_df, gsm_order)
                    return df_synced, meta_df, gsms, gsm_order, f"Automated TAR ({tar_method}; {msg})", None, detected_type
                    
    except Exception as e:
        if 'local_tar_path' in locals() and os.path.exists(local_tar_path):
            os.remove(local_tar_path)
        return None, meta_df, gsms, gsm_order, None, f"Erro no processamento do TAR: {str(e)}", detected_type
        
    return None, meta_df, gsms, gsm_order, None, "Nenhuma matriz de contagem válida encontrada.", detected_type

# ============================================================
# META-ANALYSIS ENGINE
# ============================================================

def _safe_p_from_z(z):
    z = np.asarray(z, dtype=float)
    return 2.0 * stats.norm.sf(np.abs(z))


def _run_limma_backend(expression, design_metadata, design_formula, contrast):
    """Run the bundled limma backend and return a standardized result table."""
    rscript = shutil.which("Rscript")
    backend = os.path.join(PARENT_DIR, "r", "limma_backend.R")
    if not rscript:
        raise RuntimeError(
            "limma was selected, but Rscript is unavailable. Use the publication Docker image "
            "or install R plus Bioconductor limma. Welch remains available only as an exploratory fallback."
        )
    if not os.path.exists(backend):
        raise RuntimeError("The bundled limma backend is missing.")
    with tempfile.TemporaryDirectory(prefix="lumos_limma_") as temp_dir:
        expression_path = os.path.join(temp_dir, "expression.csv")
        metadata_path = os.path.join(temp_dir, "metadata.csv")
        output_path = os.path.join(temp_dir, "results.csv")
        expression.to_csv(expression_path, index_label="Feature")
        design_metadata.to_csv(metadata_path, index_label="sample")
        completed = subprocess.run(
            [
                rscript, backend, expression_path, metadata_path, design_formula,
                contrast[2], contrast[1], output_path,
            ],
            capture_output=True, text=True, timeout=1800, check=False,
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip()
            raise RuntimeError(f"limma backend failed: {detail[-1500:]}")
        output = pd.read_csv(output_path)
    required = {"Probe_ID", "Log2FC", "SE", "PValue", "FDR"}
    if missing := required.difference(output.columns):
        raise RuntimeError("limma output is missing columns: " + ", ".join(sorted(missing)))
    return output

def _run_single_study(
    df_matrix, c_ref, c_test, mode, apply_log, mapping=None,
    study_id="Study", alpha=0.05, design_metadata=None,
    design_formula="~condition", contrast=None, continuous_engine="limma"
):
    c_ref, c_test = list(dict.fromkeys(c_ref)), list(dict.fromkeys(c_test))
    overlap = sorted(set(c_ref).intersection(c_test))
    if overlap:
        raise ValueError("Samples cannot occur in both conditions: " + ", ".join(overlap[:8]))
    if len(c_ref) < 2 or len(c_test) < 2:
        raise ValueError("At least two unique samples per condition are required.")
    if not df_matrix.columns.is_unique:
        duplicates = df_matrix.columns[df_matrix.columns.duplicated()].astype(str).unique().tolist()
        raise ValueError("Expression matrix has duplicate sample columns: " + ", ".join(duplicates[:8]))
    cols = c_ref + c_test
    missing_columns = sorted(set(cols).difference(df_matrix.columns.astype(str)))
    if missing_columns:
        raise ValueError("Selected samples are absent from the expression matrix: " + ", ".join(missing_columns[:8]))
    df_sub = df_matrix[cols].apply(pd.to_numeric, errors='coerce')
    if df_sub.empty:
        raise ValueError(f"{study_id}: empty expression matrix.")
    if mode == "RNASeq" and df_sub.isna().any().any():
        raise ValueError(
            f"{study_id}: RNA-seq matrix contains missing or non-numeric values; "
            "missing counts cannot be replaced with zero automatically."
        )
    if mode != "RNASeq":
        df_sub = df_sub.dropna(axis=0, how='any')
        if df_sub.empty:
            raise ValueError(f"{study_id}: no complete continuous-expression features remain.")

    values = df_sub.to_numpy(dtype=float)
    is_count_like = (
        np.isfinite(values).all()
        and (values >= 0).all()
        and np.allclose(values, np.rint(values), rtol=0, atol=1e-8)
    )
    mapping_applied_to_expression = False
    analysis_unit = 'Input feature'
    if mode == "RNASeq" and is_count_like and mapping is not None and not mapping.empty:
        id_to_symbol = (
            mapping.dropna(subset=['Probe_ID', 'Symbol'])
            .drop_duplicates('Probe_ID')
            .set_index('Probe_ID')['Symbol']
            .astype(str)
            .to_dict()
        )
        annotated_counts = df_sub.copy()
        original_ids = annotated_counts.index.astype(str)
        mapped_symbols = original_ids.map(id_to_symbol)
        annotated_counts['_Symbol'] = [
            symbol if pd.notna(symbol) and str(symbol) not in {'', 'nan', 'None'} else feature_id
            for symbol, feature_id in zip(mapped_symbols, original_ids)
        ]
        df_sub = annotated_counts.groupby('_Symbol', sort=False)[cols].sum()
        df_sub.index.name = 'Symbol'
        mapping_applied_to_expression = True
        analysis_unit = 'Gene-level sum across mapped count features'
    min_samples = min(len(c_ref), len(c_test))
    if mode == "RNASeq" and is_count_like:
        mask_keep = (df_sub >= 10).sum(axis=1) >= min_samples
    else:
        row_variance = df_sub.var(axis=1, ddof=1)
        mask_keep = row_variance.notna() & row_variance.gt(0)
    df_sub = df_sub[mask_keep]
    if df_sub.empty:
        raise ValueError(f"{study_id}: pre-filtering removed all genes.")

    data_vals = df_sub.values.astype(np.float32)

    high_dynamic_range = (
        np.nanmin(data_vals) >= 0 and np.nanpercentile(data_vals, 99) > 100
    )
    if mode == "Microarray":
        microarray_vals = (
            np.log2(np.clip(data_vals, 0, None) + 1.0)
            if apply_log and high_dynamic_range
            else data_vals
        )
        data_norm_vals = quantile_normalize(microarray_vals)
    else:
        if apply_log and not is_count_like and high_dynamic_range:
            data_norm_vals = np.log2(np.clip(data_vals, 0, None) + 1.0)
        else:
            data_norm_vals = data_vals

    data_norm = pd.DataFrame(data_norm_vals, columns=cols, index=df_sub.index)
    if mode == "Microarray" and mapping is not None and not mapping.empty:
        probe_to_symbol = (
            mapping.dropna(subset=['Probe_ID', 'Symbol'])
            .drop_duplicates('Probe_ID')
            .set_index('Probe_ID')['Symbol']
            .astype(str)
            .to_dict()
        )
        annotated = data_norm.copy()
        annotated['_Symbol'] = annotated.index.astype(str).map(probe_to_symbol)
        annotated = annotated[
            annotated['_Symbol'].notna()
            & ~annotated['_Symbol'].isin(['', 'nan', 'None'])
        ]
        if not annotated.empty:
            data_norm = annotated.groupby('_Symbol', sort=False)[cols].median()
            data_norm.index.name = 'Symbol'
            mapping_applied_to_expression = True
            analysis_unit = 'Gene-level median across mapped probes'
    m1 = data_norm[c_ref].values
    m2 = data_norm[c_test].values

    if design_metadata is None:
        design_metadata = pd.DataFrame(
            {"condition": ["C"] * len(c_ref) + ["T"] * len(c_test)}, index=cols
        )
        design_formula = "~condition"
        contrast = ["condition", "T", "C"]
    design_metadata = design_metadata.loc[cols].copy()
    contrast = contrast or ["condition", "T", "C"]

    if mode == "RNASeq" and is_count_like:
        if not HAS_DESEQ2:
            raise ImportError(
                "Raw integer counts were detected, but PyDESeq2 is unavailable. "
                "Install the declared dependencies instead of applying a t-test to counts."
            )
        validated_counts = validate_count_matrix(df_sub)
        dds = DeseqDataSet(
            counts=validated_counts.T,
            metadata=design_metadata,
            design=design_formula,
            refit_cooks=True,
            quiet=True
        )
        dds.deseq2()

        try:
            dds.vst()
            pca_input = orient_gene_sample_matrix(
                np.asarray(dds.layers['vst_counts']), df_sub.shape[0], cols
            )
            gc.collect()
        except Exception:
            pca_input = np.log2(df_sub.to_numpy(dtype=float) + 1.0)

        stat_res = DeseqStats(
            dds,
            contrast=contrast,
            alpha=alpha,
            cooks_filter=True,
            independent_filter=True,
            quiet=True,
        )
        stat_res.summary()
        res_df = stat_res.results_df

        res = pd.DataFrame({
            'Probe_ID': res_df.index,
            'BaseMean': res_df['baseMean'] if 'baseMean' in res_df.columns else np.nan,
            'Log2FC': res_df['log2FoldChange'],
            'SE': res_df['lfcSE'] if 'lfcSE' in res_df.columns else np.nan,
            'PValue': res_df['pvalue'],
            'FDR': res_df['padj']
        }).dropna(subset=['Log2FC', 'PValue'])
        res['AnalysisMethod'] = 'PyDESeq2 Wald test'

        data_norm = pd.DataFrame(pca_input, columns=cols, index=df_sub.index)
    elif continuous_engine == "limma":
        res = _run_limma_backend(
            data_norm, design_metadata, design_formula, contrast
        ).dropna(subset=["Log2FC", "SE", "PValue"])
        res["AnalysisMethod"] = (
            "Quantile normalization + limma empirical Bayes linear model"
            if mode == "Microarray"
            else "limma empirical Bayes linear model"
        )
    else:
        lfc = np.nanmean(m2, axis=1) - np.nanmean(m1, axis=1)
        n1, n2 = len(c_ref), len(c_test)
        v1 = np.nanvar(m1, axis=1, ddof=1)
        v2 = np.nanvar(m2, axis=1, ddof=1)
        se = np.sqrt(v1 / n1 + v2 / n2)
        pvals = stats.ttest_ind(
            m2, m1, axis=1, equal_var=False, nan_policy='omit'
        ).pvalue
        mask = np.isfinite(pvals) & np.isfinite(se) & (se > 0)
        fdr = np.full(pvals.shape, np.nan)
        if mask.any():
            fdr[mask] = fdrcorrection(pvals[mask])[1]
        res = pd.DataFrame({
            'Probe_ID': data_norm.index,
            'MeanExpression': np.nanmean(data_norm.to_numpy(dtype=float), axis=1),
            'Log2FC': lfc,
            'SE': se,
            'PValue': pvals,
            'FDR': fdr
        }).dropna(subset=['Log2FC', 'SE', 'PValue'])
        res['AnalysisMethod'] = (
            'Quantile normalization + Welch t-test (exploratory; covariates ignored)'
            if mode == 'Microarray'
            else 'Welch t-test on continuous expression (exploratory; covariates ignored)'
        )
    if mapping_applied_to_expression:
        res['Symbol'] = res['Probe_ID']
        res['AnalysisUnit'] = analysis_unit
    elif mapping is not None:
        res = res.merge(mapping, on='Probe_ID', how='left')
        res['Symbol'] = res['Symbol'].replace(['nan', 'None'], np.nan).fillna(res['Probe_ID'])
        res['AnalysisUnit'] = 'Input feature mapped after statistical testing'
    else:
        res['Symbol'] = res['Probe_ID']
        res['AnalysisUnit'] = 'Input feature'

    res['Study'] = study_id
    res.index = pd.RangeIndex(len(res))
    res.index.name = None
    return res, data_norm

def _combine_meta_effects(study_results, min_studies=2, method="Paule-Mandel"):
    if not HAS_META:
        raise ImportError("statsmodels meta-analysis is unavailable. Install statsmodels >= 0.14.")

    all_res = pd.concat(study_results, ignore_index=True)
    all_res = all_res.replace([np.inf, -np.inf], np.nan)
    all_res = all_res.dropna(subset=['Symbol', 'Log2FC', 'SE'])
    all_res = all_res[all_res['SE'] > 0]

    rows = []
    failures = []
    for gene, g in all_res.groupby('Symbol', sort=False):
        g = g.drop_duplicates(subset=['Study'], keep='first')
        if len(g) < min_studies:
            continue

        effects = g['Log2FC'].to_numpy(dtype=float)
        variances = np.square(g['SE'].to_numpy(dtype=float))
        studies = g['Study'].astype(str).tolist()

        try:
            method_re = "iterated" if method == "Paule-Mandel" else "chi2"
            meta = combine_effects(
                effects, variances,
                method_re=method_re,
                row_names=studies,
                use_t=False,
                alpha=0.05
            )

            w = 1.0 / variances
            fixed_eff = np.sum(w * effects) / np.sum(w)
            q = float(np.sum(w * (effects - fixed_eff) ** 2))
            df_q = len(effects) - 1
            i2 = max(0.0, (q - df_q) / q) * 100.0 if q > 0 and df_q > 0 else 0.0

            sf = meta.summary_frame()
            re_row = sf.loc['random effect'] if 'random effect' in sf.index else sf.iloc[-3]
            eff = float(re_row['eff'])
            se_eff = float(re_row['sd_eff'])
            ci_low = float(re_row['ci_low'])
            ci_upp = float(re_row['ci_upp'])
            z = eff / se_eff if se_eff > 0 else np.nan
            p = float(_safe_p_from_z(z)) if np.isfinite(z) else np.nan

            tau2 = float(getattr(meta, 'tau2', np.nan))

            rows.append({
                'Symbol': gene,
                'Meta_Log2FC': eff,
                'Meta_SE': se_eff,
                'CI95_Low': ci_low,
                'CI95_High': ci_upp,
                'PValue': p,
                'FDR': np.nan,
                'I2_percent': i2,
                'Tau2': tau2,
                'Q': q,
                'Q_df': df_q,
                'Q_PValue': float(stats.chi2.sf(q, df_q)) if df_q > 0 else np.nan,
                'N_Studies': len(g),
                'Studies': ', '.join(studies),
                'Estimator': method,
            })
        except Exception as exc:
            failures.append(f"{gene}: {type(exc).__name__}: {exc}")
            continue

    meta_df = pd.DataFrame(rows)
    if meta_df.empty:
        return meta_df

    valid = meta_df['PValue'].notna()
    meta_df.loc[valid, 'FDR'] = fdrcorrection(meta_df.loc[valid, 'PValue'].values)[1]
    meta_df['Status'] = 'Not Significant'
    meta_df.attrs['failures'] = failures
    return meta_df.sort_values(['FDR', 'PValue'], na_position='last').reset_index(drop=True)

def _make_meta_volcano(meta_df, fc_thr, fdr_thr=0.05):
    df = meta_df.copy()
    df['-log10(FDR)'] = -np.log10(df['FDR'].clip(lower=1e-300))
    df['Status'] = 'Not Significant'
    sig = df['FDR'] < fdr_thr
    df.loc[sig & (df['Meta_Log2FC'] >= fc_thr), 'Status'] = 'Up-regulated'
    df.loc[sig & (df['Meta_Log2FC'] <= -fc_thr), 'Status'] = 'Down-regulated'
    return px.scatter(
        df, x='Meta_Log2FC', y='-log10(FDR)', color='Status', hover_name='Symbol',
        color_discrete_map={
            'Not Significant': '#cccccc',
            'Up-regulated': 'firebrick',
            'Down-regulated': 'dodgerblue'
        }
    ).update_layout(
        template='simple_white', height=600,
        xaxis_title='Meta-analysis log2 Fold Change',
        yaxis_title='-log10 adjusted p-value'
    )

def _make_forest_plot(gene, all_study_results, meta_df):
    gene_rows = []
    for sr in all_study_results:
        x = sr[sr['Symbol'].astype(str).str.upper() == str(gene).upper()]
        if not x.empty:
            r = x.iloc[0]
            gene_rows.append((r['Study'], r['Log2FC'], r['SE']))

    m = meta_df[meta_df['Symbol'].astype(str).str.upper() == str(gene).upper()]
    if m.empty or not gene_rows:
        return None
    mr = m.iloc[0]

    labels = [x[0] for x in gene_rows] + ['Random-effects meta-analysis']
    effects = [x[1] for x in gene_rows] + [mr['Meta_Log2FC']]
    lows = [x[1] - 1.96 * x[2] for x in gene_rows] + [mr['CI95_Low']]
    highs = [x[1] + 1.96 * x[2] for x in gene_rows] + [mr['CI95_High']]
    y = list(range(len(labels)))

    fig = go.Figure()
    for yi, lo, hi, eff in zip(y, lows, highs, effects):
        fig.add_trace(go.Scatter(x=[lo, hi], y=[yi, yi], mode='lines', showlegend=False))
        fig.add_trace(go.Scatter(x=[eff], y=[yi], mode='markers', showlegend=False))
    fig.add_vline(x=0, line_dash='dash', line_color='black')
    fig.update_layout(
        template='simple_white', height=max(350, 80 * len(labels)),
        title=f'Forest plot: {gene}',
        yaxis=dict(tickmode='array', tickvals=y, ticktext=labels, autorange='reversed'),
        xaxis_title='log2 Fold Change'
    )
    return fig


# ============================================================
# APP PRINCIPAL E RENDERIZAÇÃO
# ============================================================

def _reset_analysis_state():
    keys = [
        'analysis_done', 'res', 'norm_df', 'df_diff', 'meta_results',
        'study_results', 'meta_norm_df', 'fig_v', 'fig_ma', 'fig_pca', 'fig_up',
        'fig_down', 'fig_h', 'fig_meta_v', 'fig_forest', 'genes_significativos',
        'k_res', 'g_res', 'string_df', 'tf_regulators',
        'mt_studies', 'mt_gse_ids', 'mt_results', 'mt_valid_gses',
        'mt_study_results', 'mt_norm_df', 'mt_configured',
        'single_contrasts', 'single_contrast_results',
        'single_contrast_ids', 'single_contrast_name_map',
        'active_single_contrast',
    ]
    for key in keys:
        st.session_state.pop(key, None)
    for key in list(st.session_state):
        if key.startswith(("mt_design_", "mt_ctrl_", "mt_dis_", "mt_label_")):
            st.session_state.pop(key, None)

def _study_label_cols(meta_df):
    candidatas = ['source_name_ch1', 'source', 'title', 'characteristics_ch1']
    default_cols = [c for c in meta_df.columns if any(cand in c.lower() for cand in candidatas)]
    if not default_cols and len(meta_df.columns) > 1:
        default_cols = [meta_df.columns[1]]
    return default_cols


def _unique_sample_labels(meta_df, label_columns):
    """Return one UI row per GSM with an accession-prefixed unique label."""
    return make_unique_sample_labels(meta_df, label_columns)


def _build_study_design(study, c_ref, c_test):
    meta = _unique_sample_labels(study["meta_df"], [])
    selected_accessions = list(dict.fromkeys(study["c_ref"] + study["c_test"]))
    selected_columns = c_ref + c_test
    lookup = dict(zip(selected_accessions, selected_columns))
    meta = meta[meta["Accession"].isin(selected_accessions)].copy()
    meta["_matrix_sample"] = meta["Accession"].map(lookup)
    settings = study.get("design_settings", {})
    return build_design_metadata(
        sample_ids=c_ref + c_test,
        conditions=["Control"] * len(c_ref) + ["Disease"] * len(c_test),
        source_metadata=meta,
        sample_key="_matrix_sample",
        categorical_covariates=settings.get("categorical", []),
        continuous_covariates=settings.get("continuous", []),
        batch_column=settings.get("batch"),
        pair_column=settings.get("pair"),
    )

def _render_single_dataset(mode, apply_log, gene_area, p_thr, fc_thr, max_plot):
    if 'meta_df' not in st.session_state:
        st.stop()

    df_matrix = st.session_state.get('df')
    if df_matrix is None:
        st.error("Matrix data is unavailable for this dataset. Please check the GSE ID or try a different dataset.")
        st.stop()

    meta = st.session_state['meta_df'].copy()
    label_cols = st.session_state.get('label_cols', ['Accession', 'Title'])
    meta = _unique_sample_labels(meta, label_cols)
    label_to_gsm = dict(zip(meta['display_label'], meta['Accession']))

    st.subheader("🛠️ Group Management")
    if 'groups' not in st.session_state:
        st.session_state['groups'] = {}

    c_i, c_b = st.columns([3, 1])
    def add_group_callback():
        new_g_name = st.session_state.new_group_input.strip()
        if new_g_name:
            st.session_state['groups'][new_g_name] = []
            st.session_state.new_group_input = ""
    c_i.text_input("New Group Name:", key="new_group_input")
    c_b.button("➕ Add Group", width="stretch", on_click=add_group_callback)

    cols_g = st.columns(2)
    for i, (g_name, g_samples) in enumerate(list(st.session_state['groups'].items())):
        avail = [lab for lab in meta['display_label'] if lab not in sum(st.session_state['groups'].values(), []) or lab in g_samples]
        st.session_state['groups'][g_name] = cols_g[i % 2].multiselect(
            f"**{g_name}**", avail, default=[s for s in g_samples if s in avail], key=f"s_{g_name}"
        )
        if cols_g[i % 2].button("🗑️ Remove", key=f"del_{g_name}"):
            del st.session_state['groups'][g_name]
            st.rerun()

    matrix_cols = set(df_matrix.columns.astype(str))

    if len(st.session_state['groups']) >= 2:
        group_names = list(st.session_state['groups'].keys())
        c1, c2 = st.columns(2)
        ref_g = c1.selectbox("Reference (Control):", group_names)
        test_g = c2.selectbox("Test:", [g for g in group_names if g != ref_g])
        default_contrast_label = f"{test_g} vs {ref_g}"
        contrast_key = re.sub(r"[^A-Za-z0-9_]+", "_", default_contrast_label)
        contrast_label = st.text_input(
            "Analysis label",
            value=default_contrast_label,
            key=f"single_contrast_label_{contrast_key}",
            help="Use a unique label to retain alternative models of the same comparison.",
        ).strip()

        excluded_design_cols = {"Accession", "display_label"}
        design_columns = [c for c in meta.columns if c not in excluded_design_cols]
        with st.expander("🧪 Experimental design: batch, covariates and pairing", expanded=False):
            st.caption(
                "Lumos fits these variables before the condition effect. Do not select a variable "
                "that is perfectly confounded with the comparison."
            )
            dc1, dc2 = st.columns(2)
            batch_col = dc1.selectbox("Batch variable (optional)", [None] + design_columns)
            pair_col = dc2.selectbox("Subject / pair identifier (optional)", [None] + design_columns)
            unavailable = {x for x in [batch_col, pair_col] if x}
            remaining = [c for c in design_columns if c not in unavailable]
            categorical_covariates = st.multiselect(
                "Additional categorical covariates", remaining
            )
            continuous_covariates = st.multiselect(
                "Continuous covariates (numeric)",
                [c for c in remaining if c not in categorical_covariates],
            )
            detected_for_engine = st.session_state.get('detected_type', mode)
            if detected_for_engine == "RNASeq":
                continuous_engine = "pydeseq2"
                st.info("Raw integer counts will use PyDESeq2 with the displayed design formula.")
            else:
                engine_label = st.radio(
                    "Continuous-expression engine",
                    ["limma empirical Bayes (publication)", "Welch t-test (exploratory)"],
                )
                continuous_engine = "limma" if engine_label.startswith("limma") else "welch"

        if st.button("🔥 Run Analysis", width="stretch"):
            with st.spinner("Processing Differential Expression..."):
                if not contrast_label:
                    st.error("Enter a non-empty analysis label.")
                    st.stop()
                try:
                    c_ref = resolve_sample_columns(
                        df_matrix.columns,
                        [label_to_gsm[lab] for lab in st.session_state['groups'][ref_g]],
                    )
                    c_test = resolve_sample_columns(
                        df_matrix.columns,
                        [label_to_gsm[lab] for lab in st.session_state['groups'][test_g]],
                    )
                except ValueError as exc:
                    st.error(f"Sample alignment failed: {exc}")
                    st.stop()
                if len(c_ref) < 2 or len(c_test) < 2:
                    st.error("Select at least two unambiguous samples in each group.")
                    st.stop()
                if not c_ref or not c_test:
                    st.error(f"❌ Mapping error. Columns available: {list(matrix_cols)[:5]}")
                    st.stop()

                mapping = st.session_state.get('mapping')
                try:
                    meta_for_design = meta.copy()
                    resolved_all = resolve_sample_columns(
                        df_matrix.columns, meta_for_design["Accession"].tolist()
                    )
                    lookup = dict(zip(meta_for_design["Accession"], resolved_all))
                    meta_for_design["_matrix_sample"] = meta_for_design["Accession"].map(lookup)
                    meta_for_design = meta_for_design.dropna(subset=["_matrix_sample"])
                    design_metadata, design_formula, contrast, design_name_map = build_design_metadata(
                        sample_ids=c_ref + c_test,
                        conditions=[ref_g] * len(c_ref) + [test_g] * len(c_test),
                        source_metadata=meta_for_design,
                        sample_key="_matrix_sample",
                        categorical_covariates=categorical_covariates,
                        continuous_covariates=continuous_covariates,
                        batch_column=batch_col,
                        pair_column=pair_col,
                    )
                    analytical_mode = st.session_state.get('detected_type')
                    if analytical_mode not in {'RNASeq', 'Microarray'}:
                        analytical_mode = mode
                    res, data_norm = _run_single_study(
                        df_matrix, c_ref, c_test, analytical_mode, apply_log,
                        mapping=mapping, study_id=st.session_state.get('gse_id', 'GSE'),
                        alpha=p_thr,
                        design_metadata=design_metadata,
                        design_formula=design_formula,
                        contrast=contrast,
                        continuous_engine=continuous_engine,
                    )
                except Exception as e:
                    st.error(f"❌ Analysis failed: {e}")
                    st.stop()

                st.session_state.update({
                    'res': res, 'norm_df': data_norm, 'analysis_done': True,
                    'rn': ref_g, 'tn': test_g, 'c_ref': c_ref, 'c_test': c_test,
                    'design_metadata': design_metadata, 'design_formula': design_formula,
                    'design_contrast': contrast, 'design_name_map': design_name_map,
                })
                collection = st.session_state.setdefault('single_contrasts', {})
                collection[contrast_label] = {
                    'res': res.copy(), 'norm_df': data_norm.copy(),
                    'rn': ref_g, 'tn': test_g,
                    'c_ref': list(c_ref), 'c_test': list(c_test),
                    'design_metadata': design_metadata.copy(),
                    'design_formula': design_formula,
                    'design_contrast': contrast,
                    'design_name_map': dict(design_name_map),
                }
                st.session_state['active_single_contrast'] = contrast_label
                st.session_state['single_contrast_view'] = contrast_label
                st.rerun()

    if st.session_state.get('analysis_done'):
        contrast_collection = st.session_state.get('single_contrasts', {})
        if contrast_collection:
            st.divider()
            st.subheader("🧪 Analyses saved for this dataset")
            sc1, sc2 = st.columns([4, 1])
            labels = list(contrast_collection)
            current_label = st.session_state.get('active_single_contrast', labels[-1])
            if current_label not in labels:
                current_label = labels[-1]
            active_label = sc1.selectbox(
                "Analysis shown in the figures below",
                labels,
                index=labels.index(current_label),
                key="single_contrast_view",
            )
            active = contrast_collection[active_label]
            st.session_state.update({
                'res': active['res'], 'norm_df': active['norm_df'],
                'rn': active['rn'], 'tn': active['tn'],
                'c_ref': active['c_ref'], 'c_test': active['c_test'],
                'design_metadata': active['design_metadata'],
                'design_formula': active['design_formula'],
                'design_contrast': active['design_contrast'],
                'design_name_map': active['design_name_map'],
                'active_single_contrast': active_label,
            })
            if sc2.button("🗑️ Remove analysis", key="remove_single_contrast"):
                contrast_collection.pop(active_label, None)
                if contrast_collection:
                    replacement_label = list(contrast_collection)[-1]
                    replacement = contrast_collection[replacement_label]
                    st.session_state.update({
                        'res': replacement['res'], 'norm_df': replacement['norm_df'],
                        'rn': replacement['rn'], 'tn': replacement['tn'],
                        'c_ref': replacement['c_ref'], 'c_test': replacement['c_test'],
                        'design_metadata': replacement['design_metadata'],
                        'design_formula': replacement['design_formula'],
                        'design_contrast': replacement['design_contrast'],
                        'design_name_map': replacement['design_name_map'],
                        'active_single_contrast': replacement_label,
                    })
                    st.session_state['single_contrast_view'] = replacement_label
                else:
                    for key in ('res', 'norm_df', 'analysis_done', 'active_single_contrast'):
                        st.session_state.pop(key, None)
                    st.session_state.pop('single_contrast_view', None)
                st.rerun()
            st.caption(
                "Run another comparison above to add or replace an analysis. "
                "The same label is replaced deterministically."
            )

        res, df_norm = st.session_state['res'].copy(), st.session_state['norm_df']
        if gene_area:
            res = res[res['Symbol'].str.upper().isin([g.strip().upper() for g in gene_area.split('\n') if g.strip()])]
        res['-log10(FDR)'] = -np.log10(res['FDR'].fillna(1).clip(lower=1e-300))
        res['Status'] = classify_differential_expression(
            res, 'Log2FC', 'FDR', p_thr, fc_thr
        )
        df_diff = res[res['Status'] != 'Not Significant'].sort_values('FDR')
        st.session_state['df_diff'] = df_diff
        st.session_state['genes_significativos'] = df_diff['Symbol'].unique().tolist()
        up_list = df_diff[df_diff['Status'] == 'Up-regulated']
        down_list = df_diff[df_diff['Status'] == 'Down-regulated']

        st.subheader(f"🚀 Results: {st.session_state['tn']} vs {st.session_state['rn']}")
        m1, m2, m3 = st.columns(3)
        m1.metric("Total DEGs", len(df_diff)); m2.metric("Up-Regulated 📈", len(up_list)); m3.metric("Down-Regulated 📉", len(down_list))
        t1, t2, t3, t4 = st.tabs(["📊 Volcano & PCA", "🔝 Top DEGs", "🔥 Heatmap", "📋 Data Table"])
        
        with t1:
            fig_v = px.scatter(res, x='Log2FC', y='-log10(FDR)', color='Status', hover_name='Symbol',
                color_discrete_map={'Not Significant':'#cccccc','Up-regulated':'firebrick','Down-regulated':'dodgerblue'})
            fig_v.add_vline(x=fc_thr, line_dash='dash', line_color='black'); fig_v.add_vline(x=-fc_thr, line_dash='dash', line_color='black'); fig_v.add_hline(y=-np.log10(p_thr), line_dash='dash', line_color='black')
            st.session_state['fig_v'] = fig_v.update_layout(template='simple_white', height=600, xaxis_title='log2 Fold Change', yaxis_title='-log10 adjusted p-value')
            st.plotly_chart(st.session_state['fig_v'], width="stretch")
            mean_col = 'BaseMean' if 'BaseMean' in res.columns else 'MeanExpression'
            if mean_col in res.columns:
                ma_df = res[pd.to_numeric(res[mean_col], errors='coerce').gt(0)].copy()
                ma_df['Mean abundance'] = np.log10(pd.to_numeric(ma_df[mean_col], errors='coerce'))
                fig_ma = px.scatter(
                    ma_df, x='Mean abundance', y='Log2FC', color='Status', hover_name='Symbol',
                    color_discrete_map={'Not Significant':'#cccccc','Up-regulated':'firebrick','Down-regulated':'dodgerblue'},
                    title='MA-style effect–abundance plot'
                ).update_layout(template='simple_white', height=550, xaxis_title=f'log10({mean_col})')
                fig_ma.add_hline(y=0, line_dash='dash', line_color='black')
                st.session_state['fig_ma'] = fig_ma
                st.plotly_chart(fig_ma, width="stretch")
            try:
                pca = PCA(n_components=2); pc = pca.fit_transform(df_norm.fillna(0).T)
                df_pc = pd.DataFrame(pc, columns=['PC1','PC2'], index=df_norm.columns)
                df_pc['Group'] = [st.session_state['rn'] if x in st.session_state['c_ref'] else st.session_state['tn'] for x in df_pc.index]
                variance = pca.explained_variance_ratio_ * 100
                st.session_state['fig_pca'] = px.scatter(
                    df_pc, x='PC1', y='PC2', color='Group', hover_name=df_pc.index
                ).update_layout(
                    template='simple_white',
                    xaxis_title=f'PC1 ({variance[0]:.1f}%)',
                    yaxis_title=f'PC2 ({variance[1]:.1f}%)'
                )
                st.plotly_chart(st.session_state['fig_pca'], width="stretch")
            except Exception: pass
            
        with t2:
            ca, cb = st.columns(2)
            if not up_list.empty:
                st.session_state['fig_up'] = px.bar(up_list.head(20), x='Symbol', y='Log2FC', color='Log2FC', color_continuous_scale='Reds', title='Top 20 UP'); ca.plotly_chart(st.session_state['fig_up'], width="stretch")
            if not down_list.empty:
                st.session_state['fig_down'] = px.bar(down_list.head(20), x='Symbol', y='Log2FC', color='Log2FC', color_continuous_scale='Blues_r', title='Top 20 DOWN'); cb.plotly_chart(st.session_state['fig_down'], width="stretch")
                
        with t3:
            if not df_diff.empty:
                valid_ids = [idx for idx in df_diff.head(max_plot)['Probe_ID'] if idx in df_norm.index]
                
                c_ref = st.session_state['c_ref']
                c_test = st.session_state['c_test']
                rn = st.session_state['rn']
                tn = st.session_state['tn']
                
                ordered_cols = c_ref + c_test
                ordered_cols = [c for c in ordered_cols if c in df_norm.columns]
                ordered_cols = list(dict.fromkeys(ordered_cols))
                
                df_slice = df_norm.loc[valid_ids, ordered_cols]
                df_slice = df_slice.loc[:, ~df_slice.columns.duplicated()] 
                
                h_mat = df_slice.values
                h_z = (h_mat - np.mean(h_mat, axis=1, keepdims=True)) / (np.std(h_mat, axis=1, keepdims=True) + 1e-9)
                
                y_labels = df_diff.head(len(valid_ids))['Symbol'].tolist()
                x_labels = [f"{rn} ({c})" if c in c_ref else f"{tn} ({c})" for c in df_slice.columns]

                col1, col2 = st.columns(2)
                cluster_r = col1.checkbox("Cluster Genes (Rows)", value=True, key="cluster_r")
                cluster_c = col2.checkbox("Cluster Samples (Columns)", value=False, key="cluster_c")
                st.caption("Tip: Keep 'Cluster Samples' disabled to naturally separate samples by your selected groups.")

                if cluster_r and len(valid_ids) > 1:
                    try:
                        row_link = hierarchy.linkage(distance.pdist(h_z, metric='euclidean'), method='ward')
                        row_order = hierarchy.leaves_list(row_link)
                        h_z = h_z[row_order, :]
                        y_labels = [y_labels[i] for i in row_order]
                    except Exception as e:
                        st.warning(f"Row clustering failed: {e}")

                if cluster_c and h_z.shape[1] > 1:
                    try:
                        col_link = hierarchy.linkage(distance.pdist(h_z.T, metric='euclidean'), method='ward')
                        col_order = hierarchy.leaves_list(col_link)
                        h_z = h_z[:, col_order]
                        x_labels = [x_labels[i] for i in col_order]
                    except Exception as e:
                        st.warning(f"Column clustering failed: {e}")

                st.session_state['fig_h'] = go.Figure(data=go.Heatmap(
                    z=h_z, x=x_labels, y=y_labels, colorscale='RdBu_r', zmid=0
                )).update_layout(
                    height=max(400, len(valid_ids) * 20),
                    xaxis_tickangle=-45
                )
                
                st.plotly_chart(st.session_state['fig_h'], width="stretch")
                
        with t4:
            st.dataframe(res, width="stretch")
            st.download_button('📥 Download complete results', data=res.to_csv(index=False), file_name=f"{st.session_state.get('gse_id','DE')}_complete_results.csv", mime='text/csv')
            manifest = {
                'lumos_version': LUMOS_VERSION,
                'generated_utc': datetime.utcnow().replace(microsecond=0).isoformat() + 'Z',
                'dataset': st.session_state.get('gse_id'),
                'matrix_source': st.session_state.get('matrix_source'),
                'detected_experiment_type': st.session_state.get('detected_type'),
                'comparison': {'reference': st.session_state['rn'], 'test': st.session_state['tn']},
                'samples': {'reference': st.session_state['c_ref'], 'test': st.session_state['c_test']},
                'design': {
                    'formula': st.session_state.get('design_formula', '~condition'),
                    'contrast': st.session_state.get('design_contrast'),
                    'source_to_model_column': st.session_state.get('design_name_map', {}),
                    'metadata': st.session_state.get('design_metadata', pd.DataFrame()).reset_index().to_dict('records'),
                },
                'analysis_method': sorted(res['AnalysisMethod'].dropna().astype(str).unique().tolist()) if 'AnalysisMethod' in res else [],
                'thresholds': {'FDR': p_thr, 'absolute_log2FC': fc_thr},
                'runtime': {'python': sys.version.split()[0], 'platform': platform.platform()},
                'packages': package_versions(['streamlit', 'pydeseq2', 'pandas', 'numpy', 'scipy', 'statsmodels', 'scikit-learn']),
            }
            st.download_button('📥 Download analysis manifest', data=json.dumps(manifest, indent=2), file_name=f"{st.session_state.get('gse_id','DE')}_manifest.json", mime='application/json')

        contrast_collection = st.session_state.get('single_contrasts', {})
        if contrast_collection:
            st.divider()
            st.header("🧬 Integrated view of analyses from this dataset")
            st.caption(
                "Contrasts are aligned by gene and kept as separate effect estimates. "
                "Because they may reuse samples, Lumos does not treat them as independent "
                "studies and does not calculate a pooled meta-analytic effect."
            )
            summary_rows = []
            for label, entry in contrast_collection.items():
                current = entry['res'].copy()
                current['Status'] = classify_differential_expression(
                    current, 'Log2FC', 'FDR', p_thr, fc_thr
                )
                summary_rows.append({
                    'Analysis': label,
                    'Reference samples': len(entry['c_ref']),
                    'Test samples': len(entry['c_test']),
                    'DEGs': int(current['Status'].ne('Not Significant').sum()),
                    'Up': int(current['Status'].eq('Up-regulated').sum()),
                    'Down': int(current['Status'].eq('Down-regulated').sum()),
                    'Model': entry.get('design_formula', '~condition'),
                })
            st.dataframe(pd.DataFrame(summary_rows), width="stretch", hide_index=True)

            if len(contrast_collection) < 2:
                st.info("Run at least one more comparison to create the combined gene-alignment table for APP and PG.")
            else:
                try:
                    aligned, safe_to_label = align_named_contrasts(
                        {label: entry['res'] for label, entry in contrast_collection.items()},
                        p_thr, fc_thr,
                    )
                except ValueError as exc:
                    st.error(f"Could not align the saved analyses: {exc}")
                else:
                    safe_ids = list(safe_to_label)
                    st.session_state['single_contrast_results'] = aligned
                    st.session_state['single_contrast_ids'] = safe_ids
                    st.session_state['single_contrast_name_map'] = safe_to_label
                    st.success(
                        "Combined analysis is ready. Open APP to choose the consistent signature "
                        "or the profile of an individual contrast; PG will receive the APP selection."
                    )

                    ai1, ai2, ai3, ai4 = st.tabs([
                        "📋 Gene alignment", "🔥 Effect heatmap",
                        "🔁 Pairwise effects", "📥 Export",
                    ])
                    with ai1:
                        pattern_counts = aligned['Cross_Contrast_Pattern'].value_counts()
                        mc1, mc2, mc3 = st.columns(3)
                        mc1.metric("Genes aligned", f"{len(aligned):,}")
                        mc2.metric(
                            "Significant in every analysis",
                            int(aligned['N_Significant_Contrasts'].eq(len(safe_ids)).sum()),
                        )
                        mc3.metric("Discordant", int(pattern_counts.get('Discordant', 0)))
                        st.dataframe(
                            aligned.sort_values(
                                ['N_Significant_Contrasts', 'Symbol'],
                                ascending=[False, True],
                            ),
                            width="stretch", hide_index=True, height=520,
                        )

                    with ai2:
                        heat_candidates = aligned[
                            aligned['N_Significant_Contrasts'].gt(0)
                        ].copy()
                        if heat_candidates.empty:
                            st.info("No gene passed the current thresholds in any saved analysis.")
                        else:
                            effect_columns = [f"Log2FC_{safe}" for safe in safe_ids]
                            heat_candidates['max_abs_effect'] = (
                                heat_candidates[effect_columns]
                                .apply(pd.to_numeric, errors='coerce').abs().max(axis=1)
                            )
                            heat_candidates = heat_candidates.nlargest(max_plot, 'max_abs_effect')
                            z = heat_candidates[effect_columns].apply(pd.to_numeric, errors='coerce').T
                            z.index = [safe_to_label[safe] for safe in safe_ids]
                            z.columns = heat_candidates['Symbol']
                            fig_contrast_heat = go.Figure(go.Heatmap(
                                z=z.values, x=z.columns, y=z.index,
                                colorscale='RdBu_r', zmid=0,
                                colorbar_title='log2FC',
                            )).update_layout(
                                template='simple_white', height=max(400, 70 * len(safe_ids)),
                                xaxis_title='Gene', yaxis_title='Analysis',
                            )
                            st.plotly_chart(fig_contrast_heat, width="stretch")

                    with ai3:
                        pair_cols = st.columns(2)
                        left_safe = pair_cols[0].selectbox(
                            "X-axis analysis", safe_ids,
                            format_func=lambda x: safe_to_label[x],
                            key="single_pair_left",
                        )
                        right_options = [x for x in safe_ids if x != left_safe]
                        right_safe = pair_cols[1].selectbox(
                            "Y-axis analysis", right_options,
                            format_func=lambda x: safe_to_label[x],
                            key="single_pair_right",
                        )
                        pair_df = aligned.dropna(
                            subset=[f"Log2FC_{left_safe}", f"Log2FC_{right_safe}"]
                        ).copy()
                        fig_pair = px.scatter(
                            pair_df,
                            x=f"Log2FC_{left_safe}", y=f"Log2FC_{right_safe}",
                            color='Cross_Contrast_Pattern', hover_name='Symbol',
                            title=f"{safe_to_label[right_safe]} versus {safe_to_label[left_safe]}",
                        ).update_layout(template='simple_white', height=600)
                        fig_pair.add_hline(y=0, line_dash='dash', line_color='grey')
                        fig_pair.add_vline(x=0, line_dash='dash', line_color='grey')
                        st.plotly_chart(fig_pair, width="stretch")

                    with ai4:
                        st.download_button(
                            "📥 Download combined contrast table",
                            aligned.to_csv(index=False),
                            file_name=f"{st.session_state.get('gse_id','GSE')}_multi_contrast.csv",
                            mime='text/csv',
                        )
                        contrast_manifest = {
                            'lumos_version': LUMOS_VERSION,
                            'dataset': st.session_state.get('gse_id'),
                            'analysis_scope': 'single_dataset_multi_contrast',
                            'integration_method': 'gene_alignment_without_effect_pooling',
                            'thresholds': {'FDR': p_thr, 'absolute_log2FC': fc_thr},
                            'contrast_column_map': safe_to_label,
                            'analyses': summary_rows,
                        }
                        st.download_button(
                            "📥 Download multi-contrast manifest",
                            json.dumps(contrast_manifest, indent=2),
                            file_name=f"{st.session_state.get('gse_id','GSE')}_multi_contrast_manifest.json",
                            mime='application/json',
                        )


def _render_meta_mode(mode, apply_log, p_thr, fc_thr, min_studies, meta_method, gene_area="", max_plot_meta=50):
    st.header("🧬 Multi-Dataset Meta-analysis")
    st.caption(
        "Use two or more independent datasets from the same tissue and a comparable "
        "biological contrast. Each dataset is analyzed separately before gene-level "
        "effects are combined by random-effects meta-analysis."
    )

    gse_text = st.text_area("GSE IDs (one per line or separated by commas):", placeholder="GSE186063\nGSE205748")
    gse_ids = [x.strip().upper() for x in re.split(r'[\n,;]+', gse_text) if x.strip()]
    gse_ids = list(dict.fromkeys(gse_ids))

    if st.button("📥 Fetch All Datasets", disabled=len(gse_ids) < 2, width="stretch"):
        _reset_analysis_state()
        studies = {}
        progress = st.progress(0, text="Fetching datasets...")
        for i, gse in enumerate(gse_ids):
            try:
                df, meta_df, gsms, gsm_order, source, err, detected_type = get_geo_full_data(gse, mode)
                if meta_df is None or df is None:
                    st.warning(f"⚠️ {gse}: {err or 'dataset unavailable'}")
                    continue
                id_type = detect_index_type(df.index.astype(str).tolist())
                analytical_type = detected_type if detected_type in {'RNASeq', 'Microarray'} else mode
                mapping, msg = (get_gene_mapping_rnaseq(tuple(df.index.astype(str).tolist()), id_type) if analytical_type == 'RNASeq' else get_gene_mapping_microarray(gse))
                studies[gse] = {
                    'gse_id': gse, 'df': df, 'meta_df': meta_df, 'gsm_order': gsm_order,
                    'source': source, 'detected_type': analytical_type, 'mapping': mapping,
                    'mapping_msg': msg, 'label_cols': _study_label_cols(meta_df),
                    'groups': {'Control': [], 'Disease': []}
                }
            except Exception as e:
                st.warning(f"⚠️ {gse}: {e}")
            progress.progress((i+1)/len(gse_ids), text=f"Fetched {i+1}/{len(gse_ids)}")
        progress.empty()
        st.session_state['meta_studies'] = studies
        st.session_state['meta_gse_ids'] = list(studies.keys())
        st.session_state['meta_ready'] = True
        st.rerun()

    studies = st.session_state.get('meta_studies', {})
    if not studies:
        st.info("Enter at least two GSE IDs and click 'Fetch All Datasets'.")
        return

    st.success(f"Loaded {len(studies)} dataset(s). Configure Control and Disease samples for each study below.")
    configured = {}
    for gse, study in studies.items():
        with st.expander(f"📚 {gse} — {study['detected_type']} — {study['df'].shape[0]} genes × {study['df'].shape[1]} samples", expanded=True):
            meta = study['meta_df'].copy()
            meta_cols = list(meta.columns)
            valid_defaults = [c for c in study['label_cols'] if c in meta_cols]
            
            selected_labels = st.multiselect(
                f"Select metadata columns to display in dropdowns for {gse}:", 
                meta_cols, 
                default=valid_defaults, 
                key=f"label_sel_{gse}"
            )
            
            if not selected_labels:
                selected_labels = ['Accession'] if 'Accession' in meta_cols else [meta_cols[0]]
            
            meta = _unique_sample_labels(meta, selected_labels)

            label_to_gsm = dict(zip(meta['display_label'], meta['Accession']))
            sample_labels = meta['display_label'].tolist()
            
            c1, c2 = st.columns(2)
            ctrl = c1.multiselect("Control samples", sample_labels, key=f"meta_ctrl_{gse}")
            dis = c2.multiselect("Disease/Test samples", [x for x in sample_labels if x not in ctrl], key=f"meta_dis_{gse}")
            
            st.dataframe(meta.drop(columns=['display_label'], errors='ignore'), width="stretch", height=180)
            
            configured[gse] = {
                **study,
                'c_ref': [label_to_gsm[x] for x in ctrl],
                'c_test': [label_to_gsm[x] for x in dis],
            }

    if st.button("🚀 Run Multi-Dataset Meta-analysis", width="stretch"):
        if not HAS_META:
            st.error("statsmodels meta-analysis module is unavailable. Please upgrade statsmodels (>= 0.14 recommended).")
            return
        study_results = []
        study_norms = {}
        errors = []
        progress = st.progress(0, text="Running study-level analyses...")
        valid_count = 0
        for i, (gse, study) in enumerate(configured.items()):
            if len(study['c_ref']) < 2 or len(study['c_test']) < 2:
                errors.append(f"{gse}: select at least 2 Control and 2 Disease samples.")
                continue
            try:
                c_ref = resolve_sample_columns(study['df'].columns, study['c_ref'])
                c_test = resolve_sample_columns(study['df'].columns, study['c_test'])
                design_metadata, design_formula, contrast, _ = _build_study_design(
                    study, c_ref, c_test
                )
                sr, norm = _run_single_study(
                    study['df'], c_ref, c_test,
                    study['detected_type'], apply_log,
                    mapping=study['mapping'], study_id=gse, alpha=p_thr,
                    design_metadata=design_metadata, design_formula=design_formula,
                    contrast=contrast,
                    continuous_engine=study.get('design_settings', {}).get('engine', 'limma'),
                )
                study_results.append(sr)
                study_norms[gse] = norm
                valid_count += 1
            except Exception as e:
                errors.append(f"{gse}: {e}")
            progress.progress((i+1)/len(configured), text=f"Analyzed {i+1}/{len(configured)}")
        progress.empty()

        if errors:
            for e in errors: st.warning(f"⚠️ {e}")
        if valid_count < min_studies:
            st.error(f"At least {min_studies} valid studies are required for the meta-analysis.")
            return

        try:
            meta_df = _combine_meta_effects(study_results, min_studies=min_studies, method=meta_method)
        except Exception as e:
            st.error(f"❌ Meta-analysis failed: {e}")
            return
        if meta_df.empty:
            st.error("No genes had sufficient information across the selected datasets.")
            return
        meta_failures = meta_df.attrs.get('failures', [])
        if meta_failures:
            st.warning(
                f"Meta-analysis failed for {len(meta_failures)} gene(s); no estimator substitution was performed. "
                "Successful genes remain available in the results table."
            )

        meta_df['Status'] = 'Not Significant'
        meta_df.loc[(meta_df['FDR'] < p_thr) & (meta_df['Meta_Log2FC'] >= fc_thr), 'Status'] = 'Up-regulated'
        meta_df.loc[(meta_df['FDR'] < p_thr) & (meta_df['Meta_Log2FC'] <= -fc_thr), 'Status'] = 'Down-regulated'

        st.session_state.update({
            'analysis_done': True,
            'meta_results': meta_df,
            'study_results': study_results,
            'meta_norm_df': study_norms,
            'gse_id': 'META-' + '-'.join([x.replace('GSE','') for x in configured]),
            'rn': 'Control', 'tn': 'Disease',
            'meta_mode': True,
            'meta_method': meta_method,
            'meta_configured': configured
        })
        st.rerun()

    if st.session_state.get('meta_results') is not None:
        meta_df = st.session_state['meta_results'].copy()

        target_genes = []
        if gene_area.strip():
            target_genes = [g.strip().upper() for g in gene_area.split('\n') if g.strip()]
            meta_df = meta_df[meta_df['Symbol'].str.upper().isin(target_genes)].copy()

        meta_degs = meta_df[meta_df['Status'] != 'Not Significant'].copy()
        
        st.divider(); st.header("📊 Meta-analysis Results")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Datasets", len(st.session_state.get('study_results', [])))
        m2.metric("Meta-genes", len(meta_degs))
        m3.metric("Up-regulated", int((meta_degs['Status']=='Up-regulated').sum()))
        m4.metric("Down-regulated", int((meta_degs['Status']=='Down-regulated').sum()))

        tab1, tab2, tab3, tab4, tab5 = st.tabs(["🌋 Meta-Volcano", "🌲 Forest Plot", "📋 Meta Table", "📚 Study Results", "🔥 Meta-Heatmap"])
        
        with tab1:
            fig = _make_meta_volcano(meta_df, fc_thr, p_thr)
            fig.add_vline(x=fc_thr, line_dash='dash', line_color='black'); fig.add_vline(x=-fc_thr, line_dash='dash', line_color='black'); fig.add_hline(y=-np.log10(p_thr), line_dash='dash', line_color='black')
            st.session_state['fig_meta_v'] = fig
            st.plotly_chart(fig, width="stretch")
            
        with tab2:
            options = meta_degs['Symbol'].head(100).tolist()
            if options:
                gene = st.selectbox("Gene for forest plot", options)
                figf = _make_forest_plot(gene, st.session_state['study_results'], meta_df)
                st.session_state['fig_forest'] = figf
                if figf is not None: st.plotly_chart(figf, width="stretch")
                row = meta_df[meta_df['Symbol'] == gene].iloc[0]
                st.write(f"**I²:** {row['I2_percent']:.1f}%  |  **Tau²:** {row['Tau2']:.4g}  |  **Studies:** {int(row['N_Studies'])}")
                
        with tab3:
            all_cols = list(meta_df.columns)
            default_cols = ['Symbol','Meta_Log2FC','Meta_SE','CI95_Low','CI95_High','PValue','FDR','I2_percent','Tau2','N_Studies','Studies','Status']
            valid_defaults = [c for c in default_cols if c in all_cols]

            selected_cols = st.multiselect("Select columns to display in Meta Table:", all_cols, default=valid_defaults, key="meta_cols_sel")

            st.dataframe(meta_df[selected_cols], width="stretch", height=500)
            st.download_button('📥 Download Meta-analysis CSV', data=meta_df[selected_cols].to_csv(index=False), file_name='DDEA_Meta_Analysis_results.csv', mime='text/csv')
            
        with tab4:
            for sr in st.session_state['study_results']:
                gse = sr['Study'].iloc[0]
                st.write(f"**{gse}**")
                
                sr_filtered = sr.copy()
                if target_genes:
                    sr_filtered = sr_filtered[sr_filtered['Symbol'].str.upper().isin(target_genes)]

                sr_cols = list(sr_filtered.columns)
                sr_defaults = ['Symbol','Log2FC','SE','PValue','FDR']
                sr_valid = [c for c in sr_defaults if c in sr_cols]

                sel_cols_sr = st.multiselect(f"Columns for {gse}:", sr_cols, default=sr_valid, key=f"cols_{gse}")
                st.dataframe(sr_filtered[sel_cols_sr].sort_values('FDR').head(50), width="stretch")

        with tab5:
            if not meta_degs.empty:
                top_meta_genes = meta_degs.head(max_plot_meta)['Symbol'].tolist()
                
                z_matrices = []
                for gse, df_norm_study in st.session_state['meta_norm_df'].items():
                    conf = st.session_state['meta_configured'][gse]
                    c_ref = conf['c_ref']
                    c_test = conf['c_test']
                    ordered_study_cols = c_ref + c_test

                    sr_study = next(sr for sr in st.session_state['study_results'] if sr['Study'].iloc[0] == gse)
                    mapping_dict = dict(zip(sr_study['Probe_ID'], sr_study['Symbol']))
                    
                    df_norm_mapped = df_norm_study.copy()
                    df_norm_mapped.index = df_norm_mapped.index.map(mapping_dict)
                    
                    available_genes = [g for g in top_meta_genes if g in df_norm_mapped.index and pd.notna(g)]
                    if not available_genes:
                        continue
                        
                    df_norm_mapped = df_norm_mapped.loc[available_genes].groupby(level=0).mean()
                    
                    valid_cols = [c for c in ordered_study_cols if c in df_norm_mapped.columns]
                    valid_cols = list(dict.fromkeys(valid_cols)) 
                    
                    df_slice = df_norm_mapped.loc[:, valid_cols]
                    df_slice = df_slice.loc[:, ~df_slice.columns.duplicated()]
                    
                    mat = df_slice.values
                    z_mat = (mat - np.mean(mat, axis=1, keepdims=True)) / (np.std(mat, axis=1, keepdims=True) + 1e-9)
                    
                    df_z = pd.DataFrame(z_mat, index=df_slice.index, columns=[f"{gse}_{c}" for c in df_slice.columns])
                    z_matrices.append(df_z)
                
                if z_matrices:
                    combined_df_z = pd.concat(z_matrices, axis=1, join='inner')
                    
                    if not combined_df_z.empty:
                        h_z_meta = combined_df_z.values
                        y_labels_meta = combined_df_z.index.tolist()
                        
                        x_labels_meta = []
                        for col in combined_df_z.columns:
                            gse, c = col.split('_', 1)
                            conf = st.session_state['meta_configured'][gse]
                            grp = st.session_state['rn'] if c in conf['c_ref'] else st.session_state['tn']
                            x_labels_meta.append(f"{gse} | {grp} ({c})")
                        
                        col1, col2 = st.columns(2)
                        cluster_meta_r = col1.checkbox("Cluster Genes (Rows)", value=True, key="cluster_meta_r")
                        cluster_meta_c = col2.checkbox("Cluster Samples (Columns)", value=False, key="cluster_meta_c")
                        st.caption("Tip: Keep 'Cluster Samples' disabled to naturally separate samples by Dataset and Phenotype.")
                        
                        if cluster_meta_r and h_z_meta.shape[0] > 1:
                            try:
                                row_link = hierarchy.linkage(distance.pdist(h_z_meta, metric='euclidean'), method='ward')
                                row_order = hierarchy.leaves_list(row_link)
                                h_z_meta = h_z_meta[row_order, :]
                                y_labels_meta = [y_labels_meta[i] for i in row_order]
                            except Exception as e:
                                st.warning(f"Row clustering failed: {e}")
                                
                        if cluster_meta_c and h_z_meta.shape[1] > 1:
                            try:
                                col_link = hierarchy.linkage(distance.pdist(h_z_meta.T, metric='euclidean'), method='ward')
                                col_order = hierarchy.leaves_list(col_link)
                                h_z_meta = h_z_meta[:, col_order]
                                x_labels_meta = [x_labels_meta[i] for i in col_order]
                            except Exception as e:
                                st.warning(f"Column clustering failed: {e}")
                                
                        fig_hm_meta = go.Figure(data=go.Heatmap(
                            z=h_z_meta, x=x_labels_meta, y=y_labels_meta, colorscale='RdBu_r', zmid=0
                        )).update_layout(
                            height=max(400, len(y_labels_meta) * 20),
                            xaxis_tickangle=-45
                        )
                        
                        st.plotly_chart(fig_hm_meta, width="stretch")
                    else:
                        st.warning("No common genes found across datasets to plot heatmap.")
                else:
                    st.warning("Could not build a combined expression matrix.")

def _render_multi_tissue_mode(mode, apply_log, p_thr, fc_thr, gene_area=""):
    st.header("🧬 Multi-Tissue Transcriptomic Comparison")
    st.caption(
        "Use two or more datasets representing different tissues. Each dataset is "
        "analyzed independently and results are aligned by gene symbol; effects are "
        "not pooled across biologically different tissues."
    )

    gse_text = st.text_area("GSE IDs (one per line or separated by commas):", placeholder="GSE205748\nGSE117769")
    gse_ids = [x.strip().upper() for x in re.split(r'[\n,;]+', gse_text) if x.strip()]
    gse_ids = list(dict.fromkeys(gse_ids))

    if st.button("📥 Fetch Datasets for Comparison", disabled=len(gse_ids) < 2, width="stretch"):
        _reset_analysis_state()
        studies = {}
        progress = st.progress(0, text="Fetching datasets...")
        for i, gse in enumerate(gse_ids):
            try:
                df, meta_df, gsms, gsm_order, source, err, detected_type = get_geo_full_data(gse, mode)
                if meta_df is None or df is None:
                    st.warning(f"⚠️ {gse}: {err or 'dataset unavailable'}")
                    continue
                id_type = detect_index_type(df.index.astype(str).tolist())
                analytical_type = detected_type if detected_type in {'RNASeq', 'Microarray'} else mode
                mapping, msg = (get_gene_mapping_rnaseq(tuple(df.index.astype(str).tolist()), id_type) if analytical_type == 'RNASeq' else get_gene_mapping_microarray(gse))
                studies[gse] = {
                    'gse_id': gse, 'df': df, 'meta_df': meta_df, 'gsm_order': gsm_order,
                    'source': source, 'detected_type': analytical_type, 'mapping': mapping,
                    'label_cols': _study_label_cols(meta_df)
                }
            except Exception as e:
                st.warning(f"⚠️ {gse}: {e}")
            progress.progress((i+1)/len(gse_ids), text=f"Fetched {i+1}/{len(gse_ids)}")
        progress.empty()
        st.session_state['mt_studies'] = studies
        st.session_state['mt_gse_ids'] = list(studies.keys())
        st.rerun()

    studies = st.session_state.get('mt_studies', {})
    if not studies:
        st.info("Enter at least two GSE IDs and click 'Fetch Datasets for Comparison'.")
        return

    st.success(f"Loaded {len(studies)} datasets. Configure Control and Disease samples for each.")
    configured = {}
    for gse, study in studies.items():
        with st.expander(f"📚 {gse} — {study['detected_type']} — {study['df'].shape[0]} genes", expanded=True):
            meta = study['meta_df'].copy()
            meta_cols = list(meta.columns)
            valid_defaults = [c for c in study['label_cols'] if c in meta_cols]
            selected_labels = st.multiselect(f"Select metadata columns for {gse}:", meta_cols, default=valid_defaults, key=f"mt_label_{gse}")
            if not selected_labels: selected_labels = ['Accession'] if 'Accession' in meta_cols else [meta_cols[0]]
            
            meta = _unique_sample_labels(meta, selected_labels)
            label_to_gsm = dict(zip(meta['display_label'], meta['Accession']))
            sample_labels = meta['display_label'].tolist()
            
            c1, c2 = st.columns(2)
            ctrl = c1.multiselect("Control samples", sample_labels, key=f"mt_ctrl_{gse}")
            dis = c2.multiselect("Disease/Test samples", [x for x in sample_labels if x not in ctrl], key=f"mt_dis_{gse}")
            
            configured[gse] = {
                **study,
                'c_ref': [label_to_gsm[x] for x in ctrl],
                'c_test': [label_to_gsm[x] for x in dis],
            }

    if st.button("🚀 Run Multi-Tissue Analysis", width="stretch"):
        merged_results = pd.DataFrame()
        errors = []
        valid_gses = []
        study_results = []
        study_norms = {}
        progress = st.progress(0, text="Running independent analyses...")
        
        for i, (gse, study) in enumerate(configured.items()):
            if len(study['c_ref']) < 2 or len(study['c_test']) < 2:
                errors.append(f"{gse}: select at least 2 Control and 2 Disease samples.")
                continue
            try:
                c_ref = resolve_sample_columns(study['df'].columns, study['c_ref'])
                c_test = resolve_sample_columns(study['df'].columns, study['c_test'])
                design_metadata, design_formula, contrast, _ = _build_study_design(
                    study, c_ref, c_test
                )
                sr, norm = _run_single_study(
                    study['df'], c_ref, c_test,
                    study['detected_type'], apply_log,
                    mapping=study['mapping'], study_id=gse, alpha=p_thr,
                    design_metadata=design_metadata, design_formula=design_formula,
                    contrast=contrast,
                    continuous_engine=study.get('design_settings', {}).get('engine', 'limma'),
                )
                
                sr_sub = collapse_symbol_results(sr)
                study_results.append(sr)
                study_norms[gse] = norm
                
                sr_sub['Status'] = 'Not Significant'
                sr_sub.loc[(sr_sub['FDR'] < p_thr) & (sr_sub['Log2FC'] >= fc_thr), 'Status'] = 'Up-regulated'
                sr_sub.loc[(sr_sub['FDR'] < p_thr) & (sr_sub['Log2FC'] <= -fc_thr), 'Status'] = 'Down-regulated'
                
                sr_sub = sr_sub.rename(columns={
                    'Log2FC': f'Log2FC_{gse}', 'PValue': f'PValue_{gse}', 
                    'FDR': f'FDR_{gse}', 'Status': f'Status_{gse}'
                })
                
                if merged_results.empty:
                    merged_results = sr_sub
                else:
                    merged_results = pd.merge(merged_results, sr_sub, on='Symbol', how='outer')
                
                valid_gses.append(gse)
                
            except Exception as e:
                errors.append(f"{gse}: {e}")
            progress.progress((i+1)/len(configured), text=f"Analyzed {i+1}/{len(configured)}")
        progress.empty()

        if errors:
            for e in errors: st.warning(f"⚠️ {e}")
        if len(valid_gses) < 2:
            st.error("At least 2 datasets must succeed to perform a cross-tissue comparison.")
            return

        merged_results = annotate_cross_dataset_patterns(merged_results, valid_gses)
        st.session_state.update({
            'mt_results': merged_results,
            'mt_valid_gses': valid_gses,
            'mt_study_results': study_results,
            'mt_norm_df': study_norms,
            'mt_configured': configured,
            'analysis_mode': 'multi_tissue',
        })
        st.rerun()

    if st.session_state.get('mt_results') is not None:
        df_mt = st.session_state['mt_results'].copy()
        gses = st.session_state['mt_valid_gses']
        
        target_genes = []
        if gene_area.strip():
            target_genes = [g.strip().upper() for g in gene_area.split('\n') if g.strip()]
            df_mt = df_mt[df_mt['Symbol'].str.upper().isin(target_genes)].copy()
            
        st.divider(); st.header("📊 Multi-Tissue Alignment Results")
        
        df_degs_only = df_mt[df_mt['N_Significant_Datasets'] > 0].copy()
        strict_conserved = df_degs_only['Cross_Tissue_Pattern'].isin(
            ['Conserved Up', 'Conserved Down']
        )
        dataset_specific = df_degs_only['N_Significant_Datasets'] == 1
        discordant = df_degs_only['Cross_Tissue_Pattern'] == 'Discordant'

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Datasets / tissues", len(gses))
        m2.metric("DEGs in ≥1 tissue", len(df_degs_only))
        m3.metric("Strictly conserved", int(strict_conserved.sum()))
        m4.metric("Tissue-specific", int(dataset_specific.sum()))
        if discordant.any():
            st.warning(
                f"{int(discordant.sum())} gene(s) change direction across tissues; "
                "these are labelled Discordant and are not treated as conserved."
            )

        t1, t2, t3, t4 = st.tabs([
            "📋 Cross-Tissue Table", "🌋 Pairwise Scatter",
            "🔥 Effect Heatmap", "📚 Dataset Results",
        ])

        with t1:
            pattern_options = sorted(df_degs_only['Cross_Tissue_Pattern'].dropna().unique())
            selected_patterns = st.multiselect(
                "Patterns to display", pattern_options, default=pattern_options,
                key="mt_pattern_filter",
            )
            displayed = df_degs_only[
                df_degs_only['Cross_Tissue_Pattern'].isin(selected_patterns)
            ] if selected_patterns else df_degs_only.iloc[0:0]
            st.dataframe(displayed, width="stretch", height=500)
            st.download_button(
                '📥 Download complete Multi-Tissue CSV',
                data=df_mt.to_csv(index=False),
                file_name='DDEA_MultiTissue_results.csv', mime='text/csv',
            )
            manifest = {
                'analysis': 'multi-tissue independent differential expression',
                'datasets': gses,
                'design_per_dataset': '~condition',
                'fdr_threshold': p_thr,
                'absolute_effect_threshold': fc_thr,
                'strict_conservation': 'significant with the same direction in every dataset',
                'combination_rule': 'outer join by Symbol; no pooled cross-tissue effect',
            }
            st.download_button(
                '📥 Download Multi-Tissue manifest', json.dumps(manifest, indent=2),
                'DDEA_MultiTissue_manifest.json', mime='application/json',
            )

        with t2:
            s1, s2 = st.columns(2)
            gse_x = s1.selectbox("X-axis dataset", gses, index=0, key="mt_scatter_x")
            y_options = [g for g in gses if g != gse_x]
            gse_y = s2.selectbox("Y-axis dataset", y_options, index=0, key="mt_scatter_y")
            plot_df = df_degs_only.dropna(
                subset=[f'Log2FC_{gse_x}', f'Log2FC_{gse_y}']
            ).copy()
            sx, sy = plot_df[f'Status_{gse_x}'], plot_df[f'Status_{gse_y}']
            plot_df['Pairwise_Behavior'] = 'Other / significant in another tissue'
            plot_df.loc[(sx == 'Up-regulated') & (sy == 'Up-regulated'), 'Pairwise_Behavior'] = 'Conserved Up'
            plot_df.loc[(sx == 'Down-regulated') & (sy == 'Down-regulated'), 'Pairwise_Behavior'] = 'Conserved Down'
            plot_df.loc[(sx != 'Not Significant') & (sy == 'Not Significant'), 'Pairwise_Behavior'] = f'Specific {gse_x}'
            plot_df.loc[(sx == 'Not Significant') & (sy != 'Not Significant'), 'Pairwise_Behavior'] = f'Specific {gse_y}'
            plot_df.loc[
                ((sx == 'Up-regulated') & (sy == 'Down-regulated'))
                | ((sx == 'Down-regulated') & (sy == 'Up-regulated')),
                'Pairwise_Behavior',
            ] = 'Discordant'
            fig = px.scatter(
                plot_df, x=f'Log2FC_{gse_x}', y=f'Log2FC_{gse_y}',
                color='Pairwise_Behavior', hover_name='Symbol',
                hover_data=['Cross_Tissue_Pattern', 'N_Significant_Datasets'],
                color_discrete_map={
                    'Discordant': '#8E44AD', 'Conserved Up': '#D62728',
                    'Conserved Down': '#1F77B4', f'Specific {gse_x}': '#F39C12',
                    f'Specific {gse_y}': '#16A085',
                    'Other / significant in another tissue': '#BDBDBD',
                },
                title=f"Effect-size alignment: {gse_x} vs {gse_y}",
            )
            fig.add_hline(y=0, line_dash='dash', line_color='black')
            fig.add_vline(x=0, line_dash='dash', line_color='black')
            fig.update_layout(template='simple_white', height=600)
            st.plotly_chart(fig, width="stretch")

        with t3:
            effect_cols = [f'Log2FC_{g}' for g in gses]
            fdr_cols = [f'FDR_{g}' for g in gses]
            heat = df_degs_only.copy()
            heat['_best_fdr'] = heat[fdr_cols].min(axis=1, skipna=True)
            heat = heat.sort_values(
                ['N_Significant_Datasets', '_best_fdr'], ascending=[False, True]
            ).head(60)
            if heat.empty:
                st.info("No significant genes are available for the effect heatmap.")
            else:
                matrix = heat[effect_cols].to_numpy(dtype=float)
                fig_heat = go.Figure(go.Heatmap(
                    z=matrix, x=gses, y=heat['Symbol'], colorscale='RdBu_r',
                    zmid=0, colorbar=dict(title='log2FC'),
                    hovertemplate='Gene=%{y}<br>Dataset=%{x}<br>log2FC=%{z:.3f}<extra></extra>',
                ))
                fig_heat.update_layout(
                    template='simple_white', height=max(500, len(heat) * 18),
                    title='Cross-tissue effect-size heatmap (up to 60 genes)',
                )
                st.plotly_chart(fig_heat, width="stretch")

        with t4:
            for sr in st.session_state.get('mt_study_results', []):
                gse = str(sr['Study'].iloc[0])
                sr_display = sr.copy()
                if target_genes:
                    sr_display = sr_display[
                        sr_display['Symbol'].astype(str).str.upper().isin(target_genes)
                    ]
                preferred = ['Symbol', 'Log2FC', 'SE', 'PValue', 'FDR', 'AnalysisMethod']
                columns = [c for c in preferred if c in sr_display.columns]
                st.markdown(f"#### {gse}")
                st.dataframe(
                    sr_display.sort_values('FDR', na_position='last')[columns].head(100),
                    width="stretch", hide_index=True,
                )
                st.download_button(
                    f'📥 Download complete {gse} results',
                    sr_display.to_csv(index=False), f'{gse}_DDEA_results.csv',
                    mime='text/csv', key=f'mt_download_{gse}',
                )

def run_app():
    st.title("Diagonal Differential Expression Alley 🧬")

    analysis_mode = st.sidebar.radio("Analysis mode", ["Single Dataset", "Multi-Dataset Meta-analysis", "Multi-Tissue Comparison"], index=0)
    mode = st.sidebar.radio("Experiment Type", ["Microarray", "RNASeq"], index=1)
    apply_log = st.sidebar.checkbox("Apply Log2 when needed", value=True)
    p_thr = st.sidebar.slider("FDR threshold", 0.001, 0.10, 0.05, format="%.3f")
    fc_thr = st.sidebar.slider("Min Abs Log2FC", 0.0, 10.0, 1.0, step=0.1)
    max_plot = st.sidebar.number_input("Top genes for heatmap", value=50, min_value=5)

    if analysis_mode == "Multi-Dataset Meta-analysis":
        min_studies = st.sidebar.number_input("Minimum studies per gene", min_value=2, max_value=20, value=2, step=1)
        meta_method = st.sidebar.selectbox("Random-effects estimator", ["Paule-Mandel", "DerSimonian-Laird"], index=0)
        meta_gene_area = st.sidebar.text_area("Genes Highlight/Filter (1 per line, empty for all):", key="meta_genes_input")
        _render_meta_mode(mode, apply_log, p_thr, fc_thr, int(min_studies), meta_method, meta_gene_area, max_plot_meta=max_plot)
        return
    elif analysis_mode == "Multi-Tissue Comparison":
        mt_gene_area = st.sidebar.text_area("Genes Highlight/Filter (1 per line):", key="mt_genes_input")
        _render_multi_tissue_mode(mode, apply_log, p_thr, fc_thr, mt_gene_area)
        return

    # SINGLE DATASET MODE
    with st.sidebar:
        st.header("1. GEO Input")
        gse_input = st.text_input("GSE ID:", placeholder="Ex: GSE117769")
        fetch_btn = st.button("🚀 Fetch Data", width="stretch")

    if fetch_btn and gse_input:
        with st.spinner("🚀 Retrieving data..."):
            df, meta_df, gsms, gsm_order, source, err, detected_type = get_geo_full_data(gse_input.strip(), mode)
            if meta_df is not None:
                _reset_analysis_state(); st.session_state['groups'] = {}
                st.session_state.update({'mode': mode, 'gse_id': gse_input.strip(), 'meta_df': meta_df, 'gsm_order': gsm_order, 'df': df, 'detected_type': detected_type})
                if df is not None:
                    id_type = detect_index_type(df.index.tolist())
                    analytical_type = detected_type if detected_type in {'RNASeq', 'Microarray'} else mode
                    mapping, msg = get_gene_mapping_rnaseq(tuple(df.index.astype(str).tolist()), id_type) if analytical_type == 'RNASeq' else get_gene_mapping_microarray(gse_input.strip())
                    st.session_state['detected_type'] = analytical_type
                    st.session_state.update({'mapping': mapping, 'id_type': id_type, 'matrix_source': source, 'label_cols': ['Accession'] + _study_label_cols(meta_df)})
                    st.success(f"✅ Matrix via **{source}** — {df.shape[0]} genes × {df.shape[1]} samples")
                st.rerun()
            else:
                st.error(f"❌ Failed to fetch {gse_input.strip()}"); st.warning(err)

    if 'meta_df' not in st.session_state:
        st.info("Enter a GSE ID and click Fetch Data.")
        return

    meta = st.session_state['meta_df'].copy()
    with st.expander("📊 Metadata Explorer"): st.dataframe(meta, width="stretch")
    st.sidebar.header("2. Parameters")
    meta_cols = [c for c in meta.columns if c != 'Accession']
    default_cols = _study_label_cols(meta)
    selected = st.sidebar.multiselect("Select metadata columns", meta_cols, default=default_cols)
    st.session_state['label_cols'] = ['Accession'] + selected
    gene_area = st.sidebar.text_area("Genes Highlight (1 per line):")

    _render_single_dataset(mode, apply_log, gene_area, p_thr, fc_thr, max_plot)

    if st.session_state.get('analysis_done') and not st.session_state.get('meta_mode'):
        st.sidebar.header("3. Export Report")
        if st.sidebar.button("🌐 Generate HTML Report"):
            ss = st.session_state; gse_id = ss.get('gse_id', 'N/A'); dd = ss.get('df_diff'); res = ss.get('res')
            has = lambda s: dd is not None and bool((dd['Status'] == s).any())
            methods = sorted(res['AnalysisMethod'].dropna().astype(str).unique().tolist()) if isinstance(res, pd.DataFrame) and 'AnalysisMethod' in res else []
            ctx = {'gse_id': gse_id, 'mode': ss.get('mode', mode), 'tn': ss.get('tn', ''), 'rn': ss.get('rn', ''),
                   'n_ref': len(ss.get('c_ref', [])), 'n_test': len(ss.get('c_test', [])), 'p_thr': p_thr, 'fc_thr': fc_thr, 'max_plot': max_plot, 'df_diff': dd,
                   'methods': methods, 'design_formula': ss.get('design_formula', '~condition'), 'version': LUMOS_VERSION,
                   'figures': {'fig_v': ss.get('fig_v'), 'fig_ma': ss.get('fig_ma'), 'fig_pca': ss.get('fig_pca'), 'fig_h': ss.get('fig_h'),
                               'fig_up': ss.get('fig_up') if has('Up-regulated') else None, 'fig_down': ss.get('fig_down') if has('Down-regulated') else None}}
            with st.spinner("Generating HTML report..."):
                html_str = generate_html_report(ctx)
                safe_id = re.sub(r'\W', '_', gse_id)
                st.sidebar.download_button("📥 Download HTML Report", data=html_str.encode('utf-8'), file_name=f"DDEA_Report_{safe_id}_{datetime.now().strftime('%Y%m%d')}.html", mime="text/html", on_click="ignore")

if __name__ == '__main__':
    run_app()
