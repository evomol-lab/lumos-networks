import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
import requests
import gzip
import io
import re
import gc
import tarfile
from Bio import Entrez
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.decomposition import PCA

HEADERS = {'User-Agent': 'DDEA/3.0 (Streamlit App)'}
Entrez.email = "ddea.tool@example.com"

# ============================================================
# NORMALIZAÇÃO E UTILITÁRIOS
# ============================================================

def quantile_normalize(df_values):
    if df_values.size == 0:
        return df_values
    mat = df_values.astype(np.float32)
    sorted_mat = np.sort(mat, axis=0)
    rank_mean = sorted_mat.mean(axis=1).astype(np.float32)
    del sorted_mat; gc.collect()
    indices = np.argsort(mat, axis=0)
    norm_mat = np.empty_like(mat, dtype=np.float32)
    for i in range(mat.shape[1]):
        norm_mat[indices[:, i], i] = rank_mean
    del mat, indices; gc.collect()
    return norm_mat

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

# ============================================================
# MAPEAMENTO — MICROARRAY E RNA-SEQ
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
    except Exception as e:
        return None, None

def _strip_ensembl_version(ensembl_id: str) -> str:
    return ensembl_id.split('.')[0] if '.' in ensembl_id else ensembl_id

@st.cache_data(show_spinner=False)
def get_gene_mapping_rnaseq(index_ids: tuple, id_type: str):
    if id_type == 'symbol': return pd.DataFrame({'Probe_ID': list(index_ids), 'Symbol': list(index_ids)}), "Já contém Gene Symbols."
    if id_type not in ('entrez', 'ensembl'): return None, f"Tipo de ID '{id_type}' não suportado."

    if id_type == 'entrez':
        ids_clean = [str(i).strip() for i in index_ids if str(i).strip().isdigit()]
        scope, original_to_clean = "entrezgene", {i: i for i in ids_clean}
    else:
        ids_raw = [str(i).strip() for i in index_ids if str(i).strip().startswith('ENS')]
        scope = "ensembl.gene"
        original_to_clean = {i: _strip_ensembl_version(i) for i in ids_raw}
        clean_to_original = {clean: orig for orig, clean in original_to_clean.items()}
        ids_clean = list(clean_to_original.keys())

    if not ids_clean: return None, "Nenhum ID válido."

    results = []
    for i in range(0, len(ids_clean), 1000):
        chunk = ids_clean[i:i + 1000]
        try:
            resp = requests.post("https://mygene.info/v3/query", data={"q": ",".join(chunk), "scopes": scope, "fields": "symbol", "species": "human"}, timeout=30)
            resp.raise_for_status()
            for item in resp.json():
                if item.get('notfound'): continue
                query_id = str(item.get('query', ''))
                symbol = item.get('symbol', None)
                orig_id = clean_to_original.get(query_id, query_id) if id_type == 'ensembl' else query_id
                results.append({"Probe_ID": orig_id, "Symbol": symbol})
        except: return None, "Erro na API MyGene.info"

    if not results: return None, "Sem resultados MyGene."

    mapping_df = pd.DataFrame(results).drop_duplicates('Probe_ID')
    mapping_df["Symbol"] = mapping_df["Symbol"].fillna(mapping_df["Probe_ID"])
    return mapping_df, f"{(mapping_df['Symbol'] != mapping_df['Probe_ID']).sum()}/{len(mapping_df)} IDs convertidos."

# ============================================================
# PARSE E EXTRAÇÃO DE ARQUIVOS GEO
# ============================================================

def _parse_matrix_bytes(raw_bytes, filename=""):
    try:
        content = gzip.decompress(raw_bytes) if filename.endswith('.gz') else raw_bytes
        for sep in ['\t', ',']:
            try:
                df = pd.read_csv(io.BytesIO(content), sep=sep, index_col=0, low_memory=False)
                df.index = df.index.astype(str).str.strip().str.replace('"', '')
                df.columns = [str(c).strip().replace('"', '') for c in df.columns]
                num_cols = df.select_dtypes(include=[np.number]).shape[1]
                if num_cols >= 2 and df.shape[0] >= 10:
                    return df.select_dtypes(include=[np.number])
            except: continue
        return None
    except: return None

def _parse_series_type(series_type_raw: str) -> str:
    s = series_type_raw.lower()
    if 'high throughput sequencing' in s or 'sequencing' in s: return 'RNASeq'
    if 'array' in s: return 'Microarray'
    return 'unknown'

def _try_series_matrix(gse_id):
    num = gse_id.replace("GSE", "")
    prefix = f"GSE{num[:-3]}nnn" if len(num) > 3 else "GSEnnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    try:
        r = requests.get(url, timeout=60, headers=HEADERS)
        r.raise_for_status()
        with gzip.open(io.BytesIO(r.content), 'rt') as f:
            titles, gsms, char_lines, gsm_order, series_type_raw = [], [], [], [], ""
            df = pd.DataFrame()
            for line in f:
                if line.startswith('!Series_type'): series_type_raw = line.split('\t')[1].strip().replace('"', '') if '\t' in line else ""
                if line.startswith('!Sample_title'): titles = [t.strip().replace('"', '') for t in line.split('\t')[1:]]
                if line.startswith('!Sample_geo_accession'):
                    gsms = [t.strip().replace('"', '') for t in line.split('\t')[1:]]
                    gsm_order = gsms[:]
                if line.startswith('!Sample_characteristics_ch1'): char_lines.append(line.split('\t')[1:])
                if line.startswith('ID_REF') or line.startswith('"ID_REF"'):
                    df = pd.read_csv(f, sep='\t', header=None, low_memory=True)
                    break

            detected_type = _parse_series_type(series_type_raw)
            meta_dict = {"Accession": gsms, "Title": titles}
            for row in char_lines:
                if row:
                    key = row[0].split(':')[0].strip() if ':' in row[0] else f"Info_{len(meta_dict)}"
                    meta_dict[key] = [v.split(': ')[1].strip() if ': ' in v else v.strip() for v in row][:len(gsms)]
            meta_df = pd.DataFrame(meta_dict)

            if df.empty or len(df.columns) < 2: return None, meta_df, gsms, gsm_order, detected_type
            df = df.set_index(0)
            df.index = df.index.astype(str).str.strip().str.replace('"', '')
            df.rename(columns={col: gsm_order[i] for i, col in enumerate(df.columns) if i < len(gsm_order)}, inplace=True)
            if df.select_dtypes(include=[np.number]).shape[1] < 2: return None, meta_df, gsms, gsm_order, detected_type
            return df.select_dtypes(include=[np.number]), meta_df, gsms, gsm_order, detected_type
    except: return None, None, None, [], 'unknown'

def _list_supplementary_urls(gse_id):
    num = gse_id.replace("GSE", "")
    prefix = f"GSE{num[:-3]}nnn" if len(num) > 3 else "GSEnnn"
    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/suppl/"
    try:
        r = requests.get(base_url, timeout=20, headers=HEADERS)
        if r.status_code != 200: return []
        files = re.findall(r'href="([^"]+\.(txt|tsv|csv|tar)(\.gz)?)"', r.text, re.IGNORECASE)
        return [base_url + f[0] for f in files]
    except: return []

def _score_suppl_file(url):
    name = url.lower()
    score = 0
    if any(k in name for k in ['count', 'raw', 'read', 'htseq', 'matrix', 'tar']): score += 3
    if any(k in name for k in ['norm', 'fpkm', 'rpkm', 'tpm']): score += 1
    if 'log' in name: score -= 1
    return score

def _try_extract_symbol_mapping_from_suppl(gse_id, ensembl_ids, log_cb=None):
    urls = sorted(_list_supplementary_urls(gse_id), key=_score_suppl_file, reverse=True)
    if not urls: return None, "Nenhum arquivo suplementar encontrado."
    return None, "Extração suplementar complexa pulada nesta versão para evitar instabilidade."

def _try_supplementary(gse_id, log_cb=None):
    urls = sorted(_list_supplementary_urls(gse_id), key=_score_suppl_file, reverse=True)
    for url in urls:
        try:
            if log_cb: log_cb(f"Suplementar detectado: `{url.split('/')[-1]}`")
            f_res = requests.get(url, timeout=60, headers=HEADERS)
            f_res.raise_for_status()
            
            if '.tar' in url.lower():
                if log_cb: log_cb("Extraindo arquivo TAR na memória...")
                with tarfile.open(fileobj=io.BytesIO(f_res.content), mode='r:*') as tar:
                    for member in tar.getmembers():
                        if member.isfile() and any(ext in member.name.lower() for ext in ['.txt', '.csv', '.tsv']):
                            f = tar.extractfile(member)
                            df = _parse_matrix_bytes(f.read(), member.name)
                            if df is not None and df.shape[1] >= 2: return df
            else:
                df = _parse_matrix_bytes(f_res.content, url)
                if df is not None and df.shape[1] >= 2 and df.shape[0] >= 10: return df
        except: continue
    return None

def _try_ncbi_generated(gse_id, log_cb=None):
    candidates = [
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_raw_counts_GEO.txt.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file&file={gse_id}_norm_counts_GEO.txt.gz",
    ]
    for url in candidates:
        try:
            if log_cb: log_cb(f"NCBI-generated: `{url.split('file=')[-1]}`")
            r = requests.get(url, timeout=60, headers=HEADERS, allow_redirects=True)
            if r.status_code == 200 and len(r.content) > 1000:
                df = _parse_matrix_bytes(r.content)
                if df is not None and df.shape[1] >= 2: return df
        except: continue
    return None

def _sync_suppl_columns_with_gsms(df_suppl, gsm_order):
    cols, gsm_set = list(df_suppl.columns), set(gsm_order)
    direct = [c for c in cols if c in gsm_set]
    if len(direct) >= 2: return df_suppl[direct], "interseção direta"
    if len(cols) == len(gsm_order): return df_suppl.rename(columns=dict(zip(cols, gsm_order))), "renomeação posicional"
    if len(cols) - 1 == len(gsm_order): return df_suppl.iloc[:, 1:].rename(columns=dict(zip(df_suppl.columns[1:], gsm_order))), "renomeação ignorando 1ª col"
    return df_suppl, "não sincronizado"

def get_geo_full_data(gse_id, mode, log_cb=None):
    if log_cb: log_cb("Buscando metadados e Series Matrix...")
    df_matrix, meta_df, gsms, gsm_order, detected_type = _try_series_matrix(gse_id)
    if meta_df is None: return None, None, None, [], None, "Falha ao acessar o GEO. Verifique o ID.", 'unknown'

    if mode == "Microarray" and df_matrix is not None:
        return df_matrix, meta_df, gsms, gsm_order, "Series Matrix", None, detected_type
    
    if df_matrix is not None and df_matrix.shape[1] >= 2 and mode != "RNASeq":
        return df_matrix, meta_df, gsms, gsm_order, "Series Matrix", None, detected_type

    if log_cb: log_cb("Matriz insuficiente/vazia. Buscando arquivos suplementares...")
    df_suppl = _try_supplementary(gse_id, log_cb=log_cb)
    if df_suppl is not None:
        df_synced, sync_method = _sync_suppl_columns_with_gsms(df_suppl, gsm_order)
        return df_synced, meta_df, gsms, gsm_order, f"Supplementary ({sync_method})", None, detected_type

    if log_cb: log_cb("Tentando NCBI-generated counts...")
    df_ncbi = _try_ncbi_generated(gse_id, log_cb=log_cb)
    if df_ncbi is not None:
        df_synced, sync_method = _sync_suppl_columns_with_gsms(df_ncbi, gsm_order)
        return df_synced, meta_df, gsms, gsm_order, f"NCBI-generated ({sync_method})", None, detected_type

    return None, meta_df, gsms, gsm_order, None, None, detected_type

# ============================================================
# APP PRINCIPAL E GERENCIAMENTO DE ESTADO
# ============================================================

def reset_analysis_state():
    for k in ['df', 'meta_df', 'res', 'analysis_done', 'mapping', 'mapping_msg', 'raw_gpl', 'mode', 'norm_df', 'rn', 'tn', 'gse_id', 'gsms', 'gsm_order', 'matrix_source', 'id_type', 'detected_type']:
        st.session_state.pop(k, None)
    st.session_state['groups'] = {}
    st.session_state['group_field_key'] = st.session_state.get('group_field_key', 0) + 1

def all_assigned_samples(groups: dict, exclude: str = None) -> set:
    assigned = set()
    for name, samples in groups.items():
        if name != exclude: assigned.update(samples)
    return assigned

def run_app():
    st.set_page_config(layout="wide", page_title="DDEA Analytics Master")
    st.title("Diagonal Differential Expression Alley 🧬")

    if 'groups' not in st.session_state: st.session_state['groups'] = {}
    if 'group_field_key' not in st.session_state: st.session_state['group_field_key'] = 0

    with st.sidebar:
        st.header("1. GEO Input")
        mode = st.radio("Experiment Type:", ["Microarray", "RNASeq"])
        gse_input = st.text_input("GSE ID:", placeholder="Ex: GSE117769")
        fetch_btn = st.button("🚀 Fetch Data", use_container_width=True)

        if 'meta_df' in st.session_state:
            st.divider()
            det = st.session_state.get('detected_type', 'unknown')
            if det != 'unknown':
                tipo_label = "🧬 RNA-Seq" if det == 'RNASeq' else "🔬 Microarray"
                st.caption(f"🏷️ Tecnologia detectada: **{tipo_label}**")
                if st.session_state.get('mode') and st.session_state['mode'] != det:
                    st.warning(f"⚠️ Tipo selecionado difere do detectado (**{det}**).")
            if st.session_state.get('matrix_source'): st.caption(f"📦 Fonte: **{st.session_state['matrix_source']}**")
            if st.session_state.get('id_type'): st.caption(f"🔑 Tipo de ID: **{st.session_state['id_type']}**")
            if st.session_state.get('mapping_msg'): st.caption(f"🔗 {st.session_state['mapping_msg']}")

            st.divider()
            all_cols = list(st.session_state['meta_df'].columns)
            extra_cols = [c for c in all_cols if c != 'Accession']
            extra_selected = st.multiselect("Informações no seletor:", extra_cols, default=["Title"] if "Title" in extra_cols else extra_cols[:1])
            label_cols = ['Accession'] + extra_selected

            st.divider()
            st.header("2. Parameters")
            gene_area = st.text_area("Genes Highlight (1 per line):")
            p_thr = st.slider("FDR threshold:", 0.001, 0.10, 0.05, format="%.3f")
            fc_thr = st.slider("Min Abs Log2FC:", 0.0, 10.0, 1.0, step=0.1)
            use_limma = st.checkbox("Usar modelo linear (Limma)", value=True) if mode == "Microarray" else False
            max_plot = st.number_input("Top genes p/ Heatmap:", value=50, min_value=5)
        else:
            label_cols = ['Accession', 'Title']
            gene_area = ""; p_thr = 0.05; fc_thr = 1.0; use_limma = False; max_plot = 50

    if fetch_btn and gse_input:
        if gse_input.strip() != st.session_state.get('gse_id', '') or mode != st.session_state.get('mode', ''): reset_analysis_state()
        
        log_ph, log_lines = st.empty(), []
        def log_cb(msg): log_lines.append(f"› {msg}"); log_ph.info("\n".join(log_lines))
        
        with st.spinner("🚀 Buscando dados..."):
            df, meta_df, gsms, gsm_order, source, err, detected_type = get_geo_full_data(gse_input.strip(), mode, log_cb=log_cb)
        
        log_ph.empty()
        if err: st.error(err)
        elif meta_df is None: st.error("Falha ao acessar o GEO.")
        else:
            st.session_state.update({'mode': mode, 'gse_id': gse_input.strip(), 'matrix_source': source, 'gsm_order': gsm_order, 'detected_type': detected_type})
            
            if mode == "Microarray":
                with st.spinner("🔬 Buscando mapeamento GPL..."):
                    map_res, raw_gpl = get_gene_mapping_microarray(gse_input.strip())
                st.session_state.update({'raw_gpl': raw_gpl, 'mapping': map_res, 'id_type': 'probe'})
            else:
                if df is not None and not df.empty:
                    id_type = detect_index_type(df.index.tolist())
                    st.session_state['id_type'] = id_type
                    mapping, msg = None, ""
                    if id_type == 'ensembl':
                        with st.spinner("🔍 Buscando mapeamento suplementar..."):
                            mapping, msg = _try_extract_symbol_mapping_from_suppl(gse_input.strip(), df.index.astype(str).tolist(), log_cb=log_cb)
                    if mapping is None:
                        with st.spinner(f"🔗 Mapeando IDs ({id_type}) via MyGene..."):
                            mapping, msg = get_gene_mapping_rnaseq(tuple(df.index.astype(str).tolist()), id_type)
                    st.session_state.update({'mapping': mapping, 'mapping_msg': msg})
                else:
                    st.session_state.update({'mapping': None, 'mapping_msg': "", 'id_type': 'unknown'})
                    
            st.session_state.update({'df': df, 'meta_df': meta_df, 'gsms': gsms, 'analysis_done': False})
            if source and df is not None: st.success(f"✅ Matriz via **{source}** — {df.shape[0]} genes × {df.shape[1]} amostras")
            st.rerun()

    if 'meta_df' not in st.session_state:
        st.info("Insira um GSE ID e clique em **Fetch Data**.")
        return

    meta = st.session_state['meta_df'].copy()
    with st.expander("📊 Metadata Explorer", expanded=False): st.dataframe(meta, use_container_width=True)

    def make_label(row): return " | ".join([str(row[c]) for c in label_cols if c in row.index]) if label_cols else str(row.get('Accession', ''))
    meta['display_label'] = meta.apply(make_label, axis=1)
    all_labels = meta['display_label'].tolist()
    label_to_gsm = dict(zip(meta['display_label'], meta['Accession']))

    st.subheader("🛠️ Group Management")
    groups = st.session_state['groups']
    c_i, c_b = st.columns([3, 1])
    with c_i: new_g = st.text_input("New Group Name:", key=f"gf_{st.session_state['group_field_key']}", placeholder="Control, Treatment...")
    with c_b: 
        st.write("")
        if st.button("➕ Add Group", use_container_width=True) and new_g.strip():
            if new_g.strip() not in groups: groups[new_g.strip()] = []
            st.session_state['group_field_key'] += 1; st.rerun()

    updated_groups = {}
    cols_g = st.columns(2)
    for i, (g_name, g_samples) in enumerate(list(groups.items())):
        avail = [lab for lab in all_labels if lab not in all_assigned_samples(groups, g_name)]
        with cols_g[i % 2]:
            st.markdown(f"**{g_name}**")
            updated_groups[g_name] = st.multiselect("Samples:", avail, default=[s for s in g_samples if s in avail], key=f"s_{g_name}", label_visibility='collapsed')
            if st.button("🗑️ Remove", key=f"del_{g_name}"): del st.session_state['groups'][g_name]; st.rerun()
    st.session_state['groups'] = updated_groups

    # ============================================================
    # UPLOAD MANUAL (Restaurado Integralmente)
    # ============================================================
    df_current = st.session_state.get('df')
    matrix_empty = df_current is None or (isinstance(df_current, pd.DataFrame) and df_current.empty)

    if matrix_empty:
        st.warning("📦 Nenhuma matriz de expressão encontrada automaticamente. Faça upload do arquivo de counts (.txt.gz / .tsv.gz / .tsv / .txt / .csv).")
        up_file = st.file_uploader("Upload counts file", type=["gz", "tsv", "txt", "csv"])
        if up_file:
            try:
                raw = up_file.read()
                df_up = _parse_matrix_bytes(raw)
                if df_up is None:
                    st.error("Não foi possível parsear o arquivo. Verifique se é TSV/CSV numérico.")
                else:
                    gsm_order = st.session_state.get('gsm_order', [])
                    if gsm_order:
                        df_up, sync_msg = _sync_suppl_columns_with_gsms(df_up, gsm_order)
                        st.info(f"Sincronização: {sync_msg}")

                    id_type = detect_index_type(df_up.index.tolist())
                    st.session_state['id_type'] = id_type
                    with st.spinner(f"🔗 Mapeando IDs ({id_type}) → Gene Symbol..."):
                        mapping, msg = get_gene_mapping_rnaseq(tuple(df_up.index.astype(str).tolist()), id_type)
                    st.session_state.update({'mapping': mapping, 'mapping_msg': msg, 'df': df_up, 'matrix_source': "Upload manual"})
                    st.success(f"✅ {df_up.shape[0]} genes × {df_up.shape[1]} amostras. {msg}")
                    st.rerun()
            except Exception as e: st.error(f"Erro ao ler arquivo: {e}")
        return

    # ============================================================
    # MOTOR DE ANÁLISE DIFERENCIAL
    # ============================================================
    df_matrix = st.session_state['df']
    matrix_cols = set(df_matrix.columns.astype(str))

    if len(updated_groups) >= 2:
        st.divider()
        c1, c2 = st.columns(2)
        group_names = list(updated_groups.keys())
        ref_g = c1.selectbox("Referência (Controle):", group_names)
        test_g = c2.selectbox("Teste:", [g for g in group_names if g != ref_g])

        if st.button("🔥 Run Analysis", use_container_width=True):
            with st.spinner("Analisando..."):
                c_ref = [label_to_gsm[lab] for lab in updated_groups[ref_g] if label_to_gsm.get(lab) in matrix_cols]
                c_test = [label_to_gsm[lab] for lab in updated_groups[test_g] if label_to_gsm.get(lab) in matrix_cols]

                if not c_ref or not c_test:
                    st.error("❌ Nenhuma amostra válida encontrada nos grupos.")
                    st.stop()

                df_subset = df_matrix[c_ref + c_test].copy()
                df_subset = df_subset.apply(pd.to_numeric, errors='coerce').fillna(0)
                data_vals = df_subset.values.astype(np.float32)

                if st.session_state.mode == "Microarray": 
                    data_norm_vals = quantile_normalize(data_vals)
                else: 
                    data_norm_vals = np.log2(np.clip(data_vals, 0, None) + 1.0)

                data_norm = pd.DataFrame(data_norm_vals, columns=c_ref + c_test, index=df_subset.index)
                m1, m2 = data_norm[c_ref].values, data_norm[c_test].values
                is_integer = np.all(np.equal(np.mod(data_vals, 1), 0))

                if st.session_state.mode == "RNASeq" and is_integer:
                    data_ds = df_subset[df_subset.sum(axis=1) > 0].astype(int)
                    if data_ds.empty:
                        st.error("❌ Matriz zerada após filtros."); st.stop()
                    meta_ds = pd.DataFrame({'cond': ['C']*len(c_ref) + ['T']*len(c_test)}, index=c_ref + c_test)
                    
                    try:
                        dds = DeseqDataSet(counts=data_ds.T, metadata=meta_ds, design_factors="cond", quiet=True)
                        dds.deseq2()
                        stat_res = DeseqStats(dds, contrast=["cond", "T", "C"], quiet=True)
                        stat_res.summary()
                        res_df = stat_res.results_df.copy()
                        res = pd.DataFrame({'Probe_ID': res_df.index.astype(str), 'Log2FC': res_df['log2FoldChange'].values, 'FDR': res_df['padj'].values}).dropna()
                    except Exception as e:
                        st.error(f"Erro PyDESeq2: {e}"); st.stop()
                else:
                    if st.session_state.mode == "RNASeq": st.warning("⚠️ Dados decimais detectados. Usando Teste Linear/OLS.")
                    if use_limma:
                        p_l, f_l = [], []
                        for row in range(len(data_norm)):
                            y, x = np.concatenate([m1[row], m2[row]]), sm.add_constant(np.concatenate([np.zeros(len(c_ref)), np.ones(len(c_test))]))
                            mod = sm.OLS(y, x).fit()
                            p_l.append(mod.pvalues[1]); f_l.append(mod.params[1])
                        pvals, lfc = np.array(p_l), np.array(f_l)
                    else:
                        lfc = np.nanmean(m2, axis=1) - np.nanmean(m1, axis=1)
                        pvals = stats.ttest_ind(m2, m1, axis=1, equal_var=False, nan_policy='omit').pvalue

                    mask = ~np.isnan(pvals)
                    fdr = np.full(pvals.shape, np.nan)
                    fdr[mask] = fdrcorrection(pvals[mask])[1]
                    res = pd.DataFrame({'Probe_ID': df_subset.index.astype(str), 'Log2FC': lfc, 'FDR': fdr}).dropna()

                mapping = st.session_state.get('mapping')
                if mapping is not None:
                    res = res.merge(mapping, on='Probe_ID', how='left')
                    res['Symbol'] = res['Symbol'].replace(['nan', '---', ' ', 'None'], np.nan).fillna(res['Probe_ID'])
                else: res['Symbol'] = res['Probe_ID']

                st.session_state.update({'res': res, 'norm_df': data_norm, 'analysis_done': True, 'rn': ref_g, 'tn': test_g, 'c_ref': c_ref, 'c_test': c_test})
                gc.collect()
                st.rerun()

    # ============================================================
    # RESULTADOS E VISUALIZAÇÕES
    # ============================================================
    if st.session_state.get('analysis_done'):
        res, df_norm = st.session_state['res'].copy(), st.session_state['norm_df']
        c_ref, c_test = st.session_state['c_ref'], st.session_state['c_test']

        if gene_area:
            s_list = [g.strip().upper() for g in gene_area.split('\n') if g.strip()]
            res = res[res['Symbol'].str.upper().isin(s_list)]

        res['Status'] = 'Not Significant'
        res.loc[(res['FDR'] < p_thr) & (res['Log2FC'] >= fc_thr), 'Status'] = 'Up-regulated'
        res.loc[(res['FDR'] < p_thr) & (res['Log2FC'] <= -fc_thr), 'Status'] = 'Down-regulated'
        res['-log10(FDR)'] = -np.log10(res['FDR'] + 1e-300)

        df_diff = res[res['Status'] != 'Not Significant'].sort_values('FDR')
        st.subheader(f"🚀 Analysis Results: {st.session_state['tn']} vs {st.session_state['rn']}")

        tab1, tab2, tab3 = st.tabs(["📊 Volcano & PCA", "🔥 Heatmap Top DEGs", "📋 Data Table"])

        with tab1:
            fig_v = px.scatter(res, x='Log2FC', y='-log10(FDR)', color='Status', 
                               color_discrete_map={'Not Significant': 'lightgrey', 'Up-regulated': 'red', 'Down-regulated': 'blue'},
                               hover_name='Symbol', hover_data={'FDR': ':.2e'})
            fig_v.add_vline(x=fc_thr, line_dash="dash", line_color="black", opacity=0.5)
            fig_v.add_vline(x=-fc_thr, line_dash="dash", line_color="black", opacity=0.5)
            fig_v.add_hline(y=-np.log10(p_thr), line_dash="dash", line_color="black", opacity=0.5)
            fig_v.update_layout(template="simple_white", height=600)
            st.plotly_chart(fig_v, use_container_width=True)

            st.divider()
            st.markdown("### Principal Component Analysis (PCA)")
            try:
                pca = PCA(n_components=2)
                components = pca.fit_transform(df_norm.fillna(0).T)
                pca_df = pd.DataFrame(components, columns=['PC1', 'PC2'], index=df_norm.columns)
                pca_df['Group'] = [st.session_state['rn'] if x in c_ref else st.session_state['tn'] for x in pca_df.index]
                
                fig_pca = px.scatter(pca_df, x='PC1', y='PC2', color='Group', text=pca_df.index,
                                     title=f"PCA (PC1 {pca.explained_variance_ratio_[0]*100:.1f}%, PC2 {pca.explained_variance_ratio_[1]*100:.1f}%)")
                fig_pca.update_traces(textposition='top center', marker=dict(size=12))
                fig_pca.update_layout(template="simple_white", height=500)
                st.plotly_chart(fig_pca, use_container_width=True)
            except Exception as e: st.warning(f"PCA indisponível: {e}")

        with tab2:
            if df_diff.empty: st.info("Sem DEGs para o Heatmap.")
            else:
                top_genes = df_diff.head(max_plot).copy()
                valid_ids = [idx for idx in top_genes['Probe_ID'].values if idx in df_norm.index]
                if valid_ids:
                    h_mat = df_norm.loc[valid_ids].values
                    h_z = (h_mat - np.mean(h_mat, axis=1, keepdims=True)) / (np.std(h_mat, axis=1, keepdims=True) + 1e-9)
                    fig_h = go.Figure(data=go.Heatmap(z=h_z, x=df_norm.columns, y=top_genes.iloc[:len(valid_ids)]['Symbol'].tolist(), colorscale='RdBu_r', zmid=0))
                    fig_h.update_layout(height=max(400, len(valid_ids)*20))
                    st.plotly_chart(fig_h, use_container_width=True)

        with tab3:
            if df_diff.empty: st.info("Nenhum resultado significativo.")
            else: st.dataframe(df_diff[['Symbol', 'Log2FC', 'FDR', 'Status']].rename(columns={'Symbol': 'Gene Symbol'}), use_container_width=True)
        csv = df_diff[['Symbol', 'Log2FC', 'FDR', 'Status']].to_csv(index=False)
        st.download_button(
                label="📥 Baixar Resultados (CSV) para Análise de Redes",
                data=csv,
                file_name=f"DEGs_{st.session_state['tn']}_vs_{st.session_state['rn']}.csv",
                mime="text/csv",
            )
    # ----------------------------------------------------------
    # DEBUG
    # ----------------------------------------------------------
    st.divider()
    with st.expander("🕵️ Debug", expanded=False):
        c_db1, c_db2 = st.columns(2)
        with c_db1:
            if st.session_state.get('mode') == "Microarray" and st.session_state.get('raw_gpl') is not None:
                st.subheader("Raw GPL")
                st.dataframe(st.session_state['raw_gpl'], use_container_width=True)
            elif st.session_state.get('mode') == "RNASeq":
                mapping = st.session_state.get('mapping')
                if mapping is not None:
                    st.subheader(f"Mapping ({st.session_state.get('id_type', '?')})")
                    st.dataframe(mapping.head(50), use_container_width=True)
                else: st.info("Sem mapeamento.")
        with c_db2:
            if st.session_state.get('res') is not None:
                st.subheader("Merge Result (Top 20)")
                st.dataframe(st.session_state['res'].head(20), use_container_width=True)
            df_m = st.session_state.get('df')
            if df_m is not None:
                st.subheader("Matriz Bruta (Top 15 cols)")
                st.write(list(df_m.columns[:15]))
                st.caption(f"Shape: {df_m.shape}")

if __name__ == '__main__':
    run_app()
