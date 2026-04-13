import streamlit as st
import pandas as pd
import polars as pl  # ATUALIZADO
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
import os
from datetime import datetime
from fpdf import FPDF
from Bio import Entrez
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
import itertools
from sklearn.decomposition import PCA

# Impede erros de importação se as bibliotecas não estiverem presentes
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    HAS_DESEQ2 = True
except ImportError:
    HAS_DESEQ2 = False

HEADERS = {'User-Agent': 'DDEA/4.0 (Streamlit App; Academic Version)'}
Entrez.email = "ddea.tool@example.com"

# 1. Encontrar o diretório base do projeto
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 2. Caminho para o logo
logo_path = os.path.join(BASE_DIR, "assets", "logos", "DDEA.png")

with st.sidebar:
    if os.path.exists(logo_path):
        st.image(logo_path, use_container_width=True)
    else:
        st.error(f"Erro: Logo não encontrado em {logo_path}")

# ============================================================
# PDF REPORT GENERATION CLASS
# ============================================================

class PDF(FPDF):
    def __init__(self, orientation='P', unit='mm', format='A4'):
        super().__init__(orientation, unit, format)
        self.set_auto_page_break(auto=True, margin=15)
        self.alias_nb_pages() 
        self.slytherins_logo_path = 'logo_ddea_pdf.png' 
        self.font_name = 'Helvetica'
        
        self.font_paths = {
            'Regular': 'fonts/DejaVuSans.ttf',
            'Bold': 'fonts/DejaVuSans-Bold.ttf',
            'Italic': 'fonts/DejaVuSans-Oblique.ttf',
            'BoldItalic': 'fonts/DejaVuSans-BoldOblique.ttf'
        }
        
        try:
            if os.path.exists(self.font_paths['Regular']):
                self.add_font('DejaVu', '', self.font_paths['Regular'])
                if os.path.exists(self.font_paths['Bold']):
                    self.add_font('DejaVu', 'B', self.font_paths['Bold'])
                if os.path.exists(self.font_paths['Italic']):
                    self.add_font('DejaVu', 'I', self.font_paths['Italic'])
                self.font_name = 'DejaVu'
            else:
                self.font_name = 'Helvetica'
        except Exception: 
            self.font_name = 'Helvetica' 
        
        self.set_font(self.font_name, '', 10) 

    def _standardize_text(self, text):
        text = str(text) 
        replacements = { "’": "'", "‘": "'", "“": '"', "”": '"', "–": "-", "—": "--"}
        for unicode_char, ascii_char in replacements.items():
            text = text.replace(unicode_char, ascii_char)
        return text

    def header(self):
        current_font_family = self.font_family; current_font_style = self.font_style; current_font_size = self.font_size_pt
        self.set_font(self.font_name, 'B', 10 if self.page_no() == 1 else 8)
        if self.page_no() == 1:
            try:
                if os.path.exists(self.slytherins_logo_path): self.image(self.slytherins_logo_path, x=10, y=8, w=30); self.ln(5) 
            except: pass
            self.set_font(self.font_name, 'B', 18)
            self.cell(0, 10, self._standardize_text('DDEA Analysis Report'), 0, 1, 'C')
            self.set_font(self.font_name, '', 10)
            description_text = ('Diagonal Differential Expression Alley (DDEA)\n'
                                'Developed by the EvoMol-Lab (github.com/evomol-lab).\n'
                                'BioME, UFRN, Brazil (bioinfo.imd.ufrn.br).')
            self.multi_cell(0, 5, self._standardize_text(description_text), 0, 'C'); self.ln(5)
            self.set_font(self.font_name, 'I', 9)
            self.cell(0, 8, self._standardize_text(f'Report Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'), 0, 1, 'C'); self.ln(10)
        else: 
            self.set_font(self.font_name, 'I', 8)
            self.cell(0, 10, self._standardize_text('DDEA Report'), 0, 0, 'L')
            self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', 0, 0, 'R'); self.ln(10) 
        self.set_font(current_font_family, style=current_font_style, size=current_font_size)

    def footer(self):
        self.set_y(-15); self.set_font(self.font_name, 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', 0, 0, 'C')

    def chapter_title(self, title):
        self.set_font(self.font_name, 'B', 14); self.ln(10) 
        self.cell(0, 10, self._standardize_text(title), 0, 1, 'L'); self.ln(4) 

    def section_title(self, title): 
        self.set_font(self.font_name, 'B', 12); self.ln(5)
        self.cell(0, 8, self._standardize_text(title), 0, 1, 'L')
        
    def body_text(self, text): 
        self.set_font(self.font_name, '', 9) 
        self.multi_cell(0, 5, self._standardize_text(text)); self.ln(2)

    def add_metric(self, label, value):
        self.set_font(self.font_name, 'B', 10)
        self.cell(60, 7, self._standardize_text(label) + ":", ln=0) 
        self.set_font(self.font_name, '', 10)
        self.multi_cell(0, 7, self._standardize_text(value), ln=1); self.ln(1)

    def add_plotly_figure_to_pdf(self, fig, title, caption=None, fig_width_mm=170):
        if fig is None: return
        self.section_title(title) 
        if caption: self.body_text(caption) 
        try:
            img_bytes = fig.to_image(format="png", scale=1.5) 
            img_file = io.BytesIO(img_bytes)
            page_width = self.w - 2 * self.l_margin 
            img_render_width = min(fig_width_mm, page_width) 
            x_pos = (self.w - img_render_width) / 2
            self.image(img_file, x=x_pos, w=img_render_width); self.ln(5) 
        except Exception as e:
            self.set_font(self.font_name, 'I', 9); self.set_text_color(255, 0, 0) 
            self.multi_cell(0, 5, self._standardize_text(f"Error embedding plot '{title}': {e}"))
            self.set_text_color(0, 0, 0); self.ln()

def generate_pdf_report(report_elements_list):
    pdf = PDF()
    pdf.add_page()
    pdf.current_part = 0 
    for element in report_elements_list:
        part = element.get("part")
        if pdf.current_part != part: 
            if part == 1: pdf.chapter_title("Differential Expression Overview")
            elif part == 2: pdf.add_page(); pdf.chapter_title("Top Differentially Expressed Genes")
            elif part == 3: pdf.add_page(); pdf.chapter_title("Clustering & Variation")
            pdf.current_part = part
        if element["type"] == "metric":
            pdf.add_metric(element["label"], element["value"])
        elif element["type"] == "plot": 
            pdf.add_plotly_figure_to_pdf(element["fig"], element["title"], element.get("caption"))
    pdf_output_object = pdf.output(dest='S')
    return bytes(pdf_output_object) if not isinstance(pdf_output_object, str) else pdf_output_object.encode('latin-1')

# ============================================================
# NORMALIZAÇÃO E UTILITÁRIOS (ATUALIZADO COM POLARS/GC)
# ============================================================

def quantile_normalize(df_values):
    if df_values.size == 0: return df_values
    # ATUALIZAÇÃO: Usa float32 para economizar RAM
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
# PARSE E EXTRAÇÃO GEO (ATUALIZADO COM POLARS)
# ============================================================

def _parse_matrix_bytes(raw_bytes, filename=""):
    try:
        content = gzip.decompress(raw_bytes) if filename.endswith('.gz') else raw_bytes
        # ATUALIZAÇÃO: Usa Polars para carregar dados grandes com eficiência
        try:
            df_pl = pl.read_csv(io.BytesIO(content), infer_schema_length=0, ignore_errors=True)
            df = df_pl.to_pandas()
            df = df.set_index(df.columns[0])
            df.index = df.index.astype(str).str.strip().str.replace('"', '')
            del df_pl, content; gc.collect()
            
            # Filtro Acadêmico
            cols_to_drop = ['Length', 'Chr', 'Start', 'End', 'Strand', 'Geneid', 'gene_name']
            df = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
            
            if 'ReadCount' in df.columns and len(df.columns) == 1:
                df = df.rename(columns={'ReadCount': filename.split('.')[0].split('_')[0]})

            return df.select_dtypes(include=[np.number]).astype(np.float32)
        except: 
            return None
    except: return None

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

            if df.empty or len(df.columns) < 2: return None, meta_df, gsms, gsm_order, det_type
            df = df.set_index(0)
            df.index = df.index.astype(str).str.strip().str.replace('"', '')
            df.rename(columns={col: gsm_order[i] for i, col in enumerate(df.columns) if i < len(gsm_order)}, inplace=True)
            return df.select_dtypes(include=[np.number]).astype(np.float32), meta_df, gsms, gsm_order, det_type
    except: return None, None, None, [], 'unknown'

def _sync_suppl_columns_with_gsms(df_suppl, gsm_order):
    cols, gsm_set = list(df_suppl.columns), set(gsm_order)
    direct = [c for c in cols if c in gsm_set]
    if len(direct) >= 2: return df_suppl[direct], "interseção"
    if len(cols) == len(gsm_order): return df_suppl.rename(columns=dict(zip(cols, gsm_order))), "posicional"
    return df_suppl, "não sincronizado"

def get_geo_full_data(gse_id, mode, log_cb=None):
    if log_cb: log_cb("Buscando metadados...")
    df_matrix, meta_df, gsms, gsm_order, detected_type = _try_series_matrix(gse_id)
    if meta_df is None: return None, None, None, [], None, "Falha GEO.", 'unknown'
    if df_matrix is not None and df_matrix.shape[1] >= 2 and (mode != "RNASeq" or df_matrix.max().max() > 50):
        return df_matrix, meta_df, gsms, gsm_order, "Series Matrix", None, detected_type

    num = gse_id.replace("GSE", "")
    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{f'GSE{num[:-3]}nnn' if len(num) > 3 else 'GSEnnn'}/{gse_id}/suppl/"
    try:
        r = requests.get(base_url, timeout=20, headers=HEADERS)
        if r.status_code == 200:
            urls = [base_url + f[0] for f in re.findall(r'href="([^"]+\.(txt|tsv|csv|tar)(\.gz)?)"', r.text, re.IGNORECASE)]
            for u in sorted(urls, key=lambda x: 3 if 'count' in x.lower() or 'raw' in x.lower() else 0, reverse=True):
                df_suppl = _parse_matrix_bytes(requests.get(u, timeout=60, headers=HEADERS).content, u)
                if df_suppl is not None and df_suppl.shape[1] >= 2:
                    df_synced, msg = _sync_suppl_columns_with_gsms(df_suppl, gsm_order)
                    return df_synced, meta_df, gsms, gsm_order, f"Suppl ({msg})", None, detected_type
    except: pass
    return None, meta_df, gsms, gsm_order, None, None, detected_type

# ============================================================
# APP PRINCIPAL (LÓGICA DE DDEA INTEGRADA)
# ============================================================

def run_app():
    st.set_page_config(layout="wide", page_title="DDEA Analytics Master")
    
    with st.sidebar:
        st.title("DDEA Analytics")
        st.divider()
    
    st.title("Diagonal Differential Expression Alley 🧬")

    if 'groups' not in st.session_state: st.session_state['groups'] = {}

    with st.sidebar:
        st.header("1. GEO Input")
        mode = st.radio("Experiment Type:", ["Microarray", "RNASeq"])
        gse_input = st.text_input("GSE ID:", placeholder="Ex: GSE117769")
        fetch_btn = st.button("🚀 Fetch Data", use_container_width=True)

        if 'meta_df' in st.session_state:
            st.divider()
            meta_df_cols = st.session_state['meta_df'].columns.tolist()
            candidatas = ['source_name_ch1', 'source', 'title', 'characteristics_ch1']
            default_cols = [c for c in meta_df_cols if any(cand in c.lower() for cand in candidatas)]
            if not default_cols and len(meta_df_cols) > 1:
                default_cols = [meta_df_cols[1]]

            st.header("2. Parameters")
            extra_selected = st.multiselect("Informações no seletor:", [c for c in meta_df_cols if c != 'Accession'], default=default_cols)
            label_cols = ['Accession'] + extra_selected

            st.divider()
            apply_log = st.checkbox("Aplicar Log2", value=True)
            gene_area = st.text_area("Genes Highlight (1 per line):")
            p_thr = st.slider("FDR threshold (P-adj):", 0.001, 0.10, 0.05, format="%.3f")
            fc_thr = st.slider("Min Abs Log2FC:", 0.0, 10.0, 1.0, step=0.1)
            max_plot = st.number_input("Top genes p/ Heatmap:", value=50, min_value=5)

            if st.session_state.get('analysis_done'):
                if st.button("📄 Generate PDF Report", use_container_width=True):
                    elements = [
                        {"type": "metric", "label": "GSE ID", "value": st.session_state.get('gse_id', 'N/A'), "part": 1},
                        {"type": "metric", "label": "Comparison", "value": f"{st.session_state['tn']} vs {st.session_state['rn']}", "part": 1},
                        {"type": "metric", "label": "Parameters Used", "value": f"FDR < {p_thr}, |Log2FC| > {fc_thr}", "part": 1},
                        {"type": "metric", "label": "Total DEGs Found", "value": str(len(st.session_state.get('df_diff', []))), "part": 1},
                        {"type": "plot", "title": "Volcano Plot", "fig": st.session_state.get('fig_v'), "part": 1},
                        {"type": "plot", "title": "PCA Analysis", "fig": st.session_state.get('fig_pca'), "part": 3},
                        {"type": "plot", "title": "Expression Heatmap", "fig": st.session_state.get('fig_h'), "part": 3}
                    ]
                    with st.spinner("Gerando PDF..."):
                        pdf_bytes = generate_pdf_report(elements)
                        st.download_button(label="📥 Download PDF", data=pdf_bytes, file_name="DDEA_Report.pdf", mime="application/pdf")

    if fetch_btn and gse_input:
        with st.spinner("🚀 Buscando dados..."):
            df, meta_df, gsms, gsm_order, source, err, detected_type = get_geo_full_data(gse_input.strip(), mode)
            if meta_df is not None:
                st.session_state.update({'mode': mode, 'gse_id': gse_input.strip(), 'meta_df': meta_df, 'gsm_order': gsm_order, 'df': df, 'detected_type': detected_type})
                if df is not None:
                    id_type = detect_index_type(df.index.tolist())
                    mapping, _ = get_gene_mapping_rnaseq(tuple(df.index.astype(str).tolist()), id_type) if mode == "RNASeq" else get_gene_mapping_microarray(gse_input.strip())
                    st.session_state.update({'mapping': mapping, 'id_type': id_type, 'matrix_source': source})
                st.rerun()

    if 'meta_df' not in st.session_state: return

    meta = st.session_state['meta_df'].copy()
    meta['display_label'] = meta.apply(lambda row: " | ".join([str(row[c]) for c in label_cols if c in row.index]), axis=1)
    label_to_gsm = dict(zip(meta['display_label'], meta['Accession']))

    st.subheader("🛠️ Group Management")
    c_i, c_b = st.columns([3, 1])
    def add_group_callback():
        if st.session_state.new_group_input.strip():
            st.session_state['groups'][st.session_state.new_group_input.strip()] = []
            st.session_state.new_group_input = ""
    c_i.text_input("New Group Name:", key="new_group_input")
    c_b.button("➕ Add Group", use_container_width=True, on_click=add_group_callback)

    cols_g = st.columns(2)
    for i, (g_name, g_samples) in enumerate(list(st.session_state['groups'].items())):
        avail = [lab for lab in meta['display_label'] if lab not in sum(st.session_state['groups'].values(), []) or lab in g_samples]
        st.session_state['groups'][g_name] = cols_g[i % 2].multiselect(f"**{g_name}**", avail, default=g_samples, key=f"s_{g_name}")
        if cols_g[i % 2].button("🗑️ Remove", key=f"del_{g_name}"): del st.session_state['groups'][g_name]; st.rerun()

    if st.session_state.get('df') is None:
        st.warning("Matriz não encontrada. Faça upload manual.")
        return

    if len(st.session_state['groups']) >= 2:
        group_names = list(st.session_state['groups'].keys())
        c1, c2 = st.columns(2)
        ref_g = c1.selectbox("Controle:", group_names)
        test_g = c2.selectbox("Teste:", [g for g in group_names if g != ref_g])

        if st.button("🔥 Run Analysis", use_container_width=True):
            with st.spinner("Processando..."):
                def find_col(gsm_id, columns):
                    for col in columns:
                        if gsm_id in col: return col
                    return None
                
                c_ref = [find_col(label_to_gsm[lab], st.session_state['df'].columns) for lab in st.session_state['groups'][ref_g] if find_col(label_to_gsm[lab], st.session_state['df'].columns)]
                c_test = [find_col(label_to_gsm[lab], st.session_state['df'].columns) for lab in st.session_state['groups'][test_g] if find_col(label_to_gsm[lab], st.session_state['df'].columns)]
                
                # ATUALIZAÇÃO: Uso de Float32 e limpeza agressiva
                df_sub = st.session_state['df'][c_ref + c_test].astype(np.float32).fillna(0)
                df_sub = df_sub[df_sub.sum(axis=1) > 0]
                
                data_norm_vals = quantile_normalize(df_sub.values)
                data_norm = pd.DataFrame(data_norm_vals, columns=c_ref + c_test, index=df_sub.index)
                
                # Teste T Simplificado (Rápido para RAM limitada)
                lfc = np.nanmean(data_norm[c_test].values, axis=1) - np.nanmean(data_norm[c_ref].values, axis=1)
                pvals = stats.ttest_ind(data_norm[c_test], data_norm[c_ref], axis=1).pvalue
                fdr = fdrcorrection(np.nan_to_num(pvals, nan=1.0))[1]
                
                res = pd.DataFrame({'Probe_ID': df_sub.index, 'Log2FC': lfc, 'FDR': fdr}).dropna()
                mapping = st.session_state.get('mapping')
                if mapping is not None:
                    res = res.merge(mapping, on='Probe_ID', how='left')
                    res['Symbol'] = res['Symbol'].fillna(res['Probe_ID'])
                else: res['Symbol'] = res['Probe_ID']

                st.session_state.update({'res': res, 'norm_df': data_norm, 'analysis_done': True, 'rn': ref_g, 'tn': test_g, 'c_ref': c_ref, 'c_test': c_test})
                st.rerun()

    if st.session_state.get('analysis_done'):
        res = st.session_state['res'].copy()
        res['Status'] = 'Not Significant'
        res.loc[(res['FDR'] < p_thr) & (res['Log2FC'] >= fc_thr), 'Status'] = 'Up-regulated'
        res.loc[(res['FDR'] < p_thr) & (res['Log2FC'] <= -fc_thr), 'Status'] = 'Down-regulated'
        
        df_diff = res[res['Status'] != 'Not Significant'].sort_values('FDR')
        st.session_state['df_diff'] = df_diff

        st.subheader(f"Resultados: {st.session_state['tn']} vs {st.session_state['rn']}")
        m1, m2, m3 = st.columns(3)
        m1.metric("DEGs", len(df_diff))
        m2.metric("Up 📈", len(df_diff[df_diff['Status'] == 'Up-regulated']))
        m3.metric("Down 📉", len(df_diff[df_diff['Status'] == 'Down-regulated']))

        t1, t2 = st.tabs(["📊 Gráficos", "📋 Tabela"])
        with t1:
            st.plotly_chart(px.scatter(res, x='Log2FC', y=-np.log10(res['FDR']+1e-300), color='Status', hover_name='Symbol', template="simple_white"), use_container_width=True)
        with t2:
            st.dataframe(df_diff[['Symbol', 'Log2FC', 'FDR', 'Status']], use_container_width=True)

if __name__ == '__main__':
    run_app()
