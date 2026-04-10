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
                # Registra Regular
                self.add_font('DejaVu', '', self.font_paths['Regular'])
                # Registra Bold (B)
                if os.path.exists(self.font_paths['Bold']):
                    self.add_font('DejaVu', 'B', self.font_paths['Bold'])
                # Registra Italic (I) - RESOLVE O ERRO ATUAL
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
# NORMALIZAÇÃO E UTILITÁRIOS
# ============================================================

def quantile_normalize(df_values):
    if df_values.size == 0: return df_values
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
# PARSE E EXTRAÇÃO GEO
# ============================================================

def _parse_matrix_bytes(raw_bytes, filename=""):
    try:
        content = gzip.decompress(raw_bytes) if filename.endswith('.gz') else raw_bytes
        for sep in ['\t', ',']:
            try:
                df = pd.read_csv(io.BytesIO(content), sep=sep, index_col=0, low_memory=False)
                df.index = df.index.astype(str).str.strip().str.replace('"', '')
                
                # Filtro Acadêmico: Remove colunas técnicas do featureCounts
                cols_to_drop = ['Length', 'Chr', 'Start', 'End', 'Strand', 'Geneid', 'gene_name']
                df = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
                
                if 'ReadCount' in df.columns and len(df.columns) == 1:
                    df = df.rename(columns={'ReadCount': filename.split('.')[0].split('_')[0]})

                if df.select_dtypes(include=[np.number]).shape[1] >= 1:
                    return df.select_dtypes(include=[np.number])
            except: continue
    except: return None

def _try_series_matrix(gse_id):
    num = gse_id.replace("GSE", "")
    prefix = f"GSE{num[:-3]}nnn" if len(num) > 3 else "GSEnnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    try:
        r = requests.get(url, timeout=60, headers=HEADERS)
        if r.status_code != 200: return None, None, None, [], 'unknown'
        with gzip.open(io.BytesIO(r.content), 'rt') as f:
            titles, gsms, char_lines, gsm_order, stype = [], [], [], [], ""
            df = pd.DataFrame()
            for line in f:
                if line.startswith('!Series_type'): stype = line.lower()
                elif line.startswith('!Sample_title'): titles = [t.strip().replace('"', '') for t in line.split('\t')[1:]]
                elif line.startswith('!Sample_geo_accession'):
                    gsms = [t.strip().replace('"', '') for t in line.split('\t')[1:]]
                    gsm_order = gsms[:]
                elif line.startswith('!Sample_characteristics_ch1'): char_lines.append(line.split('\t')[1:])
                elif line.startswith('ID_REF') or line.startswith('"ID_REF"'):
                    df = pd.read_csv(f, sep='\t', header=None, low_memory=True)
                    break

            det_type = 'RNASeq' if 'sequencing' in stype else 'Microarray' if 'array' in stype else 'unknown'
            meta_dict = {"Accession": gsms, "Title": titles}
            for row in char_lines:
                if row: meta_dict[row[0].split(':')[0].strip() if ':' in row[0] else "Info"] = [v.split(': ')[1].strip() if ': ' in v else v.strip() for v in row][:len(gsms)]
            meta_df = pd.DataFrame(meta_dict)

            if df.empty or len(df.columns) < 2: return None, meta_df, gsms, gsm_order, det_type
            df = df.set_index(0)
            df.index = df.index.astype(str).str.strip().str.replace('"', '')
            df.rename(columns={col: gsm_order[i] for i, col in enumerate(df.columns) if i < len(gsm_order)}, inplace=True)
            return df.select_dtypes(include=[np.number]), meta_df, gsms, gsm_order, det_type
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
# APP PRINCIPAL
# ============================================================

def run_app():
    st.set_page_config(layout="wide", page_title="DDEA Analytics Master")
    with st.sidebar:
        if os.path.exists("logo_ddea_streamlit.png"):
            st.image("logo_ddea_streamlit.png", use_container_width=True)
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
            det = st.session_state.get('detected_type', 'unknown')
            
            if det != 'unknown':
                st.caption(f"🏷️ Detectado no GEO: **{det}**")
            else:
                st.caption(f"🏷️ Análise configurada para: **{mode}**")
            
            st.divider()
            extra_cols = [c for c in st.session_state['meta_df'].columns if c != 'Accession']
            extra_selected = st.multiselect("Informações no seletor:", extra_cols, default=["Title"] if "Title" in extra_cols else extra_cols[:1])
            label_cols = ['Accession'] + extra_selected

            st.divider()
            st.header("2. Parameters")
            apply_log = st.checkbox("Aplicar Log2 (Desmarque se os dados já estiverem em log)", value=True)
            gene_area = st.text_area("Genes Highlight (1 per line):")
            p_thr = st.slider("FDR threshold (P-adj):", 0.001, 0.10, 0.05, format="%.3f")
            fc_thr = st.slider("Min Abs Log2FC:", 0.0, 10.0, 1.0, step=0.1)
            use_limma = st.checkbox("Usar modelo linear", value=True) if mode == "Microarray" else False
            max_plot = st.number_input("Top genes p/ Heatmap:", value=50, min_value=5)

            st.divider()
            st.header("3. Export Report")
            if st.session_state.get('analysis_done'):
                if st.button("📄 Generate PDF Report", use_container_width=True):
                    elements = [
                        {"type": "metric", "label": "Comparison", "value": f"{st.session_state['tn']} vs {st.session_state['rn']}", "part": 1},
                        {"type": "metric", "label": "DEGs Count", "value": str(len(st.session_state.get('df_diff', []))), "part": 1},
                        {"type": "plot", "title": "Volcano Plot", "fig": st.session_state.get('fig_v'), "part": 1},
                        {"type": "plot", "title": "Top Up-regulated Genes", "fig": st.session_state.get('fig_up'), "part": 2},
                        {"type": "plot", "title": "Top Down-regulated Genes", "fig": st.session_state.get('fig_down'), "part": 2},
                        {"type": "plot", "title": "PCA Analysis", "fig": st.session_state.get('fig_pca'), "part": 3},
                        {"type": "plot", "title": "Expression Heatmap", "fig": st.session_state.get('fig_h'), "part": 3}
                    ]
                    pdf_bytes = generate_pdf_report(elements)
                    st.download_button("📥 Download PDF Report", pdf_bytes, "DDEA_Report.pdf", "application/pdf", use_container_width=True)
            else:
                st.info("Execute a análise para habilitar o relatório.")
        else:
            label_cols = ['Accession', 'Title']; gene_area = ""; p_thr = 0.05; fc_thr = 1.0; use_limma = False; max_plot = 50; apply_log = True

    # ---- FETCH LOGIC ----
    if fetch_btn and gse_input:
        with st.spinner("🚀 Buscando dados..."):
            df, meta_df, gsms, gsm_order, source, err, detected_type = get_geo_full_data(gse_input.strip(), mode)
            if meta_df is not None:
                keys_to_reset = ['analysis_done', 'res', 'norm_df', 'df_diff', 'groups', 
                                 'fig_v', 'fig_pca', 'fig_up', 'fig_down', 'fig_h']
                for key in keys_to_reset:
                    if key in st.session_state:
                        if key == 'groups':
                            st.session_state[key] = {}
                        else:
                            del st.session_state[key]

                st.session_state.update({
                    'mode': mode, 
                    'gse_id': gse_input.strip(), 
                    'meta_df': meta_df, 
                    'gsm_order': gsm_order, 
                    'df': df,
                    'detected_type': detected_type
                })
                
                if df is not None:
                    id_type = detect_index_type(df.index.tolist())
                    mapping, msg = get_gene_mapping_rnaseq(tuple(df.index.astype(str).tolist()), id_type) if mode == "RNASeq" else get_gene_mapping_microarray(gse_input.strip())
                    st.session_state.update({'mapping': mapping, 'id_type': id_type, 'matrix_source': source})
                    st.success(f"✅ Matriz via **{source}** — {df.shape[0]} genes × {df.shape[1]} amostras")
                
                st.rerun()

    if 'meta_df' not in st.session_state: return

    # ---- METADATA & GROUPS ----
    meta = st.session_state['meta_df'].copy()
    with st.expander("📊 Metadata Explorer"): st.dataframe(meta, use_container_width=True)
    meta['display_label'] = meta.apply(lambda row: " | ".join([str(row[c]) for c in label_cols if c in row.index]), axis=1)
    label_to_gsm = dict(zip(meta['display_label'], meta['Accession']))

    st.subheader("🛠️ Group Management")
    c_i, c_b = st.columns([3, 1])
    def add_group_callback():
        new_g_name = st.session_state.new_group_input.strip()
        if new_g_name:
            st.session_state['groups'][new_g_name] = []
            st.session_state.new_group_input = ""

    c_i.text_input("New Group Name:", key="new_group_input")
    c_b.button("➕ Add Group", use_container_width=True, on_click=add_group_callback)

    cols_g = st.columns(2)
    for i, (g_name, g_samples) in enumerate(list(st.session_state['groups'].items())):
        avail = [lab for lab in meta['display_label'] if lab not in sum(st.session_state['groups'].values(), []) or lab in g_samples]
        st.session_state['groups'][g_name] = cols_g[i % 2].multiselect(f"**{g_name}**", avail, default=[s for s in g_samples if s in avail], key=f"s_{g_name}")
        if cols_g[i % 2].button("🗑️ Remove", key=f"del_{g_name}"):
            del st.session_state['groups'][g_name]
            st.rerun()

    # ---- UPLOAD MANUAL MULTI-FILE ----
    if st.session_state.get('df') is None:
        st.warning("📦 Matriz não encontrada. Faça upload de UM ou MAIS arquivos de contagem.")
        up_files = st.file_uploader("Upload counts files", accept_multiple_files=True)
        if up_files:
            combined_df = pd.DataFrame()
            for uf in up_files:
                dft = _parse_matrix_bytes(uf.read(), uf.name)
                if dft is not None: combined_df = dft if combined_df.empty else combined_df.join(dft, how='outer')
            if not combined_df.empty:
                st.session_state['df'] = combined_df.fillna(0)
                st.rerun()
        return

    # ---- ANALYSIS ENGINE ----
    df_matrix = st.session_state['df']
    matrix_cols = set(df_matrix.columns.astype(str))

    if len(st.session_state['groups']) >= 2:
        st.divider()
        group_names = list(st.session_state['groups'].keys())
        
        c1, c2 = st.columns(2)
        ref_g = c1.selectbox("Referência (Controle):", group_names)
        test_g = c2.selectbox("Teste:", [g for g in group_names if g != ref_g])

        if st.button("🔥 Run Analysis", use_container_width=True):
            with st.spinner("Processando Estatística Acadêmica..."):
                
                # BUSCA FLEXÍVEL DE COLUNAS
                def find_col(gsm_id, columns):
                    for col in columns:
                        if gsm_id in col: return col
                    return None

                c_ref = [c for c in [find_col(label_to_gsm[lab], matrix_cols) for lab in st.session_state['groups'][ref_g]] if c]
                c_test = [c for c in [find_col(label_to_gsm[lab], matrix_cols) for lab in st.session_state['groups'][test_g]] if c]

                if not c_ref or not c_test:
                    st.error(f"❌ Erro de Mapeamento. Colunas disponíveis na matriz: {list(matrix_cols)[:5]}")
                    st.stop()

                df_sub = df_matrix[c_ref + c_test].apply(pd.to_numeric, errors='coerce').fillna(0)
                # Remove genes com soma zero (ruído)
                df_sub = df_sub[df_sub.sum(axis=1) > 0]
                data_vals = df_sub.values.astype(np.float32)

                # NORMALIZAÇÃO
                if st.session_state.mode == "Microarray": 
                    data_norm_vals = quantile_normalize(data_vals)
                else: 
                    if apply_log and np.nanmax(data_vals) > 50:
                        data_norm_vals = np.log2(np.clip(data_vals, 0, None) + 1.0)
                    else:
                        data_norm_vals = data_vals
                        
                data_norm = pd.DataFrame(data_norm_vals, columns=c_ref + c_test, index=df_sub.index)
                m1, m2 = data_norm[c_ref].values, data_norm[c_test].values

                # TESTES ESTATÍSTICOS
                is_int = np.all(np.equal(np.mod(data_vals, 1), 0))
                if st.session_state.mode == "RNASeq" and is_int and HAS_DESEQ2:
                    meta_ds = pd.DataFrame({'cond': ['C']*len(c_ref) + ['T']*len(c_test)}, index=c_ref + c_test)
                    dds = DeseqDataSet(counts=df_sub.T.astype(int), metadata=meta_ds, design_factors="cond", quiet=True)
                    dds.deseq2()
                    stat_res = DeseqStats(dds, contrast=["cond", "T", "C"], quiet=True)
                    stat_res.summary()
                    res_df = stat_res.results_df
                    res = pd.DataFrame({'Probe_ID': res_df.index, 'Log2FC': res_df['log2FoldChange'], 'FDR': res_df['padj']}).dropna()
                else:
                    lfc = np.nanmean(m2, axis=1) - np.nanmean(m1, axis=1)
                    pvals = stats.ttest_ind(m2, m1, axis=1, equal_var=False, nan_policy='omit').pvalue
                    mask = ~np.isnan(pvals)
                    fdr = np.full(pvals.shape, np.nan)
                    fdr[mask] = fdrcorrection(pvals[mask])[1]
                    res = pd.DataFrame({'Probe_ID': df_sub.index, 'Log2FC': lfc, 'FDR': fdr}).dropna()

                # ANOTAÇÃO DE GENES
                mapping = st.session_state.get('mapping')
                if mapping is not None:
                    res = res.merge(mapping, on='Probe_ID', how='left')
                    res['Symbol'] = res['Symbol'].replace(['nan', 'None'], np.nan).fillna(res['Probe_ID'])
                else: res['Symbol'] = res['Probe_ID']

                st.session_state.update({'res': res, 'norm_df': data_norm, 'analysis_done': True, 'rn': ref_g, 'tn': test_g, 'c_ref': c_ref, 'c_test': c_test})
                st.rerun()

    # ---- RESULTS DASHBOARD ----
    if st.session_state.get('analysis_done'):
        res, df_norm = st.session_state['res'].copy(), st.session_state['norm_df']
        
        if gene_area:
            res = res[res['Symbol'].str.upper().isin([g.strip().upper() for g in gene_area.split('\n') if g.strip()])]

        res['Status'] = 'Not Significant'
        res.loc[(res['FDR'] < p_thr) & (res['Log2FC'] >= fc_thr), 'Status'] = 'Up-regulated'
        res.loc[(res['FDR'] < p_thr) & (res['Log2FC'] <= -fc_thr), 'Status'] = 'Down-regulated'
        res['-log10(FDR)'] = -np.log10(res['FDR'] + 1e-300)

        df_diff = res[res['Status'] != 'Not Significant'].sort_values('FDR')
        st.session_state['df_diff'] = df_diff # Salva para o PDF
        up_list = df_diff[df_diff['Status'] == 'Up-regulated']
        down_list = df_diff[df_diff['Status'] == 'Down-regulated']

        st.subheader(f"🚀 Results: {st.session_state['tn']} vs {st.session_state['rn']}")
        m1, m2, m3 = st.columns(3)
        m1.metric("Total DEGs (FDR < {})".format(p_thr), len(df_diff))
        m2.metric("Up-Regulated 📈", len(up_list))
        m3.metric("Down-Regulated 📉", len(down_list))

        t1, t2, t3, t4 = st.tabs(["📊 Volcano & PCA", "🔝 Top DEGs", "🔥 Heatmap", "📋 Data Table"])

        with t1:
            fig_v = px.scatter(res, x='Log2FC', y='-log10(FDR)', color='Status', color_discrete_map={'Not Significant': 'lightgrey', 'Up-regulated': 'red', 'Down-regulated': 'blue'}, hover_name='Symbol')
            fig_v.add_vline(x=fc_thr, line_dash="dash"); fig_v.add_vline(x=-fc_thr, line_dash="dash")
            fig_v.add_hline(y=-np.log10(p_thr), line_dash="dash")
            st.session_state['fig_v'] = fig_v.update_layout(template="simple_white", height=500)
            st.plotly_chart(st.session_state['fig_v'], use_container_width=True)
            
            try:
                pca = PCA(n_components=2)
                pc = pca.fit_transform(df_norm.fillna(0).T)
                df_pc = pd.DataFrame(pc, columns=['PC1', 'PC2'], index=df_norm.columns)
                df_pc['Group'] = [st.session_state['rn'] if x in st.session_state['c_ref'] else st.session_state['tn'] for x in df_pc.index]
                # SALVANDO o PCA
                st.session_state['fig_pca'] = px.scatter(df_pc, x='PC1', y='PC2', color='Group', text=df_pc.index).update_traces(textposition='top center')
                st.plotly_chart(st.session_state['fig_pca'], use_container_width=True)
            except: pass

        with t2:
            ca, cb = st.columns(2)
            if not up_list.empty:
                # SALVANDO Top Up
                st.session_state['fig_up'] = px.bar(up_list.head(20), x='Symbol', y='Log2FC', color='Log2FC', color_continuous_scale='Reds', title="Top 20 UP")
                ca.plotly_chart(st.session_state['fig_up'], use_container_width=True)
            if not down_list.empty:
                # SALVANDO Top Down
                st.session_state['fig_down'] = px.bar(down_list.head(20), x='Symbol', y='Log2FC', color='Log2FC', color_continuous_scale='Blues_r', title="Top 20 DOWN")
                cb.plotly_chart(st.session_state['fig_down'], use_container_width=True)

        with t3:
            if not df_diff.empty:
                valid_ids = [idx for idx in df_diff.head(max_plot)['Probe_ID'] if idx in df_norm.index]
                h_mat = df_norm.loc[valid_ids].values
                h_z = (h_mat - np.mean(h_mat, axis=1, keepdims=True)) / (np.std(h_mat, axis=1, keepdims=True) + 1e-9)
                # SALVANDO Heatmap
                st.session_state['fig_h'] = go.Figure(data=go.Heatmap(z=h_z, x=df_norm.columns, y=df_diff.head(len(valid_ids))['Symbol'].tolist(), colorscale='RdBu_r', zmid=0)).update_layout(height=max(400, len(valid_ids)*20))
                st.plotly_chart(st.session_state['fig_h'], use_container_width=True)

        with t4:
            st.dataframe(df_diff[['Symbol', 'Log2FC', 'FDR', 'Status']], use_container_width=True)
            st.download_button("📥 Baixar CSV", df_diff[['Symbol', 'Log2FC', 'FDR', 'Status']].to_csv(index=False), f"DEGs.csv", "text/csv")

if __name__ == '__main__':
    run_app()
