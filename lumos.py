#!/usr/bin/env python3
"""Lumos Networks: from transcriptomes to regulatory networks.

Engine and command line of the Lumos Networks suite (EvoMol-Lab, BioME/UFRN):

    dataset    read a GEO series or uploaded tables (expression values + sample metadata)
    de         differential expression between two groups (DDEA): PyDESeq2 for raw counts,
               limma (R) or Welch for other values, with batch, covariates and paired designs
    integrate  align contrasts of one dataset, meta-analyse independent datasets, or compare
               tissues gene by gene
    enrich     over-representation analysis with an explicit background (APP)

STRING and PrioriGraph (PG) networks, propagation, sensitivity and HTML reports are built on
demand by the web interface (LumosGUI.py) from the files these commands write.

    python lumos.py dataset --geo GSE117769 --outdir run/data
    python lumos.py de --dataset run/data --groups groups.json --reference Healthy --test RA --outdir run/de
    python lumos.py enrich --source run/de --genes-up up.txt --genes-down down.txt --outdir run/app

Every command writes summary.json, manifest.json (when applicable) and status.json in --outdir
and prints '@@progress <fraction> <message>' lines that the web interface follows.
"""
import argparse
import gzip
import importlib.util
import io
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback
import zipfile
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from lumos_core import (__version__, align_named_contrasts, annotate_cross_dataset_patterns,
                        build_design_metadata, classify_differential_expression, collapse_symbol_results,
                        conserved_signature, network_sensitivity_analysis, orient_gene_sample_matrix,
                        package_versions, prepare_ora_inputs, quantile_normalize, random_walk_with_restart,
                        sanitize_design_name, validate_count_matrix)

APP_DIR = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
LIMMA_BACKEND = os.path.join(APP_DIR, 'r', 'limma_backend.R')
CACHE_DIR = os.environ.get('LUMOS_CACHE_DIR') or os.path.join(tempfile.gettempdir(), 'lumos_cache')
USER_AGENT = f'LumosNetworks/{__version__} (EvoMol-Lab, BioME/UFRN; academic use)'
GEO_FTP = 'https://ftp.ncbi.nlm.nih.gov/geo'
MYGENE = 'https://mygene.info/v3/query'
ENSEMBL_REST = 'https://rest.ensembl.org'
ENRICHR = 'https://maayanlab.cloud/Enrichr'
SPEEDRICHR = 'https://maayanlab.cloud/speedrichr/api'
STRING_API = 'https://version-12-0.string-db.org/api'

DEFAULT_PATHWAY_LIBRARIES = ('KEGG_2021_Human', 'GO_Biological_Process_2023')
DEFAULT_TF_LIBRARIES = ('TRRUST_Transcription_Factors_2019', 'JASPAR_PWM_Human_2025')
TF_LIBRARIES = ('TRRUST_Transcription_Factors_2019', 'JASPAR_PWM_Human_2025', 'ChEA_2022',
                'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X')
# Whole words only: a plain substring test would also drop "Proliferation" or "Migration".
NON_HUMAN_TERM = re.compile(r'\b(?:mouse|mus musculus|rat|murine)\b', re.I)

# Ensembl gene ID prefixes -> NCBI taxonomy (used for gene-symbol look-ups and STRING)
ENSEMBL_TAXA = {'ENSG': 9606, 'ENSMUSG': 10090, 'ENSRNOG': 10116, 'ENSDARG': 7955,
                'ENSGALG': 9031, 'ENSSSCG': 9823, 'ENSBTAG': 9913, 'ENSMMUG': 9544,
                'ENSCAFG': 9615, 'ENSOARG': 9940, 'ENSECAG': 9796}

PCA_GENES = 500          # most variable genes used for the PCA (as DESeq2's plotPCA)
HEATMAP_POOL = 300       # genes (lowest p-values) whose values are sent for the heatmap
TOP_TERMS = 300          # enrichment terms kept per library in summary.json
MAX_DOWNLOAD_MB = float(os.environ.get('LUMOS_MAX_DOWNLOAD_MB', '2048'))

ID_PATTERNS = (
    ('ensembl', re.compile(r'^ENS[A-Z]*G\d{11}(\.\d+)?(_PAR_Y)?$')),
    ('probe', re.compile(r'^(AFFX[-_].+|.+_(s_|x_|a_|i_|f_|g_|r_)?at|.+_st|ILMN_\d+|A_\d+_P\d+|'
                         r'[A-Z]{1,3}_?\d{6,}(\.\d+)?)$')),
    ('entrez', re.compile(r'^\d+$')),
    ('symbol', re.compile(r'^[A-Za-z][A-Za-z0-9\-\.]{0,19}$')),
)
ANNOTATION_COLUMN = re.compile(
    r'^(chr|chrom|chromosome|start|end|stop|strand|length|gene_?length|effective_?length|width|gc|'
    r'(gene_?)?biotype|gene_?type|description|transcript_?ids?|entrez(_?id)?|tpm|fpkm|rpkm)$', re.I)
SYMBOL_COLUMN = re.compile(r'^(gene_?name|gene_?symbol|symbol|hgnc_symbol|external_gene_name|'
                           r'mgi_symbol|genesymbol)$', re.I)
COUNT_COLUMN = re.compile(r'(numreads|expected_count|count|reads)', re.I)
NON_GENE_ROW = re.compile(r'^(__|N_(unmapped|multimapping|noFeature|ambiguous)$)')
TABLE_FILE = re.compile(r'\.(txt|tsv|csv|tab)(\.gz)?$', re.I)


class LumosError(Exception):
    """An error explained to the user (no traceback)."""


def progress(fraction, message):
    print(f'@@progress {fraction:.3f} {message}', flush=True)


def say(message):
    print(message, flush=True)


# ---------------------------------------------------------------------------- small helpers

_SESSION = None


def session():
    """HTTP session with retries and back-off for busy public services (429, 5xx, timeouts)."""
    global _SESSION
    if _SESSION is None:
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
        retry = Retry(total=3, connect=3, read=2, status=3, backoff_factor=0.6,
                      status_forcelist=(429, 500, 502, 503, 504), allowed_methods=frozenset({'GET', 'POST'}),
                      raise_on_status=False)
        _SESSION = requests.Session()
        _SESSION.mount('https://', HTTPAdapter(max_retries=retry, pool_connections=8, pool_maxsize=8))
        _SESSION.headers['User-Agent'] = USER_AGENT
    return _SESSION


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    tmp = f'{path}.tmp'
    with open(tmp, 'w') as fh:
        json.dump(obj, fh, separators=(',', ':'), allow_nan=False)
    os.replace(tmp, path)


def fnum(x, sig=4):
    """A float rounded to `sig` significant digits, or None (JSON has no NaN)."""
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return float(f'{x:.{sig}g}') if np.isfinite(x) else None


def num_list(values, sig=4):
    return [fnum(x, sig) for x in np.asarray(values, dtype=float)]


def is_number(text):
    if text.upper() in ('NA', 'NAN', 'NULL', ''):
        return True
    try:
        float(text)
        return True
    except ValueError:
        return False


def file_stem(name):
    base = os.path.basename(name)
    base = re.sub(r'\.gz$', '', base, flags=re.I)
    base = re.sub(r'\.(txt|tsv|csv|tab)$', '', base, flags=re.I)
    return re.sub(r'(\.counts?|_ReadsPerGene\.out|\.genes\.results|\.quant(\.sf)?)$', '', base, flags=re.I)


def benjamini_hochberg(p):
    """BH-adjusted p-values; NaN stays NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    q = p[ok]
    n = q.size
    order = np.argsort(q)
    ranked = q[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


def strip_version(gene_id):
    return re.sub(r'\.\d+(_PAR_Y)?$', '', gene_id)


def detect_id_type(ids):
    """'ensembl', 'probe', 'entrez', 'symbol' or 'unknown' from a sample of identifiers."""
    ids = [str(i).strip() for i in ids if str(i).strip()]
    if not ids:
        return 'unknown'
    step = max(1, len(ids) // 500)
    sample = ids[::step][:500]
    scores = {name: sum(1 for s in sample if rx.match(s)) for name, rx in ID_PATTERNS}
    best = max(scores, key=scores.get)
    return best if scores[best] >= 0.3 * len(sample) else 'unknown'


def ensembl_taxid(ids):
    for gene_id in ids[:50]:
        m = re.match(r'^(ENS[A-Z]*G)\d', str(gene_id))
        if m and m.group(1) in ENSEMBL_TAXA:
            return ENSEMBL_TAXA[m.group(1)]
    return None


def value_kind(values):
    """'counts' (non-negative integers), 'log' (log-scale values) or 'continuous'."""
    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        raise LumosError('The expression table has no numeric values.')
    if v.size > 200000:
        v = np.random.default_rng(0).choice(v, 200000, replace=False)
    if v.min() < 0 or v.max() < 50:
        return 'log'
    return 'counts' if np.allclose(v, np.round(v)) else 'continuous'


VALUE_KIND_TEXT = {'counts': 'raw counts', 'log': 'log-scale values', 'continuous': 'non-integer values (e.g. TPM, FPKM, intensities)'}
ID_TYPE_TEXT = {'ensembl': 'Ensembl gene IDs', 'probe': 'microarray probe IDs', 'entrez': 'NCBI Gene (Entrez) IDs',
                'symbol': 'gene symbols', 'unknown': 'unrecognised identifiers'}


def download(url, label, span=(0.0, 1.0), dest=None, limit_mb=None, timeout=120):
    """Download `url` (to `dest` or to memory), reporting progress inside `span`."""
    try:
        r = session().get(url, stream=True, timeout=timeout)
    except Exception as exc:
        raise LumosError(f'Could not reach {url} ({exc}).')
    if r.status_code != 200:
        raise LumosError(f'Download failed (HTTP {r.status_code}): {url}')
    total = int(r.headers.get('content-length') or 0)
    if limit_mb and total > limit_mb * 2 ** 20:
        r.close()
        raise LumosError(f'{os.path.basename(url)} has {total / 2 ** 20:.0f} MB, more than the '
                         f'{limit_mb:.0f} MB limit (LUMOS_MAX_DOWNLOAD_MB).')
    fh = open(dest, 'wb') if dest else io.BytesIO()
    done, last = 0, 0.0
    try:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
            done += len(chunk)
            if time.time() - last > 1.0:
                last = time.time()
                frac = span[0] + (span[1] - span[0]) * (done / total if total else 0.0)
                size = f'{done / 2 ** 20:.0f} MB' + (f' of {total / 2 ** 20:.0f} MB' if total else '')
                progress(frac, f'{label} ({size})')
    finally:
        if dest:
            fh.close()
    return dest if dest else fh.getvalue()


# ---------------------------------------------------------------------------- reading tables

def _decode(data, name):
    if data[:2] == b'\x1f\x8b':
        data = gzip.decompress(data)
    return data.decode('utf-8', errors='replace')


def _sniff(lines):
    """Separator and split function for a delimited text table."""
    first = lines[0]
    counts = {'\t': first.count('\t'), ',': first.count(','), ';': first.count(';')}
    sep = max(counts, key=counts.get)
    if counts[sep] == 0:
        return r'\s+', lambda line: [f.strip().strip('"') for f in line.split()]
    return sep, lambda line: [f.strip().strip('"') for f in line.split(sep)]


def read_table(data, name):
    """Expression values (genes x columns) and an optional gene-symbol column from a text table.

    Handles CSV/TSV/whitespace tables, gzip, R-style headers (one name fewer than the
    columns), headerless files (HTSeq), featureCounts comment lines and annotation columns.
    """
    text = _decode(data, name)
    lines = [ln for ln in text.splitlines() if ln.strip()]
    start = 0
    while start < len(lines) and lines[start].lstrip().startswith(('#', '!')):
        start += 1
    lines = lines[start:]
    if len(lines) < 2:
        raise LumosError(f'{name}: not a table (fewer than two lines).')
    sep, split = _sniff(lines)
    head, second = split(lines[0]), split(lines[1])
    has_header = any(not is_number(x) for x in head[1:] if x != '') or len(head) < len(second)
    if has_header:
        if len(head) == len(second) - 1:
            head = ['gene_id'] + head
        names = [h if h else f'column_{i}' for i, h in enumerate(head)]
        body = lines[1:]
    else:
        names = ['gene_id'] + [file_stem(name) if len(second) == 2 else f'{file_stem(name)}_{i}'
                               for i in range(1, len(second))]
        body = lines
    try:
        df = pd.read_csv(io.StringIO('\n'.join(body)), sep=sep, header=None, names=names, dtype=str,
                         index_col=False, quotechar='"', engine='python' if sep == r'\s+' else 'c')
    except Exception as exc:
        raise LumosError(f'{name}: could not read the table ({exc}).')
    ids = df.iloc[:, 0].astype(str).str.strip()
    numeric, symbols = {}, None
    for col in df.columns[1:]:
        values = pd.to_numeric(df[col], errors='coerce')
        filled = values.notna().mean()
        label = str(col).strip()
        if SYMBOL_COLUMN.match(label) and filled < 0.5:
            symbols = pd.Series(df[col].astype(str).str.strip().values, index=ids.values)
            continue
        if ANNOTATION_COLUMN.match(label):
            continue
        if filled >= 0.9:
            numeric[label] = values.values
    if not numeric:
        raise LumosError(f'{name}: no numeric sample columns were found.')
    values = pd.DataFrame(numeric, index=ids.values)
    values.index.name = 'gene_id'
    keep = (values.index != '') & (values.index != 'nan') & ~values.index.str.match(NON_GENE_ROW.pattern)
    values = values[keep]
    if symbols is not None:
        symbols = symbols[keep]
    return values, symbols


def single_sample(values, sample):
    """Reduce a per-sample file (salmon, RSEM, STAR...) to its read-count column."""
    if values.shape[1] == 1:
        col = values.columns[0]
    else:
        counts = [c for c in values.columns if COUNT_COLUMN.search(str(c))]
        col = counts[0] if counts else values.sum().idxmax()
    s = values[col].copy()
    s.name = sample
    return s


def looks_single_sample(values):
    cols = [str(c).lower() for c in values.columns]
    quant = ('numreads', 'expected_count', 'tpm', 'fpkm', 'effectivelength', 'unstranded')
    return values.shape[1] == 1 or any(any(q in c for q in quant) for c in cols)


def collapse_duplicates(values, notes):
    if values.index.has_duplicates:
        n = int(values.index.duplicated().sum())
        kind = value_kind(values.values)
        values = values.groupby(level=0, sort=False).sum(min_count=1) if kind == 'counts' \
            else values.groupby(level=0, sort=False).mean()
        notes.append(f'{n} duplicated gene identifiers were merged ({"summed" if kind == "counts" else "averaged"}).')
    return values


def read_expression_files(paths, notes):
    """One expression matrix from one or more uploaded files (matrices or one file per sample)."""
    pieces, symbols = [], None
    for i, path in enumerate(paths):
        name = os.path.basename(path)
        progress(0.1 + 0.6 * i / max(1, len(paths)), f'Reading {name}')
        with open(path, 'rb') as fh:
            values, sym = read_table(fh.read(), name)
        if symbols is None and sym is not None:
            symbols = sym
        if len(paths) > 1 and looks_single_sample(values):
            gsm = re.search(r'GSM\d+', name)
            pieces.append(single_sample(values, gsm.group(0) if gsm else file_stem(name)).to_frame())
        else:
            pieces.append(values)
    return join_pieces(pieces, notes), symbols


def join_pieces(pieces, notes):
    if len(pieces) == 1:
        return pieces[0]
    seen = set()
    for p in pieces:
        p.index = p.index.astype(str)
        dup = [c for c in p.columns if c in seen]
        if dup:
            p.columns = [f'{c}.{len(seen)}' if c in seen else c for c in p.columns]
        seen.update(p.columns)
    pieces = [p[~p.index.duplicated()] for p in pieces]
    joined = pd.concat(pieces, axis=1, join='outer')
    missing = int(joined.isna().any(axis=1).sum())
    if missing:
        notes.append(f'{missing} genes were missing from some files; their missing values were set to 0.')
    return joined.fillna(0)


def read_metadata(path):
    with open(path, 'rb') as fh:
        text = _decode(fh.read(), path)
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise LumosError('The sample table needs a header line and one line per sample.')
    sep, split = _sniff(lines)
    head = split(lines[0])
    rows = [split(ln) for ln in lines[1:]]
    width = max(len(head), max(len(r) for r in rows))
    head = head + [f'column_{i}' for i in range(len(head), width)]
    head = [h if h else f'column_{i}' for i, h in enumerate(head)]
    rows = [r + [''] * (width - len(r)) for r in rows]
    return pd.DataFrame(rows, columns=head, dtype=str)


# ---------------------------------------------------------------------------- GEO

def geo_url(acc, *parts):
    acc = acc.upper()
    kind = {'GSE': 'series', 'GPL': 'platforms', 'GSM': 'samples'}[acc[:3]]
    num = acc[3:]
    stub = acc[:3] + (num[:-3] if len(num) > 3 else '') + 'nnn'
    return '/'.join([GEO_FTP, kind, stub, acc] + list(parts))


def list_geo_dir(url):
    try:
        r = session().get(url + '/', timeout=30)
    except Exception:
        return []
    if r.status_code != 200:
        return []
    names = re.findall(r'href="([^"?]+)"', r.text)
    return [n for n in names if not n.startswith(('http', '/', '#'))]


def parse_series_matrix(fh):
    """Series information, sample metadata and expression table of a GEO series matrix."""
    info = {'title': '', 'types': [], 'platforms': [], 'taxid': None, 'summary': ''}
    fields, chars, table = {}, {}, []
    in_table = False
    for raw in fh:
        line = raw.rstrip('\r\n')
        if in_table:
            if line.startswith('!series_matrix_table_end'):
                break
            table.append(line)
            continue
        if line.startswith('!series_matrix_table_begin'):
            in_table = True
            continue
        if not line.startswith('!'):
            continue
        key, _, rest = line.partition('\t')
        values = [v.strip().strip('"').strip() for v in rest.split('\t')] if rest else []
        if key == '!Series_title' and values:
            info['title'] = values[0]
        elif key == '!Series_type':
            info['types'] += values
        elif key == '!Series_platform_id':
            info['platforms'] += values
        elif key == '!Series_platform_taxid' and values and not info['taxid']:
            info['taxid'] = values[0]
        elif key == '!Series_summary':
            info['summary'] += ' '.join(values) + ' '
        elif key == '!Sample_geo_accession':
            fields['Accession'] = values
        elif key == '!Sample_title':
            fields['Title'] = values
        elif key == '!Sample_source_name_ch1':
            fields['Source'] = values
        elif key == '!Sample_platform_id':
            fields['Platform'] = values
        elif key == '!Sample_characteristics_ch1':
            for i, v in enumerate(values):
                if not v:
                    continue
                k, sep, val = v.partition(':')
                k = k.strip() if sep and 0 < len(k.strip()) <= 60 else 'characteristics'
                val = val.strip() if sep and k != 'characteristics' else v
                col = chars.setdefault(k, {})
                col[i] = f'{col[i]}; {val}' if i in col else val
    n = len(fields.get('Accession', []))
    meta = pd.DataFrame({'Accession': fields.get('Accession', [])})
    for col in ('Title', 'Source'):
        if col in fields and len(fields[col]) == n:
            meta[col] = fields[col]
    for k, col in chars.items():
        meta[k] = [col.get(i, '') for i in range(n)]
    if 'Platform' in fields and len(set(fields['Platform'])) > 1:
        meta['Platform'] = fields['Platform']
    expr = None
    if len(table) > 1:
        expr = pd.read_csv(io.StringIO('\n'.join(table)), sep='\t', index_col=0, dtype=str)
        expr.index = expr.index.astype(str).str.strip().str.strip('"')
        expr.columns = [str(c).strip().strip('"') for c in expr.columns]
        expr = expr.apply(pd.to_numeric, errors='coerce')
        expr = expr.loc[:, expr.notna().mean() > 0.5]
        if expr.shape[1] < 2 or expr.shape[0] < 10:
            expr = None
    info['taxid'] = int(info['taxid']) if info['taxid'] and str(info['taxid']).isdigit() else None
    return info, meta, expr


def geo_series(acc, notes):
    """Sample metadata (and the series-matrix table, if any) of a GEO series."""
    progress(0.02, f'Looking up {acc} on GEO')
    names = sorted(n for n in list_geo_dir(geo_url(acc, 'matrix')) if n.endswith('_series_matrix.txt.gz'))
    if not names:
        raise LumosError(f'{acc} was not found on GEO, or it has no series matrix yet.')
    infos, metas, tables = [], [], []
    for i, name in enumerate(names):
        data = download(geo_url(acc, 'matrix', name), f'Downloading {name}',
                        span=(0.03 + 0.12 * i / len(names), 0.03 + 0.12 * (i + 1) / len(names)))
        info, meta, expr = parse_series_matrix(gzip.open(io.BytesIO(data), 'rt', errors='replace'))
        infos.append(info)
        metas.append(meta)
        tables.append(expr)
    info = infos[0]
    info['platforms'] = sorted({p for i in infos for p in i['platforms']})
    info['types'] = sorted({t for i in infos for t in i['types']})
    meta = pd.concat(metas, ignore_index=True).fillna('')
    table = None
    found = [t for t in tables if t is not None]
    if found:
        table = max(found, key=lambda t: t.shape[1])
        if len(found) > 1:
            notes.append(f'{acc} has {len(names)} series matrix files (one per platform); the expression '
                         f'table with most samples was used ({table.shape[1]} samples).')
    say(f'{acc}: {info["title"]}')
    say(f'{len(meta)} samples; platform(s) {", ".join(info["platforms"])}; {"; ".join(info["types"])}')
    return info, meta, table


def platform_of(info):
    types = ' '.join(info['types']).lower()
    if 'high throughput sequencing' in types:
        return 'rnaseq'
    if 'array' in types:
        return 'microarray'
    return 'unknown'


def rank_supplementary(names):
    """Text tables in a GEO supplementary folder, most likely count matrices first."""
    scored = []
    for n in names:
        low = n.lower()
        if not TABLE_FILE.search(low) or low == 'filelist.txt':
            continue
        score = 0
        score += 4 if 'count' in low else 0
        score += 2 if 'raw' in low else 0
        score += 1 if any(k in low for k in ('matrix', 'expression', 'gene')) else 0
        score -= 2 if any(k in low for k in ('tpm', 'fpkm', 'rpkm', 'cpm', 'norm')) else 0
        score -= 3 if any(k in low for k in ('transcript', 'isoform', 'exon')) else 0
        score -= 6 if any(k in low for k in ('annot', 'metadata', 'readme', 'sample_info', 'barcode')) else 0
        scored.append((score, n))
    return [n for s, n in sorted(scored, key=lambda x: -x[0])]


def geo_raw_tar(acc, name, outdir, notes):
    """Join the per-sample files of a GEO *_RAW.tar archive into one matrix."""
    path = os.path.join(outdir, 'input', name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    download(geo_url(acc, 'suppl', name), f'Downloading {name}', span=(0.4, 0.75), dest=path,
             limit_mb=MAX_DOWNLOAD_MB, timeout=600)
    pieces = []
    try:
        with tarfile.open(path, 'r') as tar:
            members = [m for m in tar.getmembers() if m.isfile() and TABLE_FILE.search(m.name)]
            for i, m in enumerate(members):
                if i % 5 == 0:
                    progress(0.75 + 0.15 * i / max(1, len(members)), f'Reading sample files ({i + 1} of {len(members)})')
                fh = tar.extractfile(m)
                if fh is None:
                    continue
                fname = os.path.basename(m.name)
                try:
                    values, _ = read_table(fh.read(), fname)
                except LumosError as exc:
                    notes.append(str(exc))
                    continue
                gsm = re.search(r'GSM\d+', fname)
                pieces.append(single_sample(values, gsm.group(0) if gsm else file_stem(fname)).to_frame())
    finally:
        os.remove(path)
    if len(pieces) < 2:
        return None
    return join_pieces(pieces, notes)


def geo_expression(acc, info, meta, table, platform, outdir, notes):
    """Expression matrix of a GEO series: series matrix, supplementary matrix or RAW.tar."""
    if platform == 'microarray' and table is not None:
        return table, None, 'GEO series matrix'
    supp = list_geo_dir(geo_url(acc, 'suppl'))
    best = None
    for i, name in enumerate(rank_supplementary(supp)[:6]):
        try:
            data = download(geo_url(acc, 'suppl', name), f'Downloading {name}', span=(0.15, 0.6),
                            limit_mb=MAX_DOWNLOAD_MB, timeout=300)
            values, symbols = read_table(data, name)
        except LumosError as exc:
            say(f'Skipped {name}: {exc}')
            continue
        if values.shape[1] < 2:
            continue
        matched = sum(1 for c in match_samples(meta, list(values.columns))[0] if c)
        say(f'{name}: {values.shape[0]} genes x {values.shape[1]} columns, {matched} matched to samples')
        if matched >= 0.5 * len(meta):
            return values, symbols, f'GEO supplementary file {name}'
        if best is None or matched > best[0]:
            best = (matched, values, symbols, f'GEO supplementary file {name}')
    if table is not None and value_kind(table.values) == 'counts':
        return table, None, 'GEO series matrix'
    tar = next((n for n in supp if n.upper().endswith('_RAW.TAR')), None)
    if tar:
        try:
            values = geo_raw_tar(acc, tar, outdir, notes)
            if values is not None:
                return values, None, f'GEO per-sample files ({tar})'
        except LumosError as exc:
            notes.append(str(exc))
    if best is not None:
        notes.append('Fewer than half of the samples could be matched to the columns of the expression file.')
        return best[1], best[2], best[3]
    if table is not None:
        return table, None, 'GEO series matrix'
    return None, None, None


# ---------------------------------------------------------------------------- samples

def _norm(text):
    return re.sub(r'[^a-z0-9]', '', str(text).lower())


def match_samples(meta, columns):
    """For each metadata row, the expression column holding that sample (or None), and how it matched."""
    cols = [str(c) for c in columns]
    n = len(meta)
    out, how = [None] * n, [None] * n
    used = set()
    ids = [str(v) for v in meta.iloc[:, 0]] if n else []
    titles = [str(v) for v in meta['Title']] if 'Title' in meta.columns else [''] * n

    def assign(i, col, method):
        out[i], how[i] = col, method
        used.add(col)

    exact = {c: c for c in cols}
    norm_cols = {}
    for c in cols:
        norm_cols.setdefault(_norm(c), []).append(c)
    for i, sid in enumerate(ids):                                   # 1. the identifier itself
        if sid in exact and sid not in used:
            assign(i, sid, 'identifier')
    for i, sid in enumerate(ids):                                   # 2. identifier inside a column name
        if out[i] or not sid:
            continue
        hits = [c for c in cols if c not in used and re.search(rf'(^|[^A-Za-z0-9]){re.escape(sid)}([^0-9]|$)', c)]
        if len(hits) == 1:
            assign(i, hits[0], 'identifier')
    for i, title in enumerate(titles):                              # 3. the title (or its first part)
        if out[i] or not title:
            continue
        hits = [c for c in cols if c not in used and (
            c == title or _norm(c) == _norm(title) or
            (len(c) >= 3 and title.startswith(c) and not title[len(c):len(c) + 1].isalnum()))]
        if len(hits) == 1:
            assign(i, hits[0], 'title')
    for i, sid in enumerate(ids):                                   # 4. normalised identifier
        if out[i]:
            continue
        hits = [c for c in norm_cols.get(_norm(sid), []) if c not in used]
        if len(hits) == 1:
            assign(i, hits[0], 'identifier')
    if not any(out) and len(cols) == n:                             # 5. same number: by position
        for i, c in enumerate(cols):
            assign(i, c, 'position')
    return out, how


def guess_group(name):
    """Sample name without its replicate number: ctrl_1 -> ctrl."""
    g = re.sub(r'[\s_.\-]*(rep|replicate|r|s|sample)?[\s_.\-]*\d+$', '', str(name), flags=re.I)
    return g or str(name)


def sample_table(meta, columns, notes, find_id_column=False):
    """Rows of the sample table: one per sample, with the expression column that holds it.

    With `find_id_column` (uploaded sample tables) the column that names the samples best is
    moved first and becomes the sample identifier; GEO tables keep the accession first.
    """
    columns = [str(c) for c in columns]
    if meta is None or meta.empty:
        meta = pd.DataFrame({'Sample': columns, 'Name without replicate number': [guess_group(c) for c in columns]})
        matched, how = columns[:], ['identifier'] * len(columns)
    else:
        meta = meta.fillna('').astype(str)
        if columns and find_id_column:
            exact, normed = set(columns), {_norm(c) for c in columns}
            best = max(range(meta.shape[1]), key=lambda j: sum(
                1 for v in meta.iloc[:, j] if v in exact or _norm(v) in normed))
            if best:
                meta = meta[[meta.columns[best]] + [c for k, c in enumerate(meta.columns) if k != best]]
        matched, how = match_samples(meta, columns) if columns else ([None] * len(meta), [None] * len(meta))
    rows, seen = [], set()
    for i in range(len(meta)):
        sid = str(meta.iloc[i, 0]) or f'sample_{i + 1}'
        while sid in seen:
            sid += '*'
        seen.add(sid)
        rows.append({'id': sid, 'column': matched[i], 'values': [str(v) for v in meta.iloc[i]]})
    used = {m for m in matched if m}
    extra = [c for c in columns if c not in used]
    for c in extra:
        rows.append({'id': c if c not in seen else f'{c}*', 'column': c, 'extra': True,
                     'values': [c] + [''] * (meta.shape[1] - 1)})
    if columns:
        n_ok = sum(1 for m in matched if m)
        methods = sorted({h for h in how if h})
        notes.append(f'{n_ok} of {len(meta)} samples were matched to expression columns'
                     + (f' (by {", ".join(methods)})' if methods else '') + '.')
        if 'position' in methods:
            notes.append('Samples were matched to columns by their order only: check the "Expression column" '
                         'of the sample table before defining groups.')
        if extra:
            notes.append(f'{len(extra)} expression columns have no metadata; they are listed at the end of the '
                         'sample table.')
    return {'columns': list(meta.columns), 'rows': rows}


# ---------------------------------------------------------------------------- dataset command

def cmd_dataset(a):
    out = a.outdir
    notes = []
    info = {'title': '', 'platforms': [], 'taxid': None, 'types': []}
    meta, symbols, source = None, None, ''
    accession = None
    if a.geo:
        accession = a.geo.strip().upper()
        if not re.match(r'^GSE\d{1,8}$', accession):
            raise LumosError('Use a GEO series accession such as GSE117769.')
        info, meta, table = geo_series(accession, notes)
        platform = a.platform if a.platform != 'auto' else platform_of(info)
        if a.expression:
            values, symbols = read_expression_files(a.expression, notes)
            source = 'uploaded files'
        else:
            values, symbols, source = geo_expression(accession, info, meta, table, platform, out, notes)
        name = accession
    else:
        if not a.expression:
            raise LumosError('Give a GEO accession or at least one expression file.')
        values, symbols = read_expression_files(a.expression, notes)
        if a.metadata:
            meta = read_metadata(a.metadata)
        platform = a.platform if a.platform != 'auto' else 'unknown'
        source = 'uploaded files'
        name = a.name or file_stem(a.expression[0])
    progress(0.9, 'Checking the expression table')
    summary = {'kind': 'dataset', 'version': __version__, 'name': name, 'accession': accession,
               'title': info['title'], 'platform': platform, 'gpl': info['platforms'][0] if info['platforms'] else None,
               'taxid': info['taxid'], 'source': source, 'needs_expression': values is None, 'notes': notes,
               'files': {}}
    if values is None:
        notes.append('No usable expression table was found for this series. Upload its count files (one '
                     'matrix or one file per sample) together with the accession.')
        summary.update({'n_genes': 0, 'n_columns': 0, 'id_type': None, 'value_kind': None,
                        'samples': sample_table(meta, [], [])})
        write_json(os.path.join(out, 'summary.json'), summary)
        return
    values.index = values.index.astype(str)
    values = collapse_duplicates(values, notes)
    kind = value_kind(values.values)
    id_type = detect_id_type(list(values.index))
    taxid = info['taxid'] or (ensembl_taxid(list(values.index)) if id_type == 'ensembl' else None)
    samples = sample_table(meta, list(values.columns), notes, find_id_column=not accession)
    say(f'Expression table: {values.shape[0]} genes x {values.shape[1]} columns, {VALUE_KIND_TEXT[kind]}, '
        f'{ID_TYPE_TEXT[id_type]}.')
    if kind != 'counts' and platform == 'rnaseq':
        notes.append(f'The values are {VALUE_KIND_TEXT[kind]}, not raw counts: DESeq2 needs raw counts, so '
                     'the Welch t-test on log2 values will be used.')
    progress(0.93, 'Saving the expression table')
    values.to_csv(os.path.join(out, 'expression.tsv.gz'), sep='\t',
                  compression={'method': 'gzip', 'compresslevel': 3}, float_format='%.6g')
    if symbols is not None:
        symbols = symbols[~symbols.index.duplicated()]
        symbols = symbols[symbols.str.len() > 0]
        symbols.rename('symbol').rename_axis('id').to_csv(os.path.join(out, 'gene_symbols_file.tsv'), sep='\t')
    meta_out = pd.DataFrame([r['values'] for r in samples['rows']], columns=samples['columns'])
    meta_out['Expression column'] = [r['column'] or '' for r in samples['rows']]
    meta_out.to_csv(os.path.join(out, 'samples.csv'), index=False)
    summary.update({'n_genes': int(values.shape[0]), 'n_columns': int(values.shape[1]), 'id_type': id_type,
                    'value_kind': kind, 'taxid': taxid, 'samples': samples,
                    'gene_examples': [str(g) for g in values.index[:5]],
                    'files': {'expression': 'expression.tsv.gz', 'samples': 'samples.csv'}})
    write_json(os.path.join(out, 'summary.json'), summary)
    progress(1.0, 'Done')


# ---------------------------------------------------------------------------- gene symbols

def _mygene(ids, scope, species):
    found = {}
    for i in range(0, len(ids), 1000):
        batch = ids[i:i + 1000]
        progress(0.55 + 0.15 * i / max(1, len(ids)), f'Looking up gene symbols (MyGene.info, {i + len(batch)} of {len(ids)})')
        for attempt in range(3):
            try:
                r = session().post(MYGENE, data={'q': ','.join(batch), 'scopes': scope, 'fields': 'symbol',
                                                 'species': species}, timeout=60)
                r.raise_for_status()
                for hit in r.json():
                    q, sym = str(hit.get('query', '')), hit.get('symbol')
                    if not hit.get('notfound') and sym and q not in found:
                        found[q] = str(sym)
                break
            except Exception as exc:
                if attempt == 2:
                    raise LumosError(f'MyGene.info did not answer ({exc}).')
                time.sleep(2 * (attempt + 1))
    return found


def _ensembl_rest(ids):
    found = {}
    for i in range(0, len(ids), 1000):
        batch = ids[i:i + 1000]
        progress(0.55 + 0.15 * i / max(1, len(ids)), f'Looking up gene symbols (Ensembl REST, {i + len(batch)} of {len(ids)})')
        r = session().post(f'{ENSEMBL_REST}/lookup/id', json={'ids': batch}, timeout=180,
                           headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
        r.raise_for_status()
        for gene_id, rec in r.json().items():
            if rec and rec.get('display_name'):
                found[gene_id] = rec['display_name']
    return found


def _platform_table(text_lines):
    rows, header, inside = [], None, False
    for line in text_lines:
        if line.startswith('!platform_table_begin'):
            inside = True
            continue
        if line.startswith('!platform_table_end'):
            break
        if inside:
            if header is None:
                header = line.rstrip('\r\n').split('\t')
            else:
                rows.append(line.rstrip('\r\n').split('\t'))
    return header, rows


NONCODING_SYMBOL = re.compile(r'^(MIR\d|LOC\d|LINC\d|SNOR[AD]\d|SNHG\d)|-(AS|IT|OT)\d*$')


def pick_symbol(value):
    """One symbol from a multi-gene annotation ('MIR4640///DDR1' -> DDR1): the first that is
    not a microRNA, LOC/LINC, snoRNA host or antisense transcript, else the first."""
    names = [v.strip() for v in value.split('///') if v.strip() and v.strip() not in ('---', 'NA')]
    if not names:
        return '', False
    coding = [n for n in names if not NONCODING_SYMBOL.search(n)]
    return (coding or names)[0], len(set(names)) > 1


def platform_symbols(gpl):
    """Probe -> gene symbol from the GEO annotation of a microarray platform."""
    header, rows = None, []
    for url in (geo_url(gpl, 'annot', f'{gpl}.annot.gz'), geo_url(gpl, 'soft', f'{gpl}_family.soft.gz')):
        try:
            r = session().get(url, stream=True, timeout=120)
            if r.status_code != 200:
                continue
            r.raw.decode_content = True
            progress(0.6, f'Reading the {gpl} platform annotation')
            header, rows = _platform_table(gzip.open(r.raw, 'rt', errors='replace'))
            r.close()
            if header:
                break
        except Exception as exc:
            say(f'{url}: {exc}')
    if not header:
        raise LumosError(f'The annotation of platform {gpl} could not be read.')
    lower = [h.strip().lower() for h in header]
    for name in ('gene symbol', 'gene_symbol', 'symbol', 'genesymbol', 'ilmn_gene', 'gene_assignment'):
        if name in lower:
            col = lower.index(name)
            break
    else:
        raise LumosError(f'The annotation of platform {gpl} has no gene-symbol column.')
    out, ambiguous = {}, 0
    for row in rows:
        if len(row) <= col or not row[0]:
            continue
        value = row[col].strip().strip('"')
        if lower[col] == 'gene_assignment':                 # "NM_001 // SYMBOL // description /// ..."
            value = '///'.join(p.split('//')[1].strip() for p in value.split('///') if p.count('//') >= 1)
        symbol, multi = pick_symbol(value)
        if symbol:
            out[row[0].strip().strip('"')] = symbol
            ambiguous += multi
    say(f'{gpl}: {len(out)} probes with a gene symbol, {ambiguous} of them annotated to several genes.')
    return out, ambiguous


def gene_symbols(ids, dataset, ds_dir):
    """Gene symbol of each identifier (cached in the dataset folder), the source and a note."""
    id_type = dataset.get('id_type')
    if id_type == 'symbol':
        return {i: i for i in ids}, 'the identifiers are gene symbols', None
    cache_path = os.path.join(ds_dir, 'gene_symbols.tsv')
    cache, source, note = {}, None, None
    file_symbols = os.path.join(ds_dir, 'gene_symbols_file.tsv')
    if os.path.exists(file_symbols):
        cache.update(pd.read_csv(file_symbols, sep='\t', dtype=str, index_col=0)['symbol'].dropna().to_dict())
        source = 'the gene-name column of the expression file'
    if os.path.exists(cache_path):
        c = pd.read_csv(cache_path, sep='\t', dtype=str, index_col=0, keep_default_na=False)
        cache.update({k: v for k, v in c['symbol'].to_dict().items()})
        meta_path = cache_path + '.json'
        if os.path.exists(meta_path):
            meta = read_json(meta_path)
            source, note = meta.get('source', source), meta.get('note')
    missing = [i for i in ids if i not in cache]
    if missing:
        found = {}
        if dataset.get('platform') == 'microarray' and dataset.get('gpl'):
            found, ambiguous = platform_symbols(dataset['gpl'])
            source = f'the GEO annotation of {dataset["gpl"]}'
            if ambiguous:
                note = (f'{ambiguous} probes of {dataset["gpl"]} are annotated to several genes; the first '
                        'symbol that is not a microRNA, LOC/LINC or antisense transcript was used.')
        elif id_type in ('ensembl', 'entrez'):
            query = {strip_version(i) if id_type == 'ensembl' else i: i for i in missing}
            species = dataset.get('taxid') or 'all'
            try:
                hits = _mygene(list(query), 'ensembl.gene' if id_type == 'ensembl' else 'entrezgene', species)
                source = 'MyGene.info'
            except LumosError as exc:
                say(str(exc))
                hits = {}
            if id_type == 'ensembl' and len(hits) < 0.5 * len(query):
                try:
                    hits.update(_ensembl_rest([q for q in query if q not in hits]))
                    source = 'MyGene.info and Ensembl REST' if source else 'Ensembl REST'
                except Exception as exc:
                    say(f'Ensembl REST did not answer ({exc}).')
            found = {query[q]: s for q, s in hits.items() if q in query}
        else:
            return {}, None, None
        for i in missing:
            cache[i] = found.get(i, '')
        pd.Series(cache, name='symbol').rename_axis('id').to_csv(cache_path + '.tmp', sep='\t')
        os.replace(cache_path + '.tmp', cache_path)
        write_json(cache_path + '.json', {'source': source, 'note': note})
    return {i: cache.get(i, '') for i in ids}, source, note


# ---------------------------------------------------------------------------- differential expression

def pca(matrix, n_genes=PCA_GENES):
    """PCA of the samples (columns) on the most variable genes (rows)."""
    m = np.asarray(matrix, dtype=float)
    m = m[np.all(np.isfinite(m), axis=1)]
    if m.shape[0] < 3 or m.shape[1] < 3:
        return None
    var = m.var(axis=1)
    top = m[np.argsort(var)[::-1][:n_genes]]
    centred = (top - top.mean(axis=1, keepdims=True)).T
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    coords = u * s
    explained = s ** 2 / np.sum(s ** 2)
    k = min(3, coords.shape[1])
    return {'pcs': [num_list(coords[:, i]) for i in range(k)], 'variance': num_list(explained[:k]),
            'n_genes': int(top.shape[0])}


_LIMMA = None


def limma_available():
    """(available, version) of R/Bioconductor limma, checked once per process."""
    global _LIMMA
    if _LIMMA is None:
        rscript, version = shutil.which('Rscript'), None
        if rscript and os.path.exists(LIMMA_BACKEND):
            try:
                r = subprocess.run([rscript, '-e', 'suppressPackageStartupMessages(library(limma)); '
                                    'cat(as.character(packageVersion("limma")))'],
                                   capture_output=True, text=True, timeout=120)
                version = r.stdout.strip() if r.returncode == 0 else None
            except Exception:
                version = None
        _LIMMA = (bool(version), version)
    return _LIMMA


def design_frame(dataset, sample_ids):
    """Sample metadata (one row per selected sample, unique column names) for the design."""
    cols, seen = [], {}
    for c in dataset['samples']['columns']:
        seen[c] = seen.get(c, 0) + 1
        cols.append(c if seen[c] == 1 else f'{c} ({seen[c]})')
    wanted = set(sample_ids)
    rows = [r for r in dataset['samples']['rows'] if r['id'] in wanted]
    frame = pd.DataFrame([r['values'] for r in rows], columns=cols)
    frame['_sample'] = [r['id'] for r in rows]
    return frame


def run_deseq2(counts, design, formula, contrast, fdr, cpus):
    """PyDESeq2 Wald test of the condition coefficient and VST values for the plots."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    try:
        counts = validate_count_matrix(counts)
    except ValueError as exc:
        raise LumosError(f'DESeq2: {exc}')
    progress(0.25, 'DESeq2: size factors, dispersions and fold changes')
    dds = DeseqDataSet(counts=counts.T, metadata=design, design=formula, refit_cooks=True, quiet=False,
                       n_cpus=cpus)
    dds.deseq2()
    progress(0.6, 'DESeq2: Wald tests')
    st = DeseqStats(dds, contrast=contrast, alpha=fdr, cooks_filter=True, independent_filter=True, quiet=False,
                    n_cpus=cpus)
    st.summary()
    res = st.results_df
    progress(0.7, 'Variance-stabilising transformation for the plots')
    try:
        dds.vst(use_design=False)
        vst = orient_gene_sample_matrix(np.asarray(dds.layers['vst_counts']), counts.shape[0], list(counts.columns))
        viz = pd.DataFrame(vst, index=counts.index, columns=counts.columns)
        transform = 'variance-stabilising transformation (VST)'
    except Exception as exc:
        say(f'VST failed ({exc}); using log2 of the normalised counts.')
        normed = orient_gene_sample_matrix(np.asarray(dds.layers['normed_counts']), counts.shape[0],
                                           list(counts.columns))
        viz = pd.DataFrame(np.log2(normed + 1), index=counts.index, columns=counts.columns)
        transform = 'log2(normalised counts + 1)'
    table = pd.DataFrame({'Mean': res['baseMean'].values, 'Log2FC': res['log2FoldChange'].values,
                          'SE': res['lfcSE'].values, 'Stat': res['stat'].values,
                          'PValue': res['pvalue'].values, 'FDR': res['padj'].values}, index=counts.index)
    return table, viz, transform


def run_limma(values, design, formula, contrast, workdir):
    """Bioconductor limma (lmFit + eBayes, trend and robust) through r/limma_backend.R."""
    ok, version = limma_available()
    if not ok:
        raise LumosError('limma needs R (Rscript) with the Bioconductor limma package. Install them (the '
                         'Docker/Apptainer images include them) or choose the Welch t-test.')
    tmp = os.path.join(workdir, 'limma_tmp')
    os.makedirs(tmp, exist_ok=True)
    paths = {k: os.path.join(tmp, f'{k}.csv') for k in ('expression', 'metadata', 'results')}
    values.to_csv(paths['expression'], index_label='Feature')
    design.to_csv(paths['metadata'], index_label='sample')
    progress(0.4, f'limma {version}: linear models and empirical Bayes moderation')
    r = subprocess.run([shutil.which('Rscript'), LIMMA_BACKEND, paths['expression'], paths['metadata'], formula,
                        contrast[2], contrast[1], paths['results']], capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        raise LumosError('limma failed: ' + (r.stderr.strip() or r.stdout.strip())[-1500:])
    out = pd.read_csv(paths['results'], dtype={'Probe_ID': str})
    shutil.rmtree(tmp, ignore_errors=True)
    out = out.set_index('Probe_ID').reindex(values.index)
    return pd.DataFrame({'Mean': out['MeanExpression'], 'Log2FC': out['Log2FC'], 'SE': out['SE'],
                         'Stat': out['Moderated_t'], 'PValue': out['PValue'], 'FDR': out['FDR'],
                         'B': out['B_statistic']}, index=values.index), version


def run_welch(values, ref_ids, test_ids):
    from scipy import stats
    ref, test = values[ref_ids].values, values[test_ids].values
    with np.errstate(all='ignore'):
        lfc = np.nanmean(test, axis=1) - np.nanmean(ref, axis=1)
        se = np.sqrt(np.nanvar(ref, axis=1, ddof=1) / ref.shape[1] + np.nanvar(test, axis=1, ddof=1) / test.shape[1])
        t, p = stats.ttest_ind(test, ref, axis=1, equal_var=False, nan_policy='omit')
    t, p = np.asarray(t, dtype=float), np.asarray(p, dtype=float)
    return pd.DataFrame({'Mean': np.nanmean(values.values, axis=1), 'Log2FC': lfc, 'SE': se, 'Stat': t,
                         'PValue': p, 'FDR': benjamini_hochberg(p)}, index=values.index)


def status_of(lfc, fdr, fdr_cut, lfc_cut):
    lfc, fdr = np.asarray(lfc, dtype=float), np.asarray(fdr, dtype=float)
    sig = np.isfinite(fdr) & (fdr < fdr_cut)
    out = np.full(lfc.shape, 'Not significant', dtype=object)
    out[sig & (lfc >= lfc_cut)] = 'Up-regulated'
    out[sig & (lfc <= -lfc_cut)] = 'Down-regulated'
    return out


def runtime_info(extra=()):
    info = {'python': sys.version.split()[0], 'platform': platform.platform(),
            'packages': dict(package_versions(['flask', 'numpy', 'pandas', 'scipy', 'pydeseq2', 'statsmodels',
                                               'networkx', 'plotly', 'requests', *extra]))}
    return info


def cmd_de(a):
    out, ds_dir = a.outdir, a.dataset
    dataset = read_json(os.path.join(ds_dir, 'summary.json'))
    groups = read_json(a.groups)
    if a.reference not in groups or a.test not in groups:
        raise LumosError('The reference and test groups must both be defined.')
    if a.reference == a.test:
        raise LumosError('Choose two different groups.')
    col_of = {r['id']: r['column'] for r in dataset['samples']['rows']}
    ref_ids = list(dict.fromkeys(s for s in groups[a.reference] if col_of.get(s)))
    test_ids = list(dict.fromkeys(s for s in groups[a.test] if col_of.get(s)))
    if set(ref_ids) & set(test_ids):
        raise LumosError('A sample cannot belong to both groups.')
    if len(ref_ids) < 2 or len(test_ids) < 2:
        raise LumosError('Each group needs at least two samples with expression values (replicates).')
    label = (a.label or f'{a.test} vs {a.reference}').strip()[:80]
    notes = []
    progress(0.03, 'Reading the expression table')
    expr = pd.read_csv(os.path.join(ds_dir, 'expression.tsv.gz'), sep='\t', index_col=0)
    expr.index = expr.index.astype(str)
    sub = expr[[col_of[s] for s in ref_ids + test_ids]].apply(pd.to_numeric, errors='coerce')
    sub.columns = ref_ids + test_ids
    kind = value_kind(sub.values)
    lib = sub.fillna(0).sum(axis=0) if kind == 'counts' else None
    say(f'{label}: {a.test} ({len(test_ids)} samples) vs {a.reference} ({len(ref_ids)} samples); '
        f'{sub.shape[0]} genes, {VALUE_KIND_TEXT[kind]}.')

    # ---- experimental design (lumos_core.build_design_metadata validates rank and pairs)
    try:
        design, formula, contrast, name_map = build_design_metadata(
            sample_ids=ref_ids + test_ids, conditions=['reference'] * len(ref_ids) + ['test'] * len(test_ids),
            source_metadata=design_frame(dataset, ref_ids + test_ids), sample_key='_sample',
            categorical_covariates=a.covariate, continuous_covariates=a.continuous,
            batch_column=a.batch or None, pair_column=a.pair or None)
    except ValueError as exc:
        raise LumosError(f'Experimental design: {exc}')
    has_terms = formula != '~condition'
    say(f'Model: {formula} (condition: test = {a.test}, reference = {a.reference})')

    # ---- statistical engine
    has_deseq2 = importlib.util.find_spec('pydeseq2') is not None
    engine = a.engine
    if engine == 'auto':
        if kind == 'counts':
            if not has_deseq2:
                raise LumosError('Raw counts were detected but PyDESeq2 is not installed (pip install pydeseq2).')
            engine = 'deseq2'
        elif limma_available()[0]:
            engine = 'limma'
        else:
            engine = 'welch'
            notes.append('R/limma is not installed, so the Welch t-test was used as an exploratory fallback '
                         '(it cannot adjust for design terms). The Docker and Apptainer images include limma.')
    if engine == 'deseq2':
        if not has_deseq2:
            raise LumosError('PyDESeq2 is not installed (pip install pydeseq2).')
        if kind != 'counts':
            raise LumosError(f'DESeq2 needs raw integer counts, but these are {VALUE_KIND_TEXT[kind]}. '
                             'Use limma (or the Welch t-test).')
    if engine == 'welch' and has_terms:
        notes.append('The Welch t-test ignores the design terms (batch, covariates, pairs); use DESeq2 or limma '
                     'to adjust for them.')

    # ---- gene symbols are needed before testing when the analysis unit is the gene symbol
    symbols, symbol_source, symbol_note = {}, None, None
    to_symbols = a.unit == 'symbols' and dataset.get('id_type') != 'symbol'
    if to_symbols:
        progress(0.06, 'Gene symbols (analysis unit)')
        symbols, symbol_source, symbol_note = gene_symbols(list(sub.index), dataset, ds_dir)
        if not any(symbols.values()):
            raise LumosError('No gene symbols could be retrieved; analyse the features as provided instead.')
        if symbol_note:
            notes.append(symbol_note)

    # ---- filtering, transformation and test
    progress(0.1, 'Filtering low-expression genes')
    min_n = min(len(ref_ids), len(test_ids))
    unit_text = 'features as provided'
    if kind == 'counts':
        if sub.isna().any().any():
            raise LumosError('The count table has missing values for some selected samples; missing counts are '
                             'not replaced by zero automatically.')
        if to_symbols:
            sub = sub.groupby([symbols.get(i) or i for i in sub.index], sort=False).sum()
            unit_text = 'gene symbols (counts of features with the same symbol summed)'
        ok = (sub >= a.min_count).sum(axis=1) >= min_n
        rule = f'count >= {a.min_count:g} in at least {min_n} samples'
    else:
        before = len(sub)
        sub = sub.dropna(axis=0, how='any')
        if len(sub) < before:
            notes.append(f'{before - len(sub)} features with missing values were removed.')
        ok = sub.var(axis=1) > 0
        rule = 'complete, non-constant values'
    removed = int((~ok).sum())
    sub = sub[ok]
    if sub.empty:
        raise LumosError('The filter removed every gene; check the minimum count and the data scale.')
    say(f'Filter ({rule}): {sub.shape[0]} genes kept, {removed} removed.')

    if engine == 'deseq2':
        table, viz, transform = run_deseq2(sub, design, formula, contrast, a.fdr, a.cpus)
        method_text = 'DESeq2 (PyDESeq2) Wald test'
        limma_version = None
    else:
        progress(0.2, 'Transforming the values')
        if kind == 'counts' and a.transform != 'none':
            values, transform = np.log2(sub / sub.sum(axis=0) * 1e6 + 1), 'log2(CPM + 1)'
        elif kind == 'continuous' and a.transform != 'none':
            values, transform = np.log2(sub.clip(lower=0) + 1), 'log2(x + 1)'
        elif kind == 'log' and a.transform == 'log2':
            values, transform = np.log2(sub.clip(lower=0) + 1), 'log2(x + 1)'
        else:
            values, transform = sub.astype(float), 'none (values used as given)'
        normalization = a.normalization
        if normalization == 'auto':
            normalization = 'quantile' if dataset.get('platform') == 'microarray' else 'none'
        if normalization == 'quantile':
            values = pd.DataFrame(quantile_normalize(values.values), index=values.index, columns=values.columns)
            transform += ', quantile normalisation'
        if to_symbols and kind != 'counts':
            mapped = [symbols.get(i) or '' for i in values.index]
            values = values[[m != '' for m in mapped]]
            values = values.groupby([m for m in mapped if m], sort=False).median()
            unit_text = 'gene symbols (median of the features with the same symbol)'
            values = values[values.var(axis=1) > 0]
        limma_version = None
        if engine == 'limma':
            table, limma_version = run_limma(values, design, formula, contrast, out)
            method_text = 'limma empirical Bayes (trend, robust)'
        else:
            progress(0.3, 'Welch t-tests')
            table = run_welch(values, ref_ids, test_ids)
            method_text = 'Welch t-test (exploratory), Benjamini-Hochberg FDR'
        viz = values

    # ---- gene symbols as labels
    ids = list(table.index)
    if to_symbols:
        named = set(v for v in symbols.values() if v)
        labels_ = {i: (i if i in named else '') for i in ids}
        n_mapped = sum(1 for i in ids if labels_[i])
        symbols, id_type = labels_, 'symbol'
    else:
        id_type = dataset.get('id_type')
        if a.symbols:
            progress(0.75, 'Gene symbols')
            try:
                symbols, symbol_source, symbol_note = gene_symbols(ids, dataset, ds_dir)
                if symbol_note:
                    notes.append(symbol_note)
            except LumosError as exc:
                notes.append(f'Gene symbols could not be retrieved: {exc}')
            except Exception as exc:
                notes.append(f'Gene symbols could not be retrieved ({exc}).')
        n_mapped = sum(1 for i in ids if symbols.get(i))
        if symbol_source and id_type != 'symbol':
            say(f'Gene symbols: {n_mapped} of {len(ids)} identifiers ({symbol_source}).')
    table.insert(0, 'Symbol', [symbols.get(i) or '' for i in ids])
    table.insert(0, 'Gene_ID', ids)
    table['Status'] = status_of(table['Log2FC'], table['FDR'], a.fdr, a.lfc)
    table = table.sort_values('PValue', na_position='last')
    n_nan = int(table['FDR'].isna().sum())
    if n_nan and engine == 'deseq2':
        notes.append(f'{n_nan} genes have no adjusted p-value (DESeq2 independent filtering or count outliers); '
                     'they are not plotted in the volcano plot.')

    progress(0.85, 'Writing the results')
    table.to_csv(os.path.join(out, 'de_results.csv'), index=False, float_format='%.6g')
    viz.loc[table['Gene_ID']].to_csv(os.path.join(out, 'normalized_expression.tsv.gz'), sep='\t',
                                     float_format='%.5g', compression={'method': 'gzip', 'compresslevel': 3})
    pool = table[table['FDR'].notna()].head(HEATMAP_POOL)['Gene_ID']
    pc = pca(viz.values)
    rows = {r['id']: r for r in dataset['samples']['rows']}
    title_col = dataset['samples']['columns'].index('Title') if 'Title' in dataset['samples']['columns'] else None
    sample_labels = [rows[s]['values'][title_col] if title_col is not None else s for s in ref_ids + test_ids]
    up = int((table['Status'] == 'Up-regulated').sum())
    down = int((table['Status'] == 'Down-regulated').sum())
    say(f'{up} up-regulated and {down} down-regulated genes (FDR < {a.fdr:g}, |log2FC| >= {a.lfc:g}).')
    design_info = {'formula': formula,
                   'terms': {k: v for k, v in name_map.items()}, 'batch': a.batch or None, 'pair': a.pair or None,
                   'covariates': list(a.covariate), 'continuous': list(a.continuous)}
    files = {'results': 'de_results.csv', 'normalized': 'normalized_expression.tsv.gz',
             'manifest': 'manifest.json', 'zip': f'lumos_{dataset["name"]}_{safe_name(label)}_DDEA.zip'}
    summary = {
        'kind': 'de', 'version': __version__, 'label': label, 'name': dataset['name'],
        'title': dataset.get('title', ''), 'dataset_job': os.path.basename(os.path.normpath(ds_dir)),
        'comparison': {'reference': a.reference, 'test': a.test, 'reference_samples': ref_ids,
                       'test_samples': test_ids, 'sample_labels': sample_labels},
        'parameters': {'engine': engine, 'method': engine, 'fdr': a.fdr, 'lfc': a.lfc, 'min_count': a.min_count,
                       'transform': a.transform, 'normalization': a.normalization, 'unit': a.unit,
                       'symbols': bool(a.symbols)},
        'design': design_info, 'method_text': method_text, 'transform_text': transform, 'unit_text': unit_text,
        'value_kind': kind, 'mean_scale': 'log' if engine == 'deseq2' else 'linear',
        'platform': dataset.get('platform'), 'id_type': id_type, 'taxid': dataset.get('taxid'),
        'symbols': {'source': symbol_source, 'mapped': n_mapped, 'total': len(ids),
                    'available': bool(n_mapped) and id_type != 'symbol'},
        'counts': {'tested': int(len(table)), 'filtered_out': removed, 'up': up, 'down': down,
                   'no_fdr': n_nan, 'rule': rule},
        'genes': {'id': [str(i) for i in table['Gene_ID']], 'symbol': list(table['Symbol']),
                  'lfc': num_list(table['Log2FC']), 'fdr': num_list(table['FDR']), 'p': num_list(table['PValue']),
                  'mean': num_list(table['Mean'])},
        'pca': pc, 'library_sizes': num_list(lib.values, 6) if lib is not None else None,
        'heatmap': {'genes': list(pool), 'values': [num_list(viz.loc[g].values, 3) for g in pool]},
        'notes': notes, 'files': files,
    }
    manifest = {
        'lumos_version': __version__, 'generated_utc': utc_now(), 'analysis': label,
        'dataset': {k: dataset.get(k) for k in ('name', 'accession', 'title', 'source', 'platform', 'gpl', 'taxid',
                                                'id_type', 'value_kind')},
        'comparison': {'reference': a.reference, 'test': a.test, 'reference_samples': ref_ids,
                       'test_samples': test_ids, 'expression_columns': {s: col_of[s] for s in ref_ids + test_ids}},
        'design': {'formula': formula, 'contrast': ['condition', a.test, a.reference],
                   'source_to_model_column': name_map,
                   'metadata': design.assign(condition=design['condition'].map(
                       {'reference': a.reference, 'test': a.test})).reset_index().to_dict('records')},
        'method': method_text, 'analysis_unit': unit_text, 'filter': rule, 'values_for_plots': transform,
        'thresholds': {'FDR': a.fdr, 'absolute_log2FC': a.lfc},
        'gene_symbols': {'source': symbol_source, 'mapped': n_mapped, 'total': len(ids)},
        'results': {'tested': int(len(table)), 'up': up, 'down': down, 'without_FDR': n_nan},
        'runtime': runtime_info(), 'R_limma': limma_version, 'notes': notes,
    }
    write_json(os.path.join(out, 'summary.json'), summary)
    with open(os.path.join(out, 'manifest.json'), 'w') as fh:
        json.dump(manifest, fh, indent=2, default=str)
    with open(os.path.join(out, 'groups.json'), 'w') as fh:
        json.dump(groups, fh, indent=1)
    make_zip(out, files['zip'], ['de_results.csv', 'normalized_expression.tsv.gz', 'manifest.json', 'groups.json',
                                 'summary.json', 'run.log'])
    progress(1.0, 'Done')


def utc_now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def safe_name(text):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(text)).strip('_')[:60] or 'analysis'


def make_zip(folder, name, files):
    with zipfile.ZipFile(os.path.join(folder, name), 'w', zipfile.ZIP_DEFLATED) as z:
        for f in files:
            path = os.path.join(folder, f)
            if os.path.exists(path):
                z.write(path, f)


# ---------------------------------------------------------------------------- integration of analyses

INTEGRATION_TEXT = {
    'contrasts': 'contrasts of one dataset aligned by gene (no pooled effect)',
    'meta': 'random-effects meta-analysis of independent datasets',
    'tissues': 'independent datasets (tissues) aligned by gene (no pooled effect)',
}


def load_analysis(path):
    s = read_json(os.path.join(path, 'summary.json'))
    if s.get('kind') != 'de':
        raise LumosError(f'{path} is not a differential-expression analysis.')
    res = pd.read_csv(os.path.join(path, 'de_results.csv'), dtype={'Gene_ID': str, 'Symbol': str})
    res['Symbol'] = res['Symbol'].fillna('').astype(str).str.strip()
    res['Symbol'] = res['Symbol'].where(res['Symbol'] != '', res['Gene_ID']).str.upper()
    return s, res


def combine_meta_effects(studies, min_studies=2, estimator='Paule-Mandel'):
    """Gene-level random-effects meta-analysis of log2 fold changes (statsmodels combine_effects)."""
    from scipy import stats
    from statsmodels.stats.meta_analysis import combine_effects
    frames = []
    for sid, res in studies.items():
        f = res[['Symbol', 'Log2FC', 'SE', 'PValue']].copy()
        f = f.replace([np.inf, -np.inf], np.nan).dropna(subset=['Log2FC', 'SE'])
        f = f[f['SE'] > 0].sort_values('PValue').drop_duplicates('Symbol')     # best feature per gene
        f['Study'] = sid
        frames.append(f)
    allr = pd.concat(frames, ignore_index=True)
    groups = [(g, x) for g, x in allr.groupby('Symbol', sort=False) if len(x) >= min_studies]
    method_re = 'iterated' if estimator == 'Paule-Mandel' else 'chi2'
    rows, failures = [], 0
    for k, (gene, g) in enumerate(groups):
        if k % 2000 == 0:
            progress(0.2 + 0.6 * k / max(1, len(groups)), f'Meta-analysis: {k} of {len(groups)} genes')
        eff = g['Log2FC'].to_numpy(dtype=float)
        var = np.square(g['SE'].to_numpy(dtype=float))
        try:
            res = combine_effects(eff, var, method_re=method_re, row_names=list(g['Study']), use_t=False, alpha=0.05)
            sf = res.summary_frame()
            re_row = sf.loc['random effect'] if 'random effect' in sf.index else sf.iloc[-3]
            w = 1.0 / var
            fixed = np.sum(w * eff) / np.sum(w)
            q = float(np.sum(w * (eff - fixed) ** 2))
            dfq = len(eff) - 1
            est, se_ = float(re_row['eff']), float(re_row['sd_eff'])
            z = est / se_ if se_ > 0 else np.nan
            rows.append({'Symbol': gene, 'Meta_Log2FC': est, 'Meta_SE': se_, 'CI95_Low': float(re_row['ci_low']),
                         'CI95_High': float(re_row['ci_upp']),
                         'PValue': float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                         'I2_percent': max(0.0, (q - dfq) / q) * 100.0 if q > 0 and dfq > 0 else 0.0,
                         'Tau2': float(getattr(res, 'tau2', np.nan)), 'Q': q, 'Q_df': dfq,
                         'Q_PValue': float(stats.chi2.sf(q, dfq)) if dfq > 0 else np.nan,
                         'N_Studies': int(len(g)), 'Studies': ', '.join(g['Study']), 'Estimator': estimator})
        except Exception:
            failures += 1
    meta = pd.DataFrame(rows)
    if meta.empty:
        return meta, failures
    meta['FDR'] = benjamini_hochberg(meta['PValue'].to_numpy())
    return meta.sort_values(['FDR', 'PValue'], na_position='last').reset_index(drop=True), failures


def meta_expression_heatmap(analyses, profiles, genes):
    """Row z-scores (within each study) of the top meta-analysis genes across all studies' samples."""
    blocks, samples, groups, studies = [], [], [], []
    for (s, res, path), p in zip(analyses, profiles):
        norm = pd.read_csv(os.path.join(path, 'normalized_expression.tsv.gz'), sep='\t', index_col=0)
        norm.index = norm.index.astype(str)
        key = dict(zip(res['Gene_ID'].astype(str), res['Symbol']))
        norm = norm[norm.index.isin(key)]
        norm = norm.groupby(norm.index.map(key)).mean()
        norm = norm.reindex([g for g in genes if g in norm.index])
        cmp_ = s['comparison']
        cols = [c for c in cmp_['reference_samples'] + cmp_['test_samples'] if c in norm.columns]
        m = norm[cols].to_numpy(dtype=float)
        z = (m - np.nanmean(m, axis=1, keepdims=True)) / (np.nanstd(m, axis=1, keepdims=True) + 1e-9)
        blocks.append(pd.DataFrame(z, index=norm.index, columns=[f'{p["id"]}|{c}' for c in cols]))
        samples += cols
        groups += [cmp_['reference'] if c in cmp_['reference_samples'] else cmp_['test'] for c in cols]
        studies += [p['label']] * len(cols)
    joined = pd.concat(blocks, axis=1, join='inner')
    if joined.empty:
        return None
    return {'genes': list(joined.index), 'samples': samples, 'groups': groups, 'studies': studies,
            'z': [num_list(r, 3) for r in joined.to_numpy()]}


def cmd_integrate(a):
    out = a.outdir
    if len(a.analyses) < 2:
        raise LumosError('Choose at least two analyses.')
    progress(0.05, 'Reading the analyses')
    analyses = []
    for path in a.analyses:
        s, res = load_analysis(path)
        analyses.append((s, res, path))
    datasets = [(s['dataset_job'], s['name']) for s, _, _ in analyses]
    if a.mode == 'contrasts' and len(set(datasets)) > 1:
        raise LumosError('Contrasts must come from the same dataset; use the meta-analysis or the multi-tissue '
                         'comparison for different datasets.')
    if a.mode in ('meta', 'tissues') and len(set(datasets)) < len(datasets):
        raise LumosError('Each analysis must come from a different, independent dataset (two of them share one). '
                         'Compare contrasts of one dataset with the "contrasts" integration.')
    labels = []
    for s, _, _ in analyses:
        base = s['label'] if a.mode == 'contrasts' else s['name']
        if a.mode != 'contrasts' and sum(1 for x, _, _ in analyses if x['name'] == s['name']) > 1:
            base = f'{s["name"]} {s["label"]}'
        lab, k = base, 2
        while lab in labels:
            lab, k = f'{base} ({k})', k + 1
        labels.append(lab)
    ids, used = [], set()
    for lab in labels:
        sid, k = sanitize_design_name(lab), 2
        while sid in used:
            sid, k = f'{sanitize_design_name(lab)}_{k}', k + 1
        used.add(sid)
        ids.append(sid)
    profiles = [{'id': i, 'label': lab, 'job': os.path.basename(os.path.normpath(p)), 'dataset': s['name'],
                 'comparison': f'{s["comparison"]["test"]} vs {s["comparison"]["reference"]}',
                 'method': s.get('method_text'), 'tested': s['counts']['tested']}
                for i, lab, (s, _, p) in zip(ids, labels, analyses)]
    notes, summary_table, extra = [], {}, {}
    progress(0.15, f'Integrating {len(analyses)} analyses')
    if a.mode == 'meta':
        studies = {i: res for i, (_, res, _) in zip(ids, analyses)}
        missing_se = [lab for lab, (_, res, _) in zip(labels, analyses) if res['SE'].isna().all()]
        if missing_se:
            raise LumosError('These analyses have no standard errors (re-run them with this version): '
                             + ', '.join(missing_se))
        meta, failures = combine_meta_effects(studies, a.min_studies, a.estimator)
        if meta.empty:
            raise LumosError('No gene had enough information across the selected datasets.')
        if failures:
            notes.append(f'The meta-analysis failed for {failures} genes (no substitute estimator was used).')
        meta['Status'] = status_of(meta['Meta_Log2FC'], meta['FDR'], a.fdr, a.lfc)
        for i, (_, res, _) in zip(ids, analyses):
            best = res.sort_values('PValue').drop_duplicates('Symbol').set_index('Symbol')
            for col, src in (('Log2FC', 'Log2FC'), ('SE', 'SE'), ('FDR', 'FDR')):
                meta[f'{col}_{i}'] = meta['Symbol'].map(best[src])
        table = meta
        summary_table = {
            'symbol': list(meta['Symbol']), 'lfc': num_list(meta['Meta_Log2FC']), 'se': num_list(meta['Meta_SE']),
            'ci_low': num_list(meta['CI95_Low']), 'ci_high': num_list(meta['CI95_High']),
            'p': num_list(meta['PValue']), 'fdr': num_list(meta['FDR']), 'i2': num_list(meta['I2_percent'], 3),
            'tau2': num_list(meta['Tau2']), 'q_p': num_list(meta['Q_PValue']), 'n': [int(x) for x in meta['N_Studies']],
            'profiles': {i: {'lfc': num_list(meta[f'Log2FC_{i}']), 'se': num_list(meta[f'SE_{i}']),
                             'fdr': num_list(meta[f'FDR_{i}'])} for i in ids}}
        top = list(meta.loc[meta['Status'] != 'Not significant', 'Symbol'].head(50))
        if top:
            progress(0.85, 'Expression heatmap of the top genes')
            try:
                extra['heatmap'] = meta_expression_heatmap(analyses, profiles, top)
            except Exception as exc:
                notes.append(f'The combined expression heatmap could not be built ({exc}).')
        i2 = meta['I2_percent'].dropna()
        extra['heterogeneity'] = {'median_i2': fnum(i2.median()), 'mean_i2': fnum(i2.mean()),
                                  'high_i2': int((i2 > 50).sum())}
        counts = meta['Status'].value_counts().to_dict()
    else:
        if a.mode == 'contrasts':
            try:
                table, safe_to_label = align_named_contrasts(
                    {lab: res for lab, (_, res, _) in zip(labels, analyses)}, a.fdr, a.lfc)
            except ValueError as exc:
                raise LumosError(str(exc))
            ids = list(safe_to_label)
            for p, i in zip(profiles, ids):
                p['id'] = i
            pattern_col, nsig_col = 'Cross_Contrast_Pattern', 'N_Significant_Contrasts'
        else:
            merged = pd.DataFrame()
            for i, (_, res, _) in zip(ids, analyses):
                c = collapse_symbol_results(res)
                c['Status'] = classify_differential_expression(c, 'Log2FC', 'FDR', a.fdr, a.lfc)
                c = c.rename(columns={'Log2FC': f'Log2FC_{i}', 'PValue': f'PValue_{i}', 'FDR': f'FDR_{i}',
                                      'Status': f'Status_{i}'})
                merged = c if merged.empty else pd.merge(merged, c, on='Symbol', how='outer')
            table = annotate_cross_dataset_patterns(merged, ids)
            for i, lab in zip(ids, labels):
                table['Cross_Tissue_Pattern'] = table['Cross_Tissue_Pattern'].replace(
                    {f'Specific {i} (Up)': f'Specific {lab} (Up)', f'Specific {i} (Down)': f'Specific {lab} (Down)'})
            table['Analysis_Scope'] = 'multi_dataset_alignment'
            pattern_col, nsig_col = 'Cross_Tissue_Pattern', 'N_Significant_Datasets'
        table = table.sort_values([nsig_col, 'Symbol'], ascending=[False, True]).reset_index(drop=True)
        code = {'Up-regulated': 1, 'Down-regulated': -1}
        summary_table = {
            'symbol': list(table['Symbol']), 'pattern': list(table[pattern_col].fillna('')),
            'n_sig': [int(x) for x in table[nsig_col].fillna(0)],
            'profiles': {i: {'lfc': num_list(table[f'Log2FC_{i}']), 'fdr': num_list(table[f'FDR_{i}']),
                             'status': [code.get(v, 0) for v in table[f'Status_{i}']]} for i in ids}}
        counts = table[pattern_col].value_counts().to_dict()
    table.to_csv(os.path.join(out, 'integration.csv'), index=False, float_format='%.6g')
    say(f'{INTEGRATION_TEXT[a.mode].capitalize()}: {len(table)} genes.')
    files = {'table': 'integration.csv', 'manifest': 'manifest.json', 'zip': f'lumos_integration_{a.mode}.zip'}
    summary = {'kind': 'integrate', 'version': __version__, 'mode': a.mode, 'mode_text': INTEGRATION_TEXT[a.mode],
               'profiles': profiles, 'thresholds': {'fdr': a.fdr, 'lfc': a.lfc},
               'estimator': a.estimator if a.mode == 'meta' else None,
               'min_studies': a.min_studies if a.mode == 'meta' else None,
               'taxid': next((s.get('taxid') for s, _, _ in analyses if s.get('taxid')), None),
               'n_genes': int(len(table)), 'pattern_counts': {str(k): int(v) for k, v in counts.items()},
               'table': summary_table, 'notes': notes, 'files': files, **extra}
    manifest = {'lumos_version': __version__, 'generated_utc': utc_now(), 'integration': a.mode,
                'description': INTEGRATION_TEXT[a.mode],
                'analyses': [{**p, 'manifest': read_json_safe(os.path.join(path, 'manifest.json'))}
                             for p, (_, _, path) in zip(profiles, analyses)],
                'thresholds': {'FDR': a.fdr, 'absolute_log2FC': a.lfc},
                'rules': {'contrasts': 'outer join by gene symbol; no pooled effect (contrasts may share samples)',
                          'tissues': 'outer join by gene symbol; strict conservation = significant with the same '
                                     'direction in every dataset; no pooled effect',
                          'meta': f'random-effects ({a.estimator}) combination of gene-level log2 fold changes and '
                                  f'standard errors; genes in at least {a.min_studies} studies; BH FDR'}[a.mode],
                'runtime': runtime_info()}
    write_json(os.path.join(out, 'summary.json'), summary)
    with open(os.path.join(out, 'manifest.json'), 'w') as fh:
        json.dump(manifest, fh, indent=2, default=str)
    make_zip(out, files['zip'], ['integration.csv', 'manifest.json', 'summary.json', 'run.log'])
    progress(1.0, 'Done')


def read_json_safe(path):
    try:
        return read_json(path)
    except Exception:
        return None


# ---------------------------------------------------------------------------- enrichment (APP)

def tf_symbol(term):
    """TF name of an Enrichr regulator term: 'TP53 human' -> TP53, 'RELA 24523406 ChIP-Seq...' -> RELA."""
    return re.split(r'[\s_]', str(term).strip())[0].upper()


def _cloudscraper():
    """Fallback client when Enrichr answers 403 to cloud hosts."""
    try:
        import cloudscraper
    except ImportError:
        return None
    return cloudscraper.create_scraper()


def library_sets(name):
    """Gene sets of an Enrichr library (cached on disk for 30 days)."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, f'enrichr_{safe_name(name)}.txt')
    if not os.path.exists(path) or time.time() - os.path.getmtime(path) > 30 * 86400:
        say(f'Downloading the {name} library from Enrichr')
        r = session().get(f'{ENRICHR}/geneSetLibrary', params={'mode': 'text', 'libraryName': name}, timeout=300)
        if r.status_code == 403:
            s = _cloudscraper()
            if s:
                r = s.get(f'{ENRICHR}/geneSetLibrary', params={'mode': 'text', 'libraryName': name}, timeout=300)
        r.raise_for_status()
        with open(path + '.tmp', 'w') as fh:
            fh.write(r.text)
        os.replace(path + '.tmp', path)
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            genes = {g.split(',')[0].strip().upper() for g in parts[2:] if g.strip()}
            if genes:
                sets[parts[0]] = genes
    if not sets:
        raise LumosError(f'Enrichr returned an empty library for {name}.')
    return sets


class ORA:
    """Over-representation analysis against an explicit background: Speedrichr first, then a local
    one-sided hypergeometric test on the same Enrichr library (Benjamini-Hochberg within library)."""

    def __init__(self, background):
        self.background = sorted({str(g).strip().upper() for g in background if str(g).strip()})
        self.speed = True
        self.bg_id = None
        self.list_ids = {}
        self.diagnostics = []

    def _speed_post(self, path, **kw):
        r = session().post(f'{SPEEDRICHR}/{path}', timeout=(10, 60), **kw)
        if r.status_code == 403:
            s = _cloudscraper()
            if s:
                r = s.post(f'{SPEEDRICHR}/{path}', timeout=(10, 60), **kw)
        r.raise_for_status()
        return r.json()

    def _speedrichr(self, query, library):
        if self.bg_id is None:
            self.bg_id = self._speed_post('addbackground', data={'background': '\n'.join(self.background)})['backgroundid']
        key = tuple(query)
        if key not in self.list_ids:
            self.list_ids[key] = self._speed_post('addList', files={
                'list': (None, '\n'.join(query)), 'description': (None, 'Lumos Networks')})['userListId']
        data = self._speed_post('backgroundenrich', data={'userListId': self.list_ids[key],
                                                          'backgroundid': self.bg_id, 'backgroundType': library})
        rows = data.get(library, [])
        return pd.DataFrame([{'Term': r[1], 'PValue': r[2], 'OddsRatio': r[3], 'CombinedScore': r[4],
                              'Genes': [str(g).upper() for g in r[5]], 'FDR': r[6]} for r in rows])

    def _local(self, query, library):
        from scipy import stats
        sets = library_sets(library)
        universe = set(self.background) & set().union(*sets.values())
        q = set(query) & universe
        if not q or len(universe) < 2:
            return pd.DataFrame()
        n_pop, n_q, rows = len(universe), len(q), []
        for term, genes in sets.items():
            tb = universe & genes
            if not tb:
                continue
            ov = sorted(q & tb)
            k, m = len(ov), len(tb)
            p = float(stats.hypergeom.sf(k - 1, n_pop, m, n_q)) if k else 1.0
            odds = (k + 0.5) * (n_pop - m - n_q + k + 0.5) / ((n_q - k + 0.5) * (m - k + 0.5))
            rows.append({'Term': term, 'PValue': p, 'OddsRatio': odds, 'CombinedScore': np.nan, 'Genes': ov})
        df = pd.DataFrame(rows)
        if df.empty:
            return df
        df['FDR'] = benjamini_hochberg(df['PValue'].to_numpy())
        return df[df['Genes'].map(bool)]

    def run(self, genes, library):
        query, background, excluded = prepare_ora_inputs(genes, self.background)
        df, engine = pd.DataFrame(), None
        if self.speed:
            try:
                df, engine = self._speedrichr(query, library), 'Speedrichr (explicit background)'
            except Exception as exc:
                self.speed = False
                self.diagnostics.append(f'Speedrichr unavailable ({type(exc).__name__}: {exc}); the local '
                                        'hypergeometric test on the Enrichr libraries was used.')
                say(self.diagnostics[-1])
        if engine is None:
            df = self._local(query, library)
            engine = 'local hypergeometric test (Enrichr library, explicit background)'
        if not df.empty:
            df['Overlap'] = df['Genes'].map(len)
        return df, engine, {'query': len(query), 'background': len(background), 'excluded': len(excluded)}


def source_evidence(src_dir, focus, fdr, lfc):
    """Evidence table (Symbol, Effect, FDR), measured-gene universe and a description of a gene source."""
    s = read_json(os.path.join(src_dir, 'summary.json'))
    if s['kind'] == 'de':
        _, res = load_analysis(src_dir)
        c = collapse_symbol_results(res)
        ev = c.rename(columns={'Log2FC': 'Effect'})[['Symbol', 'Effect', 'FDR']]
        return ev, list(c['Symbol']), {'kind': 'de', 'label': s['label'], 'dataset': s['name'],
                                       'comparison': s['comparison']}, s
    if s['kind'] != 'integrate':
        raise LumosError('The gene source must be a differential-expression analysis or an integration.')
    table = pd.read_csv(os.path.join(src_dir, 'integration.csv'))
    table['Symbol'] = table['Symbol'].astype(str).str.upper()
    ids = [p['id'] for p in s['profiles']]
    names = {p['id']: p['label'] for p in s['profiles']}
    if s['mode'] == 'meta':
        ev = table.rename(columns={'Meta_Log2FC': 'Effect'})[['Symbol', 'Effect', 'FDR']].dropna()
        return ev, list(ev['Symbol']), {'kind': 'meta', 'label': 'meta-analysis effect'}, s
    if focus.startswith('profile:'):
        pid = focus.split(':', 1)[1]
        if pid not in ids:
            raise LumosError(f'Unknown profile {pid}.')
        ev = table[['Symbol', f'Log2FC_{pid}', f'FDR_{pid}']].rename(
            columns={f'Log2FC_{pid}': 'Effect', f'FDR_{pid}': 'FDR'}).dropna()
        return ev, list(ev['Symbol']), {'kind': 'profile', 'label': f'profile {names[pid]}'}, s
    cols = [c for i in ids for c in (f'Log2FC_{i}', f'FDR_{i}')]
    eligible = table.dropna(subset=cols)
    try:
        conserved, _ = conserved_signature(table, ids, fdr, lfc)
    except ValueError as exc:
        raise LumosError(str(exc))
    ev = conserved[['Symbol', 'Effect', 'FDR']]
    word = 'dataset' if s['mode'] == 'tissues' else 'analysis'
    return ev, list(eligible['Symbol']), {'kind': 'consistent',
                                          'label': f'consistent signature (same direction in every {word})'}, s


def terms_rows(df):
    return [{'term': r.Term, 'direction': r.Direction, 'p': fnum(r.PValue), 'fdr': fnum(r.FDR),
             'odds': fnum(r.OddsRatio), 'genes': list(r.Genes)} for r in df.itertuples()]


def cmd_enrich(a):
    out = a.outdir

    def read_genes(path):
        if not path:
            return []
        with open(path) as fh:
            return list(dict.fromkeys(g.strip().upper() for g in fh.read().split() if g.strip()))

    up, down = read_genes(a.genes_up), read_genes(a.genes_down)
    if len(set(up) | set(down)) < 3:
        raise LumosError('Select at least three genes for the enrichment analysis.')
    progress(0.03, 'Reading the gene source')
    evidence, universe, source, src_summary = source_evidence(a.source, a.focus, a.fdr, a.lfc)
    custom = read_genes(a.background)
    background = custom or universe
    background_text = 'custom list' if custom else f'all genes tested in the source ({source["label"]})'
    ora = ORA(background)
    lists = ({'Up-regulated': up, 'Down-regulated': down} if a.strategy == 'separate'
             else {'Combined': list(dict.fromkeys(up + down))})
    lists = {k: v for k, v in lists.items() if v}
    libraries = [l for l in a.libraries.split(',') if l]
    tf_libraries = [l for l in a.tf_libraries.split(',') if l]
    all_libs = [(l, 'pathway') for l in libraries] + [(l, 'tf') for l in tf_libraries]
    results, regulators, sizes = [], [], {}
    for i, (lib, role) in enumerate(all_libs):
        progress(0.1 + 0.8 * i / len(all_libs), f'Enrichment: {lib}')
        frames, engines = [], set()
        for direction, genes in lists.items():
            try:
                df, engine, size = ora.run(genes, lib)
            except ValueError as exc:
                say(f'{lib} ({direction}): {exc}')
                continue
            except Exception as exc:
                raise LumosError(f'{lib}: {exc}')
            sizes[direction] = size
            engines.add(engine)
            if not df.empty:
                frames.append(df.assign(Direction=direction))
        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=['Term', 'PValue', 'OddsRatio', 'CombinedScore', 'Genes', 'FDR', 'Overlap', 'Direction'])
        if not a.keep_nonhuman and not df.empty:
            before = len(df)
            df = df[~df['Term'].astype(str).str.contains(NON_HUMAN_TERM)]
            if before - len(df):
                say(f'{lib}: {before - len(df)} mouse/rat terms removed.')
        df = df.sort_values(['FDR', 'PValue'], na_position='last').reset_index(drop=True)
        n_sig = int((df['FDR'] < 0.05).sum()) if not df.empty else 0
        say(f'{lib}: {len(df)} terms, {n_sig} with FDR < 0.05 ({", ".join(sorted(engines)) or "no query"}).')
        fname = f'enrichment_{safe_name(lib)}.csv'
        if not df.empty:
            df.assign(Genes=df['Genes'].apply(';'.join), Library=lib, Engine='; '.join(sorted(engines))).to_csv(
                os.path.join(out, fname), index=False, float_format='%.6g')
        entry = {'name': lib, 'role': role, 'engine': '; '.join(sorted(engines)), 'n_terms': int(len(df)),
                 'n_significant': n_sig, 'file': fname if not df.empty else None, 'terms': terms_rows(df.head(TOP_TERMS))}
        if role == 'tf' and not df.empty:
            for r in df.itertuples():
                regulators.append({'tf': tf_symbol(r.Term), 'term': r.Term, 'library': lib, 'direction': r.Direction,
                                   'p': fnum(r.PValue), 'fdr': fnum(r.FDR), 'targets': list(r.Genes)})
        results.append(entry)
    regulators.sort(key=lambda r: (r['fdr'] if r['fdr'] is not None else 1.0))
    if regulators:
        pd.DataFrame([{**r, 'targets': ';'.join(r['targets'])} for r in regulators]).rename(columns={
            'tf': 'TF', 'term': 'Term', 'library': 'Library', 'direction': 'Direction', 'p': 'PValue', 'fdr': 'FDR',
            'targets': 'Targets'}).to_csv(os.path.join(out, 'regulators.csv'), index=False, float_format='%.6g')
    evidence.to_csv(os.path.join(out, 'evidence.csv'), index=False, float_format='%.6g')
    with open(os.path.join(out, 'background.txt'), 'w') as fh:
        fh.write('\n'.join(ora.background))
    ev = evidence.drop_duplicates('Symbol').set_index('Symbol')
    genes = [{'symbol': g, 'direction': 'up' if g in set(up) else 'down',
              'effect': fnum(ev['Effect'].get(g)), 'fdr': fnum(ev['FDR'].get(g))}
             for g in dict.fromkeys(up + down)]
    src_job = os.path.basename(os.path.normpath(a.source))
    files = {'zip': 'lumos_APP.zip', 'manifest': 'manifest.json', 'regulators': 'regulators.csv' if regulators else None,
             'evidence': 'evidence.csv', 'background': 'background.txt'}
    summary = {'kind': 'enrich', 'version': __version__, 'source': {**source, 'job': src_job, 'focus': a.focus},
               'name': src_summary.get('name') or source['label'], 'taxid': src_summary.get('taxid'),
               'strategy': a.strategy, 'background': {'source': background_text, 'size': len(ora.background)},
               'query_sizes': sizes, 'thresholds': {'fdr': a.fdr, 'lfc': a.lfc}, 'genes': genes,
               'libraries': results, 'regulators': regulators, 'tf_libraries': tf_libraries,
               'notes': ora.diagnostics, 'files': files}
    manifest = {'lumos_version': __version__, 'generated_utc': utc_now(), 'source': summary['source'],
                'strategy': a.strategy, 'thresholds': {'FDR': a.fdr, 'absolute_log2FC': a.lfc},
                'up_regulated_genes': up, 'down_regulated_genes': down,
                'background': {'source': background_text, 'size': len(ora.background)},
                'libraries': [{'name': e['name'], 'role': e['role'], 'engine': e['engine'], 'terms': e['n_terms']}
                              for e in results],
                'non_human_terms_removed': not a.keep_nonhuman, 'diagnostics': ora.diagnostics,
                'runtime': runtime_info(('cloudscraper',))}
    write_json(os.path.join(out, 'summary.json'), summary)
    with open(os.path.join(out, 'manifest.json'), 'w') as fh:
        json.dump(manifest, fh, indent=2, default=str)
    make_zip(out, files['zip'], [e['file'] for e in results if e['file']] +
             ['regulators.csv', 'evidence.csv', 'background.txt', 'manifest.json', 'summary.json', 'run.log'])
    progress(1.0, 'Done')


# ---------------------------------------------------------------------------- STRING (on demand)

def string_network(genes, species=9606, score=400, network_type='functional'):
    """STRING v12 associations among `genes` (nodes keep the submitted names)."""
    genes = list(dict.fromkeys(g for g in genes if g))
    if len(genes) < 2:
        raise LumosError('At least two genes are needed for a STRING network.')
    s = session()
    try:
        r = s.post(f'{STRING_API}/json/get_string_ids', timeout=60,
                   data={'identifiers': '\r'.join(genes), 'species': species, 'limit': 1, 'echo_query': 1,
                         'caller_identity': 'lumos-networks'})
        r.raise_for_status()
        mapping = {m['stringId']: m['queryItem'] for m in r.json()}
        if len(mapping) < 2:
            return {'nodes': [], 'edges': [], 'unmapped': genes, 'mapped': len(mapping)}
        time.sleep(1)       # STRING asks API clients to pause between calls
        r = s.post(f'{STRING_API}/json/network', timeout=90,
                   data={'identifiers': '\r'.join(mapping), 'species': species, 'required_score': score,
                         'network_type': network_type, 'caller_identity': 'lumos-networks'})
        r.raise_for_status()
        rows = r.json()
    except Exception as exc:
        raise LumosError(f'STRING did not answer ({exc}). The service may be busy: try again in a moment.')
    edges, seen = [], set()
    for e in rows:
        a_, b_ = mapping.get(e['stringId_A'], e['preferredName_A']), mapping.get(e['stringId_B'], e['preferredName_B'])
        key = tuple(sorted((a_, b_)))
        if a_ == b_ or key in seen:
            continue
        seen.add(key)
        edges.append({'a': a_, 'b': b_, 'score': fnum(e.get('score'), 3)})
    nodes = sorted({x for e in edges for x in (e['a'], e['b'])})
    return {'nodes': nodes, 'edges': edges, 'mapped': len(mapping),
            'unmapped': [g for g in genes if g not in mapping.values()],
            'provenance': {'version': '12.0', 'species': species, 'required_score': score,
                           'network_type': network_type, 'queried_utc': utc_now()}}


# ---------------------------------------------------------------------------- PrioriGraph (on demand)

def _pg_inputs(enrich_dir, n_terms, term_fdr, tf_fdr):
    s = read_json(os.path.join(enrich_dir, 'summary.json'))
    if not s.get('regulators'):
        raise LumosError('No regulator enrichment is available: run the APP with a regulator library.')
    best = {}
    for lib in s['libraries']:
        if lib['role'] != 'pathway':
            continue
        for t in lib['terms']:
            if t['fdr'] is None or t['fdr'] > term_fdr:
                continue
            cur = best.get(t['term'])
            if cur is None:
                best[t['term']] = {'term': t['term'], 'library': lib['name'], 'fdr': t['fdr'], 'genes': set(t['genes']),
                                   'directions': {t['direction']}}
            else:
                cur['genes'] |= set(t['genes'])
                cur['directions'].add(t['direction'])
                cur['fdr'] = min(cur['fdr'], t['fdr'])
    terms = sorted(best.values(), key=lambda t: t['fdr'])[:n_terms]
    rows = [{'TF': r['tf'], 'Target': g, 'Library': r['library'], 'Direction': r['direction'], 'Term': r['term'],
             'FDR': r['fdr']} for r in s['regulators'] if r['fdr'] is not None and r['fdr'] <= tf_fdr
            for g in r['targets']]
    reg = pd.DataFrame(rows, columns=['TF', 'Target', 'Library', 'Direction', 'Term', 'FDR'])
    reg = reg.sort_values('FDR').drop_duplicates(['TF', 'Target', 'Library', 'Direction'])
    ev = pd.read_csv(os.path.join(enrich_dir, 'evidence.csv'))
    ev['Symbol'] = ev['Symbol'].astype(str).str.upper()
    status = {g['symbol']: g['direction'] for g in s['genes']}
    return s, terms, reg, ev, status


def standard_topology(reg, network_genes):
    """Directed TF -> target metrics: degrees, HITS, PageRank; betweenness on the undirected skeleton."""
    import networkx as nx
    edges = reg[reg['Target'].isin(network_genes)][['TF', 'Target']].drop_duplicates()
    if edges.empty:
        return []
    g = nx.DiGraph()
    g.add_edges_from(edges.itertuples(index=False, name=None))
    pagerank = nx.pagerank(g, alpha=0.85)
    betweenness = nx.betweenness_centrality(g.to_undirected(), normalized=True)
    try:
        hubs, authorities = nx.hits(g, max_iter=1000, tol=1e-8, normalized=True)
    except Exception:
        hubs, authorities = {}, {}
    tfs, targets = set(edges['TF']), set(edges['Target'])
    rows = [{'node': n, 'type': 'TF and target' if n in tfs and n in targets else 'Transcription factor' if n in tfs
             else 'Target gene', 'out_degree': int(g.out_degree(n)), 'in_degree': int(g.in_degree(n)),
             'hits_hub': fnum(hubs.get(n)), 'hits_authority': fnum(authorities.get(n)),
             'pagerank': fnum(pagerank.get(n)), 'betweenness': fnum(betweenness.get(n))} for n in g.nodes]
    rows.sort(key=lambda r: (-(r['hits_hub'] or 0), -(r['pagerank'] or 0), -r['out_degree']))
    return rows


def priorigraph(enrich_dir, n_terms=6, term_fdr=0.05, tf_fdr=0.05):
    """PrioriGraph: candidate regulators -> target genes <- top enriched terms, with rankings."""
    s, terms, reg, ev, status = _pg_inputs(enrich_dir, n_terms, term_fdr, tf_fdr)
    targets_of = reg.groupby('TF')['Target'].apply(set).to_dict()
    all_targets = set(reg['Target'])
    nodes, edges, network_genes, used_terms = {}, [], set(), []
    for t in terms:
        genes = t['genes'] & all_targets
        if not genes:
            continue
        pid = f'pathway:{t["term"]}'
        used_terms.append(t)
        nodes[pid] = {'id': pid, 'label': t['term'], 'kind': 'pathway', 'fdr': t['fdr'], 'library': t['library'],
                      'directions': sorted(t['directions'])}
        for gene in sorted(genes):
            network_genes.add(gene)
            nodes.setdefault(gene, {'id': gene, 'label': gene, 'kind': 'gene', 'status': status.get(gene)})
            edges.append({'from': pid, 'to': gene, 'kind': 'pathway'})
    pairs = set()
    for gene in sorted(network_genes):
        for tf, tg in targets_of.items():
            if gene in tg and (tf, gene) not in pairs:
                pairs.add((tf, gene))
                node = nodes.setdefault(tf, {'id': tf, 'label': tf, 'status': status.get(tf)})
                node['kind'] = 'tf'
                edges.append({'from': tf, 'to': gene, 'kind': 'regulation'})
    if not nodes:
        return {'nodes': [], 'edges': [], 'regulators': [], 'targets': [], 'mapping': [], 'evidence': [],
                'metrics': [], 'edge_table': [], 'stats': {}, 'genes': []}
    shown = reg[reg['Target'].isin(network_genes)]
    submitted = set(status)
    ev_map = ev.drop_duplicates('Symbol').set_index('Symbol')
    regulators = []
    for tf, grp in shown.groupby('TF'):
        regulators.append({'tf': tf, 'targets_in_network': int(grp['Target'].nunique()),
                           'targets_in_list': len(targets_of.get(tf, set()) & submitted), 'fdr': fnum(grp['FDR'].min()),
                           'libraries': ', '.join(sorted(set(grp['Library']))), 'status': status.get(tf)})
    regulators.sort(key=lambda r: (-r['targets_in_network'], r['fdr'] if r['fdr'] is not None else 1))
    targets = []
    for gene in network_genes:
        targets.append({'gene': gene, 'regulators': int(shown.loc[shown['Target'] == gene, 'TF'].nunique()),
                        'terms': sum(1 for t in used_terms if gene in t['genes']), 'status': status.get(gene),
                        'effect': fnum(ev_map['Effect'].get(gene)), 'fdr': fnum(ev_map['FDR'].get(gene))})
    targets.sort(key=lambda r: (-r['regulators'], -r['terms']))
    mapping = []
    for t in used_terms:
        tg = t['genes'] & network_genes
        counts = shown[shown['Target'].isin(tg)].groupby('TF')['Target'].nunique().sort_values(ascending=False)
        if not counts.empty:
            mapping.append({'term': t['term'], 'library': t['library'], 'regulator': counts.index[0],
                            'regulator_targets': int(counts.iloc[0]), 'genes_in_network': len(tg), 'fdr': t['fdr']})
    evidence = sorted(({'gene': g, 'effect': fnum(ev_map['Effect'].get(g)), 'fdr': fnum(ev_map['FDR'].get(g))}
                       for g in network_genes), key=lambda r: (r['fdr'] if r['fdr'] is not None else 1))
    edge_table = [{'tf': r.TF, 'target': r.Target, 'library': r.Library, 'direction': r.Direction, 'term': r.Term,
                   'fdr': fnum(r.FDR)} for r in shown.itertuples()]
    return {'nodes': list(nodes.values()), 'edges': edges, 'regulators': regulators, 'targets': targets,
            'mapping': mapping, 'evidence': evidence, 'metrics': standard_topology(reg, network_genes),
            'edge_table': edge_table, 'genes': sorted(network_genes),
            'stats': {'nodes': len(nodes), 'edges': len(edges), 'terms': len(used_terms),
                      'tfs': sum(1 for n in nodes.values() if n['kind'] == 'tf'), 'genes': len(network_genes),
                      'regulators_used': int(reg['TF'].nunique()), 'terms_screened': len(terms),
                      'tf_libraries': s.get('tf_libraries', [])}}


def propagation(enrich_dir, n_terms=6, term_fdr=0.05, tf_fdr=0.05, restart=0.7, decay=0.5):
    """Personalised random walk with restart and signed linear propagation on the TF -> target layer."""
    graph = priorigraph(enrich_dir, n_terms, term_fdr, tf_fdr)
    if not graph['genes']:
        return {'rows': []}
    _, _, reg, ev, _ = _pg_inputs(enrich_dir, n_terms, term_fdr, tf_fdr)
    effect = ev.drop_duplicates('Symbol').set_index('Symbol')['Effect'].to_dict()
    nodes = sorted(set(graph['genes']) | set(reg['TF']))
    index = {n: i for i, n in enumerate(nodes)}
    adj = np.zeros((len(nodes), len(nodes)))
    for r in reg.itertuples():
        if r.TF in index and r.Target in index:
            adj[index[r.Target], index[r.TF]] = 1.0
    e = np.array([effect.get(n, 0.0) if np.isfinite(effect.get(n, 0.0)) else 0.0 for n in nodes])
    p, iterations, converged = random_walk_with_restart(adj, np.abs(e), restart_probability=restart)
    col = adj.sum(axis=0)
    w = np.zeros_like(adj)
    nz = col > 0
    w[:, nz] = adj[:, nz] / col[nz]
    try:
        pf = np.linalg.solve(np.eye(len(nodes)) - decay * w, e)
    except np.linalg.LinAlgError:
        pf = e
    rows = [{'gene': n, 'effect': fnum(e[i]), 'rwr': fnum(p[i]), 'linear': fnum(pf[i])} for i, n in enumerate(nodes)]
    rows.sort(key=lambda r: -(r['rwr'] or 0))
    return {'rows': rows, 'converged': bool(converged), 'iterations': int(iterations), 'restart': restart,
            'decay': decay}


def sensitivity(enrich_dir, tf_fdr=0.05, fdrs=(0.01, 0.05, 0.10), effects=(0.5, 1.0, 1.5),
                restarts=(0.5, 0.7, 0.9), decays=(0.25, 0.5, 0.75), top_k=20):
    """Stability of the network rankings over a grid of selection and propagation parameters."""
    s = read_json(os.path.join(enrich_dir, 'summary.json'))
    rows = [{'TF': r['tf'], 'Target': g} for r in s.get('regulators', [])
            if r['fdr'] is not None and r['fdr'] <= tf_fdr for g in r['targets']]
    if not rows:
        raise LumosError('No regulator passes the regulator FDR threshold.')
    ev = pd.read_csv(os.path.join(enrich_dir, 'evidence.csv'))
    ev['Symbol'] = ev['Symbol'].astype(str).str.upper()
    try:
        rankings, summary = network_sensitivity_analysis(pd.DataFrame(rows), ev, fdrs, effects, restarts, decays,
                                                         top_k)
    except ValueError as exc:
        raise LumosError(str(exc))
    if summary.empty:
        return {'configurations': 0, 'rows': []}
    cols = {'Node': 'node', 'Configurations_Present': 'present', 'Configuration_Coverage': 'coverage',
            'TopK_Frequency': 'topk', 'Median_RWR_Rank': 'median_rank', 'IQR_RWR_Rank': 'iqr_rank',
            'Mean_RWR_Score': 'mean_score', 'SD_RWR_Score': 'sd_score', 'Median_Linear_Propagation': 'median_linear'}
    out = [{v: (r[k] if k == 'Node' else fnum(r[k])) for k, v in cols.items()} for _, r in summary.iterrows()]
    return {'configurations': int(rankings['Configuration'].nunique()), 'rows': out}


def cluster_order(de_dir, gene_ids, rows=True, cols=False):
    """Ward ordering of heatmap rows (genes) and columns (samples) on row z-scores."""
    from scipy.cluster import hierarchy
    from scipy.spatial import distance
    s = read_json(os.path.join(de_dir, 'summary.json'))
    idx = {g: k for k, g in enumerate(s['heatmap']['genes'])}
    m = np.array([[np.nan if v is None else v for v in s['heatmap']['values'][idx[g]]] for g in gene_ids if g in idx])
    if m.ndim != 2 or m.shape[0] == 0:
        return {'rows': [], 'cols': []}
    z = (m - np.nanmean(m, axis=1, keepdims=True)) / (np.nanstd(m, axis=1, keepdims=True) + 1e-9)
    z = np.nan_to_num(z)
    order_r = list(range(z.shape[0]))
    order_c = list(range(z.shape[1]))
    if rows and z.shape[0] > 2:
        order_r = [int(i) for i in hierarchy.leaves_list(hierarchy.linkage(distance.pdist(z), method='ward'))]
    if cols and z.shape[1] > 2:
        order_c = [int(i) for i in hierarchy.leaves_list(hierarchy.linkage(distance.pdist(z.T), method='ward'))]
    kept = [g for g in gene_ids if g in idx]
    return {'rows': [kept[i] for i in order_r], 'cols': order_c}


# ---------------------------------------------------------------------------- HTML reports (on demand)

REPORT_MODULES = {
    'DDEA': ('ddea-icon.png', 'Diagonal Differential Expression Alley (DDEA) &mdash; Lumos Networks, EvoMol-Lab, '
             'BioME, UFRN.<br>Differential expression: PyDESeq2 (Muzellec et al., 2023) or limma (Ritchie et al., '
             '2015).'),
    'Integration': ('ddea-icon.png', 'Lumos Networks integration of analyses &mdash; EvoMol-Lab, BioME, UFRN.<br>'
                    'Random-effects meta-analysis: statsmodels (Seabold &amp; Perktold, 2010).'),
    'APP': ('app-icon.png', 'Arithmancy Pathway Profiler (APP) &mdash; EvoMol-Lab, BioME, UFRN.<br>Enrichment: '
            'Enrichr/Speedrichr (Xie et al., 2021). Protein associations: STRING v12 (Szklarczyk et al., 2023).'),
    'PG': ('pg-icon.png', 'PrioriGraph (PG) &mdash; EvoMol-Lab, BioME, UFRN.<br>Network metrics: NetworkX (Hagberg '
           'et al., 2008). Regulators: TRRUST (Han et al., 2018) / JASPAR via Enrichr. Propagation: RWR (Köhler et '
           'al., 2008).'),
}
REPORT_FORMATS = ('text', 'fixed3', 'exp2', 'int')


def html_report(payload, static_dir):
    """Self-contained interactive HTML report (lumos_report) from figures and tables sent by the page."""
    import plotly.graph_objects as go
    from lumos_report import build_html_report, table_spec
    module = payload.get('module') if payload.get('module') in REPORT_MODULES else 'DDEA'
    icon, footer = REPORT_MODULES[module]
    sections = []
    for sec in payload.get('sections', [])[:30]:
        figs = []
        for f in (sec.get('figs') or [])[:6]:
            try:
                figs.append(go.Figure(f, skip_invalid=True))
            except Exception:
                continue
        table = None
        t = sec.get('table')
        if t and t.get('columns'):
            columns = [(str(c[0]), str(c[1]), c[2] if c[2] in REPORT_FORMATS else 'text') for c in t['columns']][:40]
            df = pd.DataFrame(t.get('rows', [])[:200000], columns=[c[0] for c in columns])
            table = table_spec(df, columns, t.get('sort') or columns[0][0], bool(t.get('asc', True)),
                               filter_col=t.get('filter'), badge_col=t.get('badge'),
                               badge_map=t.get('badge_map') or {})
        sections.append({'id': str(sec.get('id', 'section')), 'title': str(sec.get('title', '')),
                         'desc': str(sec.get('desc', '')), 'figs': figs, 'table': table})
    info = [(str(k), str(v)) for k, v in payload.get('info', [])][:40]
    info.append(('Lumos Networks', __version__))
    cards = [(str(c[0]), str(c[1]), str(c[2]) if len(c) > 2 else '') for c in payload.get('cards', [])][:10]
    title = str(payload.get('title') or f'{module} report')[:120]
    return build_html_report(title, os.path.join(static_dir, icon), info, cards, sections, footer,
                             str(payload.get('file_id') or module), doc_title=title)


# ---------------------------------------------------------------------------- command line

def build_parser():
    p = argparse.ArgumentParser(prog='lumos.py', description='Lumos Networks: transcriptomics to networks.')
    p.add_argument('--version', action='version', version=f'Lumos Networks {__version__}')
    sub = p.add_subparsers(dest='command', required=True)

    d = sub.add_parser('dataset', help='read a GEO series or uploaded expression tables')
    d.add_argument('--geo', help='GEO series accession (GSE...)')
    d.add_argument('--expression', nargs='*', default=[], help='expression matrix, or one file per sample')
    d.add_argument('--metadata', help='sample table (CSV/TSV, one row per sample)')
    d.add_argument('--platform', choices=['auto', 'rnaseq', 'microarray'], default='auto')
    d.add_argument('--name', help='dataset name (uploads)')
    d.add_argument('--outdir', required=True)

    e = sub.add_parser('de', help='differential expression between two groups (DDEA)')
    e.add_argument('--dataset', required=True, help='folder written by the dataset command')
    e.add_argument('--groups', required=True, help='JSON file: {"group name": [sample ids], ...}')
    e.add_argument('--reference', required=True, help='reference (control) group')
    e.add_argument('--test', required=True, help='test group; log2FC = test vs reference')
    e.add_argument('--label', help='name of the analysis (default: "test vs reference")')
    e.add_argument('--engine', '--method', dest='engine', choices=['auto', 'deseq2', 'limma', 'welch'],
                   default='auto', help='auto: DESeq2 for raw counts, limma (R) for other values')
    e.add_argument('--batch', help='batch column of the sample table')
    e.add_argument('--pair', help='subject/pair column (paired design)')
    e.add_argument('--covariate', action='append', default=[], help='categorical covariate column (repeatable)')
    e.add_argument('--continuous', action='append', default=[], help='continuous covariate column (repeatable)')
    e.add_argument('--unit', choices=['features', 'symbols'], default='features',
                   help='analyse the features as provided, or gene symbols (counts summed, other values median)')
    e.add_argument('--transform', choices=['auto', 'log2', 'none'], default='auto')
    e.add_argument('--normalization', choices=['auto', 'quantile', 'none'], default='auto')
    e.add_argument('--min-count', type=float, default=10)
    e.add_argument('--fdr', type=float, default=0.05)
    e.add_argument('--lfc', type=float, default=1.0)
    e.add_argument('--no-symbols', dest='symbols', action='store_false', help='do not look up gene symbols')
    e.add_argument('--cpus', type=int, default=None)
    e.add_argument('--outdir', required=True)

    g = sub.add_parser('integrate', help='align or meta-analyse several differential-expression analyses')
    g.add_argument('--mode', choices=['contrasts', 'meta', 'tissues'], required=True,
                   help='contrasts of one dataset, meta-analysis (same tissue) or multi-tissue comparison')
    g.add_argument('--analyses', nargs='+', required=True, help='folders written by the de command')
    g.add_argument('--fdr', type=float, default=0.05)
    g.add_argument('--lfc', type=float, default=1.0)
    g.add_argument('--estimator', choices=['Paule-Mandel', 'DerSimonian-Laird'], default='Paule-Mandel')
    g.add_argument('--min-studies', type=int, default=2)
    g.add_argument('--outdir', required=True)

    n = sub.add_parser('enrich', help='pathway and regulator enrichment with an explicit background (APP)')
    n.add_argument('--source', required=True, help='folder of a de or integrate job')
    n.add_argument('--focus', default='analysis', help='integrations: consistent | profile:<id> | meta')
    n.add_argument('--genes-up', help='text file with up-regulated gene symbols')
    n.add_argument('--genes-down', help='text file with down-regulated gene symbols')
    n.add_argument('--strategy', choices=['separate', 'combined'], default='separate',
                   help='query up- and down-regulated genes separately (default) or together')
    n.add_argument('--background', help='custom measured-gene universe (default: genes tested in the source)')
    n.add_argument('--libraries', default=','.join(DEFAULT_PATHWAY_LIBRARIES))
    n.add_argument('--tf-libraries', default=','.join(DEFAULT_TF_LIBRARIES))
    n.add_argument('--fdr', type=float, default=0.05, help='threshold of the consistent signature')
    n.add_argument('--lfc', type=float, default=1.0, help='threshold of the consistent signature')
    n.add_argument('--keep-nonhuman', action='store_true', help='keep mouse/rat terms')
    n.add_argument('--outdir', required=True)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    status = os.path.join(a.outdir, 'status.json')
    say(f'Lumos Networks {__version__}: {a.command}')
    try:
        {'dataset': cmd_dataset, 'de': cmd_de, 'integrate': cmd_integrate, 'enrich': cmd_enrich}[a.command](a)
        write_json(status, {'state': 'done'})
        return 0
    except LumosError as exc:
        say(f'Error: {exc}')
        write_json(status, {'state': 'failed', 'error': str(exc)})
        return 1
    except Exception as exc:
        traceback.print_exc()
        write_json(status, {'state': 'failed', 'error': f'Unexpected error: {exc} (see the log).'})
        return 1


if __name__ == '__main__':
    sys.exit(main())
