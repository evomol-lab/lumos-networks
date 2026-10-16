"""Offline tests of the Lumos Networks engine (no GEO, MyGene, Enrichr or STRING access).

    pip install pytest
    python -m pytest tests
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lumos  # noqa: E402


# ---------------------------------------------------------------------------- identifiers and statistics

@pytest.mark.parametrize('ids, expected', [
    (['ENSG00000141510', 'ENSG00000146648.17', 'ENSG00000136997'], 'ensembl'),
    (['ENSMUSG00000059552', 'ENSMUSG00000020122'], 'ensembl'),
    (['1007_s_at', '1053_at', 'AFFX-BioB-5_at'], 'probe'),
    (['ILMN_1343291', 'ILMN_1343295'], 'probe'),
    (['7157', '1956', '4609'], 'entrez'),
    (['TP53', 'EGFR', 'MYC', 'HLA-A'], 'symbol'),
])
def test_detect_id_type(ids, expected):
    assert lumos.detect_id_type(ids) == expected


def test_strip_version_and_taxon():
    assert lumos.strip_version('ENSG00000146648.17') == 'ENSG00000146648'
    assert lumos.strip_version('ENSG00000182378.14_PAR_Y') == 'ENSG00000182378'
    assert lumos.ensembl_taxid(['ENSMUSG00000059552']) == 10090
    assert lumos.ensembl_taxid(['ENSG00000141510']) == 9606


def test_benjamini_hochberg_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.005, np.nan])
    # R: p.adjust(c(0.01, 0.04, 0.03, 0.005), "BH") = 0.02 0.04 0.04 0.02
    q = lumos.benjamini_hochberg(p)
    assert np.allclose(q[:4], [0.02, 0.04, 0.04, 0.02])
    assert np.isnan(q[4])


def test_pick_symbol_prefers_coding_genes():
    assert lumos.pick_symbol('MIR4640///DDR1') == ('DDR1', True)
    assert lumos.pick_symbol('EXOSC7///CLEC3B') == ('EXOSC7', True)
    assert lumos.pick_symbol('PECAM1') == ('PECAM1', False)
    assert lumos.pick_symbol('LINC00473///FOXO3-AS1') == ('LINC00473', True)
    assert lumos.pick_symbol('---') == ('', False)


def test_value_kind():
    assert lumos.value_kind(np.array([[0, 10, 500], [3, 7, 1200]])) == 'counts'
    assert lumos.value_kind(np.array([[5.2, 8.1], [12.3, 3.3]])) == 'log'
    assert lumos.value_kind(np.array([[12.5, 300.25], [0.5, 80.75]])) == 'continuous'


def test_status_uses_both_thresholds():
    st = lumos.status_of([2.0, -2.0, 0.5, 3.0, np.nan], [0.01, 0.01, 0.001, 0.2, 0.01], 0.05, 1.0)
    assert list(st) == ['Up-regulated', 'Down-regulated', 'Not significant', 'Not significant', 'Not significant']


def test_enrichr_term_helpers():
    assert lumos.tf_symbol('TP53 human') == 'TP53'
    assert lumos.tf_symbol('RELA 24523406 ChIP-Seq FIBROSARCOMA Human') == 'RELA'
    assert lumos.tf_symbol('NFIC::TLX1') == 'NFIC::TLX1'
    assert lumos.NON_HUMAN_TERM.search('Tp53 mouse')
    assert lumos.NON_HUMAN_TERM.search('STAT3 ChIP-Seq Rat liver')
    for term in ('Positive Regulation Of Cell Proliferation', 'Cell Migration', 'Carbohydrate digestion and absorption'):
        assert not lumos.NON_HUMAN_TERM.search(term)


# ---------------------------------------------------------------------------- reading tables

def test_read_table_r_style_header():
    text = '"S1"\t"S2"\t"S3"\n"G1"\t1\t2\t3\n"G2"\t4\t5\t6\n'
    values, symbols = lumos.read_table(text.encode(), 'm.txt')
    assert list(values.columns) == ['S1', 'S2', 'S3']
    assert list(values.index) == ['G1', 'G2']
    assert symbols is None


def test_read_table_headerless_htseq_gzip():
    text = 'ENSG00000000003\t120\nENSG00000000005\t0\n__no_feature\t999\n__ambiguous\t12\n'
    values, _ = lumos.read_table(gzip.compress(text.encode()), 'GSM1_sample.counts.txt.gz')
    assert values.shape == (2, 1)
    assert values.iloc[0, 0] == 120


def test_read_table_featurecounts_with_names():
    text = ('# Program:featureCounts v2.0\n'
            'Geneid\tChr\tStart\tEnd\tStrand\tLength\tgene_name\ta.bam\tb.bam\n'
            'ENSG01\tchr1\t1\t100\t+\t100\tTP53\t10\t20\n'
            'ENSG02\tchr1\t1;50\t40;90\t+;+\t80\tEGFR\t30\t40\n')
    values, symbols = lumos.read_table(text.encode(), 'fc.txt')
    assert list(values.columns) == ['a.bam', 'b.bam']
    assert symbols['ENSG02'] == 'EGFR'


def test_single_sample_prefers_count_column():
    df = pd.DataFrame({'Length': [1, 2], 'TPM': [5.0, 6.0], 'NumReads': [10, 20]}, index=['G1', 'G2'])
    s = lumos.single_sample(df, 'GSM9')
    assert s.name == 'GSM9' and list(s) == [10, 20]


# ---------------------------------------------------------------------------- samples

def test_match_samples_by_accession_title_and_position():
    meta = pd.DataFrame({'Accession': ['GSM1', 'GSM2', 'GSM3'],
                         'Title': ['C13_L001: donor C13', 'C137_L004: donor C137', 'X']})
    cols, how = lumos.match_samples(meta, ['C137_L004', 'C13_L001', 'GSM3_raw'])
    assert cols == ['C13_L001', 'C137_L004', 'GSM3_raw']
    assert how == ['title', 'title', 'identifier']
    cols, how = lumos.match_samples(meta, ['a', 'b', 'c'])
    assert cols == ['a', 'b', 'c'] and set(how) == {'position'}


def test_sample_table_without_metadata_guesses_groups():
    notes = []
    t = lumos.sample_table(None, ['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'], notes)
    assert t['columns'] == ['Sample', 'Name without replicate number']
    assert [r['values'][1] for r in t['rows']] == ['ctrl', 'ctrl', 'treat', 'treat']


def test_sample_table_lists_unmatched_columns():
    meta = pd.DataFrame({'Accession': ['GSM1', 'GSM2']})
    notes = []
    t = lumos.sample_table(meta, ['GSM1', 'GSM2', 'extra'], notes)
    assert t['rows'][-1]['id'] == 'extra' and t['rows'][-1].get('extra')


# ---------------------------------------------------------------------------- whole runs on synthetic data

def synthetic_counts(seed=1, genes=400, per_group=4):
    rng = np.random.default_rng(seed)
    base = rng.gamma(2.0, 150.0, size=genes)
    fold = np.ones(genes)
    fold[:20] = 6.0        # up in the test group
    fold[20:35] = 1 / 6    # down in the test group
    cols = [f'ctrl_{k}' for k in range(1, per_group + 1)] + [f'treat_{k}' for k in range(1, per_group + 1)]
    mat = np.column_stack([rng.negative_binomial(10, 10 / (10 + base * (fold if c.startswith('treat') else 1)))
                           for c in cols])
    return pd.DataFrame(mat, index=[f'GENE{i:04d}' for i in range(genes)], columns=cols)


def run_dataset(tmp_path, frame, sep='\t', name='counts.tsv'):
    path = tmp_path / name
    frame.to_csv(path, sep=sep)
    out = tmp_path / 'data'
    assert lumos.main(['dataset', '--expression', str(path), '--outdir', str(out)]) == 0
    return out, json.loads((out / 'summary.json').read_text())


def run_de(tmp_path, data_dir, engine, extra=(), name=None, groups=None):
    gpath = tmp_path / 'groups.json'
    gpath.write_text(json.dumps(groups or {'ctrl': [f'ctrl_{k}' for k in range(1, 5)],
                                           'treat': [f'treat_{k}' for k in range(1, 5)]}))
    out = tmp_path / (name or f'de_{engine}')
    code = lumos.main(['de', '--dataset', str(data_dir), '--groups', str(gpath), '--reference', 'ctrl',
                       '--test', 'treat', '--engine', engine, '--no-symbols', '--outdir', str(out), *extra])
    assert code == 0, (out / 'status.json').read_text()
    return out, json.loads((out / 'summary.json').read_text())


@pytest.mark.parametrize('engine', ['deseq2', 'welch'])
def test_differential_expression_finds_planted_genes(tmp_path, engine):
    if engine == 'deseq2':
        pytest.importorskip('pydeseq2')
    data_dir, ds = run_dataset(tmp_path, synthetic_counts())
    assert ds['value_kind'] == 'counts' and ds['id_type'] == 'symbol'
    assert [r['column'] for r in ds['samples']['rows']][:2] == ['ctrl_1', 'ctrl_2']
    out, s = run_de(tmp_path, data_dir, engine)
    res = pd.read_csv(out / 'de_results.csv')
    up = set(res.loc[res['Status'] == 'Up-regulated', 'Gene_ID'])
    down = set(res.loc[res['Status'] == 'Down-regulated', 'Gene_ID'])
    planted_up = {f'GENE{i:04d}' for i in range(20)}
    planted_down = {f'GENE{i:04d}' for i in range(20, 35)}
    assert len(up & planted_up) >= 15 and len(down & planted_down) >= 10
    assert len(up - planted_up) <= 3 and len(down - planted_down) <= 3
    assert res['SE'].notna().mean() > 0.9                       # standard errors for meta-analyses
    assert s['pca'] and len(s['pca']['pcs'][0]) == 8
    assert (out / s['files']['zip']).exists()
    manifest = json.loads((out / 'manifest.json').read_text())
    assert manifest['design']['formula'] == '~condition'
    assert manifest['design']['contrast'] == ['condition', 'treat', 'ctrl']


def design_dataset(tmp_path):
    """Synthetic counts with a sample table holding a batch and a subject column."""
    frame = synthetic_counts(seed=4)
    path = tmp_path / 'counts.tsv'
    frame.to_csv(path, sep='\t')
    meta = pd.DataFrame({'sample': frame.columns, 'batch': ['b1', 'b2'] * 4,
                         'subject': ['s1', 's2', 's3', 's4'] * 2, 'age': [30, 41, 52, 63, 34, 45, 56, 67]})
    meta.to_csv(tmp_path / 'samples.csv', index=False)
    out = tmp_path / 'data'
    assert lumos.main(['dataset', '--expression', str(path), '--metadata', str(tmp_path / 'samples.csv'),
                       '--outdir', str(out)]) == 0
    return out


@pytest.mark.parametrize('extra, formula', [
    (['--batch', 'batch'], '~batch + condition'),
    (['--pair', 'subject'], '~subject + condition'),
    (['--batch', 'batch', '--continuous', 'age'], '~batch + age + condition'),
])
def test_design_terms_reach_deseq2(tmp_path, extra, formula):
    pytest.importorskip('pydeseq2')
    out, s = run_de(tmp_path, design_dataset(tmp_path), 'deseq2', extra)
    assert s['design']['formula'] == formula
    assert s['counts']['up'] >= 10


def test_confounded_design_is_refused(tmp_path):
    data_dir = design_dataset(tmp_path)
    gpath = tmp_path / 'groups.json'
    gpath.write_text(json.dumps({'ctrl': ['ctrl_1', 'ctrl_3'], 'treat': ['ctrl_2', 'ctrl_4']}))   # = batch b1 vs b2
    out = tmp_path / 'de'
    assert lumos.main(['de', '--dataset', str(data_dir), '--groups', str(gpath), '--reference', 'ctrl', '--test',
                       'treat', '--batch', 'batch', '--engine', 'welch', '--no-symbols', '--outdir', str(out)]) == 1
    assert 'full rank' in json.loads((out / 'status.json').read_text())['error']


def test_continuous_values_use_limma_or_welch(tmp_path):
    frame = np.log2(synthetic_counts(seed=2) + 1)
    data_dir, ds = run_dataset(tmp_path, frame, sep=',', name='log.csv')
    assert ds['value_kind'] == 'log'
    _, s = run_de(tmp_path, data_dir, 'auto')
    assert s['parameters']['engine'] == ('limma' if lumos.limma_available()[0] else 'welch')
    assert s['counts']['up'] >= 10


@pytest.mark.skipif(not lumos.limma_available()[0], reason='R/limma is not installed')
def test_limma_with_batch(tmp_path):
    frame = np.log2(synthetic_counts(seed=5) + 1)
    path = tmp_path / 'log.tsv'
    frame.to_csv(path, sep='\t')
    pd.DataFrame({'sample': frame.columns, 'batch': ['b1', 'b2'] * 4}).to_csv(tmp_path / 's.csv', index=False)
    data_dir = tmp_path / 'data'
    assert lumos.main(['dataset', '--expression', str(path), '--metadata', str(tmp_path / 's.csv'),
                       '--outdir', str(data_dir)]) == 0
    out, s = run_de(tmp_path, data_dir, 'limma', ['--batch', 'batch'])
    assert s['method_text'].startswith('limma') and s['design']['formula'] == '~batch + condition'
    assert s['counts']['up'] >= 10


def test_deseq2_refuses_log_values(tmp_path):
    pytest.importorskip('pydeseq2')
    data_dir, _ = run_dataset(tmp_path, np.log2(synthetic_counts(seed=3) + 1))
    gpath = tmp_path / 'groups.json'
    gpath.write_text(json.dumps({'a': ['ctrl_1', 'ctrl_2'], 'b': ['treat_1', 'treat_2']}))
    out = tmp_path / 'de'
    assert lumos.main(['de', '--dataset', str(data_dir), '--groups', str(gpath), '--reference', 'a',
                       '--test', 'b', '--engine', 'deseq2', '--no-symbols', '--outdir', str(out)]) == 1
    assert 'raw integer counts' in json.loads((out / 'status.json').read_text())['error']


def test_integration_modes(tmp_path):
    data_dir, _ = run_dataset(tmp_path, synthetic_counts())
    first, _ = run_de(tmp_path, data_dir, 'welch', name='de1')
    pytest.importorskip('pydeseq2')
    second, _ = run_de(tmp_path, data_dir, 'deseq2', ['--label', 'alternative model'], name='de2')
    out = tmp_path / 'int'
    assert lumos.main(['integrate', '--mode', 'contrasts', '--analyses', str(first), str(second),
                       '--outdir', str(out)]) == 0
    s = json.loads((out / 'summary.json').read_text())
    assert s['pattern_counts'].get('Consistent Up (all contrasts)', 0) >= 10
    ev, universe, source, _ = lumos.source_evidence(str(out), 'consistent', 0.05, 1.0)
    assert source['kind'] == 'consistent' and len(ev) >= 10 and 390 <= len(universe) <= 400
    # the same dataset cannot be meta-analysed with itself
    assert lumos.main(['integrate', '--mode', 'meta', '--analyses', str(first), str(second),
                       '--outdir', str(tmp_path / 'bad')]) == 1
    # an independent dataset (another folder and name) can
    other = tmp_path / 'other'
    other.mkdir()
    data2, _ = run_dataset(other, synthetic_counts(seed=7), name='second.tsv')
    third, _ = run_de(other, data2, 'welch', name='de3')
    for mode in ('meta', 'tissues'):
        out = tmp_path / f'int_{mode}'
        assert lumos.main(['integrate', '--mode', mode, '--analyses', str(first), str(third),
                           '--outdir', str(out)]) == 0, (out / 'status.json').read_text()
        s = json.loads((out / 'summary.json').read_text())
        if mode == 'meta':
            assert s['pattern_counts'].get('Up-regulated', 0) >= 10 and s['heatmap']
        else:
            assert s['pattern_counts'].get('Conserved Up', 0) >= 10


def enrich_folder(tmp_path):
    summary = {
        'kind': 'enrich', 'taxid': 9606, 'tf_libraries': ['TRRUST'],
        'genes': [{'symbol': 'A', 'direction': 'up'}, {'symbol': 'B', 'direction': 'down'},
                  {'symbol': 'C', 'direction': 'up'}, {'symbol': 'TF2', 'direction': 'up'}],
        'libraries': [
            {'name': 'KEGG', 'role': 'pathway', 'terms': [
                {'term': 'P1', 'direction': 'Up-regulated', 'fdr': 0.001, 'genes': ['A', 'B']},
                {'term': 'P2', 'direction': 'Up-regulated', 'fdr': 0.01, 'genes': ['B', 'C', 'TF2']}]},
            {'name': 'TRRUST', 'role': 'tf', 'terms': []}],
        'regulators': [
            {'tf': 'TF1', 'term': 'TF1 human', 'library': 'TRRUST', 'direction': 'Up-regulated', 'fdr': 0.01,
             'targets': ['A', 'B', 'C']},
            {'tf': 'TF2', 'term': 'TF2 human', 'library': 'TRRUST', 'direction': 'Up-regulated', 'fdr': 0.2,
             'targets': ['C']}],
    }
    (tmp_path / 'summary.json').write_text(json.dumps(summary))
    pd.DataFrame({'Symbol': ['A', 'B', 'C', 'TF2', 'D'], 'Effect': [2.0, -1.5, 1.2, 1.1, 0.8],
                  'FDR': [0.001, 0.01, 0.02, 0.03, 0.04]}).to_csv(tmp_path / 'evidence.csv', index=False)
    return str(tmp_path)


def test_priorigraph_network_and_rankings(tmp_path):
    folder = enrich_folder(tmp_path)
    g = lumos.priorigraph(folder, n_terms=2, term_fdr=0.05, tf_fdr=1.0)
    kinds = {n['id']: n['kind'] for n in g['nodes']}
    assert kinds['TF1'] == 'tf' and kinds['A'] == 'gene' and kinds['pathway:P1'] == 'pathway'
    assert kinds['TF2'] == 'tf'                         # a regulator that is also a selected gene
    assert g['regulators'][0]['tf'] == 'TF1' and g['regulators'][0]['targets_in_network'] == 3
    assert {t['gene'] for t in g['targets']} >= {'A', 'B', 'C'}
    assert g['mapping'][0]['regulator'] == 'TF1'
    assert {m['node'] for m in g['metrics']} >= {'TF1', 'A'} and g['edge_table'][0]['library'] == 'TRRUST'
    strict = lumos.priorigraph(folder, n_terms=2, term_fdr=0.05, tf_fdr=0.05)
    assert all(n['id'] != 'TF2' for n in strict['nodes'])
    prop = lumos.propagation(folder, 2, 0.05, 1.0, 0.7, 0.5)
    assert prop['converged'] and abs(sum(r['rwr'] for r in prop['rows']) - 1) < 1e-3
    sens = lumos.sensitivity(folder, 1.0, (0.01, 0.05), (0.5, 1.0), (0.5, 0.7), (0.25, 0.5), 3)
    assert sens['configurations'] > 0 and sens['rows'][0]['node']


def test_html_report_escapes_text():
    pytest.importorskip('plotly')
    html = lumos.html_report({'module': 'DDEA', 'title': 'T <script>', 'info': [['Dataset', '<b>x</b>']],
                              'sections': [{'id': 's', 'title': 'Figure', 'figs': [{'data': [{'type': 'bar', 'x': [1], 'y': [2]}],
                                                                                    'layout': {}}],
                                            'table': {'columns': [['g', 'Gene', 'text'], ['f', 'FDR', 'exp2']],
                                                      'rows': [['A', 0.01]], 'sort': 'f'}}]},
                             os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'static'))
    assert '<script>' not in html.split('<body>')[1].split('<footer>')[0].replace('<script type="application/json"', '')[:200]
    assert '&lt;b&gt;x&lt;/b&gt;' in html and 'data-tid=' in html
