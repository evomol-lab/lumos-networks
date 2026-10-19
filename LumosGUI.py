"""Lumos Networks - web interface (Flask).

Every calculation (reading a dataset, differential expression, integration of analyses,
enrichment) runs lumos.py in its own subprocess and job directory, so long downloads and model
fits can be followed (progress + log) and cancelled from the browser. STRING and PrioriGraph
networks, propagation, sensitivity, heatmap clustering and HTML reports are built on request from
the files of finished jobs.

    python LumosGUI.py                      # http://127.0.0.1:5070
    LUMOS_HOST=0.0.0.0 python LumosGUI.py   # listen on all interfaces (containers)
"""
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid

from flask import Flask, Response, abort, jsonify, render_template, request, send_file, send_from_directory
from werkzeug.utils import secure_filename

FROZEN = getattr(sys, 'frozen', False)
APP_ROOT = sys._MEIPASS if FROZEN else os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, template_folder=os.path.join(APP_ROOT, 'templates'),
            static_folder=os.path.join(APP_ROOT, 'static'))
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024 * 1024

JOBS_DIR = os.environ.get('LUMOS_JOBS_DIR', os.path.abspath('lumos_jobs'))
os.makedirs(JOBS_DIR, exist_ok=True)
CACHE_DIR = os.environ.get('LUMOS_CACHE_DIR', os.path.join(JOBS_DIR, '_cache'))   # Enrichr libraries
os.environ['LUMOS_CACHE_DIR'] = CACHE_DIR
RETENTION_DAYS = float(os.environ.get('LUMOS_JOB_RETENTION_DAYS', '0'))   # 0 keeps every job

TABLE_EXT = ('.csv', '.tsv', '.txt', '.tab', '.csv.gz', '.tsv.gz', '.txt.gz', '.tab.gz')
JOB_ID = re.compile(r'^[0-9a-f]{12}$')
LIBRARY = re.compile(r'^[A-Za-z0-9_.\-]{1,80}$')
GENE = re.compile(r'[^A-Za-z0-9_.\-:@/]')
PROCS = {}   # job id -> subprocess.Popen
_LIBRARIES = {'time': 0, 'names': None}
UPDATE = {}   # newer image on Docker Hub, found at start-up (LUMOS_UPDATE_CHECK=1, set in the image)
_LOCK = threading.Lock()


def job_dir(job_id, kinds=None):
    """Folder of an existing job (404 otherwise); with `kinds`, the job must be finished and of that kind."""
    if not JOB_ID.match(job_id or ''):
        abort(404)
    path = os.path.join(JOBS_DIR, job_id)
    if not os.path.isdir(path):
        abort(404)
    if kinds:
        summary = os.path.join(path, 'summary.json')
        if not os.path.exists(summary):
            abort(404)
        with open(summary) as fh:
            if json.load(fh).get('kind') not in kinds:
                abort(404)
    return path


def finished_job(job_id, kinds):
    """Folder of a finished job of one of `kinds`, or a ValueError for the form."""
    if not JOB_ID.match(job_id or ''):
        raise ValueError('The previous step was not found; run it again.')
    path = os.path.join(JOBS_DIR, job_id)
    summary = os.path.join(path, 'summary.json')
    if not os.path.exists(summary):
        raise ValueError('The previous step has not finished (or was deleted); run it again.')
    with open(summary) as fh:
        if json.load(fh).get('kind') not in kinds:
            raise ValueError('The previous step is of the wrong kind.')
    return path


def purge_old_jobs():
    if RETENTION_DAYS <= 0:
        return
    limit = time.time() - RETENTION_DAYS * 86400
    for name in os.listdir(JOBS_DIR):
        path = os.path.join(JOBS_DIR, name)
        if JOB_ID.match(name) and os.path.isdir(path) and os.path.getmtime(path) < limit:
            shutil.rmtree(path, ignore_errors=True)


def save_upload(storage, dest_dir):
    name = secure_filename(storage.filename or '')
    if not name.lower().endswith(TABLE_EXT):
        raise ValueError(f'{storage.filename}: upload text tables (.csv, .tsv, .txt, optionally .gz).')
    path = os.path.join(dest_dir, name)
    base, ext = os.path.splitext(path)
    k = 1
    while os.path.exists(path):
        path = f'{base}_{k}{ext}'
        k += 1
    storage.save(path)
    return path


def number(form, key, cast, default, lo, hi):
    raw = form.get(key, '')
    if raw in ('', None):
        return default
    try:
        value = cast(raw)
    except ValueError:
        raise ValueError(f'Invalid value for {key}: {raw}')
    if not lo <= value <= hi:
        raise ValueError(f'{key} must be between {lo} and {hi}.')
    return value


def choice(form, key, allowed, default):
    return form.get(key) if form.get(key) in allowed else default


def gene_list(text, limit=20000):
    genes = [g for g in (GENE.sub('', g)[:40] for g in re.split(r'[\s,;]+', text or '') if g) if g]
    if len(genes) > limit:     # refused rather than cut silently
        raise ValueError(f'At most {limit} genes per list ({len(genes)} were given).')
    return genes


def write_lines(path, items):
    with open(path, 'w') as fh:
        fh.write('\n'.join(items))
    return path


def dataset_arguments(form, inputs):
    args = ['--platform', choice(form, 'platform', ('auto', 'rnaseq', 'microarray'), 'auto')]
    files = [f for f in request.files.getlist('expression') if f.filename]
    paths = [save_upload(f, inputs) for f in files]
    if form.get('input_mode') == 'geo':
        acc = form.get('geo', '').strip().upper()
        if not re.match(r'^GSE\d{1,8}$', acc):
            raise ValueError('Use a GEO series accession such as GSE117769.')
        args += ['--geo', acc]
    else:
        if not paths:
            raise ValueError('Upload an expression matrix (or one file per sample).')
        meta = request.files.get('metadata')
        if meta and meta.filename:
            args += ['--metadata', save_upload(meta, inputs)]
        name = re.sub(r'[^A-Za-z0-9_.\-]', '_', form.get('name', '').strip())[:60]
        if name:
            args += ['--name', name]
    if paths:
        args += ['--expression'] + paths
    return args, {'geo': form.get('geo', '').strip().upper() or None, 'uploads': [os.path.basename(p) for p in paths]}


def de_arguments(form, inputs):
    ds = finished_job(form.get('dataset'), ('dataset',))
    try:
        groups = json.loads(form.get('groups', '{}'))
    except ValueError:
        raise ValueError('The groups could not be read.')
    if not isinstance(groups, dict) or len(groups) > 50:
        raise ValueError('Define the groups in the sample table.')
    clean = {}
    for name, ids in groups.items():
        name = str(name).strip()[:60]
        if not name or not isinstance(ids, list):
            continue
        clean[name] = [str(i)[:200] for i in ids][:5000]
    ref, test = form.get('reference', '').strip(), form.get('test', '').strip()
    if ref not in clean or test not in clean or ref == test:
        raise ValueError('Choose a reference group and a different test group.')
    with open(os.path.join(ds, 'summary.json')) as fh:
        columns = set(json.load(fh)['samples']['columns'])
    design = json.loads(form.get('design', '{}') or '{}')

    def column(value):
        value = str(value or '').strip()
        if value and value not in columns:
            raise ValueError(f'Unknown design column: {value}')
        return value

    args = ['--dataset', ds, '--groups', os.path.join(inputs, 'groups.json'), '--reference', ref, '--test', test,
            '--engine', choice(form, 'engine', ('auto', 'deseq2', 'limma', 'welch'), 'auto'),
            '--unit', choice(form, 'unit', ('features', 'symbols'), 'features'),
            '--transform', choice(form, 'transform', ('auto', 'log2', 'none'), 'auto'),
            '--normalization', choice(form, 'normalization', ('auto', 'quantile', 'none'), 'auto'),
            '--min-count', str(number(form, 'min_count', float, 10, 0, 1e6)),
            '--fdr', str(number(form, 'fdr', float, 0.05, 1e-10, 1)),
            '--lfc', str(number(form, 'lfc', float, 1.0, 0, 20))]
    label = re.sub(r'[\x00-\x1f]', '', form.get('label', '')).strip()[:80]
    if label:
        args += ['--label', label]
    if column(design.get('batch')):
        args += ['--batch', column(design['batch'])]
    if column(design.get('pair')):
        args += ['--pair', column(design['pair'])]
    for c in (design.get('covariates') or [])[:10]:
        args += ['--covariate', column(c)]
    for c in (design.get('continuous') or [])[:10]:
        args += ['--continuous', column(c)]
    if form.get('symbols') == 'false':
        args.append('--no-symbols')
    with open(os.path.join(inputs, 'groups.json'), 'w') as fh:
        json.dump(clean, fh)
    return args, {'dataset': form.get('dataset'), 'reference': ref, 'test': test, 'label': label}


def integrate_arguments(form, inputs):
    ids = [i for i in form.get('analyses', '').split(',') if i][:20]
    if len(ids) < 2:
        raise ValueError('Choose at least two analyses.')
    paths = [finished_job(i, ('de',)) for i in ids]
    mode = choice(form, 'mode', ('contrasts', 'meta', 'tissues'), 'contrasts')
    args = ['--mode', mode, '--analyses', *paths,
            '--fdr', str(number(form, 'fdr', float, 0.05, 1e-10, 1)),
            '--lfc', str(number(form, 'lfc', float, 1.0, 0, 20)),
            '--estimator', choice(form, 'estimator', ('Paule-Mandel', 'DerSimonian-Laird'), 'Paule-Mandel'),
            '--min-studies', str(number(form, 'min_studies', int, 2, 2, 20))]
    return args, {'analyses': ids, 'mode': mode}


def import_arguments(form, inputs):
    paths = [save_upload(f, inputs) for f in request.files.getlist('tables') if f.filename][:10]
    if not paths:
        raise ValueError('Upload at least one result table (.csv, .tsv or .txt).')
    return ['--tables', *paths], {'uploads': [os.path.basename(p) for p in paths]}


def enrich_arguments(form, inputs):
    source = finished_job(form.get('source'), ('de', 'integrate', 'results'))
    up, down = gene_list(form.get('genes_up')), gene_list(form.get('genes_down'))
    if len(set(up) | set(down)) < 3:
        raise ValueError('Select at least three genes.')
    libs = [l for l in form.get('libraries', '').split(',') if LIBRARY.match(l)][:8]
    if not libs:
        raise ValueError('Choose at least one pathway library.')
    tfs = [l for l in form.get('tf_libraries', '').split(',') if LIBRARY.match(l)][:4]
    focus = form.get('focus', 'analysis')
    if not re.match(r'^(analysis|consistent|meta|profile:[A-Za-z0-9_]{1,80})$', focus):
        raise ValueError('Unknown gene source focus.')
    args = ['--source', source, '--focus', focus,
            '--genes-up', write_lines(os.path.join(inputs, 'genes_up.txt'), up),
            '--genes-down', write_lines(os.path.join(inputs, 'genes_down.txt'), down),
            '--strategy', choice(form, 'strategy', ('separate', 'combined'), 'separate'),
            '--libraries', ','.join(libs), '--tf-libraries', ','.join(tfs),
            '--fdr', str(number(form, 'fdr', float, 0.05, 1e-10, 1)),
            '--lfc', str(number(form, 'lfc', float, 1.0, 0, 20))]
    background = gene_list(form.get('background'), 100000)
    if background:
        args += ['--background', write_lines(os.path.join(inputs, 'background.txt'), background)]
    if form.get('keep_nonhuman') == 'true':
        args.append('--keep-nonhuman')
    tables = [save_upload(f, inputs) for f in request.files.getlist('regulator_tables') if f.filename][:10]
    if tables:
        args += ['--regulator-tables', *tables]
    return args, {'source': form.get('source'), 'focus': focus, 'n_genes': len(set(up) | set(down))}


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/documentation')
def documentation():
    path = os.path.join(APP_ROOT, 'DOCUMENTATION.md')
    if os.path.exists(path):
        return send_file(path, mimetype='text/markdown')
    return 'Documentation not found', 404


@app.route('/api/capabilities')
def capabilities():
    import importlib.util
    import lumos
    limma, limma_version = lumos.limma_available()
    return jsonify({'version': lumos.__version__,
                    'deseq2': importlib.util.find_spec('pydeseq2') is not None,
                    'limma': limma, 'limma_version': limma_version,
                    'pathway_libraries': list(lumos.DEFAULT_PATHWAY_LIBRARIES),
                    'tf_libraries_default': list(lumos.DEFAULT_TF_LIBRARIES), 'tf_libraries': list(lumos.TF_LIBRARIES),
                    'update': UPDATE or None})


@app.route('/api/enrichr/libraries')
def enrichr_libraries():
    """Gene-set libraries offered by Enrichr (cached for a day)."""
    import lumos
    with _LOCK:
        if _LIBRARIES['names'] and time.time() - _LIBRARIES['time'] < 86400:
            return jsonify({'libraries': _LIBRARIES['names']})
    try:
        r = lumos.session().get(f'{lumos.ENRICHR}/datasetStatistics', timeout=20)
        r.raise_for_status()
        names = sorted(l['libraryName'] for l in r.json()['statistics'])
    except Exception:
        return jsonify({'libraries': [], 'error': 'Enrichr could not be reached.'})
    with _LOCK:
        _LIBRARIES.update(time=time.time(), names=names)
    return jsonify({'libraries': names})


@app.route('/api/jobs', methods=['POST'])
def submit():
    form = request.form
    action = form.get('action')
    builders = {'dataset': dataset_arguments, 'de': de_arguments, 'integrate': integrate_arguments,
                'enrich': enrich_arguments, 'import': import_arguments}
    if action not in builders:
        return jsonify({'error': 'Unknown action.'}), 400
    job_id = uuid.uuid4().hex[:12]
    path = os.path.join(JOBS_DIR, job_id)
    inputs = os.path.join(path, 'input')
    os.makedirs(inputs)
    try:
        args, info = builders[action](form, inputs)
    except ValueError as exc:
        shutil.rmtree(path, ignore_errors=True)
        return jsonify({'error': str(exc)}), 400
    with open(os.path.join(path, 'job.json'), 'w') as fh:
        json.dump({'id': job_id, 'created': time.time(), 'action': action, 'info': info, 'args': args}, fh)
    if FROZEN:
        cmd = [sys.executable, '--run-lumos', action] + args + ['--outdir', path]
    else:
        cmd = [sys.executable, os.path.join(APP_ROOT, 'lumos.py'), action] + args + ['--outdir', path]
    log = open(os.path.join(path, 'run.log'), 'w')
    env = dict(os.environ, PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8', MPLBACKEND='Agg')
    PROCS[job_id] = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=path, env=env)
    log.close()
    return jsonify({'job': job_id})


@app.route('/api/jobs/<job_id>')
def job_status(job_id):
    path = job_dir(job_id)
    proc = PROCS.get(job_id)
    status_file = os.path.join(path, 'status.json')
    state, error = 'running', None
    if os.path.exists(status_file):
        with open(status_file) as fh:
            st = json.load(fh)
        state, error = st['state'], st.get('error')
    elif proc is None or proc.poll() is not None:
        if os.path.exists(os.path.join(path, 'cancelled')):
            state = 'cancelled'
        else:
            state, error = 'failed', 'The calculation stopped unexpectedly (see the log).'
    lines, prog, message = [], 0.0, ''
    log_path = os.path.join(path, 'run.log')
    if os.path.exists(log_path):
        with open(log_path, encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('@@progress'):
                    parts = line.split(maxsplit=2)
                    prog = float(parts[1])
                    message = parts[2].strip() if len(parts) > 2 else ''
                elif 'Warning' not in line and not line.startswith('  warnings.warn'):
                    lines.append(line.rstrip('\n'))
    out = {'state': state, 'error': error, 'progress': prog, 'message': message, 'log': '\n'.join(lines[-400:])}
    with open(os.path.join(path, 'job.json')) as fh:
        out['action'] = json.load(fh)['action']
    if state == 'done':
        with open(os.path.join(path, 'summary.json')) as fh:
            out['summary'] = json.load(fh)
    return jsonify(out)


@app.route('/api/jobs/<job_id>/cancel', methods=['POST'])
def cancel(job_id):
    path = job_dir(job_id)
    proc = PROCS.get(job_id)
    if proc and proc.poll() is None:
        open(os.path.join(path, 'cancelled'), 'w').close()
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
    return jsonify({'state': 'cancelled'})


@app.route('/api/jobs/<job_id>/files/<path:filename>')
def job_file(job_id, filename):
    path = job_dir(job_id)
    if filename.startswith('input/'):
        abort(404)
    return send_from_directory(path, filename, as_attachment=request.args.get('download') == '1')


def on_demand(fn):
    """Run an on-demand calculation and turn its errors into JSON messages."""
    import lumos
    try:
        return jsonify(fn())
    except (ValueError, lumos.LumosError) as exc:
        return jsonify({'error': str(exc)}), 400


def pg_params(args):
    return (number(args, 'terms', int, 6, 1, 50), number(args, 'term_fdr', float, 0.05, 1e-10, 1),
            number(args, 'tf_fdr', float, 0.05, 1e-10, 1))


@app.route('/api/jobs/<job_id>/cluster', methods=['POST'])
def cluster(job_id):
    import lumos
    path = job_dir(job_id, ('de',))
    data = request.get_json(silent=True) or {}
    genes = [str(g)[:200] for g in data.get('genes', [])][:300]
    return on_demand(lambda: lumos.cluster_order(path, genes, bool(data.get('rows')), bool(data.get('cols'))))


@app.route('/api/jobs/<job_id>/string', methods=['POST'])
def string_network(job_id):
    import lumos
    path = job_dir(job_id, ('enrich',))
    data = request.get_json(silent=True) or {}
    genes = [str(g)[:40] for g in data.get('genes', []) if g][:500]
    score = data.get('score', 400) if data.get('score') in (150, 400, 700, 900) else 400
    network_type = data.get('network_type') if data.get('network_type') in ('functional', 'physical') else 'functional'
    with open(os.path.join(path, 'summary.json')) as fh:
        taxid = json.load(fh).get('taxid') or 9606
    return on_demand(lambda: lumos.string_network(genes, taxid, score, network_type))


@app.route('/api/jobs/<job_id>/priorigraph')
def priorigraph(job_id):
    import lumos
    path = job_dir(job_id, ('enrich',))
    return on_demand(lambda: lumos.priorigraph(path, *pg_params(request.args)))


@app.route('/api/jobs/<job_id>/propagation')
def propagation(job_id):
    import lumos
    path = job_dir(job_id, ('enrich',))
    return on_demand(lambda: lumos.propagation(path, *pg_params(request.args),
                                               number(request.args, 'restart', float, 0.7, 0.01, 0.99),
                                               number(request.args, 'decay', float, 0.5, 0.0, 0.99)))


@app.route('/api/jobs/<job_id>/sensitivity', methods=['POST'])
def sensitivity(job_id):
    import lumos
    path = job_dir(job_id, ('enrich',))
    data = request.get_json(silent=True) or {}

    def grid(key, lo, hi, default):
        values = sorted({float(v) for v in data.get(key, []) if isinstance(v, (int, float)) and lo <= v <= hi})
        return values[:8] or list(default)

    return on_demand(lambda: lumos.sensitivity(
        path, number(data, 'tf_fdr', float, 0.05, 1e-10, 1),
        grid('fdrs', 1e-6, 1, (0.01, 0.05, 0.10)), grid('effects', 0, 20, (0.5, 1.0, 1.5)),
        grid('restarts', 0.01, 0.99, (0.5, 0.7, 0.9)), grid('decays', 0, 0.99, (0.25, 0.5, 0.75)),
        number(data, 'top_k', int, 20, 1, 500)))


@app.route('/api/report', methods=['POST'])
def report():
    """Self-contained interactive HTML report of the figures and tables shown on the page."""
    import lumos
    payload = request.get_json(silent=True) or {}
    try:
        html = lumos.html_report(payload, os.path.join(APP_ROOT, 'static'))
    except Exception as exc:
        return jsonify({'error': f'The report could not be built ({exc}).'}), 400
    name = secure_filename(f'Lumos_{payload.get("module", "report")}_{payload.get("file_id", "")}.html') or 'Lumos.html'
    return Response(html, mimetype='text/html', headers={'Content-Disposition': f'attachment; filename="{name}"'})


def release(tag):
    m = re.fullmatch(r'(\d+)\.(\d+)\.(\d+)', tag or '')
    return tuple(int(x) for x in m.groups()) if m else None


def check_for_update():
    """Compare this version with the release tags of the image on Docker Hub; a newer one is logged and
    announced in the page (the container cannot replace itself: `docker pull`, or `docker run --pull always`)."""
    import lumos
    repo = os.environ.get('LUMOS_UPDATE_REPO', 'evomol/lumos-networks')
    try:
        r = lumos.session().get(f'https://hub.docker.com/v2/repositories/{repo}/tags', params={'page_size': 100},
                                timeout=10)
        r.raise_for_status()
        tags = [t['name'] for t in r.json().get('results', [])]
    except Exception as exc:
        print(f'Update check: Docker Hub could not be reached ({type(exc).__name__}).', flush=True)
        return
    newest = max((t for t in tags if release(t)), key=release, default=None)
    if newest and release(newest) > (release(lumos.__version__) or (0, 0, 0)):
        UPDATE.update(version=newest, current=lumos.__version__, command=f'docker pull docker.io/{repo}:{newest}')
        print(f'Update check: Lumos Networks {newest} is available on Docker Hub (this image is '
              f'{lumos.__version__}). Update with: {UPDATE["command"]}', flush=True)
    else:
        print(f'Update check: {lumos.__version__} is the newest release on Docker Hub ({repo}).', flush=True)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--run-lumos':
        import lumos
        sys.exit(lumos.main(sys.argv[2:]))
    purge_old_jobs()
    host = os.environ.get('LUMOS_HOST', '127.0.0.1')
    port = int(os.environ.get('LUMOS_PORT', '5070'))
    url = f'http://{"127.0.0.1" if host == "0.0.0.0" else host}:{port}/'
    print(f'Lumos Networks: open {url} in your browser (Ctrl+C to stop).')
    if host in ('0.0.0.0', '::'):
        print(f'Listening on every network interface: from other computers, open http://<address of this '
              f'machine or server>:{port}/ (with Docker, the host port given to -p).')
    print(f'Jobs are stored in {JOBS_DIR}')
    if os.environ.get('LUMOS_UPDATE_CHECK') == '1':
        threading.Thread(target=check_for_update, daemon=True).start()
    if FROZEN or '--browser' in sys.argv:
        import webbrowser
        threading.Timer(1.25, lambda: webbrowser.open_new(url)).start()
    try:
        from waitress import serve
    except ImportError:
        # Development server; debug mode stays off: Werkzeug's interactive debugger would let
        # anyone who reaches the port run code.
        app.run(host=host, port=port, debug=False, threaded=True)
        return
    # Production server, one process with several threads: every request sees the running jobs (PROCS).
    serve(app, host=host, port=port, threads=int(os.environ.get('LUMOS_THREADS', '8')),
          max_request_body_size=app.config['MAX_CONTENT_LENGTH'], ident='Lumos Networks')


if __name__ == '__main__':
    main()
