'use strict';

document.addEventListener('DOMContentLoaded', () => {
    const $ = (id) => document.getElementById(id);
    const esc = (s) => String(s ?? '').replace(/[&<>"']/g, (c) => ({
        '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    })[c]);
    const fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v)) ? '–' : Number(v).toFixed(d);
    const fmtP = (v) => (v === null || v === undefined || Number.isNaN(v)) ? '–'
        : (v !== 0 && Math.abs(v) < 1e-3 ? Number(v).toExponential(2) : Number(v).toPrecision(3));
    const int = (v) => Number(v).toLocaleString('en-US');
    const css = (name) => getComputedStyle(document.body).getPropertyValue(name).trim();
    const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
    const num = (v) => (v === null || v === undefined ? NaN : v);
    const store = {
        get(k) { try { return JSON.parse(localStorage.getItem(k)); } catch (e) { return null; } },
        set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* private mode */ } },
    };

    const GROUP_COLORS = ['--chart-1', '--chart-2', '--chart-3', '--chart-4', '--chart-5', '--chart-6'];
    const ID_NAMES = { ensembl: 'Ensembl ID', entrez: 'Entrez ID', probe: 'Probe ID', symbol: 'Gene ID', unknown: 'Original ID' };
    const ID_TYPE_TEXT = { ensembl: 'Ensembl gene IDs', entrez: 'NCBI Gene (Entrez) IDs', probe: 'microarray probes',
        symbol: 'gene symbols', unknown: 'unrecognised IDs' };
    const KIND_TEXT = { counts: 'raw counts', log: 'log-scale values', continuous: 'non-integer values' };
    const PLATFORM_TEXT = { rnaseq: 'RNA-Seq', microarray: 'Microarray', unknown: 'not specified' };
    const MODE_NAMES = { contrasts: 'Contrasts of one dataset', meta: 'Meta-analysis', tissues: 'Multi-tissue comparison' };

    const state = {
        caps: null,
        files: [],
        datasets: new Map(), analyses: new Map(), integrations: new Map(),
        dsJob: null, dataset: null,            // dataset shown in steps 1-3
        groups: [],                            // [{name, color, ids: Set}]
        deJob: null, de: null, G: null,        // analysis shown in the results
        status: null,                          // Int8Array: 1 up, -1 down, 0 not significant
        labels: 'symbol', fdr: 0.05, lfc: 1, highlight: [],
        intJob: null, int: null,
        appSource: null,                       // {value, job, focus, kind, rows, universe}
        appJob: null, app: null,
        libraries: [], tfLibraries: new Set(),
        stringData: null, stringNet: null,
        pg: null, pgNet: null, pgNodes: null, pgProp: null, pgSens: null,
        pgAutoDone: false, pgRelaxed: null,          // automatic relaxation of the regulator FDR (first view only)
    };

    // ------------------------------------------------------------------ theme
    function applyTheme(light, persist) {
        document.body.classList.toggle('light-mode', light);
        $('theme-icon-dark').style.display = light ? 'none' : 'block';
        $('theme-icon-light').style.display = light ? 'block' : 'none';
        if (persist) {
            try { localStorage.setItem('lumos-theme', light ? 'light' : 'dark'); } catch (e) { /* private mode */ }
        }
        if (state.de) renderDECharts();
        if (state.int) renderIntegrationCharts();
        if (state.app) renderAppCharts();
        if (state.stringData) drawString();
        if (state.pg) drawPG();
        if (state.groups.length) renderGroups();
    }

    let savedTheme = null;
    try { savedTheme = localStorage.getItem('lumos-theme'); } catch (e) { /* private mode */ }
    applyTheme(savedTheme ? savedTheme === 'light'
        : window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches, false);
    $('theme-toggle').addEventListener('click', () => applyTheme(!document.body.classList.contains('light-mode'), true));

    function withLightTheme(fn) {
        const was = document.body.classList.contains('light-mode');
        document.body.classList.add('light-mode');
        try { return fn(); } finally { if (!was) document.body.classList.remove('light-mode'); }
    }

    // ------------------------------------------------------------------ documentation tab
    function renderMarkdown(md) {
        // keep TeX away from the Markdown parser (underscores would become emphasis)
        const math = [];
        const keep = (m) => { math.push(m); return `@@MATH${math.length - 1}@@`; };
        md = md.replace(/\$\$[\s\S]+?\$\$/g, keep).replace(/\$[^$\n]+?\$/g, keep);
        return marked.parse(md).replace(/@@MATH(\d+)@@/g, (_, i) => esc(math[Number(i)]));
    }

    let docsLoaded = false;
    function showView(view) {
        const docs = view === 'documentation';
        $('nav-calculator').classList.toggle('active', !docs);
        $('nav-documentation').classList.toggle('active', docs);
        $('calculator-view').style.display = docs ? 'none' : 'block';
        $('documentation-view').style.display = docs ? 'block' : 'none';
        document.querySelector('.container').classList.toggle('expanded', docs || !!state.dataset);
        if (docs && !docsLoaded) {
            fetch('/documentation').then((r) => {
                if (!r.ok) throw new Error('Documentation file not found on the server.');
                return r.text();
            }).then((md) => {
                const el = $('doc-content');
                el.innerHTML = renderMarkdown(md);
                if (window.renderMathInElement) {
                    renderMathInElement(el, {
                        delimiters: [{ left: '$$', right: '$$', display: true }, { left: '$', right: '$', display: false }],
                        throwOnError: false,
                    });
                }
                docsLoaded = true;
            }).catch((e) => {
                $('doc-content').innerHTML = `<div class="doc-loading">${esc(e.message)}</div>`;
            });
        }
        if (!docs) setTimeout(resizeAll, 50);
    }
    $('nav-calculator').addEventListener('click', () => showView('calculator'));
    $('nav-documentation').addEventListener('click', () => showView('documentation'));

    // ------------------------------------------------------------------ errors and URL
    function showError(msg, el = $('error-message')) {
        el.textContent = msg || '';
        el.style.display = msg ? 'block' : 'none';
    }

    function updateUrl() {
        const p = new URLSearchParams();
        if (state.datasets.size) p.set('ds', [...state.datasets.keys()].join(','));
        if (state.analyses.size) p.set('de', [...state.analyses.keys()].join(','));
        if (state.integrations.size) p.set('int', [...state.integrations.keys()].join(','));
        if (state.appJob) p.set('app', state.appJob);
        if (state.dsJob) p.set('d', state.dsJob);
        if (state.deJob) p.set('a', state.deJob);
        const q = p.toString().replace(/%2C/g, ',');
        history.replaceState(null, '', q ? `?${q}` : window.location.pathname);
    }

    // ------------------------------------------------------------------ jobs (one subprocess each)
    class Job {
        constructor(panelId, button, errorEl, onDone) {
            this.panel = $(panelId);
            this.button = button;
            this.errorEl = errorEl;
            this.onDone = onDone;
            this.bar = this.panel.querySelector('.progress-bar');
            this.text = this.panel.querySelector('.progress-text');
            this.log = this.panel.querySelector('.log-box');
            this.id = null;
            this.timer = null;
            this.panel.querySelector('.cancel-btn').addEventListener('click', () => this.cancel());
        }

        running() { return !!this.timer; }

        async start(fd, text) {
            showError('', this.errorEl);
            const r = await fetch('/api/jobs', { method: 'POST', body: fd });
            const d = await r.json();
            if (!r.ok) throw new Error(d.error || 'The job could not be started.');
            this.follow(d.job, text);
            return d.job;
        }

        follow(id, text) {
            this.id = id;
            this.panel.style.display = 'block';
            this.bar.style.width = '0%';
            this.text.textContent = text || 'Running…';
            this.log.textContent = '';
            this.button.classList.add('loading');
            this.button.disabled = true;
            this.timer = setTimeout(() => this.poll(), 400);
        }

        async poll() {
            let d;
            try {
                const r = await fetch(`/api/jobs/${this.id}`);
                d = await r.json();
            } catch (e) {
                this.timer = setTimeout(() => this.poll(), 2000);
                return;
            }
            this.bar.style.width = `${Math.round(100 * d.progress)}%`;
            if (d.message) this.text.textContent = d.message + '…';
            const atBottom = this.log.scrollTop + this.log.clientHeight >= this.log.scrollHeight - 20;
            this.log.textContent = d.log;
            if (atBottom) this.log.scrollTop = this.log.scrollHeight;
            if (d.state === 'running') {
                this.timer = setTimeout(() => this.poll(), 1000);
                return;
            }
            this.timer = null;
            this.button.classList.remove('loading');
            this.button.disabled = false;
            if (d.state === 'done') {
                this.text.textContent = 'Done.';
                setTimeout(() => { if (!this.running()) this.panel.style.display = 'none'; }, 1200);
                this.onDone(this.id, d.summary, true);
            } else if (d.state === 'cancelled') {
                this.text.textContent = 'Cancelled.';
            } else {
                this.text.textContent = 'Failed.';
                this.log.parentElement.open = true;
                showError(d.error || 'The calculation failed; see the log.', this.errorEl);
            }
            refreshButtons();
        }

        async cancel() {
            if (!this.id || !this.running()) return;
            await fetch(`/api/jobs/${this.id}/cancel`, { method: 'POST' });
        }
    }

    // ------------------------------------------------------------------ selectable, sortable, filterable table
    class DataTable {
        constructor(table, opts) {
            this.table = table;
            this.opts = opts;
            this.rows = [];
            this.filter = '';
            this.sortKey = null;
            this.sortDir = 1;
            this.selected = new Set();
            this.limit = opts.pageSize || 300;
            this.last = null;
            this.visible = [];
            table.addEventListener('click', (e) => this.click(e));
            if (opts.selectable) {
                table.addEventListener('mousedown', (e) => { if (e.shiftKey) e.preventDefault(); });
            }
        }

        setRows(rows) {
            this.rows = rows;
            this.limit = this.opts.pageSize || 300;
            this.last = null;
            this.render();
        }

        setFilter(text) {
            this.filter = text;
            this.limit = this.opts.pageSize || 300;
            this.last = null;
            this.render();
        }

        columns() { return this.opts.columns(); }

        view() {
            let rows = this.rows;
            const terms = this.filter.trim().toLowerCase().split(/\s+/).filter(Boolean);
            if (terms.length) rows = rows.filter((r) => { const t = this.opts.text(r); return terms.every((w) => t.includes(w)); });
            if (this.opts.extraFilter) rows = rows.filter(this.opts.extraFilter);
            const col = this.sortKey && this.columns().find((c) => c.key === this.sortKey);
            if (col) {
                const val = col.sort || col.value;
                const dir = this.sortDir;
                rows = rows.slice().sort((a, b) => {
                    const x = val(a), y = val(b);
                    const xn = x === null || x === undefined || Number.isNaN(x), yn = y === null || y === undefined || Number.isNaN(y);
                    if (xn || yn) return xn - yn;
                    if (typeof x === 'number' && typeof y === 'number') return (x - y) * dir;
                    return String(x).localeCompare(String(y), undefined, { numeric: true }) * dir;
                });
            }
            return rows;
        }

        enabled(r) { return !(this.opts.isDisabled && this.opts.isDisabled(r)); }

        render() {
            const cols = this.columns();
            const rows = this.view();
            this.visible = rows;
            const key = this.opts.rowKey;
            const sel = this.opts.selectable;
            let head = '<thead><tr>';
            if (sel) {
                const enabled = rows.filter((r) => this.enabled(r));
                const n = enabled.filter((r) => this.selected.has(key(r))).length;
                head += `<th class="check"><input type="checkbox" data-all aria-label="Tick every visible row"${n && n === enabled.length ? ' checked' : ''}${n && n < enabled.length ? ' data-partial' : ''}></th>`;
            }
            for (const c of cols) {
                const arrow = this.sortKey === c.key ? `<span class="arrow">${this.sortDir > 0 ? '▲' : '▼'}</span>` : '';
                head += `<th class="sortable${c.num ? ' num' : ''}" data-key="${esc(c.key)}" title="Sort">${esc(c.label)}${arrow}</th>`;
            }
            head += '</tr></thead>';
            let body = '<tbody>';
            const shown = rows.slice(0, this.limit);
            shown.forEach((r, i) => {
                const k = key(r);
                const on = this.selected.has(k);
                const dis = !this.enabled(r);
                body += `<tr data-i="${i}" class="${on ? 'selected' : ''}${dis ? ' disabled' : ''}">`;
                if (sel) body += `<td class="check"><input type="checkbox" tabindex="-1"${on ? ' checked' : ''}${dis ? ' disabled' : ''}></td>`;
                for (const c of cols) {
                    const v = c.value(r);
                    const cell = c.render ? c.render(r) : esc(v ?? '');
                    body += `<td class="${c.num ? 'num' : ''}"${!c.render && v !== null && v !== undefined ? ` title="${esc(v)}"` : ''}>${cell}</td>`;
                }
                body += '</tr>';
            });
            if (rows.length > shown.length) {
                body += `<tr><td colspan="${cols.length + (sel ? 1 : 0)}" class="table-more">Showing ${int(shown.length)} of ${int(rows.length)} rows ` +
                    '<button type="button" class="link-btn" data-more>Show more</button></td></tr>';
            } else if (!rows.length) {
                body += `<tr><td colspan="${cols.length + (sel ? 1 : 0)}" class="table-more">${esc(this.opts.empty || 'No rows.')}</td></tr>`;
            }
            body += '</tbody>';
            this.table.innerHTML = head + body;
            const all = this.table.querySelector('input[data-all]');
            if (all && all.hasAttribute('data-partial')) all.indeterminate = true;
        }

        click(e) {
            if (e.target.closest('[data-more]')) {
                this.limit += 1000;
                this.render();
                return;
            }
            if (e.target.matches('input[data-all]')) {
                const enabled = this.visible.filter((r) => this.enabled(r));
                const allOn = enabled.length && enabled.every((r) => this.selected.has(this.opts.rowKey(r)));
                enabled.forEach((r) => allOn ? this.selected.delete(this.opts.rowKey(r)) : this.selected.add(this.opts.rowKey(r)));
                this.render();
                if (this.opts.onSelect) this.opts.onSelect();
                return;
            }
            const th = e.target.closest('th[data-key]');
            if (th) {
                const k = th.dataset.key;
                if (this.sortKey === k) this.sortDir = -this.sortDir;
                else { this.sortKey = k; this.sortDir = 1; }
                this.last = null;
                this.render();
                return;
            }
            const tr = e.target.closest('tr[data-i]');
            if (!tr || !this.opts.selectable || e.target.closest('a, button, input[type=text]')) return;
            const i = Number(tr.dataset.i);
            const row = this.visible[i];
            if (!row || !this.enabled(row)) return;
            const k = this.opts.rowKey(row);
            const turnOn = !this.selected.has(k);
            if (e.shiftKey && this.last !== null) {
                const [a, b] = [Math.min(this.last, i), Math.max(this.last, i)];
                for (let j = a; j <= b; j++) {
                    const r = this.visible[j];
                    if (r && this.enabled(r)) turnOn ? this.selected.add(this.opts.rowKey(r)) : this.selected.delete(this.opts.rowKey(r));
                }
            } else {
                turnOn ? this.selected.add(k) : this.selected.delete(k);
            }
            this.last = i;
            this.render();
            if (this.opts.onSelect) this.opts.onSelect();
        }

        // columns + rows of the current view, for CSV files and HTML reports
        exportRows(limit = Infinity) {
            const cols = this.columns();
            const rows = this.view().slice(0, limit).map((r) => cols.map((c) => {
                const v = c.export ? c.export(r) : c.value(r);
                return typeof v === 'number' && !Number.isFinite(v) ? null : v;
            }));
            return { cols, rows };
        }
    }

    function downloadBlob(name, blob) {
        const a = document.createElement('a');
        a.href = URL.createObjectURL(blob);
        a.download = name;
        document.body.appendChild(a);
        a.click();
        setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    }

    function downloadCsv(name, header, rows) {
        const q = (v) => {
            let s = v === null || v === undefined ? '' : String(v);
            if (/^[=+\-@]/.test(s) && Number.isNaN(Number(s))) s = `'${s}`;          // formula injection
            return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
        };
        const text = [header.map(q).join(',')].concat(rows.map((r) => r.map(q).join(','))).join('\n');
        downloadBlob(name, new Blob([text], { type: 'text/csv' }));
    }

    function tableCsv(dt, name) {
        const { cols, rows } = dt.exportRows();
        downloadCsv(name, cols.map((c) => c.label), rows);
    }

    // ------------------------------------------------------------------ capabilities
    fetch('/api/capabilities').then((r) => r.json()).then((caps) => {
        state.caps = caps;
        if (!caps.deseq2) {
            const opt = document.querySelector('#engine option[value="deseq2"]');
            opt.disabled = true;
            opt.textContent += ' (PyDESeq2 not installed)';
        }
        if (!caps.limma) {
            const opt = document.querySelector('#engine option[value="limma"]');
            opt.disabled = true;
            opt.textContent += ' (R/limma not installed)';
        }
        state.libraries = caps.pathway_libraries.slice();
        state.tfLibraries = new Set(caps.tf_libraries_default);
        renderLibraries();
        renderTfLibraries();
        updateMethodHint();
    }).catch(() => { /* the form still works with its defaults */ });

    // ------------------------------------------------------------------ 1. datasets
    const dropZone = $('drop-zone');
    const fileInput = $('file-input');
    const TABLE_RE = /\.(csv|tsv|txt|tab)(\.gz)?$/i;
    const inputMode = () => document.querySelector('input[name="input_mode"]:checked').value;

    document.querySelectorAll('input[name="input_mode"]').forEach((el) => el.addEventListener('change', () => {
        const geo = inputMode() === 'geo';
        $('geo-zone').classList.toggle('hidden', !geo);
        $('upload-zone').classList.toggle('hidden', geo);
    }));

    function setFiles(list) {
        const files = Array.from(list || []);
        const bad = files.filter((f) => !TABLE_RE.test(f.name));
        state.files = files.filter((f) => TABLE_RE.test(f.name));
        $('file-name').textContent = state.files.length
            ? (state.files.length === 1 ? state.files[0].name : `${state.files.length} files: ${state.files.slice(0, 3).map((f) => f.name).join(', ')}${state.files.length > 3 ? '…' : ''}`)
            : (bad.length ? 'Unsupported file: use .csv, .tsv or .txt (optionally .gz)' : 'No file selected');
    }

    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach((ev) => dropZone.addEventListener(ev, (e) => {
        e.preventDefault();
        e.stopPropagation();
    }));
    ['dragenter', 'dragover'].forEach((ev) => dropZone.addEventListener(ev, () => dropZone.classList.add('dragover')));
    ['dragleave', 'drop'].forEach((ev) => dropZone.addEventListener(ev, () => dropZone.classList.remove('dragover')));
    dropZone.addEventListener('drop', (e) => setFiles(e.dataTransfer.files));
    dropZone.addEventListener('click', () => fileInput.click());
    dropZone.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); fileInput.click(); }
    });
    fileInput.addEventListener('change', () => setFiles(fileInput.files));
    $('geo-id').addEventListener('keydown', (e) => { if (e.key === 'Enter') $('load-btn').click(); });

    const datasetJob = new Job('job-dataset', $('load-btn'), $('error-message'), (id, s, fresh) => {
        state.datasets.set(id, s);
        activateDataset(id, fresh);
    });

    $('load-btn').addEventListener('click', async () => {
        showError('');
        const fd = new FormData();
        fd.append('action', 'dataset');
        const mode = inputMode();
        fd.append('input_mode', mode);
        let text;
        if (mode === 'geo') {
            const acc = $('geo-id').value.trim().toUpperCase();
            if (!/^GSE\d+$/.test(acc)) return showError('Enter a GEO series accession such as GSE117769.');
            fd.append('geo', acc);
            fd.append('platform', $('platform').value);
            for (const f of $('geo-files').files) fd.append('expression', f);
            text = `Reading ${acc} from GEO`;
        } else {
            if (!state.files.length) return showError('Choose the expression file(s) first.');
            state.files.forEach((f) => fd.append('expression', f));
            if ($('meta-file').files[0]) fd.append('metadata', $('meta-file').files[0]);
            fd.append('name', $('ds-name').value);
            fd.append('platform', $('platform-up').value);
            text = 'Reading the tables';
        }
        try {
            await datasetJob.start(fd, text + '…');
        } catch (e) {
            showError(e.message);
        }
    });

    function tile(label, value, sub, small) {
        return `<div class="tile"><div class="tile-label">${esc(label)}</div>` +
            `<div class="tile-value${small ? ' small' : ''}">${value}</div><div class="tile-sub">${sub || ''}</div></div>`;
    }

    function renderDatasetBar() {
        const bar = $('dataset-bar');
        bar.classList.toggle('hidden', !state.datasets.size);
        bar.innerHTML = '<span class="hint">Datasets in this session:</span>' + [...state.datasets.entries()].map(([id, s]) => {
            const n = [...state.analyses.values()].filter((a) => a.dataset_job === id).length;
            return `<span class="group-chip${id === state.dsJob ? ' active' : ''}" title="${esc(s.title || s.name)}">` +
                `<button type="button" class="chip-open" data-ds="${id}">${esc(s.name)}</button>` +
                `<span class="count">${n ? `${n} analys${n === 1 ? 'is' : 'es'}` : ''}</span>` +
                `<button type="button" class="chip-remove" data-ds="${id}" title="Remove from the session" aria-label="Remove ${esc(s.name)}">×</button></span>`;
        }).join('');
    }

    $('dataset-bar').addEventListener('click', (e) => {
        const open = e.target.closest('.chip-open');
        const remove = e.target.closest('.chip-remove');
        if (open) activateDataset(open.dataset.ds, false);
        if (remove) {
            const id = remove.dataset.ds;
            const own = [...state.analyses.entries()].filter(([, a]) => a.dataset_job === id).map(([k]) => k);
            if (own.length && !window.confirm(`Remove this dataset and its ${own.length} analysis(es) from the session? The files stay on the server.`)) return;
            own.forEach(removeAnalysis);
            state.datasets.delete(id);
            if (state.dsJob === id) {
                const next = [...state.datasets.keys()].pop();
                if (next) activateDataset(next, false);
                else {
                    state.dsJob = null;
                    state.dataset = null;
                    $('sec-samples').classList.add('hidden');
                    $('sec-de').classList.add('hidden');
                }
            }
            renderDatasetBar();
            updateUrl();
        }
    });

    function sourceShort(s) {
        if (s.needs_expression) return 'Sample metadata only';
        if (/^GEO supplementary file /.test(s.source || '')) return 'GEO supplementary file';
        if (/^GEO per-sample files/.test(s.source || '')) return 'GEO per-sample files';
        return s.source || '–';
    }

    function sourceDetail(s) {
        const name = (s.source || '').replace(/^GEO supplementary file /, '').replace(/^GEO per-sample files /, '');
        return [name !== s.source ? name : '', s.taxid ? `taxon ${s.taxid}` : ''].filter(Boolean).join(' · ');
    }

    function activateDataset(id, fresh) {
        const s = state.datasets.get(id);
        if (!s) return;
        state.dsJob = id;
        state.dataset = s;
        document.querySelector('.container').classList.add('expanded');
        $('sec-samples').classList.remove('hidden');
        $('samples-dataset').textContent = `· ${s.name}`;
        const rows = s.samples.rows;
        const matched = rows.filter((r) => r.column).length;
        $('ds-tiles').innerHTML =
            tile('Dataset', esc(s.name), esc(s.title || s.source), true) +
            tile('Samples', `${matched}<span class="hint"> / ${rows.length}</span>`, 'with expression values') +
            tile('Genes', s.n_genes ? int(s.n_genes) : '–', esc(ID_TYPE_TEXT[s.id_type] || '')) +
            tile('Values', esc(KIND_TEXT[s.value_kind] || '–'), esc(PLATFORM_TEXT[s.platform] || '') + (s.gpl ? ` · ${esc(s.gpl)}` : ''), true) +
            tile('Source', esc(sourceShort(s)), esc(sourceDetail(s)), true);
        $('ds-notes').innerHTML = (s.notes || []).map((n) => `<li>${esc(n)}</li>`).join('');
        showError('');
        if (s.needs_expression) {
            $('geo-files-details').open = true;
            showError('This series has no usable expression table on GEO. Add its count files in step 1 ' +
                '("Count files for this series") and load it again.');
        }
        sampleColumnsInfo();
        $('groupby-col').innerHTML = '<option value="">One group per value of…</option>' +
            s.samples.columns.map((c, j) => ({ c, j })).filter((x) => sampleColInfo[x.j].distinct > 1 && sampleColInfo[x.j].distinct <= 12)
                .map((x) => `<option value="${x.j}">${esc(x.c)} (${sampleColInfo[x.j].distinct} values)</option>`).join('');
        loadGroups();
        sampleTable.selected.clear();
        $('sample-filter').value = '';
        sampleTable.filter = '';
        sampleTable.setRows(s.samples.rows);
        updateAssignRow();
        $('sec-de').classList.toggle('hidden', !!s.needs_expression);
        renderDesignOptions();
        updateMethodHint();
        updateSymbolsToggle();
        renderDatasetBar();
        updateUrl();
        if (fresh) $('sec-samples').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    // ------------------------------------------------------------------ 2. samples and groups
    let sampleColInfo = [];

    function sampleColumnsInfo() {
        const s = state.dataset.samples;
        sampleColInfo = s.columns.map((c, j) => {
            const values = new Set(s.rows.filter((r) => !r.extra && r.column).map((r) => r.values[j]));
            return { distinct: values.size };
        });
        s.rows.forEach((r) => { r._text = (r.values.join(' ') + ' ' + (r.column || '')).toLowerCase(); });
    }

    const groupOf = (id) => state.groups.find((g) => g.ids.has(id));
    const groupColor = (g) => css(GROUP_COLORS[g.color % GROUP_COLORS.length]);

    function sampleColumns() {
        const s = state.dataset.samples;
        const showConst = $('show-constant').checked;
        const cols = [{
            key: '__group', label: 'Group', value: (r) => (groupOf(r.id) || {}).name || '',
            render: (r) => {
                const g = groupOf(r.id);
                return g ? `<span class="group-cell"><span class="dot" style="background:${groupColor(g)}"></span>${esc(g.name)}</span>` : '';
            },
        }];
        s.columns.forEach((c, j) => {
            if (!showConst && sampleColInfo[j].distinct <= 1 && j > 0) return;
            cols.push({ key: `c${j}`, label: c, value: (r) => r.values[j] });
        });
        cols.push({
            key: '__column', label: 'Expression column', value: (r) => r.column || '',
            render: (r) => r.column ? esc(r.column) : '<span class="hint">not found</span>',
        });
        return cols;
    }

    const sampleTable = new DataTable($('sample-table'), {
        selectable: true,
        rowKey: (r) => r.id,
        isDisabled: (r) => !r.column,
        columns: sampleColumns,
        text: (r) => r._text,
        onSelect: updateAssignRow,
        pageSize: 2000,
        empty: 'No sample matches the filter.',
    });

    $('sample-filter').addEventListener('input', debounce(() => sampleTable.setFilter($('sample-filter').value), 150));
    $('show-constant').addEventListener('change', () => sampleTable.render());

    function groupsKey() { return `lumos-groups-${state.dsJob}`; }

    function saveGroups() {
        if (!state.dsJob) return;
        store.set(groupsKey(), {
            groups: state.groups.map((g) => ({ name: g.name, color: g.color, ids: [...g.ids] })),
            reference: $('ref-group').value, test: $('test-group').value,
        });
    }

    function loadGroups(fallback) {
        const saved = store.get(groupsKey()) || fallback;
        const known = new Set(state.dataset.samples.rows.map((r) => r.id));
        state.groups = [];
        if (saved && Array.isArray(saved.groups)) {
            state.groups = saved.groups.map((g, i) => ({
                name: String(g.name), color: g.color ?? i, ids: new Set((g.ids || []).filter((x) => known.has(x))),
            }));
        }
        renderGroups(saved ? saved.reference : null, saved ? saved.test : null);
    }

    function uniqueName(base) {
        let name = base;
        let k = 2;
        while (state.groups.some((g) => g.name === name)) name = `${base} ${k++}`;
        return name;
    }

    function freeColor() {
        const used = new Set(state.groups.map((g) => g.color % GROUP_COLORS.length));
        for (let i = 0; i < GROUP_COLORS.length; i++) if (!used.has(i)) return i;
        return state.groups.length;
    }

    function addGroup(name, ids) {
        const g = { name: uniqueName(name), color: freeColor(), ids: new Set(ids || []) };
        for (const id of g.ids) state.groups.forEach((o) => o.ids.delete(id));
        state.groups.push(g);
        return g;
    }

    function sampleRow(id) {
        if (!state._rowById || state._rowFor !== state.dataset) {
            state._rowById = new Map(state.dataset.samples.rows.map((r) => [r.id, r]));
            state._rowFor = state.dataset;
        }
        return state._rowById.get(id);
    }

    function renderGroups(ref, test) {
        if (!state.dataset) return;
        const bar = $('group-bar');
        bar.innerHTML = state.groups.map((g, i) => {
            const n = [...g.ids].length;
            return `<span class="group-chip" data-g="${i}"><span class="dot" style="background:${groupColor(g)}"></span>` +
                `<input type="text" value="${esc(g.name)}" aria-label="Group name" maxlength="60" style="width: ${Math.max(4, Math.min(30, g.name.length + 1))}ch">` +
                `<span class="count">${n}</span><button type="button" title="Delete group" aria-label="Delete group">×</button></span>`;
        }).join('') +
            `<button type="button" class="link-btn" id="new-group-btn">+ ${sampleTable.selected.size ? 'New group from ticked samples' : 'New group'}</button>`;
        $('new-group-btn').addEventListener('click', () => {
            addGroup(`Group ${state.groups.length + 1}`, [...sampleTable.selected]);
            sampleTable.selected.clear();
            afterGroupChange();
        });
        bar.querySelectorAll('.group-chip').forEach((chip) => {
            const g = state.groups[Number(chip.dataset.g)];
            const input = chip.querySelector('input');
            input.addEventListener('change', () => {
                const name = input.value.trim().slice(0, 60);
                if (name && !state.groups.some((o) => o !== g && o.name === name)) {
                    const wasRef = $('ref-group').value === g.name, wasTest = $('test-group').value === g.name;
                    g.name = name;
                    afterGroupChange(wasRef ? name : undefined, wasTest ? name : undefined);
                } else {
                    input.value = g.name;
                }
            });
            chip.querySelector('button').addEventListener('click', () => {
                state.groups.splice(state.groups.indexOf(g), 1);
                afterGroupChange();
            });
        });
        const opts = (sel) => state.groups.map((g) => {
            const n = [...g.ids].filter((id) => sampleRow(id) && sampleRow(id).column).length;
            return `<option value="${esc(g.name)}"${g.name === sel ? ' selected' : ''}>${esc(g.name)} (${n})</option>`;
        }).join('');
        const prevRef = ref !== undefined && ref !== null ? ref : $('ref-group').value;
        const prevTest = test !== undefined && test !== null ? test : $('test-group').value;
        let r = state.groups.some((g) => g.name === prevRef) ? prevRef : null;
        let t = state.groups.some((g) => g.name === prevTest) ? prevTest : null;
        if (!r) r = guessReference();
        if (!t || t === r) t = (state.groups.find((g) => g.name !== r) || {}).name;
        $('ref-group').innerHTML = opts(r);
        $('test-group').innerHTML = opts(t);
        updateAssignRow();
        updateCompareHint();
    }

    function guessReference() {
        const ctrl = /control|healthy|normal|wild|\bwt\b|untreated|vehicle|mock|baseline|ctrl|sham|naive|dmso|placebo|reference/i;
        const g = state.groups.find((x) => ctrl.test(x.name)) || state.groups[0];
        return g ? g.name : null;
    }

    function afterGroupChange(ref, test) {
        renderGroups(ref, test);
        sampleTable.render();
        saveGroups();
    }

    function updateAssignRow() {
        const n = sampleTable.selected.size;
        $('sample-sel-count').textContent = n ? `${n} sample${n === 1 ? '' : 's'} ticked` : 'No samples ticked';
        $('assign-row').classList.toggle('active', n > 0);
        $('assign-buttons').innerHTML = state.groups.map((g, i) =>
            `<button type="button" class="link-btn assign-btn" data-g="${i}"${n ? '' : ' disabled'}>` +
            `<span class="dot" style="background:${groupColor(g)}"></span>${esc(g.name)}</button>`).join('') +
            `<button type="button" class="link-btn" data-g="new"${n ? '' : ' disabled'}>+ new group</button>`;
        $('unassign-btn').disabled = !n;
        $('clear-sel-btn').disabled = !n;
        const nb = $('new-group-btn');
        if (nb) nb.textContent = `+ ${n ? 'New group from ticked samples' : 'New group'}`;
    }

    $('assign-buttons').addEventListener('click', (e) => {
        const b = e.target.closest('button[data-g]');
        if (!b || b.disabled) return;
        const ids = [...sampleTable.selected];
        if (b.dataset.g === 'new') {
            addGroup(`Group ${state.groups.length + 1}`, ids);
        } else {
            const g = state.groups[Number(b.dataset.g)];
            state.groups.forEach((o) => ids.forEach((id) => o.ids.delete(id)));
            ids.forEach((id) => g.ids.add(id));
        }
        sampleTable.selected.clear();
        afterGroupChange();
    });

    $('unassign-btn').addEventListener('click', () => {
        state.groups.forEach((g) => sampleTable.selected.forEach((id) => g.ids.delete(id)));
        sampleTable.selected.clear();
        afterGroupChange();
    });

    $('clear-sel-btn').addEventListener('click', () => {
        sampleTable.selected.clear();
        sampleTable.render();
        updateAssignRow();
    });

    $('groupby-col').addEventListener('change', () => {
        const j = $('groupby-col').value;
        if (j === '') return;
        const col = Number(j);
        const values = new Map();
        state.dataset.samples.rows.forEach((r) => {
            if (!r.column || r.extra) return;
            const v = r.values[col] || '(empty)';
            if (!values.has(v)) values.set(v, []);
            values.get(v).push(r.id);
        });
        if (state.groups.some((g) => g.ids.size) &&
            !window.confirm('Replace the current groups with one group per value of this column?')) {
            $('groupby-col').value = '';
            return;
        }
        state.groups = [];
        [...values.entries()].forEach(([v, ids]) => addGroup(v.slice(0, 60), ids));
        $('groupby-col').value = '';
        sampleTable.selected.clear();
        afterGroupChange(guessReference(), null);
    });

    function comparison() {
        const ref = state.groups.find((g) => g.name === $('ref-group').value);
        const test = state.groups.find((g) => g.name === $('test-group').value);
        const usable = (g) => g ? [...g.ids].filter((id) => sampleRow(id) && sampleRow(id).column) : [];
        return { ref, test, nRef: usable(ref).length, nTest: usable(test).length };
    }

    function updateCompareHint() {
        const c = comparison();
        let msg = '';
        if (state.groups.length < 2) {
            msg = 'Create at least two groups to compare.';
        } else if (!c.ref || !c.test || c.ref === c.test) {
            msg = 'Choose two different groups.';
        } else if (c.nRef < 2 || c.nTest < 2) {
            msg = 'Each group needs at least two samples (replicates).';
        } else {
            msg = `log2 fold change = ${c.test.name} vs ${c.ref.name}: positive values mean higher expression in ` +
                `${c.test.name} (${c.nTest} samples) than in ${c.ref.name} (${c.nRef} samples).`;
        }
        $('compare-hint').textContent = msg;
        $('de-label').placeholder = c.ref && c.test ? `${c.test.name} vs ${c.ref.name}` : 'Test vs Reference';
        refreshButtons();
    }

    $('ref-group').addEventListener('change', () => { saveGroups(); updateCompareHint(); });
    $('test-group').addEventListener('change', () => { saveGroups(); updateCompareHint(); });

    function refreshButtons() {
        const c = state.dataset ? comparison() : {};
        $('de-btn').disabled = deJob.running() || !state.dataset || state.dataset.needs_expression ||
            !c.ref || !c.test || c.ref === c.test || c.nRef < 2 || c.nTest < 2;
        $('app-btn').disabled = appJob.running() || !state.appSource;
        updateIntegrateHint();
    }

    // ------------------------------------------------------------------ 3. differential expression form
    function updateMethodHint() {
        const s = state.dataset;
        if (!s) return;
        let auto;
        if (s.value_kind === 'counts') auto = 'DESeq2, because the values are raw counts.';
        else if (!state.caps || state.caps.limma) auto = `limma, because the values are ${KIND_TEXT[s.value_kind] || 'not raw counts'}.`;
        else auto = 'the Welch t-test (exploratory), because R/limma is not installed on this server.';
        $('method-hint').textContent = `Automatic uses ${auto}`;
    }

    function updateSymbolsToggle() {
        const t = state.dataset ? state.dataset.id_type : null;
        const box = $('symbols');
        box.closest('.toggle-switch').classList.toggle('hidden', t === 'symbol');
        const names = { ensembl: 'Ensembl IDs', entrez: 'Entrez IDs', probe: 'probe IDs' };
        $('symbols-title').textContent = `Convert ${names[t] || 'gene identifiers'} to gene symbols`;
        box.disabled = t === 'unknown' && !(state.dataset && state.dataset.platform === 'microarray');
        $('unit').querySelector('option[value="symbols"]').disabled = t === 'symbol';
        if (t === 'symbol') $('unit').value = 'features';
    }

    const safeTerm = (name) => {
        let v = String(name).trim().replace(/[^A-Za-z0-9_]/g, '_').replace(/_+/g, '_').replace(/^_|_$/g, '') || 'variable';
        if (/^\d/.test(v)) v = `v_${v}`;
        return v;
    };

    function designColumns() {
        if (!state.dataset) return [];
        return state.dataset.samples.columns.map((c, j) => ({ c, j })).filter((x) => x.j > 0 && sampleColInfo[x.j] && sampleColInfo[x.j].distinct > 1);
    }

    function renderDesignOptions() {
        const cols = designColumns();
        const rows = state.dataset.samples.rows.filter((r) => r.column && !r.extra);
        const numeric = (j) => {
            const vals = rows.map((r) => String(r.values[j]).trim()).filter((v) => v !== '');
            return vals.length > 0 && vals.filter((v) => !Number.isNaN(Number(v))).length >= 0.9 * vals.length;
        };
        const opt = cols.map((x) => `<option value="${esc(x.c)}">${esc(x.c)}</option>`).join('');
        $('design-batch').innerHTML = `<option value="">None</option>${opt}`;
        $('design-pair').innerHTML = `<option value="">None</option>${opt}`;
        const pills = (name, list, none) => list.map((x) => `<label class="pill"><input type="checkbox" name="${name}" value="${esc(x.c)}">${esc(x.c)}</label>`).join('') ||
            `<span class="hint">${none}</span>`;
        // one value per sample identifies the samples: offered as subject/pair, not as a covariate
        $('design-categorical').innerHTML = pills('cat', cols.filter((x) => sampleColInfo[x.j].distinct < rows.length && !numeric(x.j)),
            'No categorical column varies between samples.');
        $('design-continuous').innerHTML = pills('cont', cols.filter((x) => numeric(x.j)), 'No numeric column.');
        updateFormula();
    }

    function designChoice() {
        const pick = (name) => [...document.querySelectorAll(`input[name="${name}"]:checked`)].map((el) => el.value);
        return { batch: $('design-batch').value, pair: $('design-pair').value, covariates: pick('cat'), continuous: pick('cont') };
    }

    function updateFormula() {
        const d = designChoice();
        const terms = [d.batch, d.pair, ...d.covariates, ...d.continuous].filter(Boolean).map(safeTerm);
        $('formula-preview').textContent = `~ ${[...terms, 'condition'].join(' + ')}`;
    }
    ['design-batch', 'design-pair', 'design-categorical', 'design-continuous'].forEach((id) =>
        $(id).addEventListener('change', updateFormula));

    const deJob = new Job('job-de', $('de-btn'), $('error-message'), (id, s, fresh) => addAnalysis(id, s, fresh));
    const intJob = new Job('job-int', $('int-btn'), $('int-error'), (id, s, fresh) => addIntegration(id, s, fresh));
    const appJob = new Job('job-app', $('app-btn'), $('app-error'), (id, s, fresh) => showApp(id, s, fresh));

    $('de-btn').addEventListener('click', async () => {
        const c = comparison();
        const fd = new FormData();
        fd.append('action', 'de');
        fd.append('dataset', state.dsJob);
        const groups = {};
        state.groups.forEach((g) => { groups[g.name] = [...g.ids]; });
        fd.append('groups', JSON.stringify(groups));
        fd.append('reference', c.ref.name);
        fd.append('test', c.test.name);
        fd.append('label', $('de-label').value.trim());
        for (const k of ['engine', 'unit', 'transform', 'normalization']) fd.append(k, $(k).value);
        fd.append('min_count', $('min-count').value);
        fd.append('fdr', $('fdr').value);
        fd.append('lfc', $('lfc').value);
        fd.append('design', JSON.stringify(designChoice()));
        fd.append('symbols', $('symbols').checked && !$('symbols').disabled);
        try {
            await deJob.start(fd, `${$('de-label').value.trim() || `${c.test.name} vs ${c.ref.name}`}…`);
        } catch (e) {
            showError(e.message);
        }
        refreshButtons();
    });

    // ------------------------------------------------------------------ analyses in the session
    function buildGenes(s) {
        const g = s.genes;
        const n = g.id.length;
        const arr = (a) => Float64Array.from(a, (v) => (v === null ? NaN : v));
        const G = {
            n, id: g.id, sym: g.symbol, lfc: arr(g.lfc), fdr: arr(g.fdr), p: arr(g.p), mean: arr(g.mean),
            byName: new Map(), bySym: new Map(), byId: new Map(),
        };
        for (let i = 0; i < n; i++) {
            G.byId.set(g.id[i], i);
            G.byName.set(g.id[i].toUpperCase(), i);
            if (g.symbol[i]) {
                const u = g.symbol[i].toUpperCase();
                if (!G.bySym.has(u)) G.bySym.set(u, i);
                if (!G.byName.has(u)) G.byName.set(u, i);
            }
        }
        return G;
    }

    const geneCache = new Map();
    function genesOf(id) {
        if (!geneCache.has(id)) geneCache.set(id, buildGenes(state.analyses.get(id)));
        return geneCache.get(id);
    }

    const symbolsAvailable = (s) => !!(s && s.symbols && s.symbols.available);
    const labelIn = (G, s, i) => (state.labels === 'symbol' && symbolsAvailable(s) && G.sym[i]) ? G.sym[i] : G.id[i];
    const label = (i) => labelIn(state.G, state.de, i);
    const keyIn = (G, i) => (G.sym[i] || G.id[i]).toUpperCase();          // what Enrichr, STRING and integrations use

    function labelOfSymbol(sym) {
        // nodes come back as symbols; show the original ID when the gene source is an analysis and labels = ID
        const src = state.app ? state.app.source : null;
        if (src && src.kind === 'de' && state.analyses.has(src.job)) {
            const G = genesOf(src.job);
            const i = G.bySym.get(String(sym).toUpperCase());
            if (i !== undefined) return labelIn(G, state.analyses.get(src.job), i);
        }
        return sym;
    }

    function analysisName(id) {
        const s = state.analyses.get(id);
        return s ? `${s.label} · ${s.name}` : id;
    }

    function renderAnalysisSelect() {
        const byDs = new Map();
        state.analyses.forEach((s, id) => {
            if (!byDs.has(s.name)) byDs.set(s.name, []);
            byDs.get(s.name).push([id, s]);
        });
        $('analysis-select').innerHTML = [...byDs.entries()].map(([ds, list]) =>
            `<optgroup label="${esc(ds)}">${list.map(([id, s]) =>
                `<option value="${id}"${id === state.deJob ? ' selected' : ''}>${esc(s.label)} — ${s.counts.up} up, ${s.counts.down} down</option>`).join('')}</optgroup>`).join('');
    }
    $('analysis-select').addEventListener('change', () => activateAnalysis($('analysis-select').value, false));

    function addAnalysis(id, s, fresh) {
        state.analyses.set(id, s);
        geneCache.delete(id);
        renderDatasetBar();
        activateAnalysis(id, fresh);
        intTable.setRows([...state.analyses.keys()]);
        renderAppSources();
    }

    function removeAnalysis(id) {
        state.analyses.delete(id);
        geneCache.delete(id);
        intTable.selected.delete(id);
        intTable.setRows([...state.analyses.keys()]);
        if (state.deJob === id) {
            const next = [...state.analyses.keys()].pop();
            if (next) activateAnalysis(next, false);
            else { state.deJob = null; state.de = null; $('results').style.display = 'none'; }
        }
        renderAppSources();
    }

    function activateAnalysis(id, fresh) {
        const s = state.analyses.get(id);
        if (!s) return;
        state.deJob = id;
        state.de = s;
        state.G = genesOf(id);
        state.fdr = s.parameters.fdr;
        state.lfc = s.parameters.lfc;
        $('live-fdr').value = state.fdr;
        $('live-lfc').value = state.lfc;
        state.highlight = [];
        $('highlight').value = '';
        const symOk = symbolsAvailable(s);
        $('label-switch').classList.toggle('hidden', !symOk);
        $('id-label-name').textContent = ID_NAMES[s.id_type] || 'Original ID';
        if (!symOk) state.labels = 'id';
        else if (fresh) state.labels = 'symbol';
        document.querySelector(`input[name="labels"][value="${state.labels}"]`).checked = true;
        const cmp = s.comparison;
        $('results-title').textContent = `${s.name} · ${s.label}` + (s.label !== `${cmp.test} vs ${cmp.reference}` ? ` (${cmp.test} vs ${cmp.reference})` : '');
        $('zip-btn').href = `/api/jobs/${id}/files/${encodeURIComponent(s.files.zip)}?download=1`;
        $('manifest-btn').href = `/api/jobs/${id}/files/manifest.json?download=1`;
        $('results').style.display = 'block';
        document.querySelector('.container').classList.add('expanded');
        if (state.dsJob !== s.dataset_job && state.datasets.has(s.dataset_job)) activateDataset(s.dataset_job, false);
        computeStatus();
        renderAnalysisSelect();
        renderDE();
        if (state.appSource && state.appSource.kind === 'de' && state.appSource.job === id) renderAppGenes(true);
        renderFiles();
        updateUrl();
        refreshButtons();
        if (fresh) {
            switchTab('ddea');
            $('results').scrollIntoView({ behavior: 'smooth', block: 'start' });
        }
    }

    document.querySelectorAll('input[name="labels"]').forEach((el) => el.addEventListener('change', () => {
        state.labels = el.value;
        renderDE();
        if (state.appSource) renderAppGenes(false);
        if (state.app) renderAppCharts();
        if (state.stringData) drawString();
        if (state.pg) { drawPG(); renderPGView(); }
    }));

    function computeStatus() {
        const G = state.G;
        const st = new Int8Array(G.n);
        for (let i = 0; i < G.n; i++) {
            if (G.fdr[i] < state.fdr) {
                if (G.lfc[i] >= state.lfc) st[i] = 1;
                else if (G.lfc[i] <= -state.lfc) st[i] = -1;
            }
        }
        state.status = st;
    }

    const onThresholds = debounce(() => {
        const f = Number($('live-fdr').value), l = Number($('live-lfc').value);
        if (!(f > 0 && f <= 1) || !(l >= 0)) return;
        state.fdr = f;
        state.lfc = l;
        computeStatus();
        renderDE();
        if (state.appSource && state.appSource.kind === 'de' && state.appSource.job === state.deJob) {
            $('app-fdr').value = f;
            $('app-lfc').value = l;
            buildAppSource(false);
        }
    }, 300);
    $('live-fdr').addEventListener('input', onThresholds);
    $('live-lfc').addEventListener('input', onThresholds);
    $('heat-n').addEventListener('input', debounce(() => renderHeatmap(), 300));
    $('cluster-rows').addEventListener('change', () => renderHeatmap());
    $('cluster-cols').addEventListener('change', () => renderHeatmap());
    $('highlight').addEventListener('input', debounce(() => {
        state.highlight = $('highlight').value.split(/[\s,;]+/).filter(Boolean);
        renderVolcano();
    }, 300));

    function toggleHighlight(name) {
        const u = name.toUpperCase();
        const k = state.highlight.findIndex((h) => h.toUpperCase() === u);
        if (k >= 0) state.highlight.splice(k, 1);
        else state.highlight.push(name);
        $('highlight').value = state.highlight.join(' ');
        renderVolcano();
    }

    // ------------------------------------------------------------------ chart helpers
    const PLOT_CONFIG = {
        displaylogo: false, responsive: true, displayModeBar: 'hover',
        modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d', 'toggleSpikelines',
            'hoverClosestCartesian', 'hoverCompareCartesian'],
        toImageButtonOptions: { format: 'png', scale: 2 },
    };

    function axis(title, extra = {}) {
        return Object.assign({
            title: { text: title, standoff: 8 },
            gridcolor: css('--chart-grid'),
            linecolor: css('--chart-axis'),
            zerolinecolor: css('--chart-axis'),
            tickcolor: css('--chart-axis'),
            ticks: 'outside',
            ticklen: 4,
            automargin: true,
        }, extra);
    }

    function layout(xTitle, yTitle, extra = {}) {
        const l = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { family: 'Inter, system-ui, sans-serif', size: 12, color: css('--chart-text') },
            margin: { l: 56, r: 14, t: 12, b: 44 },
            xaxis: axis(xTitle, extra.xaxis || {}),
            yaxis: axis(yTitle, extra.yaxis || {}),
            hoverlabel: { font: { family: 'Inter, system-ui, sans-serif' } },
            showlegend: false,
            bargap: 0.2,
        };
        const rest = Object.assign({}, extra);
        delete rest.xaxis;
        delete rest.yaxis;
        return Object.assign(l, rest);
    }

    function legendTop() {
        return { showlegend: true, legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom', bgcolor: 'rgba(0,0,0,0)' },
            margin: { l: 56, r: 14, t: 34, b: 44 } };
    }

    function emptyFig(text) {
        return { data: [], layout: layout('', '', {
            xaxis: { visible: false }, yaxis: { visible: false },
            annotations: [{ text, showarrow: false, xref: 'paper', yref: 'paper', x: 0.5, y: 0.5, font: { color: css('--text-muted') } }],
        }) };
    }

    const plotted = new Set();
    function show(id, fig) {
        const el = $(id);
        if (!el) return;
        if (!window.Plotly) {
            el.textContent = 'Chart library unavailable. Check /static/vendor/plotly.min.js.';
            console.error('Lumos: Plotly did not load for', id);
            return;
        }
        try {
            if (fig.height) el.style.height = `${fig.height}px`;
            Promise.resolve(Plotly.react(el, fig.data, fig.layout, PLOT_CONFIG)).catch((e) => {
                console.error('Lumos chart failed:', id, e);
                el.textContent = `Chart rendering failed: ${e.message}`;
            });
            plotted.add(id);
        } catch (e) {
            console.error('Lumos chart failed:', id, e);
            el.textContent = `Chart rendering failed: ${e.message}`;
        }
    }

    function resizeAll() {
        for (const id of plotted) {
            const el = $(id);
            if (el && el.offsetParent !== null) Plotly.Plots.resize(el);
        }
        for (const net of [state.stringNet, state.pgNet]) if (net) net.redraw();
    }

    function uniqueLabels(indices) {
        const seen = new Map();
        indices.forEach((i) => { const l = label(i); seen.set(l, (seen.get(l) || 0) + 1); });
        return indices.map((i) => (seen.get(label(i)) > 1 ? `${label(i)} (${state.G.id[i]})` : label(i)));
    }

    // ------------------------------------------------------------------ DDEA results
    function renderDE() {
        const s = state.de;
        const st = state.status;
        let up = 0, down = 0;
        for (let i = 0; i < st.length; i++) { if (st[i] > 0) up++; else if (st[i] < 0) down++; }
        const c = s.counts;
        const design = s.design && s.design.formula && s.design.formula !== '~condition' ? s.design.formula : '~ condition';
        $('de-tiles').innerHTML =
            tile('Genes tested', int(c.tested), `${int(c.filtered_out)} removed by the filter`) +
            tile('Up-regulated', `<span class="status-up">${int(up)}</span>`, `higher in ${esc(s.comparison.test)}`) +
            tile('Down-regulated', `<span class="status-down">${int(down)}</span>`, `lower in ${esc(s.comparison.test)}`) +
            tile('Method', esc(s.method_text), `model ${esc(design)} · ${esc(s.unit_text || '')}`, true) +
            (s.symbols.source && s.id_type !== 'symbol'
                ? tile('Gene symbols', `${Math.round(100 * s.symbols.mapped / Math.max(1, s.symbols.total))}%`,
                    `${int(s.symbols.mapped)} of ${int(s.symbols.total)} from ${esc(s.symbols.source)}`)
                : '');
        $('de-notes').innerHTML = (s.notes || []).map((n) => `<li>${esc(n)}</li>`).join('');
        try { renderDECharts(); }
        catch (e) { console.error('Lumos DDEA visualization error:', e); showError(`Visualization error: ${e.message}`, $('error-message')); }
        degTable.render();
    }

    function renderDECharts() {
        renderVolcano();
        show('chart-ma', figMA());
        show('chart-pca', figPCA());
        show('chart-up', figTop(1));
        show('chart-down', figTop(-1));
        renderHeatmap();
        $('lib-card').classList.toggle('hidden', !state.de.library_sizes);
        if (state.de.library_sizes) show('chart-lib', figLib());
    }

    function statusTraces(xOf, yOf, hover, small) {
        const G = state.G, st = state.status;
        const lists = { 0: [], '-1': [], 1: [] };
        for (let i = 0; i < G.n; i++) {
            const x = xOf(i), y = yOf(i);
            if (Number.isFinite(x) && Number.isFinite(y)) lists[st[i]].push(i);
        }
        const tr = (list, name, color, size) => ({
            type: 'scattergl', mode: 'markers', name, x: list.map(xOf), y: list.map(yOf), customdata: list,
            text: list.map(label), hovertemplate: hover, marker: { color, size, line: { width: 0 } },
        });
        return [tr(lists[0], 'Not significant', css('--chart-muted'), 4),
            tr(lists[-1], `Down (${int(lists[-1].length)})`, css('--chart-1'), small ? 4 : 5),
            tr(lists[1], `Up (${int(lists[1].length)})`, css('--chart-2'), small ? 4 : 5)];
    }

    function figVolcano() {
        const G = state.G;
        const y = (i) => -Math.log10(Math.max(G.fdr[i], 1e-300));
        const data = statusTraces((i) => G.lfc[i], y, '<b>%{text}</b><br>log2FC %{x:.2f}<br>−log10 FDR %{y:.2f}<extra></extra>');
        const hl = state.highlight.map((h) => G.byName.get(h.toUpperCase())).filter((i) => i !== undefined && Number.isFinite(G.fdr[i]));
        if (hl.length) {
            data.push({
                type: 'scatter', mode: 'markers+text', name: 'Highlighted', x: hl.map((i) => G.lfc[i]), y: hl.map(y),
                text: hl.map(label), textposition: 'top center', customdata: hl, hoverinfo: 'skip',
                textfont: { color: css('--chart-ink'), size: 11 },
                marker: { size: 10, color: 'rgba(0,0,0,0)', line: { width: 2, color: css('--chart-ink') } }, showlegend: false,
            });
        }
        const ink = css('--chart-axis');
        const shapes = [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: -Math.log10(state.fdr), y1: -Math.log10(state.fdr), line: { color: ink, width: 1, dash: 'dash' } }];
        if (state.lfc > 0) {
            for (const v of [-state.lfc, state.lfc]) shapes.push({ type: 'line', yref: 'paper', y0: 0, y1: 1, x0: v, x1: v, line: { color: ink, width: 1, dash: 'dash' } });
        }
        const cmp = state.de.comparison;
        return { data, layout: layout(`log2 fold change (${cmp.test} vs ${cmp.reference})`, '−log10 adjusted p-value',
            Object.assign(legendTop(), { shapes, hovermode: 'closest' })) };
    }

    function renderVolcano() {
        show('chart-volcano', figVolcano());
        $('volcano-sub').textContent = `Dashed lines: FDR = ${state.fdr}` + (state.lfc > 0 ? ` and |log2FC| = ${state.lfc}` : '') +
            '. Click a point to label it.' + (state.de.counts.no_fdr ? ` ${int(state.de.counts.no_fdr)} genes without an adjusted p-value are not shown.` : '');
        const el = $('chart-volcano');
        if (!el._lumosClick) {
            el.on('plotly_click', (ev) => {
                const pt = ev.points && ev.points[0];
                if (pt && pt.customdata !== undefined) toggleHighlight(label(pt.customdata));
            });
            el._lumosClick = true;
        }
    }

    function figMA() {
        const G = state.G;
        const logx = state.de.mean_scale === 'log';
        $('ma-sub').textContent = logx ? 'DESeq2 base mean (log scale) against the fold change.'
            : `Mean of the ${state.de.transform_text} values against the fold change.`;
        const data = statusTraces((i) => (logx && !(G.mean[i] > 0) ? NaN : G.mean[i]), (i) => G.lfc[i],
            '<b>%{text}</b><br>mean %{x:.3g}<br>log2FC %{y:.2f}<extra></extra>', true);
        return { data, layout: layout(logx ? 'mean of normalised counts' : 'mean expression', 'log2 fold change',
            { xaxis: logx ? { type: 'log' } : {}, hovermode: 'closest' }) };
    }

    function figPCA() {
        const s = state.de;
        const p = s.pca;
        if (!p) return emptyFig('Too few genes or samples for a PCA.');
        const cmp = s.comparison;
        const ids = cmp.reference_samples.concat(cmp.test_samples);
        const titles = cmp.sample_labels || ids;
        const nRef = cmp.reference_samples.length;
        const tr = (name, from, to, color) => ({
            type: 'scatter', mode: 'markers', name, x: p.pcs[0].slice(from, to), y: p.pcs[1].slice(from, to),
            text: ids.slice(from, to).map((id, k) => `${id}<br>${titles[from + k]}`),
            hovertemplate: '%{text}<extra></extra>', marker: { size: 9, color, line: { width: 1, color: css('--panel-bg') } },
        });
        $('pca-sub').textContent = `${p.n_genes} most variable genes, ${s.transform_text}.`;
        return { data: [tr(cmp.reference, 0, nRef, css('--chart-1')), tr(cmp.test, nRef, ids.length, css('--chart-2'))],
            layout: layout(`PC1 (${(100 * p.variance[0]).toFixed(1)}%)`, `PC2 (${(100 * p.variance[1]).toFixed(1)}%)`, legendTop()) };
    }

    function topIndices(sign, n) {
        const out = [];
        const st = state.status;
        for (let i = 0; i < st.length && out.length < n; i++) if (st[i] === sign) out.push(i);
        return out;
    }

    function figTop(sign) {
        const G = state.G;
        const idx = topIndices(sign, 20);
        if (!idx.length) return Object.assign(emptyFig('No genes with these thresholds.'), { height: 150 });
        return {
            height: Math.max(150, 19 * idx.length + 70),
            data: [{
                type: 'bar', orientation: 'h', x: idx.map((i) => G.lfc[i]), y: uniqueLabels(idx),
                marker: { color: css(sign > 0 ? '--chart-2' : '--chart-1') }, customdata: idx.map((i) => fmtP(G.fdr[i])),
                hovertemplate: '<b>%{y}</b><br>log2FC %{x:.2f}<br>FDR %{customdata}<extra></extra>',
            }],
            layout: layout('log2 fold change', '', { yaxis: { autorange: 'reversed', ticks: '', type: 'category' }, margin: { l: 90, r: 14, t: 8, b: 44 } }),
        };
    }

    function heatmapRows() {
        const s = state.de, G = state.G;
        const n = Math.max(5, Math.min(300, Number($('heat-n').value) || 50));
        const rows = [];
        s.heatmap.genes.forEach((gid, k) => {
            const i = G.byId.get(gid);
            if (i !== undefined && state.status[i] !== 0 && rows.length < n) rows.push({ i, gid, v: s.heatmap.values[k] });
        });
        return rows;
    }

    const clusterCache = new Map();
    let heatOrder = null;
    async function renderHeatmap() {
        if (!state.de) return;
        const rows = heatmapRows();
        const wantRows = $('cluster-rows').checked, wantCols = $('cluster-cols').checked;
        heatOrder = null;
        if (rows.length > 2 && (wantRows || wantCols)) {
            const key = `${state.deJob}|${wantRows}|${wantCols}|${rows.map((r) => r.gid).join(',')}`;
            if (!clusterCache.has(key)) {
                try {
                    const r = await fetch(`/api/jobs/${state.deJob}/cluster`, {
                        method: 'POST', headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ genes: rows.map((x) => x.gid), rows: wantRows, cols: wantCols }),
                    });
                    clusterCache.set(key, r.ok ? await r.json() : null);
                } catch (e) {
                    clusterCache.set(key, null);
                }
            }
            heatOrder = clusterCache.get(key);
        }
        show('chart-heatmap', figHeatmap(rows, heatOrder));
    }

    function figHeatmap(rows, order) {
        const s = state.de;
        if (!rows.length) {
            $('heatmap-sub').textContent = '';
            return Object.assign(emptyFig('No differentially expressed genes with these thresholds.'), { height: 200 });
        }
        const clusteredRows = order && order.rows && order.rows.length === rows.length && $('cluster-rows').checked;
        if (clusteredRows) {
            const pos = new Map(order.rows.map((g, k) => [g, k]));
            rows = rows.slice().sort((a, b) => pos.get(a.gid) - pos.get(b.gid));
        }
        const cmp = s.comparison;
        let colIdx = cmp.reference_samples.concat(cmp.test_samples).map((_, k) => k);
        const clusteredCols = order && order.cols && order.cols.length === colIdx.length && $('cluster-cols').checked;
        if (clusteredCols) colIdx = order.cols;
        const z = rows.map((r) => {
            const vals = r.v.map((x) => (x === null ? NaN : x));
            const ok = vals.filter(Number.isFinite);
            const m = ok.reduce((a, b) => a + b, 0) / ok.length;
            const sd = Math.sqrt(ok.reduce((a, b) => a + (b - m) ** 2, 0) / ok.length) || 1;
            return colIdx.map((k) => (Number.isFinite(vals[k]) ? (vals[k] - m) / sd : null));
        });
        const allSamples = cmp.reference_samples.concat(cmp.test_samples);
        const nRef = cmp.reference_samples.length;
        const samples = colIdx.map((k) => allSamples[k]);
        const groupName = colIdx.map((k) => (k < nRef ? cmp.reference : cmp.test));
        const names = uniqueLabels(rows.map((r) => r.i));
        const light = document.body.classList.contains('light-mode');
        const extra = {
            xaxis: { showticklabels: samples.length <= 60, tickangle: -60, ticks: '', type: 'category', showgrid: false },
            yaxis: { autorange: 'reversed', ticks: '', type: 'category', showgrid: false },
            margin: { l: 100, r: 14, t: 26, b: samples.length <= 60 ? 90 : 20 },
        };
        if (!clusteredCols) {
            extra.shapes = [{ type: 'line', xref: 'x', yref: 'paper', x0: nRef - 0.5, x1: nRef - 0.5, y0: 0, y1: 1, line: { color: css('--chart-ink'), width: 2 } }];
            extra.annotations = [
                { xref: 'x', yref: 'paper', x: (nRef - 1) / 2, y: 1.0, yanchor: 'bottom', text: esc(cmp.reference), showarrow: false, font: { color: css('--chart-text') } },
                { xref: 'x', yref: 'paper', x: nRef + (cmp.test_samples.length - 1) / 2, y: 1.0, yanchor: 'bottom', text: esc(cmp.test), showarrow: false, font: { color: css('--chart-text') } },
            ];
        }
        $('heatmap-sub').textContent = `${rows.length} differentially expressed genes with the lowest p-values; ` +
            `row z-scores of the ${s.transform_text} values` + (clusteredRows || clusteredCols ? '; Ward clustering.' : '.');
        return {
            height: Math.max(260, 16 * rows.length + 110),
            data: [{
                type: 'heatmap', z, x: samples, y: names, zmid: 0, zmin: -2.5, zmax: 2.5,
                colorscale: [[0, css('--chart-1')], [0.5, light ? '#f4f4f2' : '#22252c'], [1, css('--chart-2')]],
                customdata: rows.map(() => groupName),
                hovertemplate: '<b>%{y}</b><br>%{x} (%{customdata})<br>z = %{z:.2f}<extra></extra>',
                colorbar: { title: { text: 'z-score', side: 'right' }, thickness: 12, outlinewidth: 0 },
            }],
            layout: layout('', '', extra),
        };
    }

    function figLib() {
        const s = state.de;
        const cmp = s.comparison;
        const samples = cmp.reference_samples.concat(cmp.test_samples);
        const nRef = cmp.reference_samples.length;
        return { data: [{
            type: 'bar', x: samples, y: s.library_sizes,
            marker: { color: samples.map((_, k) => css(k < nRef ? '--chart-1' : '--chart-2')) },
            customdata: samples.map((_, k) => (k < nRef ? cmp.reference : cmp.test)),
            hovertemplate: '%{x} (%{customdata})<br>%{y:,.0f} reads<extra></extra>',
        }], layout: layout('', 'reads', { xaxis: { showticklabels: samples.length <= 60, tickangle: -60, type: 'category' } }) };
    }

    // gene table of the analysis shown
    const degTable = new DataTable($('deg-table'), {
        columns: () => {
            const G = state.G, s = state.de;
            const st = (i) => state.status[i];
            const cols = [{ key: 'label', label: 'Gene', value: label }];
            if (symbolsAvailable(s)) {
                cols.push(state.labels === 'symbol'
                    ? { key: 'id', label: ID_NAMES[s.id_type] || 'ID', value: (i) => G.id[i] }
                    : { key: 'sym', label: 'Gene symbol', value: (i) => G.sym[i] || '' });
            }
            cols.push({
                key: 'status', label: 'Status', value: st, export: (i) => (st(i) > 0 ? 'Up-regulated' : st(i) < 0 ? 'Down-regulated' : 'Not significant'),
                render: (i) => (st(i) > 0 ? '<span class="status-up">Up</span>' : st(i) < 0 ? '<span class="status-down">Down</span>' : '<span class="hint">–</span>'),
            });
            cols.push({ key: 'lfc', label: 'log2FC', num: true, value: (i) => G.lfc[i], render: (i) => fmt(G.lfc[i]) });
            cols.push({ key: 'fdr', label: 'FDR', num: true, value: (i) => G.fdr[i], render: (i) => fmtP(G.fdr[i]) });
            cols.push({ key: 'p', label: 'p-value', num: true, value: (i) => G.p[i], render: (i) => fmtP(G.p[i]) });
            cols.push({ key: 'mean', label: s.mean_scale === 'log' ? 'Base mean' : 'Mean', num: true, value: (i) => G.mean[i], render: (i) => fmt(G.mean[i], 1) });
            return cols;
        },
        rowKey: (i) => i,
        text: (i) => `${state.G.id[i]} ${state.G.sym[i] || ''}`.toLowerCase(),
        empty: 'No genes with these thresholds.',
    });
    const degRender = degTable.render.bind(degTable);
    degTable.render = () => {
        if (!state.G) return;
        const all = $('deg-all').checked;
        const rows = [];
        for (let i = 0; i < state.G.n; i++) if (all || state.status[i] !== 0) rows.push(i);
        degTable.rows = rows;
        degRender();
    };
    $('deg-all').addEventListener('change', () => { degTable.limit = 300; degTable.render(); });
    $('deg-filter').addEventListener('input', debounce(() => degTable.setFilter($('deg-filter').value), 200));
    $('deg-export').addEventListener('click', () => tableCsv(degTable, `lumos_${state.de.name}_${state.de.label}_genes.csv`.replace(/[^A-Za-z0-9_.-]+/g, '_')));

    // ------------------------------------------------------------------ tabs
    function switchTab(name) {
        document.querySelectorAll('.tab-btn').forEach((b) => b.classList.toggle('active', b.dataset.tab === name));
        document.querySelectorAll('.tab-pane').forEach((p) => p.classList.toggle('active', p.id === `pane-${name}`));
        if (name === 'pg' && state.app && !state.pg) loadPG();
        if (name === 'app') loadLibraryList();
        setTimeout(resizeAll, 30);
        if (name === 'pg' && state.pgNet) setTimeout(() => state.pgNet.fit(), 60);
    }
    $('de-to-app').addEventListener('click', () => {
        if (!state.deJob || !state.analyses.has(state.deJob)) {
            return showError('Finish a DDEA analysis before opening APP.');
        }
        renderAppSources();
        const value = `de|${state.deJob}|analysis`;
        const selector = $('app-source');
        if ([...selector.options].some((o) => o.value === value)) {
            selector.value = value;
            buildAppSource(true);
        }
        switchTab('app');
        if (!state.appSource || !state.appSource.rows.length) {
            showError('No differentially expressed genes meet the APP filters. Review FDR and |log2FC|, or test another comparison.', $('app-error'));
        } else {
            showError('', $('app-error'));
        }
    });
    const activeTab = () => document.querySelector('.tab-btn.active').dataset.tab;
    document.querySelectorAll('.tab-btn').forEach((b) => b.addEventListener('click', () => {
        if (!b.disabled) switchTab(b.dataset.tab);
    }));

    // ------------------------------------------------------------------ Integrate
    const intTable = new DataTable($('int-table'), {
        selectable: true,
        rowKey: (id) => id,
        columns: () => [
            { key: 'label', label: 'Analysis', value: (id) => state.analyses.get(id).label },
            { key: 'ds', label: 'Dataset', value: (id) => state.analyses.get(id).name },
            { key: 'cmp', label: 'Comparison', value: (id) => { const c = state.analyses.get(id).comparison; return `${c.test} vs ${c.reference}`; } },
            { key: 'method', label: 'Method', value: (id) => state.analyses.get(id).method_text },
            { key: 'up', label: 'Up', num: true, value: (id) => state.analyses.get(id).counts.up },
            { key: 'down', label: 'Down', num: true, value: (id) => state.analyses.get(id).counts.down },
            { key: 'tested', label: 'Genes tested', num: true, value: (id) => state.analyses.get(id).counts.tested },
        ],
        text: (id) => analysisName(id).toLowerCase(),
        onSelect: updateIntegrateHint,
        empty: 'Run differential expression analyses first.',
    });

    const intMode = () => document.querySelector('input[name="int_mode"]:checked').value;
    document.querySelectorAll('input[name="int_mode"]').forEach((el) => el.addEventListener('change', updateIntegrateHint));

    function updateIntegrateHint() {
        const mode = intMode();
        document.querySelectorAll('.int-meta').forEach((el) => el.classList.toggle('hidden', mode !== 'meta'));
        const ids = [...intTable.selected].filter((id) => state.analyses.has(id));
        const distinct = new Set(ids.map((id) => state.analyses.get(id).dataset_job)).size;
        let msg = '', ok = ids.length >= 2;
        if (ids.length < 2) msg = 'Tick at least two analyses.';
        else if (mode === 'contrasts' && distinct > 1) { msg = 'Contrasts must come from the same dataset; for different datasets choose the meta-analysis or the multi-tissue comparison.'; ok = false; }
        else if (mode !== 'contrasts' && distinct < ids.length) { msg = 'Each analysis must come from a different, independent dataset; for contrasts of one dataset choose "Contrasts of one dataset".'; ok = false; }
        else if (mode === 'meta') msg = 'Use datasets of the same tissue with comparable contrasts; the effects are combined with a random-effects model.';
        else if (mode === 'tissues') msg = 'Each tissue keeps its own effect; genes are classified as conserved, shared, specific or discordant.';
        else msg = 'The contrasts may share samples, so their effects are aligned but not pooled.';
        $('int-hint').textContent = msg;
        $('int-btn').disabled = !ok || intJob.running();
    }

    $('int-btn').addEventListener('click', async () => {
        const fd = new FormData();
        fd.append('action', 'integrate');
        fd.append('analyses', [...intTable.selected].join(','));
        fd.append('mode', intMode());
        fd.append('fdr', $('int-fdr').value);
        fd.append('lfc', $('int-lfc').value);
        fd.append('estimator', $('int-estimator').value);
        fd.append('min_studies', $('int-min').value);
        try {
            await intJob.start(fd, `${MODE_NAMES[intMode()]}…`);
        } catch (e) {
            showError(e.message, $('int-error'));
        }
    });

    function addIntegration(id, s, fresh) {
        state.integrations.set(id, s);
        activateIntegration(id, fresh);
        renderAppSources();
        updateUrl();
    }

    const metaStatus = (k) => {
        const T = state.int.table, th = state.int.thresholds;
        const f = num(T.fdr[k]), l = num(T.lfc[k]);
        return f < th.fdr ? (l >= th.lfc ? 1 : l <= -th.lfc ? -1 : 0) : 0;
    };

    function activateIntegration(id, fresh) {
        const s = state.integrations.get(id);
        state.intJob = id;
        state.int = s;
        const T = s.table;
        $('int-results').classList.remove('hidden');
        $('int-title').textContent = `${MODE_NAMES[s.mode]}: ${s.profiles.map((p) => p.label).join(' · ')}`;
        const meta = s.mode === 'meta';
        document.querySelectorAll('.int-align').forEach((el) => el.classList.toggle('hidden', meta));
        document.querySelectorAll('.int-meta-results').forEach((el) => el.classList.toggle('hidden', !meta));
        const pc = s.pattern_counts;
        if (meta) {
            const h = s.heterogeneity || {};
            $('int-tiles').innerHTML = tile('Datasets', s.profiles.length, esc(s.estimator)) +
                tile('Genes analysed', int(s.n_genes), `in at least ${s.min_studies} studies`) +
                tile('Up-regulated', `<span class="status-up">${int(pc['Up-regulated'] || 0)}</span>`, `FDR < ${s.thresholds.fdr}, |log2FC| ≥ ${s.thresholds.lfc}`) +
                tile('Down-regulated', `<span class="status-down">${int(pc['Down-regulated'] || 0)}</span>`, 'random-effects estimate') +
                tile('Median I²', h.median_i2 === null || h.median_i2 === undefined ? '–' : `${fmt(h.median_i2, 1)}%`, `${int(h.high_i2 || 0)} genes with I² > 50%`);
            const top = T.symbol.map((_, k) => k).filter((k) => metaStatus(k) !== 0).slice(0, 100);
            $('forest-gene').innerHTML = (top.length ? top : T.symbol.map((_, k) => k).slice(0, 100))
                .map((k) => `<option value="${k}">${esc(T.symbol[k])}</option>`).join('');
        } else {
            const n = s.profiles.length;
            const all = T.n_sig.filter((x) => x === n).length;
            const same = (pc['Consistent Up (all contrasts)'] || pc['Conserved Up'] || 0) + (pc['Consistent Down (all contrasts)'] || pc['Conserved Down'] || 0);
            const specific = Object.entries(pc).filter(([k]) => k.startsWith('Specific')).reduce((a, [, v]) => a + v, 0);
            $('int-tiles').innerHTML = tile(s.mode === 'contrasts' ? 'Analyses' : 'Datasets', n, esc(s.mode_text)) +
                tile('Genes aligned', int(s.n_genes), `${int(T.n_sig.filter((x) => x > 0).length)} significant in at least one`) +
                tile(s.mode === 'contrasts' ? 'Significant in every analysis' : 'Significant in every dataset', int(all), `${int(same)} with the same direction`) +
                tile('Specific', int(specific), 'significant in one only') +
                tile('Discordant', int(pc.Discordant || 0), 'opposite directions');
            $('pair-x').innerHTML = s.profiles.map((p, k) => `<option value="${esc(p.id)}"${k === 0 ? ' selected' : ''}>${esc(p.label)}</option>`).join('');
            $('pair-y').innerHTML = s.profiles.map((p, k) => `<option value="${esc(p.id)}"${k === 1 ? ' selected' : ''}>${esc(p.label)}</option>`).join('');
        }
        $('int-notes').innerHTML = (s.notes || []).map((x) => `<li>${esc(x)}</li>`).join('');
        const patterns = meta ? ['Up-regulated', 'Down-regulated', 'Not significant'] : Object.keys(pc);
        $('int-pattern').innerHTML = '<option value="">All genes</option><option value="__sig">Significant somewhere</option>' +
            patterns.map((p) => `<option value="${esc(p)}">${esc(p)}</option>`).join('');
        $('int-pattern').value = '__sig';
        $('int-files').innerHTML = `<a class="link-btn" href="/api/jobs/${id}/files/integration.csv?download=1">Complete table (.csv)</a>` +
            `<a class="link-btn" href="/api/jobs/${id}/files/manifest.json?download=1">Manifest (.json)</a>`;
        intGenes.sortKey = null;
        intGenes.setRows(T.symbol.map((_, k) => k));
        renderIntegrationCharts();
        renderFiles();
        if (fresh) $('int-results').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    function renderIntegrationCharts() {
        const s = state.int;
        if (!s) return;
        if (s.mode === 'meta') {
            show('chart-meta-volcano', figMetaVolcano());
            renderForest();
            show('chart-i2', figI2());
            show('chart-meta-heat', figMetaHeat());
            const el = $('chart-meta-volcano');
            if (!el._lumosClick) {
                el.on('plotly_click', (ev) => {
                    const pt = ev.points && ev.points[0];
                    if (pt && pt.customdata !== undefined) {
                        const k = pt.customdata;
                        if (![...$('forest-gene').options].some((o) => Number(o.value) === k)) {
                            $('forest-gene').insertAdjacentHTML('afterbegin', `<option value="${k}">${esc(state.int.table.symbol[k])}</option>`);
                        }
                        $('forest-gene').value = k;
                        renderForest();
                    }
                });
                el._lumosClick = true;
            }
        } else {
            show('chart-patterns', figPatterns());
            show('chart-pair', figPair());
            show('chart-effects', figEffects());
        }
    }
    $('pair-x').addEventListener('change', () => show('chart-pair', figPair()));
    $('pair-y').addEventListener('change', () => show('chart-pair', figPair()));
    $('forest-gene').addEventListener('change', renderForest);

    function figPatterns() {
        const pc = Object.entries(state.int.pattern_counts).filter(([k]) => k !== 'Not significant').sort((a, b) => b[1] - a[1]);
        if (!pc.length) return emptyFig('No gene passed the thresholds in any analysis.');
        const color = (k) => css(/Up/.test(k) && /Consistent|Conserved|Shared/.test(k) ? '--chart-2' : /Down/.test(k) && /Consistent|Conserved|Shared/.test(k) ? '--chart-1'
            : k === 'Discordant' ? '--chart-5' : '--chart-4');
        return { height: Math.max(220, 26 * pc.length + 70), data: [{
            type: 'bar', orientation: 'h', x: pc.map(([, v]) => v), y: pc.map(([k]) => k), text: pc.map(([, v]) => String(v)), textposition: 'auto',
            marker: { color: pc.map(([k]) => color(k)) }, hovertemplate: '%{y}: %{x} genes<extra></extra>',
        }], layout: layout('genes', '', { yaxis: { autorange: 'reversed', ticks: '', type: 'category' }, margin: { l: 10, r: 14, t: 8, b: 44 } }) };
    }

    function figPair() {
        const s = state.int, T = s.table;
        const x = $('pair-x').value, y = $('pair-y').value;
        if (!T.profiles[x] || !T.profiles[y] || x === y) return emptyFig('Choose two different analyses.');
        const px = T.profiles[x], py = T.profiles[y];
        const lab = (id) => s.profiles.find((p) => p.id === id).label;
        const groups = new Map();
        const add = (name, k) => { if (!groups.has(name)) groups.set(name, []); groups.get(name).push(k); };
        T.symbol.forEach((_, k) => {
            if (!Number.isFinite(num(px.lfc[k])) || !Number.isFinite(num(py.lfc[k]))) return;
            const a = px.status[k], b = py.status[k];
            if (a && b) add(a === b ? (a > 0 ? 'Conserved up' : 'Conserved down') : 'Discordant', k);
            else if (a) add(`Specific ${lab(x)}`, k);
            else if (b) add(`Specific ${lab(y)}`, k);
            else add('Not significant', k);
        });
        const colors = { 'Conserved up': '--chart-2', 'Conserved down': '--chart-1', Discordant: '--chart-5',
            [`Specific ${lab(x)}`]: '--chart-4', [`Specific ${lab(y)}`]: '--chart-3', 'Not significant': '--chart-muted' };
        const order = ['Not significant', `Specific ${lab(x)}`, `Specific ${lab(y)}`, 'Discordant', 'Conserved down', 'Conserved up'];
        const data = order.filter((g) => groups.has(g)).map((g) => {
            const ks = groups.get(g);
            return { type: 'scattergl', mode: 'markers', name: `${g} (${ks.length})`, x: ks.map((k) => px.lfc[k]), y: ks.map((k) => py.lfc[k]),
                text: ks.map((k) => T.symbol[k]), hovertemplate: '<b>%{text}</b><br>%{x:.2f} / %{y:.2f}<extra></extra>',
                marker: { color: css(colors[g]), size: g === 'Not significant' ? 3 : 6 } };
        });
        const ink = css('--chart-axis');
        return { data, layout: layout(`log2FC · ${lab(x)}`, `log2FC · ${lab(y)}`, Object.assign(legendTop(), {
            hovermode: 'closest',
            shapes: [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: ink, width: 1, dash: 'dot' } },
                { type: 'line', yref: 'paper', y0: 0, y1: 1, x0: 0, x1: 0, line: { color: ink, width: 1, dash: 'dot' } }],
        })) };
    }

    function figEffects() {
        const s = state.int, T = s.table;
        const ids = s.profiles.map((p) => p.id);
        const cand = T.symbol.map((_, k) => k).filter((k) => T.n_sig[k] > 0)
            .map((k) => [k, Math.max(...ids.map((i) => Math.abs(num(T.profiles[i].lfc[k])) || 0))])
            .sort((a, b) => b[1] - a[1]).slice(0, 50).map(([k]) => k);
        if (!cand.length) return emptyFig('No gene passed the thresholds in any analysis.');
        const light = document.body.classList.contains('light-mode');
        const maxAbs = Math.max(1, ...cand.flatMap((k) => ids.map((i) => Math.abs(num(T.profiles[i].lfc[k])) || 0)));
        return { height: Math.max(240, 46 * ids.length + 160), data: [{
            type: 'heatmap', x: cand.map((k) => T.symbol[k]), y: s.profiles.map((p) => p.label),
            z: ids.map((i) => cand.map((k) => T.profiles[i].lfc[k])), zmid: 0, zmin: -maxAbs, zmax: maxAbs,
            colorscale: [[0, css('--chart-1')], [0.5, light ? '#f4f4f2' : '#22252c'], [1, css('--chart-2')]],
            hovertemplate: '<b>%{x}</b><br>%{y}<br>log2FC %{z:.2f}<extra></extra>', colorbar: { title: { text: 'log2FC', side: 'right' }, thickness: 12, outlinewidth: 0 },
        }], layout: layout('', '', { xaxis: { tickangle: -60, type: 'category', showgrid: false, ticks: '' }, yaxis: { type: 'category', showgrid: false, ticks: '' }, margin: { l: 140, r: 14, t: 10, b: 90 } }) };
    }

    function figMetaVolcano() {
        const T = state.int.table, th = state.int.thresholds;
        const lists = { 0: [], 1: [], '-1': [] };
        T.symbol.forEach((_, k) => { if (Number.isFinite(num(T.lfc[k])) && Number.isFinite(num(T.fdr[k]))) lists[metaStatus(k)].push(k); });
        const y = (k) => -Math.log10(Math.max(T.fdr[k], 1e-300));
        const tr = (ks, name, color, size) => ({ type: 'scattergl', mode: 'markers', name, x: ks.map((k) => T.lfc[k]), y: ks.map(y), customdata: ks,
            text: ks.map((k) => `${T.symbol[k]} · I² ${fmt(T.i2[k], 0)}% · ${T.n[k]} studies`), hovertemplate: '<b>%{text}</b><br>log2FC %{x:.2f}<br>−log10 FDR %{y:.2f}<extra></extra>',
            marker: { color: css(color), size, line: { width: 0 } } });
        const ink = css('--chart-axis');
        const shapes = [{ type: 'line', xref: 'paper', x0: 0, x1: 1, y0: -Math.log10(th.fdr), y1: -Math.log10(th.fdr), line: { color: ink, width: 1, dash: 'dash' } }];
        if (th.lfc > 0) for (const v of [-th.lfc, th.lfc]) shapes.push({ type: 'line', yref: 'paper', y0: 0, y1: 1, x0: v, x1: v, line: { color: ink, width: 1, dash: 'dash' } });
        return { data: [tr(lists[0], 'Not significant', '--chart-muted', 4), tr(lists[-1], `Down (${lists[-1].length})`, '--chart-1', 6), tr(lists[1], `Up (${lists[1].length})`, '--chart-2', 6)],
            layout: layout('meta-analysis log2 fold change', '−log10 adjusted p-value', Object.assign(legendTop(), { shapes, hovermode: 'closest' })) };
    }

    function figForest(k) {
        const s = state.int, T = s.table;
        if (k === undefined || k === null || Number.isNaN(k) || T.symbol[k] === undefined) return emptyFig('Choose a gene.');
        const rows = s.profiles.map((p) => ({ label: p.label, eff: num(T.profiles[p.id].lfc[k]), se: num(T.profiles[p.id].se[k]) }))
            .filter((r) => Number.isFinite(r.eff) && Number.isFinite(r.se));
        rows.forEach((r) => { r.lo = r.eff - 1.96 * r.se; r.hi = r.eff + 1.96 * r.se; });
        rows.push({ label: 'Random-effects estimate', eff: num(T.lfc[k]), lo: num(T.ci_low[k]), hi: num(T.ci_high[k]), meta: true });
        $('forest-sub').textContent = `${T.symbol[k]}: I² ${fmt(T.i2[k], 1)}%, τ² ${fmtP(T.tau2[k])}, FDR ${fmtP(T.fdr[k])}; 95% confidence intervals.`;
        return { height: Math.max(240, 50 * rows.length + 90), data: [{
            type: 'scatter', mode: 'markers', x: rows.map((r) => r.eff), y: rows.map((r) => r.label),
            error_x: { type: 'data', symmetric: false, array: rows.map((r) => r.hi - r.eff), arrayminus: rows.map((r) => r.eff - r.lo), color: css('--chart-text'), thickness: 1.5, width: 4 },
            marker: { size: rows.map((r) => (r.meta ? 14 : 9)), symbol: rows.map((r) => (r.meta ? 'diamond' : 'square')), color: rows.map((r) => css(r.meta ? '--accent-color' : '--chart-1')) },
            hovertemplate: '%{y}<br>log2FC %{x:.2f}<extra></extra>',
        }], layout: layout('log2 fold change', '', { yaxis: { autorange: 'reversed', type: 'category', ticks: '' }, margin: { l: 170, r: 14, t: 10, b: 44 },
            shapes: [{ type: 'line', yref: 'paper', y0: 0, y1: 1, x0: 0, x1: 0, line: { color: css('--chart-axis'), width: 1, dash: 'dash' } }] }) };
    }

    function renderForest() { show('chart-forest', figForest(Number($('forest-gene').value))); }

    function figI2() {
        const v = state.int.table.i2.filter((x) => x !== null);
        return { data: [{ type: 'histogram', x: v, xbins: { start: 0, end: 100, size: 5 }, marker: { color: css('--chart-3') }, hovertemplate: 'I² %{x}%: %{y} genes<extra></extra>' }],
            layout: layout('I² (%)', 'genes', { bargap: 0.05 }) };
    }

    function figMetaHeat() {
        const h = state.int.heatmap;
        if (!h) return Object.assign(emptyFig('No meta-analysis gene passed the thresholds.'), { height: 200 });
        const light = document.body.classList.contains('light-mode');
        const x = h.samples.map((smp, k) => `${h.studies[k]} · ${smp}`);
        return { height: Math.max(260, 16 * h.genes.length + 120), data: [{
            type: 'heatmap', z: h.z, x, y: h.genes, zmid: 0, zmin: -2.5, zmax: 2.5,
            customdata: h.genes.map(() => h.groups),
            colorscale: [[0, css('--chart-1')], [0.5, light ? '#f4f4f2' : '#22252c'], [1, css('--chart-2')]],
            hovertemplate: '<b>%{y}</b><br>%{x} (%{customdata})<br>z = %{z:.2f}<extra></extra>', colorbar: { title: { text: 'z-score', side: 'right' }, thickness: 12, outlinewidth: 0 },
        }], layout: layout('', '', { xaxis: { showticklabels: false, type: 'category', showgrid: false, ticks: '' }, yaxis: { autorange: 'reversed', type: 'category', showgrid: false, ticks: '' }, margin: { l: 90, r: 14, t: 10, b: 20 } }) };
    }

    // gene table of the integration
    const intGenes = new DataTable($('int-genes'), {
        columns: () => {
            const s = state.int;
            if (!s) return [];
            const T = s.table;
            if (s.mode === 'meta') {
                return [
                    { key: 'g', label: 'Gene', value: (k) => T.symbol[k] },
                    { key: 'st', label: 'Status', value: (k) => metaStatus(k), export: (k) => ['Down-regulated', 'Not significant', 'Up-regulated'][metaStatus(k) + 1],
                        render: (k) => (metaStatus(k) > 0 ? '<span class="status-up">Up</span>' : metaStatus(k) < 0 ? '<span class="status-down">Down</span>' : '<span class="hint">–</span>') },
                    { key: 'lfc', label: 'Meta log2FC', num: true, value: (k) => num(T.lfc[k]), render: (k) => fmt(T.lfc[k]) },
                    { key: 'ci', label: '95% CI', value: (k) => `${fmt(T.ci_low[k])} to ${fmt(T.ci_high[k])}` },
                    { key: 'fdr', label: 'FDR', num: true, value: (k) => num(T.fdr[k]), render: (k) => fmtP(T.fdr[k]) },
                    { key: 'i2', label: 'I² (%)', num: true, value: (k) => num(T.i2[k]), render: (k) => fmt(T.i2[k], 1) },
                    { key: 'tau', label: 'τ²', num: true, value: (k) => num(T.tau2[k]), render: (k) => fmtP(T.tau2[k]) },
                    { key: 'n', label: 'Studies', num: true, value: (k) => T.n[k] },
                ];
            }
            const cols = [
                { key: 'g', label: 'Gene', value: (k) => T.symbol[k] },
                { key: 'pat', label: 'Pattern', value: (k) => T.pattern[k] },
                { key: 'n', label: 'Significant in', num: true, value: (k) => T.n_sig[k] },
            ];
            s.profiles.forEach((p) => {
                cols.push({ key: `l_${p.id}`, label: `log2FC · ${p.label}`, num: true, value: (k) => num(T.profiles[p.id].lfc[k]), render: (k) => fmt(T.profiles[p.id].lfc[k]) });
                cols.push({ key: `f_${p.id}`, label: `FDR · ${p.label}`, num: true, value: (k) => num(T.profiles[p.id].fdr[k]), render: (k) => fmtP(T.profiles[p.id].fdr[k]) });
            });
            return cols;
        },
        rowKey: (k) => k,
        text: (k) => String(state.int.table.symbol[k]).toLowerCase(),
        extraFilter: (k) => {
            const v = $('int-pattern').value;
            if (!v) return true;
            const s = state.int;
            if (s.mode === 'meta') {
                if (v === '__sig') return metaStatus(k) !== 0;
                return ['Down-regulated', 'Not significant', 'Up-regulated'][metaStatus(k) + 1] === v;
            }
            return v === '__sig' ? s.table.n_sig[k] > 0 : s.table.pattern[k] === v;
        },
        empty: 'No gene in this selection.',
    });
    $('int-filter').addEventListener('input', debounce(() => intGenes.setFilter($('int-filter').value), 200));
    $('int-pattern').addEventListener('change', () => { intGenes.limit = 300; intGenes.render(); });

    // ------------------------------------------------------------------ APP: gene source
    function renderAppSources() {
        const sel = $('app-source');
        const prev = sel.value;
        let html = '';
        if (state.analyses.size) {
            html += '<optgroup label="Differential-expression analyses">' + [...state.analyses.entries()].map(([id, s]) =>
                `<option value="de|${id}|analysis">${esc(s.label)} · ${esc(s.name)}</option>`).join('') + '</optgroup>';
        }
        state.integrations.forEach((s, id) => {
            const names = s.profiles.map((p) => p.label).join(' · ');
            html += `<optgroup label="${esc(MODE_NAMES[s.mode])}: ${esc(names)}">`;
            if (s.mode === 'meta') html += `<option value="int|${id}|meta">Meta-analysis effect</option>`;
            else {
                html += `<option value="int|${id}|consistent">Consistent signature (same direction in every ${s.mode === 'tissues' ? 'dataset' : 'analysis'})</option>`;
                html += s.profiles.map((p) => `<option value="int|${id}|profile:${esc(p.id)}">Profile: ${esc(p.label)}</option>`).join('');
            }
            html += '</optgroup>';
        });
        sel.innerHTML = html;
        if ([...sel.options].some((o) => o.value === prev)) sel.value = prev;
        else if (state.deJob) sel.value = `de|${state.deJob}|analysis`;
        if (sel.value !== (state.appSource && state.appSource.value)) buildAppSource(true);
    }
    $('app-source').addEventListener('change', () => buildAppSource(true));
    $('app-fdr').addEventListener('input', debounce(() => buildAppSource(false), 300));
    $('app-lfc').addEventListener('input', debounce(() => buildAppSource(false), 300));

    // rows of the gene source: {key, label, id, sym, lfc, fdr}
    function buildAppSource(newSource) {
        const value = $('app-source').value;
        if (!value) { state.appSource = null; refreshButtons(); return; }
        const [kind, job, focus] = value.split('|');
        if (newSource) {
            const th = kind === 'de' ? (job === state.deJob ? { fdr: state.fdr, lfc: state.lfc } : state.analyses.get(job).parameters)
                : state.integrations.get(job).thresholds;
            $('app-fdr').value = th.fdr;
            $('app-lfc').value = th.lfc;
        }
        const fdr = Number($('app-fdr').value) || 0.05, lfc = Number($('app-lfc').value) || 0;
        const rows = [];
        let universe = 0;
        if (kind === 'de') {
            const s = state.analyses.get(job), G = genesOf(job);
            const keys = new Set();
            for (let i = 0; i < G.n; i++) {
                keys.add(keyIn(G, i));
                if (G.fdr[i] < fdr && Math.abs(G.lfc[i]) >= lfc) {
                    rows.push({ key: keyIn(G, i), label: labelIn(G, s, i), id: G.id[i], sym: G.sym[i], lfc: G.lfc[i], fdr: G.fdr[i] });
                }
            }
            universe = keys.size;
        } else {
            const s = state.integrations.get(job), T = s.table;
            if (focus === 'meta') {
                T.symbol.forEach((g, k) => {
                    if (Number.isFinite(num(T.fdr[k]))) universe++;
                    if (num(T.fdr[k]) < fdr && Math.abs(num(T.lfc[k])) >= lfc) rows.push({ key: g, label: g, lfc: T.lfc[k], fdr: T.fdr[k] });
                });
            } else if (focus.startsWith('profile:')) {
                const P = T.profiles[focus.slice(8)];
                T.symbol.forEach((g, k) => {
                    if (Number.isFinite(num(P.fdr[k]))) universe++;
                    if (num(P.fdr[k]) < fdr && Math.abs(num(P.lfc[k])) >= lfc) rows.push({ key: g, label: g, lfc: P.lfc[k], fdr: P.fdr[k] });
                });
            } else {
                const ids = s.profiles.map((p) => p.id);
                T.symbol.forEach((g, k) => {
                    const l = ids.map((i) => num(T.profiles[i].lfc[k])), f = ids.map((i) => num(T.profiles[i].fdr[k]));
                    if (l.concat(f).every(Number.isFinite)) universe++;
                    if (f.every((x) => x < fdr) && l.every((x) => Math.abs(x) >= lfc) && (l.every((x) => x > 0) || l.every((x) => x < 0))) {
                        rows.push({ key: g, label: g, lfc: l.reduce((a, b) => a + b, 0) / l.length, fdr: Math.max(...f) });
                    }
                });
            }
        }
        rows.sort((a, b) => a.fdr - b.fdr);
        state.appSource = { value, kind, job, focus: kind === 'de' ? 'analysis' : focus, rows, universe };
        $('bg-source-label').textContent = `Genes tested in the source (${int(universe)})`;
        renderAppGenes(true);
        refreshButtons();
    }

    const appTable = new DataTable($('app-table'), {
        selectable: true,
        columns: () => {
            const src = state.appSource;
            const cols = [{ key: 'label', label: 'Gene', value: (r) => r.label }];
            if (src && src.kind === 'de' && symbolsAvailable(state.analyses.get(src.job))) {
                cols.push(state.labels === 'symbol' ? { key: 'id', label: 'ID', value: (r) => r.id } : { key: 'sym', label: 'Gene symbol', value: (r) => r.sym || '' });
            }
            cols.push({ key: 'dir', label: 'Direction', value: (r) => Math.sign(r.lfc), render: (r) => (r.lfc > 0 ? '<span class="status-up">Up</span>' : '<span class="status-down">Down</span>') });
            cols.push({ key: 'lfc', label: src && src.focus === 'consistent' ? 'Mean log2FC' : 'log2FC', num: true, value: (r) => r.lfc, render: (r) => fmt(r.lfc) });
            cols.push({ key: 'fdr', label: src && src.focus === 'consistent' ? 'Largest FDR' : 'FDR', num: true, value: (r) => r.fdr, render: (r) => fmtP(r.fdr) });
            return cols;
        },
        rowKey: (r) => r.key,
        text: (r) => `${r.label} ${r.key} ${r.id || ''}`.toLowerCase(),
        onSelect: updateAppCount,
        pageSize: 500,
        empty: 'No gene passes these thresholds in this source.',
    });

    function renderAppGenes(resetSelection) {
        const src = state.appSource;
        if (!src) return;
        if (src.kind === 'de') {
            const s = state.analyses.get(src.job), G = genesOf(src.job);
            src.rows.forEach((r) => { const i = G.byId.get(r.id); if (i !== undefined) r.label = labelIn(G, s, i); });
        }
        if (resetSelection) appTable.selected = new Set(src.rows.map((r) => r.key));
        appTable.setRows(src.rows);
        updateAppCount();
    }

    function selectedGenes() {
        const keys = appTable.selected;
        const rows = state.appSource ? state.appSource.rows.filter((r) => keys.has(r.key)) : [];
        return { up: [...new Set(rows.filter((r) => r.lfc > 0).map((r) => r.key))], down: [...new Set(rows.filter((r) => r.lfc < 0).map((r) => r.key))] };
    }

    function updateAppCount() {
        const { up, down } = selectedGenes();
        $('app-sel-count').textContent = `${up.length + down.length} of ${appTable.rows.length} genes ticked (${up.length} up, ${down.length} down)`;
    }

    document.querySelectorAll('[data-pick]').forEach((b) => b.addEventListener('click', () => {
        const pick = b.dataset.pick;
        const rows = appTable.view();
        if (pick === 'none') appTable.selected.clear();
        else appTable.selected = new Set(rows.filter((r) => pick === 'all' || (pick === 'up' ? r.lfc > 0 : r.lfc < 0)).map((r) => r.key));
        appTable.render();
        updateAppCount();
    }));
    $('app-filter').addEventListener('input', debounce(() => appTable.setFilter($('app-filter').value), 200));
    document.querySelectorAll('input[name="app_bg"]').forEach((el) => el.addEventListener('change', () => {
        $('app-background').classList.toggle('hidden', document.querySelector('input[name="app_bg"]:checked').value !== 'custom');
    }));

    function renderLibraries() {
        $('lib-list').innerHTML = state.libraries.map((l, k) =>
            `<span class="pill">${esc(l)}<button type="button" data-k="${k}" aria-label="Remove ${esc(l)}">×</button></span>`).join('') ||
            '<span class="hint">No library: add one below.</span>';
    }
    $('lib-list').addEventListener('click', (e) => {
        const b = e.target.closest('button[data-k]');
        if (!b) return;
        state.libraries.splice(Number(b.dataset.k), 1);
        renderLibraries();
    });

    function renderTfLibraries() {
        const all = state.caps ? state.caps.tf_libraries : [];
        $('tf-libs').innerHTML = all.map((l) => `<label class="pill"><input type="checkbox" value="${esc(l)}"${state.tfLibraries.has(l) ? ' checked' : ''}>${esc(l)}</label>`).join('');
    }
    $('tf-libs').addEventListener('change', (e) => {
        if (e.target.checked) state.tfLibraries.add(e.target.value); else state.tfLibraries.delete(e.target.value);
    });

    let librariesLoaded = false;
    function loadLibraryList() {
        if (librariesLoaded) return;
        librariesLoaded = true;
        fetch('/api/enrichr/libraries').then((r) => r.json()).then((d) => {
            const sel = $('lib-add');
            if (!d.libraries.length) {
                sel.innerHTML = `<option value="">${esc(d.error || 'Library list unavailable')}</option>`;
                librariesLoaded = false;
                return;
            }
            const pathway = /KEGG|GO_|Reactome|WikiPathway|MSigDB|BioCarta|Panther|HumanCyc|BioPlanet|Elsevier|Disease|DisGeNET|Cell_Type|Tissue|HPO|OMIM|Human_Phenotype/i;
            const first = d.libraries.filter((l) => pathway.test(l));
            const rest = d.libraries.filter((l) => !pathway.test(l));
            sel.innerHTML = '<option value="">Add an Enrichr library…</option>' +
                `<optgroup label="Pathways, ontologies, diseases">${first.map((l) => `<option>${esc(l)}</option>`).join('')}</optgroup>` +
                `<optgroup label="Other libraries">${rest.map((l) => `<option>${esc(l)}</option>`).join('')}</optgroup>`;
        }).catch(() => { librariesLoaded = false; });
    }
    $('lib-add').addEventListener('change', () => {
        const v = $('lib-add').value;
        if (v && !state.libraries.includes(v) && state.libraries.length < 8) state.libraries.push(v);
        $('lib-add').value = '';
        renderLibraries();
    });

    $('app-btn').addEventListener('click', async () => {
        const src = state.appSource;
        const { up, down } = selectedGenes();
        if (up.length + down.length < 3) return showError('Tick at least three genes in the table.', $('app-error'));
        if (!state.libraries.length) return showError('Add at least one pathway library.', $('app-error'));
        const custom = document.querySelector('input[name="app_bg"]:checked').value === 'custom';
        if (custom && $('app-background').value.trim().split(/\s+/).length < 50) return showError('The custom background needs the measured genes (at least 50).', $('app-error'));
        const fd = new FormData();
        fd.append('action', 'enrich');
        fd.append('source', src.job);
        fd.append('focus', src.focus);
        fd.append('genes_up', up.join('\n'));
        fd.append('genes_down', down.join('\n'));
        fd.append('strategy', document.querySelector('input[name="app_strategy"]:checked').value);
        if (custom) fd.append('background', $('app-background').value);
        fd.append('libraries', state.libraries.join(','));
        fd.append('tf_libraries', [...state.tfLibraries].join(','));
        fd.append('fdr', $('app-fdr').value);
        fd.append('lfc', $('app-lfc').value);
        fd.append('keep_nonhuman', $('keep-nonhuman').checked);
        try {
            await appJob.start(fd, `Enrichment: ${up.length} up and ${down.length} down genes…`);
        } catch (e) {
            showError(e.message, $('app-error'));
        }
        refreshButtons();
    });

    // ------------------------------------------------------------------ APP results
    function showApp(id, s, fresh) {
        state.appJob = id;
        state.app = s;
        state.stringData = null;
        if (state.stringNet) { state.stringNet.destroy(); state.stringNet = null; }
        $('string-msg').style.display = 'flex';
        $('string-msg').textContent = 'Build the network to see the associations among the ticked genes.';
        resetPG();
        $('app-results').classList.remove('hidden');
        $('tab-pg-btn').disabled = !(s.regulators && s.regulators.length);
        const q = s.query_sizes || {};
        const sizes = Object.entries(q).map(([d, v]) => `${d.replace('-regulated', '')}: ${v.query}`).join(' · ');
        const tiles = [tile('Gene source', esc(s.source.label), s.source.dataset ? esc(s.source.dataset) : '', true),
            tile('Genes queried', int(s.genes.length), esc(sizes) + (s.strategy === 'separate' ? ' (separately)' : ' (combined)')),
            tile('Background', int(s.background.size), esc(s.background.source))];
        s.libraries.forEach((l) => tiles.push(tile(l.name, int(l.n_significant), `of ${int(l.n_terms)} terms with FDR < 0.05 · ${esc(l.engine.replace(/ \(.*\)/, ''))}`)));
        $('app-tiles').innerHTML = tiles.join('');
        $('app-notes').innerHTML = (s.notes || []).map((n) => `<li>${esc(n)}</li>`).join('');
        $('tf-csv').href = `/api/jobs/${id}/files/regulators.csv?download=1`;
        $('tf-csv').classList.toggle('hidden', !s.files.regulators);
        renderAppCharts();
        renderFiles();
        updateUrl();
        if (fresh) $('app-results').scrollIntoView({ behavior: 'smooth', block: 'start' });
    }

    const DIR = { 'Up-regulated': ['↑', '--chart-2'], 'Down-regulated': ['↓', '--chart-1'], Combined: ['•', '--chart-3'] };
    const geneLabelList = (genes) => genes.map(labelOfSymbol).join(', ');

    function figLibrary(l) {
        const terms = l.terms.slice(0, 15);
        if (!terms.length) return emptyFig('No terms.');
        return { height: 460, data: [{
            type: 'bar', orientation: 'h',
            x: terms.map((t) => -Math.log10(Math.max(t.fdr ?? 1, 1e-300))),
            y: terms.map((t, k) => `${(DIR[t.direction] || ['•'])[0]} ${t.term.length > 52 ? t.term.slice(0, 50) + '…' : t.term}${'​'.repeat(k)}`),
            text: terms.map((t) => String(t.genes.length)), textposition: 'auto',
            customdata: terms.map((t) => [t.term, fmtP(t.fdr), t.direction, geneLabelList(t.genes.slice(0, 12)) + (t.genes.length > 12 ? '…' : '')]),
            hovertemplate: '<b>%{customdata[0]}</b><br>%{customdata[2]} · FDR %{customdata[1]}<br>%{customdata[3]}<extra></extra>',
            marker: { color: terms.map((t) => css((DIR[t.direction] || ['', '--chart-3'])[1])) },
        }], layout: layout('−log10 adjusted p-value', '', {
            yaxis: { autorange: 'reversed', ticks: '', type: 'category' }, margin: { l: 10, r: 14, t: 8, b: 44 },
            shapes: [{ type: 'line', yref: 'paper', y0: 0, y1: 1, x0: -Math.log10(0.05), x1: -Math.log10(0.05), line: { color: css('--chart-axis'), width: 1, dash: 'dash' } }],
        }) };
    }

    function bestRegulators(n = 20) {
        const best = new Map();
        (state.app.regulators || []).forEach((r) => { if (!best.has(r.tf) || (r.fdr ?? 1) < (best.get(r.tf).fdr ?? 1)) best.set(r.tf, r); });
        return [...best.values()].sort((a, b) => (a.fdr ?? 1) - (b.fdr ?? 1)).slice(0, n);
    }

    function figRegulators() {
        const regs = bestRegulators(20);
        if (!regs.length) return emptyFig('No regulators.');
        const libs = [...new Set(regs.map((r) => r.library))];
        return { height: 460, data: libs.map((lib, k) => {
            const rs = regs.filter((r) => r.library === lib);
            return { type: 'bar', orientation: 'h', name: lib, x: rs.map((r) => -Math.log10(Math.max(r.fdr ?? 1, 1e-300))), y: rs.map((r) => labelOfSymbol(r.tf)),
                text: rs.map((r) => `${(DIR[r.direction] || ['•'])[0]} ${r.targets.length}`), textposition: 'auto',
                customdata: rs.map((r) => [r.term, fmtP(r.fdr), r.direction, geneLabelList(r.targets.slice(0, 15)) + (r.targets.length > 15 ? '…' : '')]),
                hovertemplate: '<b>%{y}</b> (%{customdata[0]})<br>%{customdata[2]} · FDR %{customdata[1]}<br>targets: %{customdata[3]}<extra></extra>',
                marker: { color: css(['--chart-6', '--chart-5', '--chart-4', '--chart-3'][k % 4]) } };
        }), layout: layout('−log10 adjusted p-value', '', Object.assign(legendTop(), {
            yaxis: { autorange: 'reversed', ticks: '', type: 'category', categoryorder: 'array', categoryarray: regs.map((r) => labelOfSymbol(r.tf)) },
            margin: { l: 90, r: 14, t: 34, b: 44 }, barmode: 'overlay',
            shapes: [{ type: 'line', yref: 'paper', y0: 0, y1: 1, x0: -Math.log10(0.05), x1: -Math.log10(0.05), line: { color: css('--chart-axis'), width: 1, dash: 'dash' } }],
        })) };
    }

    function renderAppCharts() {
        const s = state.app;
        const box = $('lib-charts');
        const libs = s.libraries.filter((l) => l.role === 'pathway');
        if (box.dataset.job !== state.appJob) {
            box.innerHTML = libs.map((l, k) => `
                <div class="chart-card">
                    <div class="card-head"><div><h3>${esc(l.name)}</h3>
                    <p class="chart-sub">Top 15 terms by adjusted p-value (↑ up-regulated list, ↓ down-regulated list); bar labels give the number of genes in the term.</p></div>
                    ${l.file ? `<a class="link-btn" href="/api/jobs/${state.appJob}/files/${encodeURIComponent(l.file)}?download=1">.csv</a>` : ''}</div>
                    <div class="plot" id="chart-lib-${k}" style="height: 460px;"></div>
                    <details class="explain" style="margin: 8px 0 10px;"><summary>Table (top 100 terms)</summary>
                        <div class="table-wrap scroll" style="margin-top: 8px; max-height: 360px;"><table class="data-table wrap" id="table-lib-${k}"></table></div>
                    </details>
                </div>`).join('');
            box.dataset.job = state.appJob;
        }
        libs.forEach((l, k) => {
            show(`chart-lib-${k}`, figLibrary(l));
            $(`table-lib-${k}`).innerHTML = '<thead><tr><th>Term</th><th>List</th><th class="num">Genes</th><th class="num">p-value</th><th class="num">FDR</th><th>Genes in the term</th></tr></thead><tbody>' +
                l.terms.slice(0, 100).map((t) => `<tr><td>${esc(t.term)}</td><td>${esc(t.direction)}</td><td class="num">${t.genes.length}</td><td class="num">${fmtP(t.p)}</td>` +
                    `<td class="num">${fmtP(t.fdr)}</td><td>${esc(geneLabelList(t.genes))}</td></tr>`).join('') + '</tbody>';
        });
        const hasTf = s.regulators && s.regulators.length;
        $('tf-card').classList.toggle('hidden', !hasTf);
        if (hasTf) {
            $('tf-sub').textContent = `${s.tf_libraries.join(' and ')}: candidate regulators whose annotated targets are over-represented among the selected genes ` +
                '(best term per factor, top 20; ↑/↓ list and number of targets). Target-set enrichment, not proof of binding or causality.';
            show('chart-tf', figRegulators());
        }
    }

    // ------------------------------------------------------------------ STRING network
    function topPathwayTerms(n) {
        const terms = [];
        state.app.libraries.filter((l) => l.role === 'pathway').forEach((l) => l.terms.slice(0, n).forEach((t) => terms.push(t)));
        return terms.sort((a, b) => (a.fdr ?? 1) - (b.fdr ?? 1)).slice(0, n);
    }

    $('string-btn').addEventListener('click', async () => {
        const max = Math.max(10, Math.min(500, Number($('string-max').value) || 150));
        const genes = [...state.app.genes].sort((a, b) => Math.abs(b.effect ?? 0) - Math.abs(a.effect ?? 0)).slice(0, max).map((g) => g.symbol);
        if (genes.length < 2) return showError('At least two genes are needed.', $('app-error'));
        $('string-btn').disabled = true;
        $('string-msg').style.display = 'flex';
        $('string-msg').textContent = `Asking STRING about ${genes.length} genes…`;
        try {
            const r = await fetch(`/api/jobs/${state.appJob}/string`, {
                method: 'POST', headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ genes, score: Number($('string-score').value), network_type: $('string-type').value }),
            });
            const d = await r.json();
            if (!r.ok) throw new Error(d.error || 'STRING failed.');
            state.stringData = d;
            drawString();
        } catch (e) {
            $('string-msg').textContent = e.message;
        }
        $('string-btn').disabled = false;
    });

    function geneDirection(sym) {
        const g = state.app && state.app.genes.find((x) => x.symbol === String(sym).toUpperCase());
        return g ? (g.direction === 'up' ? 1 : -1) : 0;
    }

    function netColors() {
        return {
            text: css('--text-main'), muted: css('--chart-muted'), up: css('--chart-2'), down: css('--chart-1'),
            tf: css('--chart-6'), pathway: css('--accent-color'), hit: css('--chart-3'),
            neutral: document.body.classList.contains('light-mode') ? '#9aa0aa' : '#6b7280',
        };
    }

    function drawString() {
        const d = state.stringData;
        if (state.stringNet) { state.stringNet.destroy(); state.stringNet = null; }
        if (!d.edges.length) {
            $('string-msg').style.display = 'flex';
            $('string-msg').textContent = 'STRING found no associations among these genes at this score.';
            $('string-map').innerHTML = '';
            return;
        }
        $('string-msg').style.display = 'none';
        const c = netColors();
        const terms = topPathwayTerms(25);
        const degree = new Map();
        d.edges.forEach((e) => { degree.set(e.a, (degree.get(e.a) || 0) + 1); degree.set(e.b, (degree.get(e.b) || 0) + 1); });
        const nodes = d.nodes.map((n) => {
            const st = geneDirection(n);
            const col = st > 0 ? c.up : st < 0 ? c.down : c.neutral;
            return { id: n, label: labelOfSymbol(n), shape: 'dot', size: 8 + Math.min(16, 2 * (degree.get(n) || 0)),
                color: { background: col, border: col, highlight: { background: col, border: c.text } },
                font: { color: c.text, size: 13, face: 'Inter' }, title: `${labelOfSymbol(n)} · ${degree.get(n) || 0} associations` };
        });
        const edges = d.edges.map((e) => {
            const shared = terms.filter((t) => t.genes.includes(e.a.toUpperCase()) && t.genes.includes(e.b.toUpperCase())).map((t) => t.term);
            return { from: e.a, to: e.b, width: 0.5 + 3 * Math.max(0, (e.score || 0.4) - 0.15),
                color: { color: c.muted, highlight: c.pathway, opacity: 0.8 },
                title: `${labelOfSymbol(e.a)} – ${labelOfSymbol(e.b)} (STRING score ${fmt(e.score)})\n` +
                    (shared.length ? `Shared terms: ${[...new Set(shared)].join(' | ')}` : 'No shared top term') };
        });
        state.stringNet = new vis.Network($('string-net'), { nodes, edges }, {
            physics: { barnesHut: { gravitationalConstant: -4000, springLength: 110 }, stabilization: { iterations: 250 } },
            interaction: { hover: true, tooltipDelay: 120 },
            nodes: { scaling: { label: { drawThreshold: 1 } } },
            edges: { smooth: false },
        });
        state.stringNet.once('stabilizationIterationsDone', () => state.stringNet.setOptions({ physics: false }));
        const inNet = new Set(d.nodes.map((n) => n.toUpperCase()));
        const rows = terms.map((t) => ({ t, genes: t.genes.filter((g) => inNet.has(g)) })).filter((r) => r.genes.length);
        $('string-map').innerHTML = '<thead><tr><th>Term</th><th>List</th><th class="num">Genes</th><th>Genes in the network</th><th class="num">FDR</th></tr></thead><tbody>' +
            rows.map((r) => `<tr><td>${esc(r.t.term)}</td><td>${esc(r.t.direction)}</td><td class="num">${r.genes.length}</td><td>${esc(geneLabelList(r.genes))}</td><td class="num">${fmtP(r.t.fdr)}</td></tr>`).join('') +
            '</tbody>';
        $('string-info').textContent = `STRING ${d.provenance.version}, ${d.provenance.network_type} network, score ≥ ${d.provenance.required_score / 1000}: ` +
            `${d.mapped} genes mapped${d.unmapped.length ? `, ${d.unmapped.length} not recognised` : ''}. Hover over an edge for the top terms both genes share. ` +
            'Edges are prior-knowledge associations, not inferred from this dataset.';
    }

    // ------------------------------------------------------------------ PrioriGraph
    function resetPG() {
        state.pg = null;
        state.pgProp = null;
        state.pgSens = null;
        state.pgAutoDone = false;
        state.pgRelaxed = null;
        $('pg-fdr').value = '0.05';
        $('pg-notice').classList.add('hidden');
        if (state.pgNet) { state.pgNet.destroy(); state.pgNet = null; }
        $('pg-msg').style.display = 'flex';
        $('pg-msg').textContent = 'Run the APP enrichment first.';
        $('pg-table').innerHTML = '';
        $('pg-tiles').innerHTML = '';
    }

    const pgQuery = () => new URLSearchParams({ terms: $('pg-n').value || 6, term_fdr: $('pg-term-fdr').value, tf_fdr: $('pg-fdr').value });

    async function fetchPG() {
        const r = await fetch(`/api/jobs/${state.appJob}/priorigraph?${pgQuery()}`);
        const d = await r.json();
        if (!r.ok) throw new Error(d.error || 'PrioriGraph failed.');
        return d;
    }

    async function loadPG() {
        if (!state.appJob) return;
        $('pg-msg').style.display = 'flex';
        $('pg-msg').textContent = 'Building the network…';
        try {
            let d = await fetchPG();
            if (!state.pgAutoDone && !d.nodes.length && d.reason === 'no_regulators') {
                // First view: no regulator reaches the default FDR. Show the least relaxed threshold that
                // gives a network, and say so (exploratory), instead of an empty page.
                const from = Number($('pg-fdr').value);
                const diag = d.diagnostics;
                for (const t of diag.thresholds.filter((x, k) => x > from && diag.regulator_counts[k] > 0)) {
                    $('pg-fdr').value = String(t);
                    const next = await fetchPG();
                    if (next.nodes.length) { d = next; state.pgRelaxed = { from, to: t }; break; }
                }
                if (!state.pgRelaxed) $('pg-fdr').value = String(from);
            }
            state.pgAutoDone = true;
            state.pg = d;
            state.pgProp = null;
            state.pgSens = null;
            renderPGNotice();
            drawPG();
            renderPGView();
        } catch (e) {
            $('pg-msg').textContent = e.message;
        }
    }
    $('pg-n').addEventListener('input', debounce(loadPG, 400));
    $('pg-term-fdr').addEventListener('change', loadPG);
    $('pg-fdr').addEventListener('change', () => { state.pgRelaxed = null; loadPG(); });

    const fmtT = (t) => (t >= 1 ? 'any' : String(t));

    // option labels with how many terms / regulators pass each threshold; notice explaining an empty network
    function renderPGNotice() {
        const d = state.pg, diag = d.diagnostics;
        if (!diag) return;
        const count = (list, t) => list[diag.thresholds.findIndex((x) => Math.abs(x - t) < 1e-9)];
        [...$('pg-fdr').options].forEach((o) => {
            const n = count(diag.regulator_counts, Number(o.value));
            o.textContent = `${fmtT(Number(o.value))} · ${n ?? '?'} regulator${n === 1 ? '' : 's'}`;
        });
        [...$('pg-term-fdr').options].forEach((o) => {
            const n = count(diag.term_counts, Number(o.value));
            o.textContent = `${fmtT(Number(o.value))} · ${n ?? '?'} term${n === 1 ? '' : 's'}`;
        });
        const tf = Number($('pg-fdr').value), terms = Number($('pg-term-fdr').value);
        const best = diag.best_regulator;
        const bestText = best ? ` (best: ${esc(labelOfSymbol(best.tf))}, FDR ${fmtP(best.fdr)}, ${esc(best.library)})` : '';
        const tfButtons = diag.thresholds.map((t, k) => ({ t, n: diag.regulator_counts[k] })).filter((x) => x.t > tf && x.n > 0)
            .map((x) => `<button type="button" class="link-btn" data-pg-set="pg-fdr:${x.t}">${x.t >= 1 ? 'All regulators' : `Regulators with FDR ≤ ${x.t}`} (${x.n})</button>`).join('');
        const termButtons = diag.thresholds.map((t, k) => ({ t, n: diag.term_counts[k] })).filter((x) => x.t > terms && x.n > 0)
            .map((x) => `<button type="button" class="link-btn" data-pg-set="pg-term-fdr:${x.t}">${x.t >= 1 ? 'All terms' : `Terms with FDR ≤ ${x.t}`} (${x.n})</button>`).join('');
        const strict = 0.05, nStrict = count(diag.regulator_counts, strict);
        let html = '';
        if (d.nodes.length && tf > strict) {
            html = '<strong>Exploratory network.</strong> ' + (nStrict
                ? `${nStrict} candidate regulator${nStrict === 1 ? ' has' : 's have'} FDR ≤ ${strict}; the network also includes `
                : `No candidate regulator reaches FDR ≤ ${strict} for this gene list${bestText}; the network includes `) +
                `regulators up to FDR ${fmtT(tf)} (${diag.regulators_passing}). Treat them as weak evidence.` +
                (state.pgRelaxed ? ' This threshold was chosen automatically because the default one gave an empty network.' : '') +
                `<br><button type="button" class="link-btn" data-pg-set="pg-fdr:${strict}">Back to FDR ≤ ${strict}</button>`;
        } else if (!d.nodes.length && d.reason === 'no_regulators') {
            html = `<strong>Empty network:</strong> no candidate regulator has FDR ≤ ${fmtT(tf)}${bestText}. ` +
                'Regulator enrichment needs larger gene lists; allow a higher regulator FDR (exploratory) or add regulator libraries in APP.' +
                `<br>${tfButtons}`;
        } else if (!d.nodes.length && d.reason === 'no_terms') {
            html = `<strong>Empty network:</strong> no enriched term has FDR ≤ ${fmtT(terms)}` +
                (diag.best_term_fdr !== null ? ` (best term FDR ${fmtP(diag.best_term_fdr)})` : '') + '.<br>' + termButtons;
        } else if (!d.nodes.length) {
            html = `<strong>Empty network:</strong> none of the genes in the top ${esc($('pg-n').value)} terms is a target of the ` +
                `${diag.regulators_passing} regulators with FDR ≤ ${fmtT(tf)}.<br>` +
                '<button type="button" class="link-btn" data-pg-set="pg-n:25">Use 25 terms</button>' + tfButtons;
        }
        $('pg-notice').innerHTML = html;
        $('pg-notice').classList.toggle('hidden', !html);
    }

    $('pg-notice').addEventListener('click', (e) => {
        const b = e.target.closest('button[data-pg-set]');
        if (!b) return;
        const [id, value] = b.dataset.pgSet.split(':');
        $(id).value = value;
        if (id === 'pg-fdr') state.pgRelaxed = null;
        loadPG();
    });
    $('pg-search').addEventListener('input', debounce(() => drawPG(true), 250));
    $('pg-fit').addEventListener('click', () => { if (state.pgNet) state.pgNet.fit({ animation: true }); });

    function drawPG(searchOnly) {
        const d = state.pg;
        if (!d) return;
        const c = netColors();
        const q = $('pg-search').value.trim().toUpperCase();
        const hit = (n) => q && (n.label.toUpperCase().includes(q) || (n.kind !== 'pathway' && labelOfSymbol(n.id).toUpperCase().includes(q)));
        const nodes = d.nodes.map((n) => {
            const h = hit(n);
            const stCol = n.status === 'up' ? c.up : n.status === 'down' ? c.down : c.neutral;
            if (n.kind === 'pathway') {
                const col = h ? c.hit : c.pathway;
                return { id: n.id, label: n.label.length > 32 ? n.label.slice(0, 30) + '…' : n.label, shape: 'hexagon', size: h ? 26 : 20,
                    color: { background: col, border: col }, font: { color: c.text, size: 14, face: 'Inter' },
                    title: `${n.label}\n${n.library} · FDR ${fmtP(n.fdr)} · ${n.directions.join(', ')}` };
            }
            if (n.kind === 'tf') {
                const col = h ? c.hit : c.tf;
                return { id: n.id, label: labelOfSymbol(n.id), shape: 'square', size: h ? 18 : 12,
                    color: { background: col, border: n.status ? stCol : col }, borderWidth: n.status ? 3 : 1,
                    font: { color: c.text, size: 15, face: 'Inter' },
                    title: `${labelOfSymbol(n.id)}: candidate regulator${n.status ? ` (itself ${n.status}-regulated)` : ''}` };
            }
            const col = h ? c.hit : stCol;
            return { id: n.id, label: labelOfSymbol(n.id), shape: 'dot', size: h ? 16 : 9,
                color: { background: col, border: col }, font: { color: c.text, size: 14, face: 'Inter' },
                title: `${labelOfSymbol(n.id)}${n.status ? ` (${n.status}-regulated)` : ''}` };
        });
        if (searchOnly && state.pgNet) {
            state.pgNodes.update(nodes);
            const first = d.nodes.find(hit);
            if (first) state.pgNet.focus(first.id, { scale: 1.1, animation: true });
            return;
        }
        const edges = d.edges.map((e) => (e.kind === 'regulation'
            ? { from: e.from, to: e.to, arrows: { to: { enabled: true, scaleFactor: 0.5 } }, color: { color: c.muted, opacity: 0.9 }, width: 1 }
            : { from: e.from, to: e.to, dashes: true, color: { color: c.pathway, opacity: 0.45 }, width: 1 }));
        if (state.pgNet) state.pgNet.destroy();
        $('pg-msg').style.display = d.nodes.length ? 'none' : 'flex';
        if (!d.nodes.length) {
            $('pg-msg').textContent = 'Empty network: see the explanation above.';
            state.pgNet = null;
            const diag = d.diagnostics || {};
            $('pg-tiles').innerHTML = tile('Enriched terms passing', int(diag.terms_passing ?? 0), `FDR ≤ ${fmtT(Number($('pg-term-fdr').value))}`) +
                tile('Regulators passing', int(diag.regulators_passing ?? 0), `FDR ≤ ${fmtT(Number($('pg-fdr').value))}`);
            return;
        }
        state.pgNodes = new vis.DataSet(nodes);
        state.pgNet = new vis.Network($('pg-net'), { nodes: state.pgNodes, edges: new vis.DataSet(edges) }, {
            physics: { barnesHut: { gravitationalConstant: -6000, springLength: 130, avoidOverlap: 0.2 }, stabilization: { iterations: 300 } },
            interaction: { hover: true, tooltipDelay: 120, navigationButtons: false },
            nodes: { scaling: { label: { drawThreshold: 1 } } },
            edges: { smooth: false },
        });
        state.pgNet.once('stabilizationIterationsDone', () => state.pgNet.setOptions({ physics: false }));
        const s = d.stats;
        $('pg-tiles').innerHTML = tile('Enriched terms', int(s.terms), `of ${int(s.terms_screened)} screened`) +
            tile('Candidate regulators', int(s.tfs), `of ${int(s.regulators_used)} passing the FDR`) +
            tile('Target genes', int(s.genes), 'selected genes in those terms') +
            tile('Edges', int(s.edges), esc(s.tf_libraries.join(', ')), false);
    }

    const statusText = (s) => (s === 'up' ? '<span class="status-up">up</span>' : s === 'down' ? '<span class="status-down">down</span>' : '<span class="hint">–</span>');
    const pgView = () => document.querySelector('input[name="pg_view"]:checked').value;
    const numCol = (key, lab, digits, kind = 'f') => ({ key, label: lab, num: true, value: (r) => num(r[key]),
        render: (r) => (kind === 'p' ? fmtP(r[key]) : fmt(r[key], digits)) });
    const PG_VIEWS = {
        regulators: { desc: 'Candidate regulators ranked by their unique target genes in the network (out-degree).', rows: () => state.pg.regulators, cols: () => [
            { key: 'tf', label: 'Regulator', value: (r) => labelOfSymbol(r.tf) }, numCol('targets_in_network', 'Targets in network', 0),
            numCol('targets_in_list', 'Targets among selected genes', 0), numCol('fdr', 'Best FDR', 0, 'p'),
            { key: 'libraries', label: 'Evidence', value: (r) => r.libraries }, { key: 'status', label: 'Itself DE', value: (r) => r.status || '', render: (r) => statusText(r.status) }] },
        targets: { desc: 'Target genes ranked by the number of regulators (in-degree) and enriched terms.', rows: () => state.pg.targets, cols: () => [
            { key: 'gene', label: 'Gene', value: (r) => labelOfSymbol(r.gene) }, { key: 'status', label: 'Status', value: (r) => r.status || '', render: (r) => statusText(r.status) },
            numCol('regulators', 'Regulators', 0), numCol('terms', 'Terms', 0), numCol('effect', 'Effect (log2FC)', 2), numCol('fdr', 'FDR', 0, 'p')] },
        mapping: { desc: 'Main candidate regulator of each enriched term: the factor with most targets among the term genes in the network.', rows: () => state.pg.mapping, cols: () => [
            { key: 'term', label: 'Term', value: (r) => r.term }, { key: 'library', label: 'Library', value: (r) => r.library },
            { key: 'regulator', label: 'Top regulator', value: (r) => labelOfSymbol(r.regulator) }, numCol('regulator_targets', 'Its targets in the term', 0),
            numCol('genes_in_network', 'Term genes in network', 0), numCol('fdr', 'Term FDR', 0, 'p')] },
        evidence: { desc: 'Expression evidence of the network genes in the gene source.', rows: () => state.pg.evidence, cols: () => [
            { key: 'gene', label: 'Gene', value: (r) => labelOfSymbol(r.gene) }, numCol('effect', 'Effect (log2FC)', 2), numCol('fdr', 'FDR', 0, 'p')] },
        metrics: { desc: 'HITS and PageRank on the directed regulator → target layer; betweenness on its undirected skeleton.', rows: () => state.pg.metrics, cols: () => [
            { key: 'node', label: 'Node', value: (r) => labelOfSymbol(r.node) }, { key: 'type', label: 'Type', value: (r) => r.type },
            numCol('out_degree', 'Out-degree', 0), numCol('in_degree', 'In-degree', 0), numCol('hits_hub', 'HITS hub', 4),
            numCol('hits_authority', 'HITS authority', 4), numCol('pagerank', 'PageRank', 4), numCol('betweenness', 'Betweenness', 4)] },
        propagation: { desc: 'Personalised random walk with restart seeded with |log2FC| (probability flows from regulators to targets and restarts at the seeds); signed linear propagation (I − λW)⁻¹E. Heuristic scores, not causal effects.',
            rows: () => (state.pgProp ? state.pgProp.rows : []), cols: () => [
                { key: 'gene', label: 'Node', value: (r) => labelOfSymbol(r.gene) }, numCol('effect', 'log2FC (seed)', 2), numCol('rwr', 'RWR score', 4), numCol('linear', 'Linear propagation', 3)] },
        sensitivity: { desc: 'Stability of the RWR rankings over the parameter grid: how often each node is in the top k and how much its rank varies. A sensitivity analysis, not a bootstrap or validation.',
            rows: () => (state.pgSens ? state.pgSens.rows : []), cols: () => [
                { key: 'node', label: 'Node', value: (r) => labelOfSymbol(r.node) }, numCol('topk', 'Top-k frequency', 2), numCol('coverage', 'Configuration coverage', 2),
                numCol('median_rank', 'Median rank', 1), numCol('iqr_rank', 'Rank IQR', 1), numCol('mean_score', 'Mean RWR score', 4), numCol('median_linear', 'Median linear propagation', 3)] },
    };

    const pgTable = new DataTable($('pg-table'), {
        columns: () => (state.pg ? PG_VIEWS[pgView()].cols() : []),
        rowKey: (r) => JSON.stringify(r),
        text: (r) => Object.values(r).join(' ').toLowerCase(),
        empty: 'Nothing to show yet.',
    });

    function renderPGView() {
        const v = pgView();
        $('pg-view-desc').textContent = PG_VIEWS[v].desc;
        $('pg-prop-controls').classList.toggle('hidden', v !== 'propagation');
        $('pg-sens-controls').classList.toggle('hidden', v !== 'sensitivity');
        if (!state.pg) { pgTable.setRows([]); return; }
        if (v === 'propagation' && !state.pgProp) loadPropagation();
        pgTable.sortKey = null;
        pgTable.setRows(PG_VIEWS[v].rows());
        if (v === 'sensitivity' && state.pgSens) $('pg-view-desc').textContent += ` ${state.pgSens.configurations} configurations evaluated.`;
        if (v === 'propagation' && state.pgProp) $('pg-view-desc').textContent += ` Converged: ${state.pgProp.converged ? 'yes' : 'no'} (${state.pgProp.iterations} iterations).`;
    }
    document.querySelectorAll('input[name="pg_view"]').forEach((el) => el.addEventListener('change', renderPGView));

    async function loadPropagation() {
        if (!state.appJob) return;
        const p = pgQuery();
        p.set('restart', $('pg-restart').value);
        p.set('decay', $('pg-decay').value);
        try {
            const r = await fetch(`/api/jobs/${state.appJob}/propagation?${p}`);
            const d = await r.json();
            if (!r.ok) throw new Error(d.error);
            state.pgProp = d;
        } catch (e) {
            state.pgProp = { rows: [], converged: false, iterations: 0 };
        }
        if (pgView() === 'propagation') renderPGView();
    }
    const onPropSlider = debounce(loadPropagation, 300);
    ['pg-restart', 'pg-decay'].forEach((id) => $(id).addEventListener('input', () => {
        $(`${id}-val`).textContent = Number($(id).value).toFixed(2);
        onPropSlider();
    }));

    const GRID = { fdrs: [[0.01, 0.025, 0.05, 0.1], [0.01, 0.05, 0.1]], effects: [[0, 0.5, 1, 1.5, 2], [0.5, 1, 1.5]],
        restarts: [[0.3, 0.5, 0.7, 0.9], [0.5, 0.7, 0.9]], decays: [[0.25, 0.5, 0.75], [0.25, 0.5, 0.75]] };
    document.querySelectorAll('[data-grid]').forEach((el) => {
        const [values, defaults] = GRID[el.dataset.grid];
        el.innerHTML = values.map((v) => `<label class="pill"><input type="checkbox" value="${v}"${defaults.includes(v) ? ' checked' : ''}>${v}</label>`).join('');
    });

    $('pg-sens-btn').addEventListener('click', async () => {
        const grid = {};
        document.querySelectorAll('[data-grid]').forEach((el) => {
            grid[el.dataset.grid] = [...el.querySelectorAll('input:checked')].map((x) => Number(x.value));
        });
        if (Object.values(grid).some((v) => !v.length)) {
            $('pg-view-desc').textContent = 'Choose at least one value in every row of the grid.';
            return;
        }
        $('pg-sens-btn').disabled = true;
        $('pg-view-desc').textContent = 'Evaluating the configurations…';
        try {
            const r = await fetch(`/api/jobs/${state.appJob}/sensitivity`, {
                method: 'POST', headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(Object.assign(grid, { tf_fdr: Number($('pg-fdr').value), top_k: Number($('pg-topk').value) || 20 })),
            });
            const d = await r.json();
            if (!r.ok) throw new Error(d.error);
            state.pgSens = d;
        } catch (e) {
            state.pgSens = { rows: [], configurations: 0 };
        }
        $('pg-sens-btn').disabled = false;
        renderPGView();
    });

    $('pg-export').addEventListener('click', () => { if (state.pg) tableCsv(pgTable, `lumos_PG_${pgView()}.csv`); });
    $('pg-edges').addEventListener('click', () => {
        if (!state.pg) return;
        downloadCsv('lumos_PG_edges.csv', ['Regulator', 'Target', 'Library', 'Direction', 'Term', 'FDR'],
            state.pg.edge_table.map((e) => [e.tf, e.target, e.library, e.direction, e.term, e.fdr]));
    });

    // ------------------------------------------------------------------ HTML reports
    function tableSpec(dt, limit = 50000) {
        const { cols, rows } = dt.exportRows(limit);
        const fmtOf = (c, k) => {
            if (!c.num) return 'text';
            if (/FDR|p-value|τ²|RWR score|PageRank|HITS/.test(c.label)) return 'exp2';
            return rows.every((r) => r[k] === null || Number.isInteger(r[k])) ? 'int' : 'fixed3';
        };
        const fdrCol = cols.find((c) => /FDR/.test(c.label));
        return { columns: cols.map((c, k) => [c.key, c.label, fmtOf(c, k)]), rows, sort: (fdrCol || cols[0]).key, asc: true };
    }

    function reportFig(fig) {
        const f = JSON.parse(JSON.stringify(fig));
        f.layout = Object.assign({}, f.layout, { paper_bgcolor: '#ffffff', plot_bgcolor: '#ffffff', height: fig.height || f.layout.height || 480 });
        return { data: f.data, layout: f.layout };
    }

    function pgNetworkFigure() {
        if (!state.pgNet || !state.pg) return null;
        const pos = state.pgNet.getPositions();
        const c = { up: '#eb6834', down: '#2a78d6', tf: '#c43a6f', pathway: '#96700a', muted: '#9aa0aa' };
        const lines = { regulation: { x: [], y: [] }, pathway: { x: [], y: [] } };
        state.pg.edges.forEach((e) => {
            const a = pos[e.from], b = pos[e.to];
            if (!a || !b) return;
            lines[e.kind].x.push(a.x, b.x, null);
            lines[e.kind].y.push(-a.y, -b.y, null);
        });
        const groups = { pathway: [], tf: [], up: [], down: [], other: [] };
        state.pg.nodes.forEach((n) => {
            const g = n.kind === 'pathway' ? 'pathway' : n.kind === 'tf' ? 'tf' : n.status === 'up' ? 'up' : n.status === 'down' ? 'down' : 'other';
            if (pos[n.id]) groups[g].push(n);
        });
        const style = { pathway: ['hexagon', 16, c.pathway, 'Enriched term'], tf: ['square', 11, c.tf, 'Candidate regulator'], up: ['circle', 9, c.up, 'Up-regulated gene'],
            down: ['circle', 9, c.down, 'Down-regulated gene'], other: ['circle', 9, c.muted, 'Gene'] };
        const data = [
            { type: 'scatter', mode: 'lines', x: lines.pathway.x, y: lines.pathway.y, line: { color: c.pathway, width: 1, dash: 'dot' }, hoverinfo: 'skip', name: 'term membership' },
            { type: 'scatter', mode: 'lines', x: lines.regulation.x, y: lines.regulation.y, line: { color: c.muted, width: 1 }, hoverinfo: 'skip', name: 'regulator → target' },
        ];
        Object.entries(groups).forEach(([g, list]) => {
            if (!list.length) return;
            const [symbol, size, color, name] = style[g];
            data.push({ type: 'scatter', mode: 'markers+text', name, x: list.map((n) => pos[n.id].x), y: list.map((n) => -pos[n.id].y),
                text: list.map((n) => (n.kind === 'pathway' ? (n.label.length > 30 ? n.label.slice(0, 28) + '…' : n.label) : labelOfSymbol(n.id))),
                hovertext: list.map((n) => n.label), hoverinfo: 'text', textposition: 'top center', textfont: { size: 9 },
                marker: { symbol, size, color, line: { color: '#ffffff', width: 0.5 } } });
        });
        return { height: 720, data, layout: { showlegend: true, legend: { orientation: 'h' }, hovermode: 'closest', margin: { l: 10, r: 10, t: 10, b: 10 },
            xaxis: { visible: false }, yaxis: { visible: false }, font: { family: 'Inter, sans-serif', size: 11, color: '#333' } } };
    }

    function reportPayload(tab) {
        if (tab === 'ddea' || tab === 'files') {
            const s = state.de;
            const cmp = s.comparison;
            let up = 0, down = 0;
            state.status.forEach((x) => { if (x > 0) up++; else if (x < 0) down++; });
            const figs = withLightTheme(() => ({ v: figVolcano(), ma: figMA(), pca: figPCA(), u: figTop(1), d: figTop(-1),
                h: figHeatmap(heatmapRows(), heatOrder), lib: s.library_sizes ? figLib() : null }));
            const genes = tableSpec(degTable);
            return { module: 'DDEA', title: `DDEA report: ${s.name} · ${s.label}`, file_id: `${s.name}_${s.label}`,
                info: [['Dataset', `${s.name}${s.title ? ` – ${s.title}` : ''}`], ['Comparison', `${cmp.test} (${cmp.test_samples.length} samples) vs ${cmp.reference} (${cmp.reference_samples.length} samples, reference)`],
                    ['Method', s.method_text], ['Design formula', s.design.formula], ['Analysis unit', s.unit_text], ['Filter', s.counts.rule],
                    ['Values for plots', s.transform_text], ['Thresholds', `FDR < ${state.fdr}, |log2FC| ≥ ${state.lfc}`],
                    ['Gene labels', state.labels === 'symbol' && symbolsAvailable(s) ? `gene symbols (${s.symbols.source})` : 'original identifiers']],
                cards: [['Genes tested', s.counts.tested, ''], ['Up-regulated', up, 'up'], ['Down-regulated', down, 'down']],
                sections: [
                    { id: 'volcano', title: 'Volcano plot', desc: 'log2 fold change against adjusted significance; dashed lines mark the thresholds.', figs: [reportFig(figs.v)] },
                    { id: 'ma', title: 'MA plot', desc: 'Fold change against mean expression.', figs: [reportFig(figs.ma)] },
                    { id: 'pca', title: 'PCA', desc: $('pca-sub').textContent, figs: [reportFig(figs.pca)] },
                    { id: 'top', title: 'Top differentially expressed genes', desc: 'The 20 most significant up- and down-regulated genes.', figs: [reportFig(figs.u), reportFig(figs.d)] },
                    { id: 'heatmap', title: 'Heatmap', desc: $('heatmap-sub').textContent, figs: [reportFig(figs.h)] },
                    ...(figs.lib ? [{ id: 'libraries', title: 'Library sizes', desc: 'Total counts per sample.', figs: [reportFig(figs.lib)] }] : []),
                    { id: 'genes', title: 'Gene table', desc: $('deg-all').checked ? 'Every tested gene.' : 'Genes passing both thresholds.',
                        table: Object.assign(genes, { filter: 'status', badge: 'status', badge_map: { 'Up-regulated': 'up', 'Down-regulated': 'down' } }) },
                ] };
        }
        if (tab === 'int') {
            const s = state.int;
            if (!s) throw new Error('Run an integration first.');
            const k = Number($('forest-gene').value);
            const figs = withLightTheme(() => (s.mode === 'meta'
                ? [['volcano', 'Meta-analysis volcano plot', figMetaVolcano()], ['forest', `Forest plot: ${s.table.symbol[k] || ''}`, figForest(k)],
                    ['i2', 'Heterogeneity (I²)', figI2()], ['heat', 'Expression of the top genes', figMetaHeat()]]
                : [['patterns', 'Patterns', figPatterns()], ['pair', 'Pairwise effects', figPair()], ['effects', 'Effect heatmap', figEffects()]]));
            return { module: 'Integration', title: `${MODE_NAMES[s.mode]}: ${s.profiles.map((p) => p.label).join(' · ')}`, file_id: s.mode,
                info: [['Integration', s.mode_text], ['Analyses', s.profiles.map((p) => `${p.label} (${p.dataset}: ${p.comparison}; ${p.method})`).join(' | ')],
                    ['Thresholds', `FDR < ${s.thresholds.fdr}, |log2FC| ≥ ${s.thresholds.lfc}`], ...(s.mode === 'meta' ? [['Estimator', s.estimator], ['Minimum studies per gene', s.min_studies]] : [])],
                cards: Object.entries(s.pattern_counts).filter(([key]) => key !== 'Not significant').slice(0, 6).map(([key, v]) => [key, v, /Up/.test(key) ? 'up' : /Down/.test(key) ? 'down' : '']),
                sections: figs.map(([id, title, f]) => ({ id, title, figs: [reportFig(f)] })).concat([{ id: 'genes', title: 'Genes', desc: $('int-pattern').selectedOptions[0].textContent, table: tableSpec(intGenes) }]) };
        }
        if (tab === 'app') {
            const s = state.app;
            if (!s) throw new Error('Run the enrichment first.');
            const libs = s.libraries.filter((l) => l.role === 'pathway');
            const figs = withLightTheme(() => ({ libs: libs.map(figLibrary), tf: s.regulators.length ? figRegulators() : null }));
            const termTable = (l) => ({ columns: [['term', 'Term', 'text'], ['direction', 'List', 'text'], ['n', 'Genes', 'int'], ['p', 'p-value', 'exp2'], ['fdr', 'FDR', 'exp2'], ['genes', 'Genes in the term', 'text']],
                rows: l.terms.map((t) => [t.term, t.direction, t.genes.length, t.p, t.fdr, t.genes.join(', ')]), sort: 'fdr', asc: true, filter: 'direction' });
            const sections = libs.map((l, k) => ({ id: `lib-${k}`, title: l.name, desc: `${l.engine}; ${l.n_significant} of ${l.n_terms} terms with FDR < 0.05.`, figs: [reportFig(figs.libs[k])], table: termTable(l) }));
            if (figs.tf) {
                sections.push({ id: 'regulators', title: 'Candidate regulators', desc: $('tf-sub').textContent, figs: [reportFig(figs.tf)],
                    table: { columns: [['tf', 'Regulator', 'text'], ['term', 'Term', 'text'], ['library', 'Library', 'text'], ['direction', 'List', 'text'], ['fdr', 'FDR', 'exp2'], ['targets', 'Targets', 'text']],
                        rows: s.regulators.map((r) => [r.tf, r.term, r.library, r.direction, r.fdr, r.targets.join(', ')]), sort: 'fdr', asc: true, filter: 'library' } });
            }
            if (state.stringData && state.stringData.edges.length) {
                sections.push({ id: 'string', title: 'STRING associations', desc: $('string-info').textContent,
                    table: { columns: [['a', 'Gene A', 'text'], ['b', 'Gene B', 'text'], ['score', 'Score', 'fixed3']], rows: state.stringData.edges.map((e) => [labelOfSymbol(e.a), labelOfSymbol(e.b), e.score]), sort: 'score', asc: false } });
            }
            return { module: 'APP', title: `APP enrichment: ${s.source.label}`, file_id: 'APP',
                info: [['Gene source', `${s.source.label}${s.source.dataset ? ` (${s.source.dataset})` : ''}`], ['Gene lists', s.strategy === 'separate' ? 'up- and down-regulated separately' : 'combined'],
                    ['Genes', `${s.genes.filter((g) => g.direction === 'up').length} up, ${s.genes.filter((g) => g.direction === 'down').length} down`],
                    ['Background', `${s.background.source} (${s.background.size} genes)`], ['Libraries', s.libraries.map((l) => `${l.name} (${l.engine})`).join('; ')],
                    ['Thresholds', `FDR < ${s.thresholds.fdr}, |log2FC| ≥ ${s.thresholds.lfc}`]],
                cards: [['Genes queried', s.genes.length, ''], ['Background', s.background.size, ''], ...s.libraries.map((l) => [l.name, l.n_significant, ''])].slice(0, 6),
                sections };
        }
        if (tab === 'pg') {
            if (!state.pg) throw new Error('Build the PrioriGraph first.');
            const d = state.pg;
            const sections = [{ id: 'network', title: 'Integrated network', desc: 'Hexagons: enriched terms; squares: candidate regulators; circles: target genes (orange up, blue down). Dotted edges: term membership; grey: regulator → target.',
                figs: [pgNetworkFigure()].filter(Boolean).map(reportFig) }];
            const prev = pgView();
            Object.keys(PG_VIEWS).forEach((v) => {
                const rows = PG_VIEWS[v].rows();
                if (!rows.length) return;
                document.querySelector(`input[name="pg_view"][value="${v}"]`).checked = true;
                pgTable.setRows(rows);
                sections.push({ id: v, title: v[0].toUpperCase() + v.slice(1), desc: PG_VIEWS[v].desc, table: tableSpec(pgTable) });
            });
            document.querySelector(`input[name="pg_view"][value="${prev}"]`).checked = true;
            renderPGView();
            sections.push({ id: 'edges', title: 'Edges with provenance', table: { columns: [['tf', 'Regulator', 'text'], ['target', 'Target', 'text'], ['library', 'Library', 'text'], ['direction', 'List', 'text'], ['term', 'Term', 'text'], ['fdr', 'FDR', 'exp2']],
                rows: d.edge_table.map((e) => [e.tf, e.target, e.library, e.direction, e.term, e.fdr]), sort: 'fdr', asc: true, filter: 'library' } });
            return { module: 'PG', title: `PrioriGraph: ${state.app.source.label}`, file_id: 'PG',
                info: [['Gene source', state.app.source.label], ['Enriched terms', `top ${$('pg-n').value} with FDR ≤ ${$('pg-term-fdr').value}`], ['Regulators', `FDR ≤ ${$('pg-fdr').value} (${d.stats.tf_libraries.join(', ')})`],
                    ...(Number($('pg-fdr').value) > 0.05 ? [['Note', `exploratory: includes candidate regulators with FDR up to ${$('pg-fdr').value}`]] : []),
                    ...(state.pgProp ? [['Propagation', `restart ${state.pgProp.restart}, decay ${state.pgProp.decay}`]] : [])],
                cards: [['Nodes', d.stats.nodes, ''], ['Edges', d.stats.edges, ''], ['Terms', d.stats.terms, ''], ['Regulators', d.stats.tfs, ''], ['Targets', d.stats.genes, '']],
                sections };
        }
        throw new Error('Nothing to report on this tab.');
    }

    $('report-btn').addEventListener('click', async () => {
        const btn = $('report-btn');
        const tab = activeTab();
        const errorEl = tab === 'app' || tab === 'pg' ? $('app-error') : tab === 'int' ? $('int-error') : $('error-message');
        let payload;
        try {
            payload = reportPayload(tab);
        } catch (e) {
            showError(e.message, errorEl);
            return;
        }
        btn.disabled = true;
        btn.textContent = 'Building the report…';
        try {
            const r = await fetch('/api/report', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) });
            if (!r.ok) throw new Error((await r.json()).error || 'The report could not be built.');
            const name = (r.headers.get('Content-Disposition') || '').match(/filename="([^"]+)"/);
            downloadBlob(name ? name[1] : 'Lumos_report.html', await r.blob());
        } catch (e) {
            showError(e.message, errorEl);
        }
        btn.disabled = false;
        btn.textContent = 'HTML report';
    });

    // ------------------------------------------------------------------ files
    function renderFiles() {
        const rows = [];
        const add = (job, file, desc) => rows.push(`<div class="file-row"><div><a href="/api/jobs/${job}/files/${encodeURIComponent(file)}?download=1">${esc(file)}</a>` +
            `<div class="file-desc">${esc(desc)}</div></div></div>`);
        const head = (text) => rows.push(`<h3 class="section-title" style="margin: 14px 0 4px;">${esc(text)}</h3>`);
        if (state.dataset && !state.dataset.needs_expression) {
            head(`Dataset · ${state.dataset.name}`);
            add(state.dsJob, 'expression.tsv.gz', 'Expression table as read (genes × columns), duplicated IDs merged.');
            add(state.dsJob, 'samples.csv', 'Sample table with the matched expression column.');
        }
        if (state.de) {
            head(`Analysis · ${state.de.label}`);
            add(state.deJob, 'de_results.csv', 'Every tested gene: ID, symbol, mean, log2FC, SE, statistic, p-value, FDR and status (thresholds of the run).');
            add(state.deJob, 'normalized_expression.tsv.gz', `Values used for the plots (${state.de.transform_text}).`);
            add(state.deJob, 'manifest.json', 'Provenance: dataset, samples, model formula and design metadata, method, thresholds, versions.');
            add(state.deJob, state.de.files.zip, 'The files above plus the groups, summary and run log.');
        }
        if (state.int) {
            head(`Integration · ${MODE_NAMES[state.int.mode]}`);
            add(state.intJob, 'integration.csv', 'Complete gene-level table of the integration.');
            add(state.intJob, 'manifest.json', 'Analyses (with their manifests), rule, thresholds and versions.');
            add(state.intJob, state.int.files.zip, 'All files of the integration.');
        }
        if (state.app) {
            head('APP · enrichment');
            state.app.libraries.forEach((l) => { if (l.file) add(state.appJob, l.file, `Every ${l.name} term with its list, engine and overlapping genes.`); });
            if (state.app.files.regulators) add(state.appJob, 'regulators.csv', 'Candidate regulators: term, library, list, FDR and targets.');
            add(state.appJob, 'evidence.csv', 'Expression evidence (effect, FDR) of the gene source.');
            add(state.appJob, 'background.txt', 'The explicit background (universe) of the analysis.');
            add(state.appJob, 'manifest.json', 'Query genes, background, libraries, engines and versions.');
            add(state.appJob, state.app.files.zip, 'All enrichment files.');
        }
        $('files-list').innerHTML = rows.join('');
    }

    window.addEventListener('resize', debounce(resizeAll, 200));

    // ------------------------------------------------------------------ restore the session from the address bar
    async function fetchJob(id) {
        if (!/^[0-9a-f]{12}$/.test(id || '')) return null;
        try {
            const r = await fetch(`/api/jobs/${id}`);
            return r.ok ? await r.json() : null;
        } catch (e) {
            return null;
        }
    }

    (async function restore() {
        const p = new URLSearchParams(window.location.search);
        const list = (k) => (p.get(k) || '').split(',').filter(Boolean);
        for (const id of list('ds')) {
            const d = await fetchJob(id);
            if (d && d.state === 'done') state.datasets.set(id, d.summary);
            else if (d && d.state === 'running') datasetJob.follow(id, 'Reading the dataset…');
        }
        const dsActive = state.datasets.has(p.get('d')) ? p.get('d') : [...state.datasets.keys()].pop();
        if (dsActive) activateDataset(dsActive, false);
        for (const id of list('de')) {
            const d = await fetchJob(id);
            if (d && d.state === 'done') { state.analyses.set(id, d.summary); geneCache.delete(id); }
            else if (d && d.state === 'running') deJob.follow(id, 'Differential expression…');
        }
        if (state.analyses.size) {
            const a = state.analyses.has(p.get('a')) ? p.get('a') : [...state.analyses.keys()].pop();
            const s = state.analyses.get(a);
            if (state.datasets.has(s.dataset_job) && !store.get(`lumos-groups-${s.dataset_job}`)) {
                if (state.dsJob !== s.dataset_job) activateDataset(s.dataset_job, false);
                const c = s.comparison;
                loadGroups({ groups: [{ name: c.reference, color: 0, ids: c.reference_samples }, { name: c.test, color: 1, ids: c.test_samples }], reference: c.reference, test: c.test });
                sampleTable.render();
            }
            intTable.setRows([...state.analyses.keys()]);
            activateAnalysis(a, false);
        }
        for (const id of list('int')) {
            const d = await fetchJob(id);
            if (d && d.state === 'done') { state.integrations.set(id, d.summary); activateIntegration(id, false); }
            else if (d && d.state === 'running') intJob.follow(id, 'Integration…');
        }
        renderAppSources();
        const app = await fetchJob(p.get('app'));
        if (app && app.state === 'running') appJob.follow(p.get('app'), 'Enrichment…');
        if (app && app.state === 'done') {
            const src = app.summary.source;
            const value = src.kind === 'de' ? `de|${src.job}|analysis` : `int|${src.job}|${src.focus}`;
            if ([...$('app-source').options].some((o) => o.value === value)) { $('app-source').value = value; buildAppSource(true); }
            showApp(p.get('app'), app.summary, false);
        }
        renderDatasetBar();
        refreshButtons();
    }());

    refreshButtons();
});
