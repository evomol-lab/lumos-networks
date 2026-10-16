import base64
import html
import io
import json
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
from plotly.offline import get_plotlyjs

REPORT_CSS = """*{box-sizing:border-box}body{margin:0;font:14px/1.5 -apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;color:#222;background:#f5f6f8}
nav{position:fixed;top:0;left:0;right:0;height:56px;background:#1f2933;color:#fff;display:flex;align-items:center;gap:12px;padding:0 20px;z-index:100}
nav img{height:36px}nav h1{font-size:18px;margin:0;flex:1}nav small{color:#cbd2d9}
.nav-actions{display:flex;gap:8px}
.nav-btn{background:#323f4b;color:#fff;border:1px solid #486581;padding:6px 12px;border-radius:4px;cursor:pointer;font-size:13px}
.nav-btn:hover{background:#486581}
.wrap{display:flex;margin-top:56px}aside{position:sticky;top:56px;height:calc(100vh - 56px);width:220px;flex-shrink:0;padding:20px 12px;background:#fff;border-right:1px solid #e4e7eb;overflow-y:auto}
aside a{display:block;padding:6px 10px;color:#3e4c59;text-decoration:none;border-radius:4px;margin-bottom:2px}
aside a:hover{background:#eef1f4}aside a.active{font-weight:700;background:#eef1f4;color:#1f2933}aside a.muted{opacity:.45}
main{flex:1;min-width:0;padding:20px 28px}
section{background:#fff;border:1px solid #e4e7eb;border-radius:6px;padding:16px 20px;margin-bottom:20px;scroll-margin-top:70px}
h2{margin:0 0 4px;font-size:20px}
.tog{background:none;border:0;padding:0;font:inherit;color:inherit;cursor:pointer;display:flex;align-items:center;gap:8px;width:100%;text-align:left}
.chevron{display:inline-block;transition:transform .2s;font-size:12px;color:#616e7c}
section.collapsed .chevron{transform:rotate(-90deg)}section.collapsed .sbody{display:none}
.desc{margin:0 0 14px;color:#616e7c}
table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:6px 10px;border-bottom:1px solid #e4e7eb}
.kv{max-width:640px}.kv th{width:180px;color:#52606d;font-weight:600}.cards{display:flex;gap:14px;margin-top:16px;flex-wrap:wrap}
.card{flex:1;min-width:150px;padding:14px;border-radius:6px;background:#eef1f4;color:#3e4c59}.card span{display:block;font-size:28px;font-weight:700}
.card.up{background:#fde8e8;color:#c81e1e}.card.down{background:#e1effe;color:#1c64f2}.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px}
.tw{max-height:480px;overflow:auto}.dtbl thead th{position:sticky;top:0;background:#f5f6f8;user-select:none}
th.sortable{cursor:pointer}th.sortable:hover{background:#eef1f4}.sort-icon{font-size:10px;margin-left:4px;color:#7b8794}
.badge{padding:2px 8px;border-radius:10px;font-size:12px;color:#fff;background:#7b8794}.badge.up{background:#c81e1e}.badge.down{background:#1c64f2}
.tbl-controls{display:flex;gap:10px;align-items:center;margin-bottom:12px;flex-wrap:wrap}
.tbl-controls input,.tbl-controls select{padding:6px 10px;border:1px solid #cbd2d9;border-radius:4px;font-size:13px}
.tbl-controls input{flex:1;min-width:180px}
.btn-sub{background:#eef1f4;color:#3e4c59;border:1px solid #cbd2d9;padding:6px 12px;border-radius:4px;cursor:pointer;font-size:13px}.btn-sub:hover{background:#e4e7eb}
.tcount{margin-left:auto;color:#616e7c;font-size:12px}
.pagination{display:flex;gap:6px;align-items:center;justify-content:flex-end;margin-top:12px}
.pagination button{padding:4px 10px;border:1px solid #cbd2d9;background:#fff;border-radius:4px;cursor:pointer;font-size:12px}.pagination button:disabled{opacity:.4;cursor:not-allowed}
.print-only{display:none}
footer{text-align:center;color:#7b8794;font-size:12px;padding:10px 0 30px}
@media(max-width:900px){aside{display:none}.grid{grid-template-columns:1fr}main{padding:14px}nav small{display:none}}
@page{size:A4;margin:12mm}
@media print{
*{-webkit-print-color-adjust:exact;print-color-adjust:exact}
body{background:#fff}
nav{position:static;height:auto;padding:0 0 8px;margin-bottom:12px;border-bottom:2px solid #1f2933}
.nav-actions,aside,.chevron,.tbl-controls,.pagination{display:none!important}
.tog{cursor:default}
.wrap{display:block;margin-top:0}main{padding:0;width:100%}
section{border:0;padding:0}section.collapsed{display:none}
.tw{max-height:none;overflow:visible}.dtbl thead th{position:static}
.grid{grid-template-columns:1fr}
.plot,.card,.cards,tr{break-inside:avoid;page-break-inside:avoid}
.print-only{display:block;font-size:12px;color:#52606d;margin-bottom:8px}
}"""

REPORT_JS = r"""(function () {
  'use strict';
  var hasPlotly = typeof Plotly !== 'undefined';
  var plotsIn = function (root) { return Array.prototype.slice.call(root.querySelectorAll('.js-plotly-plot')); };
  var tocLink = function (sec) { return document.querySelector('aside a[href="#' + sec.id + '"]'); };

  function setCollapsed(sec, collapsed) {
    sec.classList.toggle('collapsed', collapsed);
    var btn = sec.querySelector('.tog');
    if (btn) btn.setAttribute('aria-expanded', String(!collapsed));
    var a = tocLink(sec);
    if (a) a.classList.toggle('muted', collapsed);
    if (!collapsed && hasPlotly) plotsIn(sec).forEach(function (gd) { Plotly.Plots.resize(gd); });
  }

  var sections = Array.prototype.slice.call(document.querySelectorAll('main section'));
  sections.forEach(function (sec) {
    var btn = sec.querySelector('.tog');
    if (btn) btn.addEventListener('click', function () { setCollapsed(sec, !sec.classList.contains('collapsed')); });
  });
  document.querySelectorAll('aside a').forEach(function (a) {
    a.addEventListener('click', function () {
      var sec = document.getElementById(a.getAttribute('href').slice(1));
      if (sec && sec.classList.contains('collapsed')) setCollapsed(sec, false);
    });
  });

  // Active TOC: a section is "in view" while it crosses a band just below the fixed nav.
  if ('IntersectionObserver' in window) {
    var visible = {};
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { visible[e.target.id] = e.isIntersecting; });
      var current = sections.filter(function (s) { return visible[s.id]; })[0];
      if (!current) return;
      document.querySelectorAll('aside a.active').forEach(function (a) { a.classList.remove('active'); });
      var a = tocLink(current);
      if (a) a.classList.add('active');
    }, { rootMargin: '-60px 0px -55% 0px', threshold: 0 });
    sections.forEach(function (s) { io.observe(s); });
  }

  // ---------------- Interactive tables (independent state per div.tbl[data-tid]) ----------------
  var collator = new Intl.Collator(undefined, { numeric: true, sensitivity: 'base' });
  var FMT = {
    fixed3: function (v) { return v.toFixed(3); },
    exp2: function (v) { return v.toExponential(2); },
    int: function (v) { return String(Math.round(v)); }
  };
  function esc(s) {
    return String(s == null ? '' : s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
  }

  function initTable(root) {
    var tid = root.getAttribute('data-tid');
    var el = function (s) { return document.getElementById(tid + '-' + s); };
    var dataEl = el('data');
    if (!dataEl) return null;
    var cfg = JSON.parse(dataEl.textContent), DATA = cfg.rows, COLS = cfg.cols;
    var badges = cfg.badgeMap || {}, hasOwn = Object.prototype.hasOwnProperty;
    var TEXT = [];
    COLS.forEach(function (c, i) { if (c.f === 'text' && i !== cfg.filter) TEXT.push(i); });
    var searchEl = el('search'), filterEl = el('filter'), sizeEl = el('size'), countEl = el('count');
    var tbody = el('tbody'), pagEl = el('pag'), summaryEl = el('summary'), csvEl = el('csv'), tsvEl = el('tsv');
    var ths = Array.prototype.slice.call(root.querySelectorAll('th.sortable'));
    var DEFAULTS = { search: '', filter: 'All', size: 25, page: 1, sortCol: cfg.sort, sortAsc: cfg.asc };
    var state = Object.assign({}, DEFAULTS);

    function view() {
      var q = state.search.toLowerCase(), col = state.sortCol, dir = state.sortAsc ? 1 : -1;
      return DATA.filter(function (r) {
        return (cfg.filter == null || state.filter === 'All' || r[cfg.filter] === state.filter) &&
          (!q || TEXT.some(function (i) { return r[i] != null && String(r[i]).toLowerCase().indexOf(q) !== -1; }));
      }).sort(function (a, b) {
        var x = a[col], y = b[col];
        if (x == null || y == null) return (x == null) - (y == null); // nulls last in both directions
        return dir * (typeof x === 'number' && typeof y === 'number' ? x - y : collator.compare(String(x), String(y)));
      });
    }

    function cell(r, i) {
      var v = r[i];
      if (v == null) return '';
      if (i === cfg.badge) return '<span class="badge ' + esc(hasOwn.call(badges, v) ? badges[v] : '') + '">' + esc(v) + '</span>';
      return typeof v === 'number' && FMT[COLS[i].f] ? FMT[COLS[i].f](v) : esc(v);
    }

    function render() {
      var rows = view(), n = rows.length;
      var size = state.size === 'All' ? Math.max(n, 1) : state.size;
      var pages = Math.max(1, Math.ceil(n / size));
      state.page = Math.min(Math.max(1, state.page), pages);
      var start = (state.page - 1) * size, end = Math.min(start + size, n), shown = rows.slice(start, end);

      tbody.innerHTML = shown.map(function (r) {
        return '<tr>' + COLS.map(function (c, i) { return '<td>' + cell(r, i) + '</td>'; }).join('') + '</tr>';
      }).join('');

      if (countEl) countEl.textContent = 'Showing ' + (n ? start + 1 : 0) + '–' + end + ' of ' + n + ' (filtered from ' + DATA.length + ')';
      var arrow = state.sortAsc ? '▲' : '▼';
      if (summaryEl) {
        var parts = [];
        if (cfg.filter != null) parts.push('Filter: ' + state.filter);
        if (state.search) parts.push('search: "' + state.search + '"');
        parts.push('sorted by ' + COLS[state.sortCol].h + ' ' + arrow, 'showing ' + shown.length + ' of ' + n);
        summaryEl.textContent = parts.join(' · ');
      }
      ths.forEach(function (th) {
        var active = +th.getAttribute('data-col') === state.sortCol;
        th.setAttribute('aria-sort', active ? (state.sortAsc ? 'ascending' : 'descending') : 'none');
        th.querySelector('.sort-icon').textContent = active ? arrow : '';
      });
      if (pagEl) {
        pagEl.innerHTML = pages > 1 ? '<button type="button" data-step="-1"' + (state.page <= 1 ? ' disabled' : '') + '>Prev</button><span>Page ' +
          state.page + ' of ' + pages + '</span><button type="button" data-step="1"' + (state.page >= pages ? ' disabled' : '') + '>Next</button>' : '';
      }
    }

    if (pagEl) pagEl.addEventListener('click', function (e) {
      var step = e.target.getAttribute && e.target.getAttribute('data-step');
      if (step && !e.target.disabled) { state.page += +step; render(); }
    });
    if (searchEl) searchEl.addEventListener('input', function () { state.search = searchEl.value; state.page = 1; render(); });
    if (filterEl) filterEl.addEventListener('change', function () { state.filter = filterEl.value; state.page = 1; render(); });
    if (sizeEl) sizeEl.addEventListener('change', function () {
      state.size = sizeEl.value === 'All' ? 'All' : parseInt(sizeEl.value, 10); state.page = 1; render();
    });
    ths.forEach(function (th) {
      th.addEventListener('click', function () {
        var c = +th.getAttribute('data-col');
        state.sortAsc = state.sortCol === c ? !state.sortAsc : true;
        state.sortCol = c; state.page = 1; render();
      });
      th.addEventListener('keydown', function (e) {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); th.click(); }
      });
    });

    // Export: every filtered+sorted row (all pages). Text starting with = + - @ gets a leading apostrophe (formula injection).
    function exportRows(sep, ext) {
      var guard = function (s) { s = String(s == null ? '' : s); return /^[=+\-@]/.test(s) ? "'" + s : s; };
      var quote = sep === ','
        ? function (s) { return /[",\r\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; }
        : function (s) { return s.replace(/[\t\r\n]/g, ' '); };
      var out = function (v, i) { return v == null ? '' : COLS[i].f === 'text' ? quote(guard(v)) : String(v); };
      var lines = [COLS.map(function (c) { return quote(guard(c.h)); }).join(sep)].concat(view().map(function (r) { return r.map(out).join(sep); }));
      var blob = new Blob([lines.join('\r\n') + '\r\n'], { type: (sep === ',' ? 'text/csv' : 'text/tab-separated-values') + ';charset=utf-8' });
      var a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = cfg.file + '.' + ext;
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(function () { URL.revokeObjectURL(a.href); }, 1000);
    }
    if (csvEl) csvEl.addEventListener('click', function () { exportRows(',', 'csv'); });
    if (tsvEl) tsvEl.addEventListener('click', function () { exportRows('\t', 'tsv'); });

    render();
    return function reset() {
      state = Object.assign({}, DEFAULTS);
      if (searchEl) searchEl.value = '';
      if (filterEl) filterEl.value = 'All';
      if (sizeEl) sizeEl.value = '25';
      render();
    };
  }
  var resets = Array.prototype.slice.call(document.querySelectorAll('div.tbl[data-tid]')).map(initTable).filter(Boolean);

  // ---------------- Nav buttons ----------------
  document.getElementById('btn-reset').addEventListener('click', function () {
    resets.forEach(function (fn) { fn(); });
    sections.forEach(function (s) { setCollapsed(s, false); });
    if (hasPlotly) plotsIn(document).forEach(function (gd) {
      Plotly.restyle(gd, { visible: true });
      Plotly.relayout(gd, { 'xaxis.autorange': true, 'yaxis.autorange': true });
    });
    window.scrollTo(0, 0);
  });
  // Fixed plot width while printing so figures fit A4; restored afterwards.
  // The button waits for Plotly to finish redrawing; Ctrl+P relies on beforeprint alone.
  var printWidth = function () {
    return Promise.all(hasPlotly ? plotsIn(document).map(function (gd) { return Plotly.relayout(gd, { width: 680 }); }) : []);
  };
  document.getElementById('btn-pdf').addEventListener('click', function () {
    printWidth().then(function () { requestAnimationFrame(function () { window.print(); }); });
  });
  window.addEventListener('beforeprint', printWidth);
  window.addEventListener('afterprint', function () {
    if (hasPlotly) plotsIn(document).forEach(function (gd) {
      Plotly.relayout(gd, { width: null }).then(function () { Plotly.Plots.resize(gd); });
    });
  });
})();"""

_PLOT_CFG = {'displaylogo': False, 'responsive': True, 'toImageButtonOptions': {'format': 'svg'}}
_NUM_FMTS = ('fixed3', 'exp2', 'int')
_esc = lambda v: html.escape(str(v))


def _num(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def _text(v):
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    return str(v)


def table_spec(df, columns, sort_col, sort_asc=True, filter_col=None, badge_col=None, badge_map=None):
    names = [c for c, _, _ in columns]
    idx = lambda c: names.index(c) if c in names else None
    conv = [_num if f in _NUM_FMTS else _text for _, _, f in columns]
    rows = [] if df is None else [[fn(v) for fn, v in zip(conv, r)] for r in df.reindex(columns=names).itertuples(index=False, name=None)]
    return {'cols': [{'h': str(h), 'f': f} for _, h, f in columns], 'rows': rows, 'sort': idx(sort_col) or 0, 'asc': bool(sort_asc),
            'filter': idx(filter_col), 'badge': idx(badge_col), 'badgeMap': {str(k): str(v) for k, v in (badge_map or {}).items()}}


def _table_html(tid, spec, file_id):
    cfg = dict(spec, file=f'{file_id}_{tid}')
    cols, fc = cfg['cols'], cfg['filter']
    text = [c['h'] for i, c in enumerate(cols) if c['f'] == 'text' and i != fc]
    ths = ''.join(f'<th class="sortable" data-col="{i}" aria-sort="none" tabindex="0">{_esc(c["h"])}<span class="sort-icon"></span></th>' for i, c in enumerate(cols))
    search = (f'<input type="search" id="{tid}-search" placeholder="Search {_esc(", ".join(text))}..." aria-label="Search {_esc(", ".join(text))}">') if text else ''
    filt = ''
    if fc is not None:
        vals = {r[fc] for r in cfg['rows'] if r[fc] not in (None, '')}
        opts = ''.join(f'<option>{_esc(v)}</option>' for v in [k for k in cfg['badgeMap'] if k in vals] + sorted(vals - cfg['badgeMap'].keys()))
        filt = f'<select id="{tid}-filter" aria-label="{_esc(cols[fc]["h"])} filter"><option value="All">All</option>{opts}</select>'
    data = json.dumps(cfg, allow_nan=False).replace('</', r'<\/').replace('<!--', r'<\u0021--')
    return (f'<div class="tbl" data-tid="{tid}"><p class="print-only" id="{tid}-summary"></p><div class="tbl-controls">{search}{filt}'
            f'<select id="{tid}-size" aria-label="Rows per page"><option>25</option><option>50</option><option>100</option><option>All</option></select>'
            f'<button type="button" class="btn-sub" id="{tid}-csv">Export CSV</button><button type="button" class="btn-sub" id="{tid}-tsv">Export TSV</button>'
            f'<span class="tcount" id="{tid}-count"></span></div><div class="tw"><table class="dtbl"><thead><tr>{ths}</tr></thead><tbody id="{tid}-tbody"></tbody></table></div>'
            f'<div class="pagination" id="{tid}-pag"></div><script type="application/json" id="{tid}-data">{data}</script></div>')


def _logo_img(path, alt):
    if not path or not os.path.exists(path):
        return ''
    try:
        from PIL import Image
        with Image.open(path) as im:
            im.thumbnail((10 * 72, 72))
            buf = io.BytesIO()
            im.save(buf, format='PNG', optimize=True)
            raw = buf.getvalue()
    except Exception:
        with open(path, 'rb') as fh:
            raw = fh.read()
    return f'<img src="data:image/png;base64,{base64.b64encode(raw).decode()}" alt="{_esc(alt)}">'


def build_html_report(title, logo_path, info_rows, cards, sections, footer_html, file_id, doc_title=None,
                      stats_desc='Dataset, comparison and thresholds used in this analysis.'):
    """sections: [{id, title, desc, figs: [fig|None], table: table_spec(...)|None, table_id: optional}]; plain-text fields are escaped here, footer_html is trusted markup."""
    file_id = re.sub(r'\W', '_', str(file_id or 'report'))
    plot = lambda f: f'<div class="plot">{f.to_html(full_html=False, include_plotlyjs=False, config=_PLOT_CFG)}</div>'
    kv = ''.join(f'<tr><th>{_esc(k)}</th><td>{_esc(v)}</td></tr>' for k, v in info_rows)
    cds = ''.join(f'<div class="card {_esc(c)}"><span>{_esc(v)}</span>{_esc(l)}</div>' for l, v, c in cards)
    secs = [('general-stats', 'General Statistics', stats_desc, f'<table class="kv">{kv}</table><div class="cards">{cds}</div>')]
    for s in sections:
        plots = [plot(f) for f in s.get('figs') or () if f is not None]
        body = f'<div class="grid">{"".join(plots)}</div>' if len(plots) > 1 else ''.join(plots)
        if s.get('table') is not None:
            body += _table_html(re.sub(r'[^\w-]', '_', str(s.get('table_id') or s['id'])), s['table'], file_id)
        if body:
            secs.append((re.sub(r'[^\w-]', '_', str(s['id'])), s['title'], s.get('desc', ''), body))
    toc = ''.join(f'<a href="#{i}">{_esc(t)}</a>' for i, t, _, _ in secs)
    body = ''.join(f'<section id="{i}"><h2><button type="button" class="tog" aria-expanded="true" aria-controls="{i}-body"><span class="chevron" aria-hidden="true">&#9660;</span>{_esc(t)}</button></h2>'
                   f'<div class="sbody" id="{i}-body"><p class="desc">{_esc(d)}</p>{c}</div></section>' for i, t, d, c in secs)
    btns = '<div class="nav-actions"><button type="button" class="nav-btn" id="btn-reset">&#8634; Reset</button><button type="button" class="nav-btn" id="btn-pdf">&#128424; Save as PDF</button></div>'
    return (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{_esc(doc_title or title)}</title><style>{REPORT_CSS}</style><script type="text/javascript">{get_plotlyjs()}</script></head>'
            f'<body><nav>{_logo_img(logo_path, str(title).split(" ")[0])}<h1>{_esc(title)}</h1><small>Generated {datetime.now().strftime("%Y-%m-%d %H:%M")}</small>{btns}</nav>'
            f'<div class="wrap"><aside>{toc}</aside><main>{body}<footer>{footer_html}</footer></main></div>'
            f'<script>{REPORT_JS}</script></body></html>')


def network_figure(nodes, edges, title=None, seed=42):
    """Static plotly rendering of streamlit_agraph Node/Edge objects (or dicts); labels and titles are treated as untrusted plain text."""
    import networkx as nx
    import plotly.graph_objects as go

    def attr(o, *keys):
        for k in keys:
            v = o.get(k) if isinstance(o, dict) else getattr(o, k, None)
            if v is not None:
                return v
        return None

    def plain(v, strip_tags=False):
        s = re.sub(r'<[^>]*>', '', str(v)) if strip_tags else str(v)
        return re.sub(r'\r?\n', '<br>', html.escape(s, quote=False))

    nodes = [n for n in nodes if attr(n, 'id') is not None]
    ids = [attr(n, 'id') for n in nodes]
    links = [(attr(e, 'source', 'from'), attr(e, 'to', 'target'), e) for e in edges]
    G = nx.Graph()
    G.add_nodes_from(ids)
    G.add_edges_from((s, t) for s, t, _ in links if s in G and t in G)
    pos = nx.spring_layout(G, seed=seed) if len(G) else {}
    fig = go.Figure()
    lines, arrows = {}, {}
    for s, t, e in links:
        if s in pos and t in pos:
            c = attr(e, 'color')
            c = c if isinstance(c, str) else '#999'
            xs, ys = lines.setdefault((c, bool(attr(e, 'dashed'))), ([], []))
            xs += [pos[s][0], pos[t][0], None]
            ys += [pos[s][1], pos[t][1], None]
            if attr(e, 'directed'):
                ax, ay, sz = arrows.setdefault(c, ([], [], []))
                for f, size in ((0.6, 0), (0.7, 11)):
                    ax.append(pos[s][0] + f * (pos[t][0] - pos[s][0]))
                    ay.append(pos[s][1] + f * (pos[t][1] - pos[s][1]))
                    sz.append(size)
    for (c, dashed), (xs, ys) in lines.items():
        fig.add_trace(go.Scatter(x=xs, y=ys, mode='lines', line=dict(color=c, width=1, dash='dot' if dashed else 'solid'), hoverinfo='skip'))
    # Arrowheads: each visible marker points away from the invisible one placed just before it on the same edge.
    for c, (ax, ay, sz) in arrows.items():
        fig.add_trace(go.Scatter(x=ax, y=ay, mode='markers', hoverinfo='skip',
                                 marker=dict(symbol='arrow', angleref='previous', size=sz, color=c, line=dict(width=0))))
    symbols = {'square': 'square', 'hexagon': 'hexagon', 'diamond': 'diamond', 'star': 'star', 'triangle': 'triangle-up', 'triangleDown': 'triangle-down'}
    groups = {}
    for n, i in zip(nodes, ids):
        c = attr(n, 'color')
        g = groups.setdefault(c if isinstance(c, str) else '#97C2FC', {k: [] for k in ('x', 'y', 'text', 'hovertext', 'size', 'symbol')})
        g['x'].append(pos[i][0])
        g['y'].append(pos[i][1])
        g['text'].append(plain(attr(n, 'label', 'id')))
        g['hovertext'].append(plain(attr(n, 'title', 'label', 'id'), strip_tags=True))
        g['size'].append(min(max((_num(attr(n, 'size')) or 25) * 0.5, 4), 40))
        g['symbol'].append(symbols.get(attr(n, 'shape'), 'circle'))
    for c, g in groups.items():
        fig.add_trace(go.Scatter(x=g['x'], y=g['y'], mode='markers+text', text=g['text'], textposition='top center', textfont=dict(size=9), cliponaxis=False,
                                 hovertext=g['hovertext'], hoverinfo='text',
                                 marker=dict(color=c, size=g['size'], symbol=g['symbol'], line=dict(color='#fff', width=0.5))))
    fig.update_layout(template='simple_white', height=700, showlegend=False, hovermode='closest', margin=dict(l=60, r=60, t=60 if title else 30, b=20),
                      title=html.escape(str(title), quote=False) if title else None)
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return fig
