# -*- mode: python ; coding: utf-8 -*-
# PyInstaller build of the web interface:  pyinstaller LumosGUI.spec
# The executable starts the server and opens the browser; calculations run the same
# executable with --run-lumos (see LumosGUI.main). R/limma is not bundled: install R and
# the Bioconductor limma package on the machine to use it.
from PyInstaller.utils.hooks import collect_all

datas = [('templates', 'templates'), ('static', 'static'), ('r', 'r'), ('DOCUMENTATION.md', '.')]
binaries = []
hiddenimports = ['lumos', 'lumos_core', 'lumos_report']
for package in ('pydeseq2', 'anndata', 'networkx', 'plotly', 'statsmodels'):
    d, b, h = collect_all(package)
    datas += d
    binaries += b
    hiddenimports += h

a = Analysis(
    ['LumosGUI.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter', 'streamlit'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='LumosNetworks',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
