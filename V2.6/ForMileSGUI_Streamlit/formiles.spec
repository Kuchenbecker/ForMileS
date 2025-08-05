# -*- mode: python ; coding: utf-8 -*-
import os
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

block_cipher = None

# Coleta módulos internos do RDKit
hiddenimports = collect_submodules('rdkit')

# Coleta arquivos de dados do RDKit (se necessário)
datas = collect_data_files('rdkit')

# Arquivos do seu projeto a empacotar
project_data = [
    ('config.json', '.'),
    ('parameters.json', '.'),
    ('main.py', '.'),
    ('app.py', '.')
]

datas += project_data

a = Analysis(
    ['launch.py'],               # Script principal
    pathex=[os.getcwd()],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# Aqui removemos o console (janela silenciosa)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=False,
    name='ForMileS',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False  # <-- silencioso!
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    name='ForMileS'
)
