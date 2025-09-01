# gui.spec - PyInstaller spec file for ForMileS GUI
# Salve este arquivo na mesma pasta do gui.py, main.py e parameters.json

# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

from PyInstaller.utils.hooks import collect_submodules

a = Analysis(
    ['gui.py'],
    pathex=[],
    binaries=[],
    datas=[
        ('main.py', '.'),             # Inclui main.py na raiz do bundle
        ('parameters.json', '.'),     # Inclui parameters.json na raiz
    ],
    hiddenimports=[],
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='ForMileS',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False  # Define como GUI, sem console
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='ForMileS'
)
