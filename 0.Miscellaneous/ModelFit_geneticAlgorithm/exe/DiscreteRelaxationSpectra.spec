# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['..\\DiscreteRelaxationSpectra.py'],
             pathex=['C:\\Users\\oskat\\OneDrive - Instituto Tecnologico y de Estudios Superiores de Monterrey\\Documents\\MNT_ITESM_courses\\0.Miscellaneous\\ModelFit_geneticAlgorithm\\exe'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='DiscreteRelaxationSpectra',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True , icon='..\\icon.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='DiscreteRelaxationSpectra')
