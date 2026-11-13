"""Locate installed executables even when macOS GUI/virtualenv PATH omits them."""
import os
import shutil
import subprocess
from pathlib import Path

def find_tool(name):
    if name not in ('datasets', 'minimap2', 'mafft', 'iqtree2', 'iqtree'):
        raise ValueError('Unsupported executable')
    candidates = [shutil.which(name)]
    home = Path.home()
    candidates += [
        str(home / 'bin' / name),
        str(Path('/opt/homebrew/bin') / name),
        str(Path('/usr/local/bin') / name),
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return str(Path(candidate).resolve())
    return None

def check_tools():
    results={}
    for name in ('datasets','minimap2'):
        exe=find_tool(name)
        if not exe:
            results[name]={'ok':False, 'path':None,
                           'detail':'Исполняемый файл не найден. Проверьте ~/bin и /opt/homebrew/bin.'}
            continue
        try:
            p=subprocess.run([exe, '--version'],
                             capture_output=True, text=True, timeout=20)
            lines=(p.stdout+'\n'+p.stderr).strip()
            results[name]={'ok':p.returncode==0,'path':exe,
                           'detail':lines[-400:] if lines else ('Программа вернула код '+str(p.returncode))}
        except (OSError,subprocess.TimeoutExpired) as exc:
            results[name]={'ok':False,'path':exe,'detail':str(exc)}
    return results
