#!/usr/bin/env python3
"""Static integrity only: no attempt to interpret/execute MATLAB."""
from pathlib import Path
import hashlib,re,csv
root=Path(__file__).resolve().parents[1]
core=list((root/'+pof').glob('*.m'))
assert len(core)>=10 and all(p.stat().st_size>100 for p in core)
files=list((root/'+pof').rglob('*.m')) + list((root/'io').glob('*.m')) + list((root/'tests').glob('*.m')) + list((root/'config').glob('*.m')) + [root/'run_pof.m']
for p in files:
    s=p.read_text(encoding='utf-8-sig')
    assert s.count('function')>=1,(p,'missing MATLAB function declaration')
    assert not re.search(r'\bE:\\DAT\\|\bG:\\0604\\',s),(p,'hard-coded original drive')
    assert 'http://' not in s and 'https://' not in s,(p,'network reference in executable')
assert (root/'LICENSE_PENDING.md').exists()
with (root/'docs/FIGURE_PROVENANCE.csv').open(encoding='utf-8') as f:
    rows=list(csv.DictReader(f))
assert len(rows)==24
expected={}
for line in (root/'SHA256SUMS.txt').read_text().splitlines():
    digest,name=line.split('  ',1);expected[name]=digest
actual={}
for p in sorted(root.rglob('*')):
    if p.is_file() and p.name!='SHA256SUMS.txt' and '.git' not in p.relative_to(root).parts:
        actual[p.relative_to(root).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
assert actual==expected,'SHA256 source manifest mismatch (missing/modified/extra source).'
print('PASS: source structure, no old drive paths, figure mapping, SHA256.')
print('MATLAB execution and exact-paper reproduction: NOT TESTED.')
