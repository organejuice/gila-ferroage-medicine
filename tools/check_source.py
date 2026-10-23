"""Offline publication-integrity check; no biological execution or dependencies."""
from pathlib import Path
import ast
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    failures = []
    expected = {}
    for line in (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        if line:
            value, name = line.split('  ', 1)
            if Path(name).is_absolute() or '..' in Path(name).parts:
                failures.append('Unsafe manifest path')
                continue
            expected[name] = value
            p = ROOT / name
            if not p.is_file() or digest(p) != value:
                failures.append('Missing or modified publication file: ' + name)
    changes = json.loads((ROOT / 'PUBLICATION_CHANGES.json').read_text(encoding='utf-8'))
    allowed = {r['path']: r for r in changes['modified_source_files']}
    source = json.loads((ROOT / 'code_manifest.json').read_text(encoding='utf-8'))
    for record in source:
        name = 'project/' + record['path']
        p = ROOT / name
        wanted = record['sha256']
        if name in allowed:
            if wanted != allowed[name]['original_sha256']:
                failures.append('Privacy source-baseline mismatch: ' + name)
            wanted = allowed[name]['published_sha256']
        if not p.is_file() or digest(p) != wanted:
            failures.append('Frozen-source mismatch: ' + name)
    parsed = 0
    for name in expected:
        if name.endswith('.py'):
            try:
                ast.parse((ROOT / name).read_text(encoding='utf-8-sig'), filename=name)
                parsed += 1
            except (SyntaxError, UnicodeError) as exc:
                failures.append('Python syntax/encoding error: ' + str(exc))
    result = dict(status='FAIL' if failures else 'PASS',
                  publication_files_checked=len(expected),
                  original_source_records_checked=len(source),
                  documented_privacy_modified_files=len(allowed),
                  python_files_parsed=parsed,
                  biological_analyses_executed=False,
                  failures=failures)
    print(json.dumps(result, indent=2))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
