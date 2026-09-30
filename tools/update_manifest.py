#!/usr/bin/env python3
"""Refresh the portable bundle's hash inventory, without accepting changed results."""

import hashlib
from pathlib import Path

from repository_files import retained_files

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest = ROOT / 'MANIFEST.sha256'
    entries = [f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(ROOT)}'
               for path in retained_files(ROOT) if path != manifest]
    manifest.write_text('\n'.join(entries) + '\n')
    print(f'Recorded {len(entries)} portable calculation files.')


if __name__ == '__main__':
    main()
