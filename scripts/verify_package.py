#!/usr/bin/env python3
"""Verify distributed source and synthetic input bytes, not research inputs."""
import hashlib
from pathlib import Path, PurePosixPath
root=Path(__file__).resolve().parents[1]
count=0
for line in (root/'MANIFEST.sha256').read_text().splitlines():
 digest, name=line.split('  ',1)
 rel=PurePosixPath(name)
 if rel.is_absolute() or '..' in rel.parts:
  raise ValueError('Invalid manifest path')
 path=root.joinpath(*rel.parts)
 if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
  raise ValueError(f'Missing or changed package file: {name}')
 count+=1
print(f'[PASS] Package SHA256: {count} files')
