#!/usr/bin/env python3
"""Move superseded result artifacts into a timestamped, recoverable archive."""
from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from datetime import datetime, timezone
from pathlib import Path


TARGETS = (
    Path("results/tables"),
    Path("results/figures"),
    Path("results/2026-09-19-disease-side-deterministic"),
    Path("results/2026-09-20-drug-side-execution"),
    Path("docs/reports/2026-09-20_disease-side-completion.md"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def archive(root: Path, archive_name: str) -> Path:
    archive_root = root / "results" / "archive" / archive_name
    if archive_root.exists():
        raise FileExistsError(f"Archive already exists: {archive_root}")
    archive_root.mkdir(parents=True)
    moved = []
    for relative in TARGETS:
        source = root / relative
        if not source.exists():
            continue
        destination = archive_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(destination))
        moved.append((relative, destination))

    manifest = archive_root / "ARCHIVE_MANIFEST.tsv"
    rows = []
    for relative, destination in moved:
        files = [destination] if destination.is_file() else list(destination.rglob("*"))
        for path in sorted(p for p in files if p.is_file()):
            rows.append({
                "original_path": str(relative / path.relative_to(destination)) if destination.is_dir() else str(relative),
                "archived_path": str(path.relative_to(root)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["original_path", "archived_path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    readme = archive_root / "README.md"
    readme.write_text(
        "# Archived pre-v3 disease-side results\n\n"
        f"Archived on {datetime.now(timezone.utc).isoformat()} before execution of the frozen disease-side v3.0 protocol.\n\n"
        "The archive is recoverable; no raw data, source code, frozen protocol, or environment files were removed. "
        "See `ARCHIVE_MANIFEST.tsv` for file-level checksums.\n",
        encoding="utf-8",
    )
    return archive_root


def refresh_manifest(root: Path, archive_name: str) -> Path:
    archive_root = root / "results" / "archive" / archive_name
    if not archive_root.exists():
        raise FileNotFoundError(f"Archive does not exist: {archive_root}")
    manifest = archive_root / "ARCHIVE_MANIFEST.tsv"
    rows = []
    for path in sorted(p for p in archive_root.rglob("*") if p.is_file() and p.name not in {"ARCHIVE_MANIFEST.tsv", "README.md"}):
        rows.append({
            "original_path": "see archived path; includes recoverable failed attempts",
            "archived_path": str(path.relative_to(root)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["original_path", "archived_path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--archive-name", required=True)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    path = refresh_manifest(args.root.resolve(), args.archive_name) if args.refresh else archive(args.root.resolve(), args.archive_name)
    print(path)


if __name__ == "__main__":
    main()
