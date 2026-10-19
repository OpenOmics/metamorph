#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
ABOUT:
    Removes or compresses "cruft" left behind in a metamorph output
    directory: intermediate files that are no longer needed, or
    uncompressed files that are worth shrinking once the pipeline is
    done with them. Driven entirely by a JSON manifest (config/cleanup.json)
    so the list of files to clean up can be extended without touching
    any workflow code -- see that file for the manifest format.

USAGE:
    cleanup_cruft.py --manifest config/cleanup.json --project-dir <path> [--dry-run]
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys


def log(msg):
    print(f"[cleanup_cruft] {msg}", flush=True)


def resolve(project_dir, pattern):
    """Expands a manifest glob pattern (relative to project_dir) to
    concrete, existing paths on disk. A pattern matching nothing is not
    an error -- it just means this entry doesn't apply to this run/mode.
    """
    return sorted(glob.glob(os.path.join(project_dir, pattern), recursive=True))


def do_delete(project_dir, entries, dry_run):
    for entry in entries:
        pattern = entry["pattern"]
        matches = resolve(project_dir, pattern)
        if not matches:
            log(f"delete: no files matched `{pattern}`, skipping")
            continue
        for path in matches:
            log(f"delete: removing {path}")
            if dry_run:
                continue
            if os.path.isdir(path) and not os.path.islink(path):
                shutil.rmtree(path)
            else:
                os.remove(path)


def do_compress(project_dir, entries, dry_run):
    for entry in entries:
        pattern = entry["pattern"]
        matches = resolve(project_dir, pattern)
        if not matches:
            log(f"compress: no files matched `{pattern}`, skipping")
            continue
        for path in matches:
            if path.endswith(".bz2"):
                continue
            log(f"compress: bzip2 -9 {path}")
            if dry_run:
                continue
            subprocess.run(["bzip2", "-f", "-9", path], check=True)


def do_compress_indexed(project_dir, entries, dry_run):
    for entry in entries:
        pattern = entry["pattern"]
        preset = entry.get("tabix_preset")
        extra_args = entry.get("tabix_args", [])
        matches = resolve(project_dir, pattern)
        if not matches:
            log(f"compress_indexed: no files matched `{pattern}`, skipping")
            continue
        for path in matches:
            if path.endswith(".gz"):
                continue
            log(f"compress_indexed: bgzip -l 9 {path}")
            if not dry_run:
                subprocess.run(["bgzip", "-f", "-l", "9", path], check=True)

            gz_path = f"{path}.gz"
            tabix_cmd = ["tabix"]
            if preset:
                tabix_cmd += ["-p", preset]
            tabix_cmd += extra_args
            tabix_cmd += [gz_path]
            log(f"compress_indexed: {' '.join(tabix_cmd)}")
            if dry_run:
                continue
            subprocess.run(tabix_cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Path to cleanup.json")
    parser.add_argument(
        "--project-dir", required=True,
        help="Project output directory; manifest patterns are relative to this",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Log what would happen without deleting/compressing/indexing anything",
    )
    args = parser.parse_args()

    with open(args.manifest) as fh:
        manifest = json.load(fh)

    do_delete(args.project_dir, manifest.get("delete", []), args.dry_run)
    do_compress(args.project_dir, manifest.get("compress", []), args.dry_run)
    do_compress_indexed(args.project_dir, manifest.get("compress_indexed", []), args.dry_run)

    log("done")


if __name__ == "__main__":
    main()
