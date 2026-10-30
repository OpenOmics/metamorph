#!/usr/bin/env bash
# Runs metamorph's unit test suite (python's built-in unittest, test/ directory).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

python3 -W ignore::ResourceWarning -m unittest discover -s test -p 'test_*.py' -v
