#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Rscript tools/setup-r-packages.R
Rscript tools/check-project.R
