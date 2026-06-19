#!/usr/bin/env bash
set -euo pipefail

# Native Ubuntu/Debian setup. This installs the system R provided by apt,
# then installs the repository R packages from the CRAN snapshot configured in .Rprofile.
# For a pinned R interpreter version, use Dockerfile instead.

if ! command -v sudo >/dev/null 2>&1; then
  echo "sudo is required for native installation. Use Docker if sudo is unavailable." >&2
  exit 1
fi

sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  r-base \
  r-base-dev \
  build-essential \
  ca-certificates \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev

cd "$(dirname "$0")/.."
Rscript tools/setup-r-packages.R
Rscript tools/check-project.R
