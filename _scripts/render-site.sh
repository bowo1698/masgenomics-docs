#!/bin/bash
# Full-project Quarto render for masgenomics-docs.
#
# Rebuilds _site/ and refreshes the _freeze/ tutorial caches. The task pages
# execute masbayes / masreml live, so those packages (plus CMplot) must be
# installed in the active R library — this is why the render is run here / on
# HPC and never in CI (CI renders from the committed _freeze/ cache only).
#
# Usage:
#   bash _scripts/render-site.sh            # render, reusing valid freeze caches
#   bash _scripts/render-site.sh --clean    # rm _site + _freeze first (re-exec all)
#   bash _scripts/render-site.sh <path.qmd> # render a single file (always executes)
set -euo pipefail

# Repo root = parent of this script's directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
cd "${REPO_ROOT}"

CLEAN=0
TARGET=""
case "${1:-}" in
  --clean) CLEAN=1 ;;
  "")      ;;
  *)       TARGET="$1" ;;
esac

# R packages the live-eval task tutorials need. Fail early with a clear list
# rather than dying mid-render on a missing library() call.
REQ_PKGS="masbayes masreml CMplot ggplot2 knitr rmarkdown data.table"
missing=$(Rscript -e "p<-strsplit('${REQ_PKGS}',' ')[[1]];cat(p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)])")
if [[ -n "${missing}" ]]; then
  echo "ERROR: missing R packages: ${missing}" >&2
  echo "Install before rendering (masbayes/masreml need the Rust toolchain)." >&2
  exit 1
fi

# Keep the Rust MCMC kernel from over-subscribing the allocated cores.
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-4}}"

if [[ "${CLEAN}" == "1" ]]; then
  echo "Cleaning _site/ and _freeze/ ..."
  rm -rf _site _freeze
fi

echo "host=$(hostname) | start=$(date '+%F %T') | quarto=$(quarto --version) | rayon=${RAYON_NUM_THREADS}"
if [[ -n "${TARGET}" ]]; then
  time quarto render "${TARGET}"
else
  time quarto render
fi
echo "end=$(date '+%F %T') | output -> ${REPO_ROOT}/_site"
