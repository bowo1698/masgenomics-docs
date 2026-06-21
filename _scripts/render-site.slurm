#!/bin/bash
#SBATCH --job-name=masgenomics_docs_render
#SBATCH --partition=general
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --account=a_zenger
#SBATCH -o logs/masgenomics_docs_render_%j.out
#SBATCH -e logs/masgenomics_docs_render_%j.err

# Full Quarto render of masgenomics-docs on HPC. Rebuilds _site/ and the
# _freeze/ task caches (the tutorial pages execute masbayes/masreml live).
# After the job: rsync the refreshed _freeze/ back to the laptop repo and
# commit it so CI renders from cache without needing the packages installed.
#
# Submit from the repo root so the relative logs/ path below resolves:
#   sbatch _scripts/render-site.slurm
# Override the checkout location with: sbatch --export=REPO=/path ...

module load r/4.4.0-gfbf-2023a
module load quarto              # adjust to the cluster's quarto module name

REPO=${REPO:-/scratch/user/aguswibowo/masgenomics-docs}
mkdir -p logs

echo "host=$(hostname) | start=$(date) | REPO=${REPO}"
cd "${REPO}"
bash _scripts/render-site.sh --clean
echo "end=$(date)"
