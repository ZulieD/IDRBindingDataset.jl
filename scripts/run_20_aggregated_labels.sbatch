#!/bin/bash
#SBATCH --job-name=20_aggregated_labels
#SBATCH --nodelist=node48
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=20_aggregated_labels.o%j.out

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version
# Remove `--plot_dir ...` to skip the PNG plots.
"$JULIA" scripts/20_aggregated_labels.jl \
    --max_resolution 5.0 \
    --pdb_vote fraction \
    --plot_dir data/scratch/plots/aggregated_labels
