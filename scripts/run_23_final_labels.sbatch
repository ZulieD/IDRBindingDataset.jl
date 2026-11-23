#!/bin/bash
#SBATCH --job-name=23_final_labels
#SBATCH --nodelist=node48
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=23_final_labels.o%j.out

# Reads the outputs of scripts 18, 20 and 22 and writes
# data/derived/final_labels (final_labels.csv.gz, final_labels_sources.csv.gz).

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version
"$JULIA" scripts/23_final_labels.jl
