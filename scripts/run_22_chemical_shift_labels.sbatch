#!/bin/bash
#SBATCH --job-name=22_chemical_shift_labels
#SBATCH --nodelist=node48
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=22_chemical_shift_labels.o%j.out

# Entries are read from the BMRB cache filled by script 21 (entries missing
# from the cache are downloaded, which needs internet access). Outputs go to
# data/derived/chemical_shift_labels; settings can be changed with options,
# e.g. --require pH,temperature,buffer (see --help).

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version
"$JULIA" scripts/22_chemical_shift_labels.jl
