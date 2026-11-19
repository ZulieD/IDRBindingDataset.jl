#!/bin/bash
#SBATCH --job-name=21_disprot_bmrb_mapping
#SBATCH --nodelist=node48
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=21_disprot_bmrb_mapping.o%j.out

# The node needs internet access to https://api.bmrb.io (otherwise run the
# same command on the login node). Downloads are cached in
# data/scratch/bmrb_cache: if the job stops, resubmit it and only the missing
# entries are downloaded.

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version
"$JULIA" scripts/21_disprot_bmrb_mapping.jl
