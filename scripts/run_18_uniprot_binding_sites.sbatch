#!/bin/bash
#SBATCH --job-name=18_uniprot_binding_sites
#SBATCH --nodelist=node48
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=18_uniprot_binding_sites.o%j.out

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version
# The node needs internet access to https://rest.uniprot.org
"$JULIA" scripts/18_uniprot_binding_sites.jl
