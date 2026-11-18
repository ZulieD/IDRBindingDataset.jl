#!/bin/bash
#SBATCH --job-name=17_disprot_pdb_mapping
#SBATCH --nodelist=node48
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --output=17_disprot_pdb_mapping.o%j.out

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version
"$JULIA" --threads="${SLURM_CPUS_PER_TASK}" scripts/17_disprot_pdb_mapping.jl \
    --sifts_xml_dir /alpha/database/sift_new-2025-10-27/sift_files \
    --checkpoint_dir data/scratch/checkpoints/17_disprot_pdb_mapping \
    --batch_size 50
