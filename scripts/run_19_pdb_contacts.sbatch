#!/bin/bash
#SBATCH --job-name=19_pdb_contacts
#SBATCH --nodelist=node48
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=19_pdb_contacts.o%j.out

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version
# Add `--download_missing` if the node has internet access (RCSB / EBI), and
# `--pdb_quality_path` to copy METHOD/RESOLUTION into the output.
"$JULIA" --threads="${SLURM_CPUS_PER_TASK}" scripts/19_pdb_contacts.jl \
    --mmcif_dir /alpha/database/pdb/mmcif_files \
    --sifts_xml_dir /alpha/database/sift_new-2025-10-27/sift_files \
    --pdb_quality_path /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/AlphaConformerBenchmark/scripts/csv/pdb_information_details_V2_10_2026.csv \
    --checkpoint_dir data/scratch/checkpoints/19_pdb_contacts \
    --batch_size 50
