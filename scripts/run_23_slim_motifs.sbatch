#!/bin/bash
#SBATCH --job-name=23_slim_motifs
#SBATCH --nodelist=node48
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=23_slim_motifs.o%j.out

# Reads the ELM files downloaded by script 11 (data/scratch/elm) and writes
# data/scratch/elm/slim_labels/slim_labels.csv.gz (ELM-derived, not versioned).

# Julia 1.12 installed by juliaup (`juliaup add 1.12`).
JULIA=/home/julie.daniel/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia

# Run from the repository root so that the script finds scripts/Project.toml.
cd /stockage/EQUIPES/AMIG/MEMBERS/julie.daniel/IDRBindingDataset.jl || exit 1

"$JULIA" --version

# ELM files from script 11. If one is missing or empty, run script 11 first
# (needs internet access from the compute node).
INSTANCES=data/scratch/elm/instances.tsv
ELMS_INDEX=data/scratch/elm/elms_index.tsv
if [ ! -s "$INSTANCES" ] || [ ! -s "$ELMS_INDEX" ]; then
    echo "ELM files missing: running scripts/11_download_elm.jl"
    "$JULIA" scripts/11_download_elm.jl --timeout_seconds 600 || {
        echo "ERROR: ELM download failed (no internet on this node?)." >&2
        echo "Run scripts/11_download_elm.jl on the login node, then resubmit." >&2
        exit 1
    }
fi
for file in "$INSTANCES" "$ELMS_INDEX"; do
    if [ ! -s "$file" ]; then
        echo "ERROR: $file is still missing or empty." >&2
        exit 1
    fi
done
echo "ELM files: $(wc -l < "$INSTANCES") lines in $INSTANCES, $(wc -l < "$ELMS_INDEX") lines in $ELMS_INDEX"

"$JULIA" scripts/23_slim_motifs.jl
