#!/usr/bin/env bash

# This script runs strainphlan to reconstruct a phylogeny of Lactobacillus crispatus
# using consensus markers from metagenomes and isolates together with selected
# reference genomes from GTDB and the Isala cohort. Reference genomes are selected
# based on seeds.txt files generated after sampling with SCARAP and accessed directly from
# their respective directories without copying.

set -euo pipefail

# Set directory
dir_project="../../"

threads=16
species_tag="t__SGB7045"

# Inputs from metaphlan step
dir_phylo="$dir_project/results/strain_level_analyses/isala/phylogeny"

consensus_markers="$dir_phylo/metaphlan/consensus_markers"
targeted_markers="$dir_phylo/metaphlan/targeted_markers"

# Metaphlan database
metaphlan_db="$dir_project/data/metaphlan_db/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl"

# Reference genomes (organized structure)
dir_ref="$dir_project/data/phylogeny_references"

dir_gtdb="$dir_ref/gtdb_crispatus/genomes"
dir_isala="$dir_ref/isala_crispatus/genomes"

# Representative genome selection
gtdb_seeds="$dir_project/results/public_genomes/gtdb/representatives/seeds.txt"
isala_seeds="$dir_project/results/public_genomes/isala/representatives/seeds.txt"

# Define output directory
dout_strainphlan="$dir_phylo/strainphlan"
mkdir -p "$dout_strainphlan"

# Build reference genome list dynamically
ref_list="$dout_strainphlan/reference_list.txt"
> "$ref_list"

echo "collecting GTDB reference genomes"

while read -r genome; do
    [[ -z "$genome" || "$genome" =~ ^# ]] && continue

    for f in "$dir_gtdb/${genome}"*.fna; do
        [ -f "$f" ] || continue
        echo "$f" >> "$ref_list"
    done
done < "$gtdb_seeds"

echo "collecting Isala reference genomes"

while read -r genome; do
    [[ -z "$genome" || "$genome" =~ ^# ]] && continue

    for f in "$dir_isala/${genome}"*.fna; do
        [ -f "$f" ] || continue
        echo "$f" >> "$ref_list"
    done
done < "$isala_seeds"

# Ensure that at least one reference genome was found
[ -s "$ref_list" ] || { echo "no reference genomes found"; exit 1; }

echo "running strainphlan"

# Run StrainPhlAn

strainphlan \
    -s "$consensus_markers"/*.json.bz2 \
    -d "$metaphlan_db" \
    -m "$targeted_markers/${species_tag}.fna" \
    -r $(cat "$ref_list") \
    -o "$dout_strainphlan" \
    -n "$threads" \
    -c "$species_tag" \
    --marker_in_n_samples 50 \
    --sample_with_n_markers 100 \
    --phylophlan_mode accurate

echo "finished strainphlan"
