#!/usr/bin/env bash

# This script dereplicates whole-genome sequenced isolates as a preprocessing
# step before running further strain-level analyses using SCARAP v1.0.1.

# Assumptions:
# - Gene prediction has already been performed (e.g. using Prodigal).
# - Isolate genomes follow the documented directory structure:
#   data/isolate_genomes/<participant_ID>/genes/
# - Each genes/ directory contains predicted gene files:
#     *.ffn  (nucleotide sequences)
#     *.faa  (amino acid sequences)

# Load directories

dir_isolates="../../data/isolate_genomes"
dout="../../results/isala/drep"

threads=5

mkdir -p "$dout"

# Summary of representative isolates
rep_summary="${dout}/representatives_summary.tsv"
echo -e "participant_id\tisolate_id" > "$rep_summary"

for participant_dir in "$dir_isolates"/participant_*; do
    [ -d "$participant_dir" ] || continue

    participant_id="$(basename "$participant_dir")"
    genes_dir="$participant_dir/genes"
    raw_reads_dir="$participant_dir/raw_reads"
    assembly_dir="$participant_dir/assembly"

    [ -d "$genes_dir" ] || {
        echo "Skipping $participant_id (no genes dir)"
        continue
    }

    echo "Processing $participant_id"

    # Output structure per participant
    out_participant="$dout/$participant_id"
    core_dir="$out_participant/core"
    rep_dir="$out_participant/representatives"
    assemblies_out="$out_participant/assemblies"
    reads_out="$out_participant/raw_reads"

    mkdir -p "$out_participant" "$assemblies_out" "$reads_out"

    # Step 1: identify core genes
    scarap core \
        "$genes_dir" \
        "$core_dir" \
        -t "$threads"

    core_file="$core_dir/genes.tsv"

    [ -f "$core_file" ] || {
        echo "Core gene step failed for $participant_id"
        continue
    }

    # Step 2: dereplicate isolates
    scarap sample \
        "$genes_dir" \
        "$core_file" \
        "$rep_dir" \
        -i 0.9999 \
        -t "$threads" \
        --method mean90

    seeds_file="$rep_dir/seeds.txt"

    [ -f "$seeds_file" ] || {
        echo "No representatives found for $participant_id"
        continue
    }

    # Step 3: collect representative data
    while IFS= read -r isolate_id; do

        # Record selected isolate
        echo -e "${participant_id}\t${isolate_id}" >> "$rep_summary"

        # Copy raw reads
        if [ -d "$raw_reads_dir" ]; then
            cp "$raw_reads_dir/${isolate_id}"*_R1*.fastq.gz "$reads_out/" 2>/dev/null || true
            cp "$raw_reads_dir/${isolate_id}"*_R2*.fastq.gz "$reads_out/" 2>/dev/null || true
        fi

        # Copy assembly
        if [ -d "$assembly_dir" ]; then
            find "$assembly_dir" -name "${isolate_id}*.fna.gz" -exec cp {} "$assemblies_out/" \;
        fi

    done < "$seeds_file"

    # Step 4: unzip assemblies
    for f in "$assemblies_out"/*.fna.gz; do
        [ -f "$f" ] || continue
        gunzip -k "$f"
    done

done

echo "Finished 01_sample_genomes"
