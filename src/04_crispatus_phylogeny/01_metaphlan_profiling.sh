#!/usr/bin/env bash

# This script runs MetaPhlAn on metagenomic samples and dereplicated isolate reads
# from the Isala cohort to obtain taxonomic profiles and consensus marker sequences.
# It also extracts species-specific marker genes from the MetaPhlAn database for the
# target species (Lactobacillus crispatus). These markers will be combined with sample
# markers for downstream phylogenetic reconstruction using StrainPhlAn.

set -euo pipefail

# Set project structure and load input data
dir_project="../../"
dir_metagenomes="$dir_project/data/metagenomic_data/isala"
dir_isolates="$dir_project/results/isala/drep"

threads=16

# Output directories for phylogeny workflow
dout_phylo="$dir_project/results/strain_level_analyses/isala/phylogeny"
dout_metaphlan="$dout_phylo/metaphlan"

dout_sams="$dout_metaphlan/sams"
dout_bt2="$dout_metaphlan/bowtie2"
dout_profiles="$dout_metaphlan/profiles"
dout_markers="$dout_metaphlan/consensus_markers"
dout_targeted="$dout_metaphlan/targeted_markers"

mkdir -p "$dout_sams" "$dout_bt2" "$dout_profiles" "$dout_markers" "$dout_targeted"

# Run MetaPhlan on metagenomic samples
for participant_dir in "$dir_metagenomes"/*; do
    [ -d "$participant_dir" ] || continue

    participant_id=$(basename "$participant_dir")
    reads_dir="$participant_dir/raw_reads"

    [ -d "$reads_dir" ] || continue

    for r1 in "$reads_dir"/*_R1.fastq.gz; do
        [ -f "$r1" ] || continue

        sample=$(basename "${r1%_R1.fastq.gz}")
        r2="$reads_dir/${sample}_R2.fastq.gz"

        echo "metagenome: $participant_id - $sample"

        metaphlan \
            "$r1,$r2" \
            --input_type fastq \
            -s "$dout_sams/${participant_id}_${sample}.sam.bz2" \
            --bowtie2out "$dout_bt2/${participant_id}_${sample}.bowtie2.bz2" \
            -o "$dout_profiles/${participant_id}_${sample}_profiled.tsv" \
            --nproc "$threads"
    done
done

echo "finished metaphlan on metagenomes"

# Run metaphlan on dereplicated isolate genomes reads
for participant_dir in "$dir_isolates"/participant_*; do
    [ -d "$participant_dir" ] || continue

    participant_id=$(basename "$participant_dir" | sed 's/^participant_//')
    reads_dir="$participant_dir/raw_reads"

    [ -d "$reads_dir" ] || continue

    for r1 in "$reads_dir"/*_R1*.fastq.gz; do
        [ -f "$r1" ] || continue

        sample=$(basename "${r1%_R1*.fastq.gz}")
        r2="$reads_dir/${sample}_R2.fastq.gz"

        echo "isolate: $participant_id - $sample"

        metaphlan \
            "$r1,$r2" \
            --input_type fastq \
            -s "$dout_sams/${participant_id}_${sample}_isolate.sam.bz2" \
            --bowtie2out "$dout_bt2/${participant_id}_${sample}_isolate.bowtie2.bz2" \
            -o "$dout_profiles/${participant_id}_${sample}_isolate_profiled.tsv" \
            --nproc "$threads"
    done
done

echo "finished metaphlan on isolate reads"

# Extract consensus marker sequences from all samples
for sam_file in "$dout_sams"/*.sam.bz2; do
    [ -f "$sam_file" ] || continue

    echo "extracting markers from $(basename "$sam_file")"

    sample2markers.py \
        -i "$sam_file" \
        -o "$dout_markers" \
        -n "$threads"
done

echo "finished consensus marker extraction"

# Extract markers for the target species (Lactobacillus crispatus)
# These are reference marker sequences used by strainphlan
# The species is defined using its SGB identifier (t__SGB7045)

echo "extracting targeted markers for L. crispatus (t__SGB7045)"

extract_markers.py \
    -c t__SGB7045 \
    -o "$dout_targeted"

echo "finished marker extraction for target species"
