#!/usr/bin/env bash

# This script performs inStrain profile and inStrain compare following a reference-dependent 
# approach on metagenomic sequencing data from mother-daughter pairs in the Isala cohort
# and from the cohort of France et al (2022). For Isala, an additional step is provide to 
# also profile the isolate genomes.
#
# Dependencies: bowtie2, inStrain, wget, parse_stb.py, prodigal

# Root directory
dir_project="../../"

# Reference genome directory
dir_ref="$dir_project/data/instrain/reference_genomes"
mkdir -p "$dir_ref"

# File listing references (name=url)
fin_genomes="$dir_ref/reference_genomes.txt"

# Metagenomic reads
dir_metagenomes="$dir_project/data/metagenomic_data"

# Dereplicated isolate reads (applicable for Isala only)
dir_isolates="$dir_project/results/isala/drep"

threads_bt2=8
threads_instrain=8

cohort="isala" # For non-European cohort: "non_eu_cohort"

# Outputs for metagenomes
dout_results="$dir_project/results/strain_level_analyses/${cohort}/ref_dependent/metagenomes_only"
dout_ref="$dout_results/reference"
dout_mapping="$dout_results/mapping"
dout_profiles="$dout_results/profiles"

# Outputs for isolate integration + comparison for Isala cohort
dout_results2="$dir_project/results/strain_level_analyses/isala/ref_dependent/integrative_approach"
dout_profiles_isolates="$dout_results2/profiles_isolates"
dout_compare="$dout_results2/compare"

mkdir -p "$dir_ref" "$dout_ref" "$dout_mapping" "$dout_profiles"
mkdir -p "$dout_profiles_isolates" "$dout_compare"

# Reference genome files used by bowtie2 and instrain
fout_fna="$dout_ref/reference_genomes.fna"
fout_stb="$dout_ref/reference_genomes.stb"

# Define which genomes to include in the reference set
cat > "$fin_genomes" <<EOF
L_crispatus=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/987/235/GCA_018987235.1_ASM1898723v1/GCA_018987235.1_ASM1898723v1_genomic.fna.gz
L_iners=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/160/875/GCA_000160875.1_ASM16087v1/GCA_000160875.1_ASM16087v1_genomic.fna.gz
G_vaginalis=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/042/655/GCF_001042655.1_ASM104265v1/GCF_001042655.1_ASM104265v1_genomic.fna.gz
G_leopoldii=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/293/675/GCF_003293675.1_ASM329367v1/GCF_003293675.1_ASM329367v1_genomic.fna.gz
EOF

# Download genomes only if not already present locally
while IFS='=' read -r name url; do
    if [ ! -f "$dir_ref/${name}.fna" ]; then
        wget -q -O "$dir_ref/${name}.fna.gz" "$url"
        gunzip -f "$dir_ref/${name}.fna.gz"
    fi
done < "$fin_genomes"

# Merge all references into a single fasta for mapping
cat "$dir_ref"/*.fna > "$fout_fna" 2>/dev/null

# Create scaffold-to-bin mapping for instrain
parse_stb.py --reverse -f "$dir_ref"/*.fna -o "$fout_stb"

# Build bowtie2 index once (skip if already exists)
[ -f "$fout_fna.1.bt2" ] || bowtie2-build "$fout_fna" "$fout_fna"

# Map each metagenomic sample against the reference genomes
for participant_dir in "$dir_metagenomes"/"$cohort"/*; do
    [ -d "$participant_dir" ] || continue

    participant_id="$(basename "$participant_dir")"
    reads_dir="$participant_dir/raw_reads"
    [ -d "$reads_dir" ] || continue

    for r1 in "$reads_dir"/*_R1.fastq.gz; do
        [ -f "$r1" ] || continue

        sample="$(basename "${r1%_R1.fastq.gz}")"

        bowtie2 \
            -p "$threads_bt2" \
            -x "$fout_fna" \
            -1 "$reads_dir/${sample}_R1.fastq.gz" \
            -2 "$reads_dir/${sample}_R2.fastq.gz" \
            -S "$dout_mapping/ref_genomes_vs_${participant_id}_${sample}.sam"
    done
done

# Run instrain profiling on all mapped metagenomes
for sam in "$dout_mapping"/*.sam; do
    [ -f "$sam" ] || continue

    sample="$(basename "${sam%.sam}")"

    inStrain profile \
        "$sam" \
        "$fout_fna" \
        -o "$dout_profiles/$sample" \
        --min_read_ani 0.95 \
        -p "$threads_instrain" \
        --stb "$fout_stb"
done

echo "finished metagenome profiling"

# Process isolate reads the same way, but store results separately
for participant_dir in "$dir_isolates"/participant_*; do
    [ -d "$participant_dir" ] || continue

    participant_id="$(basename "$participant_dir")"
    reads_dir="$participant_dir/raw_reads"
    [ -d "$reads_dir" ] || continue

    for r1 in "$reads_dir"/*_R1.fastq.gz; do
        [ -f "$r1" ] || continue

        sample="$(basename "${r1%_R1.fastq.gz}")"

        sam_file="$dout_mapping/ref_genomes_vs_${participant_id}_${sample}_isolate.sam"

        bowtie2 \
            -p "$threads_bt2" \
            -x "$fout_fna" \
            -1 "$reads_dir/${sample}_R1.fastq.gz" \
            -2 "$reads_dir/${sample}_R2.fastq.gz" \
            -S "$sam_file"

        inStrain profile \
            "$sam_file" \
            "$fout_fna" \
            -o "$dout_profiles_isolates/${participant_id}_${sample}_isolate" \
            --min_read_ani 0.95 \
            -p "$threads_instrain" \
            --stb "$fout_stb"
    done
done

echo "finished profiling from isolate genomes"

# Compare strain profiles across metagenomes and isolates together
# This produces the integrative strain-level similarity results
inStrain compare \
    -i "$dout_profiles"/*.IS "$dout_profiles_isolates"/*.IS \
    -s "$fout_stb" \
    -p "$threads_instrain" \
    -o "$dout_compare" \
    --store_mismatch_locations \
    --database_mode

echo "finished integrative comparison"
