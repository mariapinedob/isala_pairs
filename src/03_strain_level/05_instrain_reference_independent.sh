
#!/usr/bin/env bash

# This script implements a reference-independent inStrain workflow by mapping
# metagenomic reads and isolate reads to each representative isolate genome
# obtained after dereplication. For each isolate, reads from all individuals
# belonging to the same family (e.g. mother-daughter pairs) are mapped
# independently.

# Set project directories and input data
dir_project="../../"

# Load directories
dir_isolates="$dir_project/results/isala/drep" # Representative isolates (after dereplication)
dir_metagenomes="$dir_project/data/metagenomic_data/isala" # Metagenomic reads per participant

# Number of threads for mapping and instrain
threads_bt2=8
threads_instrain=8

# Output structure
dout_results="$dir_project/results/strain_level_analyses/isala/ref_independent"
dout_ref="$dout_results/reference"      # Per-isolate reference genomes
dout_mapping="$dout_results/mapping"    # Read mapping files (SAM)
dout_profiles="$dout_results/profiles"  # inStrain profiles 

mkdir -p "$dout_ref" "$dout_mapping" "$dout_profiles"

# Extract family ID from participant (e.g. F01M -> F01)
get_family() {
    echo "$1" | grep -oE '^F[0-9]+'
}

# Loop over each participant that has representative isolates
for participant_dir in "$dir_isolates"/participant_*; do
    [ -d "$participant_dir" ] || continue

    # Participant naming
    participant_id=$(basename "$participant_dir")        # E.g.: participant_F01M
    short_id=${participant_id#participant_}              # E.g.: F01M
    family_id=$(get_family "$short_id")                  # E.g.: F01

    # Directory with representative assemblies
    assembly_dir="$participant_dir/assemblies"
    [ -d "$assembly_dir" ] || continue

    echo "processing $short_id (family $family_id)"

    # Identify all members of the same family (e.g. F01M, F01D)
    family_members=()
    for p in "$dir_isolates"/participant_${family_id}*; do
        [ -d "$p" ] || continue
        family_members+=("$(basename "$p" | sed 's/^participant_//')")
    done

    # Loop over each representative isolate genome
    for fna in "$assembly_dir"/*.fna; do
        [ -f "$fna" ] || continue

        isolate_name=$(basename "${fna%.fna}")
        ref_prefix="${short_id}_${isolate_name}"

        # Reference files for this isolate
        ref_fna="$dout_ref/${ref_prefix}.fna"
        ref_stb="$dout_ref/${ref_prefix}.stb"
        ref_index="${ref_fna%.fna}"

        echo "  isolate $isolate_name"

        # Copy genome and make scaffold names unique (required for instrain)
        sed "s/^>/\>${ref_prefix}_/" "$fna" > "$ref_fna"

        # Create scaffold-to-bin mapping
        parse_stb.py --reverse -f "$ref_fna" -o "$ref_stb"

        # Build bowtie2 index if not already present
        [ -f "${ref_index}.1.bt2" ] || bowtie2-build "$ref_fna" "$ref_index"

        # Map all metagenomes from this family to the isolate
        for member in "${family_members[@]}"; do

            meta_dir="$dir_metagenomes/$member/raw_reads"
            [ -d "$meta_dir" ] || continue

            for r1 in "$meta_dir"/*_R1.fastq.gz; do
                [ -f "$r1" ] || continue

                sample=$(basename "${r1%_R1.fastq.gz}")
                sam_out="$dout_mapping/${ref_prefix}_vs_${member}_${sample}.sam"

                bowtie2 \
                    -p "$threads_bt2" \
                    -x "$ref_index" \
                    -1 "$meta_dir/${sample}_R1.fastq.gz" \
                    -2 "$meta_dir/${sample}_R2.fastq.gz" \
                    -S "$sam_out"
            done
        done

        # map all isolate reads from the same family
        for member in "${family_members[@]}"; do

            reads_dir="$dir_isolates/participant_${member}/raw_reads"
            [ -d "$reads_dir" ] || continue

            for r1 in "$reads_dir"/*_R1*.fastq.gz; do
                [ -f "$r1" ] || continue

                sample=$(basename "${r1%_R1*.fastq.gz}")
                sam_out="$dout_mapping/${ref_prefix}_vs_${member}_${sample}_isolate.sam"

                bowtie2 \
                    -p "$threads_bt2" \
                    -x "$ref_index" \
                    -1 "$reads_dir/${sample}_R1.fastq.gz" \
                    -2 "$reads_dir/${sample}_R2.fastq.gz" \
                    -S "$sam_out"
            done
        done

        # Run instrain profiling on all mappings generated for this isolate
        for sam in "$dout_mapping"/${ref_prefix}_vs_*.sam; do
            [ -f "$sam" ] || continue

            sample=$(basename "${sam%.sam}")

            inStrain profile \
                "$sam" \
                "$ref_fna" \
                -o "$dout_profiles/$sample" \
                --min_read_ani 0.95 \
                -p "$threads_instrain" \
                --stb "$ref_stb"
        done

    done

done

echo "finished reference-independent profiling"

