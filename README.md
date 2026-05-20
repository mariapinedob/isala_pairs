# Isala mother-daughter pairs

This repository contains scripts used to process and analyze 16S rRNA gene sequencing data, shotgun metagenomic sequencing data, and genome data derived from vaginal isolates from mother–daughter pairs in the Isala cohort. Isala is a pioneering citizen-science program that aims to better understand the female microbiome using advanced DNA technologies and enrolled over 3,345 women in its first phase. This analysis and its results are described in the paper "*Vaginal microbiome tracking finds species and strain sharing between mothers and their adult daughters*" (Pinedo-Bardales, Erreygers et al.).

The scripts are organized into four workflows, each corresponding to a different part of the analysis: quality control, species-level analysis, and strain‑level analysis.

If you want to run these scripts yourself, you will need to update the file paths to match the location of your input files. In each script, these paths are listed at the top. In addition, you will need to install the dependencies specified in each script.

## How to run?

### Step 1: Create folders for data, results and scripts (src)

``` bash
{cd isala_pairs}
mkdir data results src
```

### Step 2: Download necessary input files

Tabular data for the analysis in **workflows 01 and 02** are available on [Zenodo](https://doi.org/10.5281/zenodo.20307421). This repository also includes tables listing ENA accession run IDs for metagenomic sequencing data, 16S rRNA gene sequencing data, and isolate genomes (along with their corresponding GenBank assemblies) required for **workflow 03**. Raw sequencing data and assembled genomes, where applicable, should be downloaded from ENA using the accession numbers provided in these tables. Gene prediction should then be performed on the assembled genomes, using tools such as [Prodigal](https://github.com/hyattpd/Prodigal), to generate nucleotide (`.ffn`) and amino acid (`.faa`) gene sequences.

For raw metagenomic sequencing data from the Isala and the cohort of France et al., the following directory structure is recommended:

```         
data/
└── metagenomic_data/
    ├── isala/
    │   ├── participant_ID_01/
    │   │   ├── raw_reads/
    │   │   │   ├── sample_1_R1.fastq.gz
    │   │   │   └── sample_1_R2.fastq.gz
    │   ├── participant_ID_02/
    │   │   └── ...
    │   └── ...
    └── non_european_cohort/
        ├── participant_ID_01/
        │   ├── raw_reads/
        │   │   ├── sample_1_R1.fastq.gz
        │   │   └── sample_1_R2.fastq.gz
        ├── participant_ID_02/
        │   └── ...
        └── ...
```

For isolate genomes, the following directory structure is recommended:

```         
data/
└── isolate_genomes/
    ├── participant_ID_01/
    │   ├── raw_reads/
    │   │   ├── isolate_1_R1.fastq.gz
    │   │   └── isolate_1_R2.fastq.gz
    │   ├── assembly/
    │   │   └── isolate_1_assembly.fna.gz
    │   ├── genes/
    │   │   ├── isolate_1.fna
    │   │   ├── isolate_1.ffn
    │   │   └── isolate_1.faa
    └── ...
```

For **workflow 04**, which includes the phylogenetic reconstruction of *Lactobacillus crispatus*, genomes from the GTDB database must be retrieved from the release R226. A reference script for downloading these genomes is available [here](https://github.com/SWittouck/legen/tree/master/src/01_prepare_genomes). In addition, *L. crispatus* genomes from the Isala cohort should be included. These genomes can be downloaded from the European Nucleotide Archive (ENA) under BioProject PRJEB105013. Gene prediction should also be performed on all assembled genomes prior to downstream analyses. The following directory structure is suggested:

```         
data/
└── phylogeny_references/
    ├── gtdb_crispatus/
    │   ├── genomes/
    │   │   ├── genome1.fna
    │   │   └── ...
    │   ├── annotations/
    │   │   ├── genome1.ffn
    │   │   ├── genome1.faa
    │   │   └── ...
    │
    ├── isala_crispatus/
    │   ├── genomes/
    │   │   ├── isolate1.fna
    │   │   ├── isolate2.fna
    │   │   └── ...
    │   ├── annotations/
    │   │   ├── isolate1.ffn
    │   │   ├── isolate2.faa
    │   │   └── ...
```

### Step 3: Data analysis

The analyses include quality-controlled comparisons between sequencing modalities, species‑level similarity assessments using permutation-based approaches and strain‑level investigations based on isolate genomes and metagenomic data. Each analysis workflow is implemented independently and is described in detail in the sections below, where the corresponding scripts, inputs, and methodological choices are explained.

**Note:** For workflow 04, representative genomes sampling should be carried out for both GTDB and Isala isolate genomes. To do this, follow the procedure described in the `01_sample_isolate_genomes` script from workflow 03. After sampling, the resulting directories should follow a consistent structure, as shown below:

```         
results/
└── derep_public_genomes/
    ├── gtdb/
    │   ├── core/
    │   └── representatives/
    │       └── seeds.txt
    │
    └── isala/
        ├── core/
        └── representatives/
            └── seeds.txt
```

## Workflows

-   `01_quality_control`: **Script for comparing 16S rRNA gene sequencing and deep metagenomic shotgun sequencing (DMGS) data**

    -   This workflow includes taxonomy harmonization, merging of high-quality16S and shotgun sequencing data and visualization.
    -   Input read count tables (`tt_isala_16s_hq`, `tt_isala_shotgun_krakendb`) and the corresponding metadata can be found in the [Zenodo repository](https://doi.org/10.5281/zenodo.20307421). These tables were created based on the raw (meta)genomic reads using in-house workflows described previously ([Lebeer et al. 2023](https://www.nature.com/articles/s41564-023-01500-0), [Vander Dock et al. 2025](https://www.sciencedirect.com/science/article/pii/S2211124725009428?ssrnid=5243601&dgcid=SSRN_redirect_SD)) and stored in [tidytacos](https://github.com/LebeerLab/tidytacos) format.
    -   Quality control of 16S rRNA gene sequencing data was performed before following the procedure described in the [first paper of the Isala cohort](https://www.nature.com/articles/s41564-023-01500-0). DNA derived from 59 out of our 62 participants passed quality control and were included in this comparative analysis.

-   `02_species_level`: **Scripts for assessing species‑level similarity between mother-daughter pairs**

    -   This workflow includes (i) species‑level exploration of vaginal microbiome composition using tidytacos objects, (ii) permutation‑based analyses comparing metagenomic profiles from mothers and daughters in the Isala cohort and the cohort described by [France et al. (2022)](https://pubmed.ncbi.nlm.nih.gov/36288274/), and (iii) performance of a generalized linear model to assess the association between the maternal dominance of L. crispatus and the dominance status of their daughters.

-   `03_strain_level`: **Scripts for investigating strain sharing within and between mother-daughter pairs**

    -   This workflow includes strain-level analyses based on both unique isolate genomes per participant and metagenomic reads using the strain-level metagenome analysis tool [inStrain](https://instrain.readthedocs.io/en/latest/tutorial.html#tutorial-2-running-instrain-using-a-public-genome-database) ([Olm et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33462508/)). 
    -   inStrain analyses were performed following reference-dependent and reference-independent approaches.

-   `04_crispatus_phylogeny:` **Scripts for performing phylogenetic analysis of *Lactobacillus crispatus***

    -   This workflow performs phylogenetic reconstruction of *Lactobacillus crispatus* by integrating isolate genomes, metagenomic data, and publicly available reference genomes from GTDB and the Isala cohort. The analysis is based on marker genes and carried out using [StrainPhlAn 4](http://segatalab.cibio.unitn.it/tools/strainphlan/) ([Blanco-Miguez et al., 2023](https://pubmed.ncbi.nlm.nih.gov/36823356/); [Truong et al. 2017](https://genome.cshlp.org/content/early/2017/02/06/gr216242116)), followed by visualization of the resulting phylogeny.
