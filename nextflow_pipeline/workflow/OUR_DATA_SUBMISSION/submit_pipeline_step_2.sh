#!/usr/bin/env bash


nextflow ../main.nf --ont_reads_fq "../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.fastq" \
    --ont_reads_txt "../../../../../../scratch/bag222/data/ont_data/09-02-2022_uky_6ad_6ct/*.txt" \
    --ref "../../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../../references/Homo_sapiens_plus_new_high_cofidence.GRCh38.109.gtf" \
    --housekeeping "../../references/hg38.HouseKeepingGenes.bed" \
    --cdna_kit "PCS111" \
    --out_dir "./GRCh38_quant_mapq10_gtf_109_and_high-confidence_OUR_DATA/" \
    --is_discovery "False" \
    --bambu_track_reads "False" \
    --mapq "10" \
    --step "2" \
    --is_chm13 "False" -resume
