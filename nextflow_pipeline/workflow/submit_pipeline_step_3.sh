#!/usr/bin/env bash


nextflow ./main.nf --bambu_rds "./results/GRCh38_no_discovery_mapq10_annotation_109/bambu_prep/*.rds" \
    --ref "../references/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --annotation "../references/Homo_sapiens.GRCh38.109.gtf" \
    --out_dir "./GRCh38_no_discovery_mapq10_annotation_109/" \
    --is_discovery "False" \
    --bambu_track_reads "False" \
    --multiqc_input "./results/GRCh38_no_discovery_mapq10_annotation_109/multiQC_input/**" \
    --multiqc_config "../references/multiqc_config.yaml" \
    --fai "./results/GRCh38_no_discovery_mapq10_annotation_109/fai/*.fai" \
    --step 3 \
    --is_chm13 "False" -resume
