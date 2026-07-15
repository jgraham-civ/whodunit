#!/usr/bin/env nextflow

/*
 * Rank similarity of proteome embeddings to reference embedding
 */
process GENERATE_CSV {

    tag "top${top_n}"
    conda "${projectDir}/envs/esm-environment.yml"

    input:
    path ref_dir
    path proteome_dir
    path proteome_fasta
    path reference_fasta
    path blastp_results
    val top_n
    val model

    output:
    path "top_${top_n}_similar_proteins.csv", emit: ranked_csv

    script:
    """
    python ${projectDir}/scripts/generate_csv.py \\
        --reference-dir ${ref_dir} \\
        --proteome-dir ${proteome_dir} \\
        --proteome-fasta ${proteome_fasta} \\
        --reference-fasta ${reference_fasta} \\
        --blast-results ${blastp_results} \\
        --model ${model} \\
        --top-n ${top_n}
    """
}