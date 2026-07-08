#!/usr/bin/env nextflow

/*
 * Rank similarity of proteome embeddings to reference embedding
 */
process RANK_SIMILARITY {

    tag "top${top_n}"

    conda "${projectDir}/envs/esm-environment.yml"

    input:
    path ref_dir
    path proteome_dir
    val top_n
    val model

    output:
    path "top_${top_n}_similar_proteins.csv", emit: ranked_csv

    script:
    """
    python ${projectDir}/scripts/similarity.py \\
        --reference-dir ${ref_dir} \\
        --proteome-dir ${proteome_dir} \\
        --model ${model} \\
        --top-n ${top_n}
    """
}