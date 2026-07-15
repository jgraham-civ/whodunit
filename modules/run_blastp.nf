#!/usr/bin/env nextflow

/*
 * Run BLASTp for reference sequence against database
 */
process RUN_BLASTP {

    conda "${projectDir}/envs/esm-environment.yml"

    input:
    path ref_fasta
    path blast_db

    output:
    path "blast_results.tsv", emit: blast_results

    script:
    """
    blastp \\
        -query ${ref_fasta} \\
        -db ${blast_db}/proteome \\
        -out blast_results.tsv \\
        -outfmt "6 sseqid stitle pident ppos bitscore evalue" \\
        -evalue 1000 \\
        -max_target_seqs 10000
    """
}