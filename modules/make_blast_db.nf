#!/usr/bin/env nextflow

/*
 * Make BLAST database from proteome
 */
process MAKE_BLAST_DB {

    conda "${projectDir}/envs/esm-environment.yml"

    input:
    path proteome_fasta

    output:
    path "blast_db", emit: blast_db

    script:
    """
    mkdir -p blast_db
    makeblastdb -in ${proteome_fasta} -dbtype prot -parse_seqids -out blast_db/proteome
    """
}