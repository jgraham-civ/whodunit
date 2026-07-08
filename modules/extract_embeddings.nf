#!/usr/bin/env nextflow

/*
 * Clean FASTA headers
 */
process CLEAN_HEADERS {
    
    tag "${fasta_file}"

    input:
    path fasta_file

    output:
    path "cleaned_${fasta_file}", emit: cleaned_fasta

    script:
    """
    awk -F'|' '/^>/{print ">"\$2; next} {print}' ${fasta_file} > cleaned_${fasta_file}
    """
}

/*
 * Compute ESM embedding for reference sequence
 */
process EMBED_REFERENCE {

    conda "${projectDir}/envs/esm-environment.yml"

    input:
    path ref_fasta
    val model

    output:
    path "reference_emb", emit: ref_emb

    script:
    """
    esm-extract ${model} ${ref_fasta} reference_emb --include mean
    """
}

/*
 * Compute ESM embeddings for proteome sequences
 */
process EMBED_PROTEOME {

    conda "${projectDir}/envs/esm-environment.yml"

    input:
    path tax_fasta
    val model

    output:
    path "proteome_emb", emit: tax_emb

    script:
    """
    esm-extract ${model} ${tax_fasta} proteome_emb --include mean --toks_per_batch 1024
    """
}