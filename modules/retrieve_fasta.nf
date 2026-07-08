#!/usr/bin/env nextflow

/*
 * Retrieve UniProt sequence for a given UniProt ID
 */
process RETRIEVE_REFERENCE {

    input:
    val ref_id

    output:
    path "${ref_id}.fasta", emit: ref_fasta

    script:
    """
    curl "https://rest.uniprot.org/uniprotkb/${ref_id}.fasta" -o ${ref_id}.fasta
    """
}

/*
 * Retrieve NCBI and UniProt sequences for a given taxonomy ID
 */
process RETRIEVE_PROTEOME {

    input:
    val tax_id

    output:
    path "${tax_id}.fasta", emit: tax_fasta

    script:
    """
    curl "https://rest.uniprot.org/uniprotkb/stream?query=(taxonomy_id:${tax_id})&format=fasta" -o ${tax_id}.fasta
    """
}