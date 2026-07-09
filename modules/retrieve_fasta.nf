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

    if ! grep -q "^>" ${ref_id}.fasta; then
        echo "ERROR: downloaded file does not contain valid FASTA content:" >&2
        cat ${ref_id}.fasta >&2
        exit 1
    fi
    """
}

/*
 * Retrieve NCBI and UniProt sequences for a given taxonomy ID
 */
process RETRIEVE_PROTEOME {

    errorStrategy 'retry'
    maxRetries 3

    input:
    val tax_id

    output:
    path "${tax_id}.fasta", emit: tax_fasta

    script:
    """
    url="https://rest.uniprot.org/uniprotkb/search?query=(taxonomy_id:${tax_id})&format=fasta&size=500"
    > ${tax_id}.fasta

    while [ -n "\$url" ]; do
        headers=\$(mktemp)
        curl -s -D "\$headers" "\$url" >> ${tax_id}.fasta

        url=\$(grep -o '<[^>]*>; rel="next"' "\$headers" | sed -E 's/<(.*)>; rel="next"/\\1/')
        rm "\$headers"
    done

    if ! grep -q "^>" ${tax_id}.fasta; then
        echo "ERROR: downloaded file does not contain valid FASTA content:" >&2
        cat ${tax_id}.fasta >&2
        exit 1
    fi
    """
}