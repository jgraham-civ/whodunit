#!/usr/bin/env nextflow

// Module INCLUDE statements
include { RETRIEVE_REFERENCE } from './modules/retrieve_fasta.nf'
include { RETRIEVE_PROTEOME } from './modules/retrieve_fasta.nf'
include { CLEAN_HEADERS as CLEAN_REF_HEADERS} from './modules/extract_embeddings.nf'
include { CLEAN_HEADERS as CLEAN_TAX_HEADERS} from './modules/extract_embeddings.nf'
include { EMBED_REFERENCE } from './modules/extract_embeddings.nf'
include { EMBED_PROTEOME } from './modules/extract_embeddings.nf'
include { RANK_SIMILARITY } from './modules/rank_similarity.nf'

/*
 * Pipeline parameters
 */
params {
   ref_id: String
   tax_id: String
   model: String
   top_n: Integer
   batch_name: String
}

workflow {

    main:
    // Create input channels (if needed)

    // Retrieve Retrieve UniProt sequences
    RETRIEVE_REFERENCE(params.ref_id)
    RETRIEVE_PROTEOME(params.tax_id)

    // Clean FASTA headers
    ref_fasta_cleaned = CLEAN_REF_HEADERS(RETRIEVE_REFERENCE.out.ref_fasta)
    tax_fasta_cleaned = CLEAN_TAX_HEADERS(RETRIEVE_PROTEOME.out.tax_fasta)

    // Compute ESM embeddings
    EMBED_REFERENCE(ref_fasta_cleaned, params.model)
    EMBED_PROTEOME(tax_fasta_cleaned, params.model)

    // Perform similarity scoring & ranking
    RANK_SIMILARITY(
        EMBED_REFERENCE.out.ref_emb,   
        EMBED_PROTEOME.out.tax_emb,
        RETRIEVE_PROTEOME.out.tax_fasta,    
        params.top_n,
        params.model 
    )

    publish:
    // Declare outputs to publish
    ref_fasta = RETRIEVE_REFERENCE.out.ref_fasta
    tax_fasta = RETRIEVE_PROTEOME.out.tax_fasta
    ref_emb = EMBED_REFERENCE.out.ref_emb
    tax_emb = EMBED_PROTEOME.out.tax_emb
    ranked_csv = RANK_SIMILARITY.out.ranked_csv

}

output {
    // Configure publish targets
    ref_fasta {
        path "${params.batch_name}/sequences"
    }
    tax_fasta {
        path "${params.batch_name}/sequences"
    }
    ref_emb {
        path "${params.batch_name}/embeddings"
    }
    tax_emb {
        path "${params.batch_name}/embeddings"
    }
    ranked_csv {
        path "${params.batch_name}/hits"
    }
}