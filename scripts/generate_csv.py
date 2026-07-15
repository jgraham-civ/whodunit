import os
import re
import csv
import torch
import argparse
from Bio import Align
from Bio.Align import substitution_matrices
from Bio import SeqIO

# ---  HANDLING EMBEDDINGS --- #

def parse_layer_from_model(model_string):
    """
    Extracts the layer number from ESM model names.
    """
    # Look for _t[number]_ in the model string and extract 'number'
    match = re.search(r'_t(\d+)_', model_string)
    if not match:
        raise ValueError(f"Could not parse layer number from model name: {model_string}. "
                         f"Expected format containing '_t[layers]_' (e.g., esm2_t12_35M_UR50D).")
    return int(match.group(1))

def load_embeddings(dir_path, layer):
    """Loads all .pt files inside a specific folder and uses filenames as labels."""
    embeddings = []
    labels = []
    
    # List all files in the directory and filter for .pt files
    filenames = [f for f in os.listdir(dir_path) if f.endswith('.pt')]
    
    for f in filenames:
        # Build the complete path to the file
        file_path = os.path.join(dir_path, f)
        
        # Extract just the filename without the .pt extension for the label
        label = os.path.splitext(f)[0]
        labels.append(label)
        
        # Load tensor data onto CPU
        data = torch.load(file_path, map_location="cpu")
        mean_embedding = data["mean_representations"][layer]
        embeddings.append(mean_embedding)
        
    if not embeddings:
        raise FileNotFoundError(f"No .pt files found in directory: {dir_path}")
        
    return labels, torch.stack(embeddings)


# --- HANDLING METADATA --- #

PE_MEANINGS = {
    "1": "Experimental evidence at protein level",
    "2": "Experimental evidence at transcript level",
    "3": "Protein inferred from homology",
    "4": "Protein predicted",
    "5": "Protein uncertain",
}

def parse_fasta_headers(fasta_path):
    """Maps accession -> (description, organism, PE meaning) from raw UniProt headers. Covers every protein, BLAST hit or not."""
    headers_dict = {}
    with open(fasta_path) as f:
        for line in f:

            # Ignore non-header lines
            if not line.startswith(">"):
                continue
        
            # Look between the first and second |, e.g taking Q47KB1 from sp|Q47KB1|DYP_THEFY
            accession = line.split("|")[1]

            # Slice off the first character (">") and remove whitespace from both ends
            rest = line[1:].strip()
            
            # Extract metadata from string
            try:
                description = rest.split(" ", 1)[1].split(" OS=")[0] # Description is between first space and " OS="
                organism = rest.split(" OS=")[1].split(" OX=")[0] # Organism is between " OS=" and " OX="
                pe_code = rest.split(" PE=")[1].split(" ")[0] # Protein existence is between " PE=" and the next space
                
                # Insert all the values into the dictionary
                headers_dict[accession] = (description, organism, PE_MEANINGS.get(pe_code, ""))

            # Insert blank values if not found
            except IndexError:
                headers_dict[accession] = ("", "", "")
    return headers_dict

# --- ALIGNMENT AND SCORING --- #

def load_sequences(fasta_path):
    """Maps accession -> sequence string, for alignment."""
    sequences_dict = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Look between the first and second |, or take whole line if no | present
        accession = record.id.split("|")[1] if "|" in record.id else record.id
        # Assign accession to dictionary as the key, and the sequence to the value
        sequences_dict[accession] = str(record.seq)
    return sequences_dict

def build_aligner():
    aligner = Align.PairwiseAligner()
    # Aligner parameters are built to match BLASTp with BLOSUM62
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    return aligner

def compute_alignment_score(aligner, matrix, ref_seq, target_seq):
    """Computes pident/ppos via local alignment, matching BLAST's own definitions as closely as possible."""
    alignment = aligner.align(ref_seq, target_seq)[0]
    aligned_ref, aligned_target = alignment[0], alignment[1]

    # Start counters at 0
    aligned_positions = 0
    matches = 0
    positives = 0

    # Scoring algorithm
    for a, b in zip(aligned_ref, aligned_target):
        # Ignore non-matches
        if a == "-" or b == "-":
            continue
        # Count aligned positions
        aligned_positions += 1
        # Count matches
        if a == b:
            matches += 1
            positives += 1
        # Count positively scored substitutions in BLOSUM62 matrix
        elif matrix[a, b] > 0:
            positives += 1

    if aligned_positions == 0:
        return "", ""

    pident = (matches / aligned_positions) * 100
    ppos = (positives / aligned_positions) * 100
    return f"{pident:.3f}", f"{ppos:.3f}"

def load_blast_results(blast_tsv):
    """Maps accession -> (pident, ppos, bitscore, evalue). Only covers proteins BLAST actually reported."""
    results_dict = {}
    with open(blast_tsv) as f:
        for line in f:

            # Extract BLAST results
            sseqid, stitle, pident, ppos, bitscore, evalue = line.rstrip("\n").split("\t")
            accession = sseqid.split("|")[1] if "|" in sseqid else sseqid

            # Assign BLAST results to dict
            if accession not in results_dict:  # keep only the first (best-scoring) HSP seen
                results_dict[accession] = (pident, ppos, bitscore, evalue)
    return results_dict

# --- RUN SCRIPT --- #

def main():

    # Define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--proteome-dir", required=True)
    parser.add_argument("--proteome-fasta", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--blast-results", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--top-n", type=int, required=True)
    args = parser.parse_args()

    # Extract model layer
    layer = parse_layer_from_model(args.model)

    # Load reference embeddings from given layer
    ref_label, ref_emb = load_embeddings(args.reference_dir, layer)
    if ref_emb.shape[0] != 1:
        raise ValueError(f"Expected exactly 1 reference file in {args.reference_dir}, found {ref_emb.shape[0]}")

    # Load proteome embeddings from given layer
    target_labels, target_embs = load_embeddings(args.proteome_dir, layer)

    # Load metadata, blast results, and target sequence dictionaries
    headers_dict = parse_fasta_headers(args.proteome_fasta) # {accession_string: (description, organism, pe_meaning)}
    blast_dict = load_blast_results(args.blast_results) # {accession_string: (pident, ppos, bitscore, evalue)}
    target_sequences_dict = load_sequences(args.proteome_fasta) # {accession_string: sequence_string}

    # Load reference sequence
    ref_seq = SeqIO.read(args.reference_fasta, "fasta").seq

    # Load aligner and substitution matrix for later calculations
    aligner = build_aligner()
    matrix = substitution_matrices.load("BLOSUM62")
   
    # Compute cosine similarities of target embeddings to reference embedding
    similarities = torch.nn.functional.cosine_similarity(ref_emb, target_embs)

    # Return top_n results and keep note of their scores and indices
    top_n = min(args.top_n, len(target_labels))
    top_scores, top_indices = torch.topk(similarities, top_n)

    # Generate output CSV
    output_filename = f"top_{args.top_n}_similar_proteins.csv"
    with open(output_filename, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "Target", "Reference", "Cosine_Similarity",
                          "Pct_Identity", "Pct_Positive", "Bitscore", "Evalue",
                          "Description", "Organism", "Protein_Existence"])

        for rank, (score, idx) in enumerate(zip(top_scores, top_indices), start=1):
            target = target_labels[idx]
            description, organism, pe = headers_dict.get(target, ("", "", ""))

            target_seq = target_sequences_dict.get(target, "")
            pident, ppos = compute_alignment_score(aligner, matrix, ref_seq, target_seq) if target_seq else ("", "")

            _, _, bitscore, evalue = blast_dict.get(target, ("", "", "", ""))

            writer.writerow([rank, target, ref_label[0], f"{score.item():.4f}",
                              pident, ppos, bitscore, evalue,
                              description, organism, pe])

if __name__ == "__main__":
    main()