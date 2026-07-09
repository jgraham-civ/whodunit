import os
import re
import csv
import torch
import argparse

# ---  HANDLING EMBEDDINGS --- #

def parse_layer_from_model(model_string):
    """
    Extracts the layer number from ESM model names.
    """
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


# --- HANDLING METDATA --- #

PE_MEANINGS = {
    "1": "Experimental evidence at protein level",
    "2": "Experimental evidence at transcript level",
    "3": "Protein inferred from homology",
    "4": "Protein predicted",
    "5": "Protein uncertain",
}

def parse_fasta_headers(fasta_path):
    """Maps accession -> (description, organism, PE meaning) from raw UniProt headers."""
    headers = {}
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            accession = line.split("|")[1]
            rest = line[1:].strip()

            description = rest.split(" ", 1)[1].split(" OS=")[0]
            organism = rest.split(" OS=")[1].split(" OX=")[0]
            pe_code = rest.split(" PE=")[1].split(" ")[0]

            headers[accession] = (description, organism, PE_MEANINGS.get(pe_code, ""))
    return headers


# --- RUN SCRIPT --- #

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--proteome-dir", required=True)
    parser.add_argument("--proteome-fasta", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--top-n", type=int, required=True)
    
    args = parser.parse_args()
    layer = parse_layer_from_model(args.model)
    
    # 1. Load Reference Embedding
    ref_label, ref_emb = load_embeddings(args.reference_dir, layer)
    if ref_emb.shape[0] != 1:
        raise ValueError(f"Expected exactly 1 reference file in {args.reference_dir}, found {ref_emb.shape[0]}")
        
    # 2. Load Proteome Search Space
    target_labels, target_embs = load_embeddings(args.proteome_dir, layer)
    headers = parse_fasta_headers(args.proteome_fasta)
    
    # 3. Compute Cosine Similarity
    similarities = torch.nn.functional.cosine_similarity(ref_emb, target_embs)
    
    # 4. Get Top-N Rankings
    top_n = min(args.top_n, len(target_labels))
    top_scores, top_indices = torch.topk(similarities, top_n)
    
    # 5. Write Output
    output_filename = f"top_{args.top_n}_similar_proteins.csv"
    with open(output_filename, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Rank", "Target", "Reference", "Cosine_Similarity", "Description", "Organism", "Protein_Existence"])
        
        for rank, (score, idx) in enumerate(zip(top_scores, top_indices), start=1):
            target = target_labels[idx]
            description, organism, pe = headers.get(target, ("", "", ""))
            writer.writerow([rank, target, ref_label[0], f"{score.item():.4f}",
                              description, organism, pe])

if __name__ == "__main__":
    main()