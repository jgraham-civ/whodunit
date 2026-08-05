# whodunit:

> A Nextflow pipeline for identifying distant protein orthologs and orphaned enzymes by comparing their ESM structural embeddings.

---

## Overview

Traditional sequence alignment algorithms (like BLAST) struggle to identify functionally similar proteins when sequence identity is low. This pipeline addresses this problem by leveraging the ESM suite of protein language models (PLMs).

Much like language models in natural language processing learn grammar and meaning from large collections of text, PLMs are trained on vast numbers of protein sequences to learn the ‘grammar’ and ‘vocabulary’ that underpin their structure and function. PLMs can implicitly learn biochemical, physical, structural and functional properties from sequence alone, outputting these properties as embeddings: Numerical, vectorised representations of the input sequence.

Whilst the meaning of each embedding dimension is abstract and not directly interpretable, the similarity of embeddings at the global level should reflect functional similarity.

By calculating the **cosine similarity** between a reference protein embedding and embeddings in a taxon-wide 'search space' and returning the best results, this pipeline acts as a homology search tool.

**Inputs:** 
* UniProt ID of a reference protein
* Taxon ID of a target organism
* An ESM model selection (e.g. `esm2_t12_35M_UR50D`)
* The number of proteins to return in the ranked CSV
* User-defined batch name

**Output:**
* A ranked CSV of the target organism's proteins, scored by cosine similarity to the reference. Alongside protein accession IDs, this CSV reports protein metadata and sequence identities. The pipeline also runs protein BLAST against a local, taxon-wide BLAST database and reports the bitscores and e-values of each hit.

---

## Requirements & Dependencies

This pipeline uses **Conda** to manage its environment automatically via Nextflow. Ensure you have Conda (or Mamba) and Nextflow installed.

* **Nextflow** (>= 26.0)
* **Conda** The pipeline handles the rest by building from the included `esm-environment.yml`:

```yml
name: esm-env
    channels:
    - pytorch
    - bioconda
    - conda-forge
  dependencies:
    - python=3.10.20
    - pytorch=2.5.1
    - biopython=1.87
    - blast=2.17.0
    - pip
    - pip:
        - git+https://github.com/facebookresearch/esm.git@2b369911bb5b4b0dda914521b9475cad1656b2ac
```

*Note: Conda is enabled by default in the `nextflow.config`. Environments are cached in `${projectDir}/conda-envs`.*

---

## Usage & Profiles

You can run this pipeline by defining a custom search or by using the pre-configured profiles defined in `nextflow.config`.

**To run a custom search:**
```bash
nextflow run main.nf \
  --ref_id "YOUR_UNIPROT_ID" \
  --tax_id "YOUR_TAXON_ID" \
  --model "esm2_t12_35M_UR50D" \
  --top_n 100 \
  --batch_name "custom_run"
```

**To run a profile:**
```bash
nextflow run main.nf -profile <profile_name>
```

For proof-of-concept, we include a few profiles:

* **PET_t12:** screens *Thermobifida fusca* proteins against Poly(ethylene terephthalate) hydrolase (PETase) in *Ideonella sakaiensis*, correctly identifying the structural homolog cutinase TfCut2 as a top hit.

* **globin_t6:** screens *Glycine max* (Soybean) proteins against human myoglobin, correctly identifying leghemoglobins as structurally similar.

* **MTase_t12:** screens *Trypanosoma brucei brucei* (strain 927/4 GUTat10.1) proteins against human N-terminal N-methyltransferase 1 (NTMT1), as an extension of my MSc Project.

---

## Limitations

1. Whole-proteome embedding generation can be computationally cumbersome

2. Global (mean-pooled) cosine similarity is simplistic and may miss localized functional signal

3. Using only final-layer, mean-pooled embeddings underuses what ESM actually encodes

4. Limited input flexibility (currently allows only UniProt published sequences)

5. No statistical significance framework for cosine similarity

6. Molecular docking analysis, as done in my MSc thesis, is not currently built into the pipeline.

These limitations offer many avenues for future work.

---

## References

Altschul, S.F., Gish, W., Miller, W., Myers, E.W. and Lipman, D.J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), pp.403–410. doi:10.1016/S0022-2836(05)80360-2.

Bateman, A., Martin, M.-J., Orchard, S., Magrane, M., Adesina, A., Ahmad, S.,
Bowler-Barnett, E.H., Bye-A-Jee, H., Carpentier, D., Denny, P., Fan, J., Garmiri, P.,
Jose, L., Hussein, A., Ignatchenko, A., Insana, G., Ishtiaq, R., Joshi, V., Jyothi, D.
and Kandasaamy, S. (2024). UniProt: the Universal Protein Knowledgebase in 2025.
Nucleic Acids Research, 53(D1). doi:https://doi.org/10.1093/nar/gkae1010.

Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., Smetanin, N., Verkuil, R., Kabeli,
O., Shmueli, Y., dos Santos Costa, A., Fazel-Zarandi, M., Sercu, T., Candido, S. and
Rives, A. (2023). Evolutionary-scale prediction of atomic-level protein structure with a
language model. Science, 379(6637), pp.1123–1130.
doi:https://doi.org/10.1126/science.ade2574.
