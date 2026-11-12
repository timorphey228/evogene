# EvoGene · Gene origin reconstruction (hackathon prototype)

Evidence-aware classification of protein-coding gene origins, initially targeting *Drosophila melanogaster*.

## Sources and scientific basis

- Peng & Zhao (2024), *The origin and structural evolution of de novo genes in Drosophila*, Nature Communications: https://www.nature.com/articles/s41467-024-45028-1
- Original LiZhaoLab research code: https://github.com/LiZhaoLab/DrosophilaDenovoGene

Upstream analysis methods: ProcessingWholeGenomeAlignments/ortholog_mapping.py, ortholog_relations.py, find_last_orf.py, check_homologs.py, extract_raw_hit.py, extract_alnscore_genewise.py. Their original pipeline requires progressive Cactus/HAL and laboratory files. **This MVP does not reproduce those alignment analyses**: it independently integrates evidence that has already been obtained. The application does not copy upstream source code.

## Run

Python 3.10+.

    python -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt
    streamlit run app.py

CLI (no third-party dependencies):

    python -m evogene.cli --evidence data/demo_evidence.csv --properties /path/to/41467_2024_45028_MOESM10_ESM.xlsx --output results.json
    python -m unittest discover -s tests -v

Windows virtual environment: .venv\Scripts\activate

## Data

The case's supplementary XLSX, Figure 2 a,c,d,e,f,g,h, contains FBID, protein length, ISD, TM probabilities, RSA, GC content, SignalP, expression, tissue specificity, research de novo gene annotation and evolutionary branch age. **These are reference metadata, not direct proof of evolutionary origin.** The original workbook can be uploaded in the app or supplied to the CLI; it is not embedded into the repo.

The evidence file has columns:

    gene_id,species,paralog_detected,distant_coding_homolog,syntenic_outgroup_locus,outgroup_locus_noncoding,homolog_search_completed,rearrangement_evidence,outgroup_count,alignment_identity,notes

Boolean values are yes, no, unknown. Unknown is **not** no. Demo evidence is *SYNTHETIC* for testing despite using gene IDs from the supplement; its classifications are **not new biological findings**.

A de novo proposal requires a syntenic noncoding locus in outgroups, at least two outgroup species, and completed search for homologs without contradictory coding homologs. This is a conservative prototype heuristic, not a validated published threshold. Weight scores are **not posterior probabilities or calibrated confidences**.

Possible outputs: de_novo, duplication, ancient_divergence, other_mechanism, undetermined. Each result includes supporting evidence, conflicts, missing evidence and qualitative scores.

## Future validation

1. Obtain real DNA/protein FASTA, genome references, Cactus/HAL liftover and synteny outputs for related species.
2. Generate gene-by-gene evidence with exact references, aligner thresholds and provenance; detect hidden homology before any de novo claim.
3. Compare candidate gene trees and alignments, validate classification on independent positive and negative controls.
4. Estimate uncertainties scientifically; add phylogenetic visualization and automated sequence-analysis jobs.

Citation: Peng & Zhao 2024. Branch labels Br1/Br9 are kept as labels; not converted to absolute ages.
