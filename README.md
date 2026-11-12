# EvoGene · Case 2 (ITMO «СИНТЕЗ», 2026)

An **evidence-first comparative-genomics prototype** for reconstructing the
possible origin of protein-coding genes in *Drosophila* and relatives.

> **Scientific scope:** screening, not a validated automated origin classifier.
> **Missing BLAST/minimap2 hits and absent GFF3 CDS annotation never prove
> de novo origin or gene loss.** No made-up posterior probabilities or trees.

## What is implemented

- Upload a target gene as DNA FASTA or supply a sequence/coordinates in a target
  assembly; use available NCBI Datasets assemblies and GFF3 for any chosen subset of up to 19 comparison species (the focal species plus the 19 others from the 20-species panel).
- Map DNA with **minimap2**, retain mapping coordinates, coverage, identity,
  CDS-annotation overlap and assembly accessions.
- Screen for additional non-overlapping loci, CDS exon structure, preliminary
  protein similarity, and approximate *flank-based* synteny.
- Attempt bounded local similarity rescue **inside** an observed syntenic window,
  when the query and GFF3 are available. This is candidate DNA homology,
  not a confirmed ancestral noncoding sequence.
- Produce **independent evidence statuses** for de novo, DNA duplication,
  retroposition, relocation/transposition, and fusion. Other mechanisms,
  ancient divergence and **insufficient evidence** remain explicit possibilities.
  The *case-required* five broad conclusions are **not** yet automatically
  resolved by robust statistical hypothesis comparison.
- Survey complete **ATG-to-stop genomic ORFs** on both DNA strands (at most
  100 kb). These are **not** spliced CDS, confirmed ancestral coding states,
  frameshifts or verified translation.
- Trace each analysis stage and show errors instead of suppressing missing evidence.

## Published supplementary sequence file

The supplied `supplementary_file_unannotated_alignments.fasta` (~53.9 MB) is a
**nonstandard three-section research export**, NOT one aligned FASTA. It contains
780 gene IDs and 17,011 entries per section (51,033 total), with 20 species
codes. Sections:

- I: `raw sequence for` — Cactus-selected genomic DNA fragments;
- II: `spaln for` — spliced sequence predictions;
- III: `genewise for` — independent spliced sequence predictions.

Open `python -m streamlit run app.py`, choose `Peng & Zhao` data mode,
upload the original `.fasta`, select **any** of the indexed 780 genes, select
raw/Spaln/Genewise and any subset of 20 species, then click analyze. You can
export only the selected sequences as **ordinary, unaligned FASTA** and use a
separate appropriate alignment workflow. Do NOT upload the 53.9 MB source as
an already aligned MSA or as the sequence of one gene. The repository does not
need to contain a copy of this large dataset.

The parser was corrected to respect all three sections, preserve header
`score`/`shift` metadata, and avoid O(N×M) alignment on very long segments.
Files matching the published layout are supported; no biological truth labels
are inferred from author selection. See `CASE_FILE_GUIDE.md`.

**How many species?** The old five-species NCBI panel was a prototype-only
hard-coded implementation detail, not a scientific limit. `auto_app.py` now
uses a **selectable 20-species study panel** (up to 19 relatives for one
focal species). `app.py` studies all 20 codes already present in the
supplementary extract without redownloading genomes. Large NCBI analyses
consume substantial disk and compute; outgroup coverage matters more than
maximizing the count.

All stages have a selectable scientific annotation in the UI, including the
question, implemented method, measured output, limitations and visualization.
The source is `evogene/analysis_notes.py`.

## Visualizations

Both `auto_app.py` (automatic NCBI workflow) and `app.py` (local data/manual
case inputs) now expose interactive panels:

1. **Genome mapping blocks** on query coordinates, with per-hit annotations,
   identity, coverage and target coordinates.
2. **Base/residue-level MSA viewer** for a separately uploaded, *already aligned*
   FASTA. Raw mapping blocks must not be mislabeled as an MSA.
3. **Flank / candidate-homolog synteny tracks** only when physical coordinates
   are available.
4. **CIGAR insertion/deletion operation counts** from minimap2, labeled as alignment gaps, not ancestral mutations.
5. **Six-frame complete-ORF tracks** from input genomic DNA, clearly marked
   as candidate coding potential.
6. **Independent hypothesis evidence and uncertainty**: status, completeness,
   unassessed tests and processing failures, never forced 100%-sum scores.
7. **Newick phylogeny visualization** with original node labels/branch lengths
   on a genuinely supplied tree. EvoGene currently **does not perform ML
   tree inference or calculate UFBoot/SH-aLRT**; these must come from
   verified aligned homologs and an external phylogenetic workflow.

The UI refuses to invent absent ORFs, ancestral sequences or trees.

## Running

Python 3.11+ and the external tools **NCBI Datasets v2** and **minimap2**
are required for the automatic analysis; the app checks their availability.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run auto_app.py
```

Alternative manual workflow (any number of manually uploaded relative FASTA files, with individually assigned GFF3 where available):

```bash
python -m streamlit run app.py
python -m unittest discover -s tests -v
```

CLI: `python -m evogene.auto_cli --help`, `python -m evogene.pipeline --help`.

**Data provenance:** all downloaded species assemblies are reported with NCBI
accessions. Reproduction of Peng & Zhao (2024) requires their published genome
alignments, precise assemblies/annotations and methods; the uploaded supplementary
spreadsheet contains phenotypic/protein descriptors, **not** enough ancestral
sequences to reconstruct origin.

Case: attached "СИНТЕЗ — Кейс 2", section 3–8. Article:
https://www.nature.com/articles/s41467-024-45028-1
Original author data: https://doi.org/10.6084/m9.figshare.19989395
Code: https://github.com/LiZhaoLab/DrosophilaDenovoGene

## Validation requirements / limitations

Automated gene-tree inference, support values, careful ancestral-ORF
reconstruction, robust sensitivity searches for deep protein homologs, formal
duplication/loss reconciliation and external biological benchmarking **remain
open work**. Synthetic unit tests test software behavior, not biological accuracy.
Review `ANALYSIS_PIPELINE.md`, `MECHANISM_V3.md` and tests for detailed scope.

For a hackathon demonstration, show at least **three independently documented
genes or controlled scenarios** with input accessions, alignments, limitations,
alternative hypotheses and complete analysis traces. The bundled
`data/demo_evidence.csv` is **synthetic** and must not be presented as
scientific validation.

In the local-genome UI, each uploaded relative FASTA gets an explicit species
name and an optional GFF3 selection; there is no fixed two-genome restriction.
NCBI automatic species are limited to the study's **20-name panel**, not to a
hard count of five. Downloading many assemblies is resource-intensive and may
leave more missing annotations; select appropriate outgroups deliberately.
