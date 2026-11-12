# EvoGene · Case 2 (ITMO «СИНТЕЗ», 2026)

An **evidence-first comparative-genomics prototype** for reconstructing the
possible origin of protein-coding genes in *Drosophila* and relatives.

> **Scientific scope:** screening, not a validated automated origin classifier.
> **Missing BLAST/minimap2 hits and absent GFF3 CDS annotation never prove
> de novo origin or gene loss.** No made-up posterior probabilities or trees.

## What is implemented

- Upload a target gene as DNA FASTA or supply a sequence/coordinates in a target
  assembly; use available NCBI Datasets assemblies and GFF3 for any chosen subset of up to 20 comparison species (the focal species plus the 20 others from the 21-species panel, including Musca as an outgroup).
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
uses a **selectable 21-species panel** (20 comparative species plus the focal species; Musca is an additional outgroup). `app.py` studies all 20 codes already present in the
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
7. **Automatic DNA-locus phylogeny + Newick** (new single-workflow dashboard).
   MAFFT + IQ-TREE2 computes an unrooted ML tree, UFBoot and SH-aLRT
   (1000 replicates); MAFFT alone enables an NJ tree with site bootstrap.
   Without these tools the app can only offer explicitly marked
   **exploratory unaligned 7-mer NJ clustering without branch support**.
   A DNA locus tree is NOT by itself a proven ortholog/gene-family tree,
   and absent reliable matches means no invented tree. Existing manual
   Newick upload remains supported.

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

Validated gene-family phylogeny for every candidate, calibrated Cactus/Genewise
ancestral-coding classification, robust sensitivity searches for deep protein homologs, formal
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


## User scenario: one or more candidate genes plus comparative sequences

There are now **three distinct supported input scenarios**, not one ambiguous
"FASTA" input:

1. **Direct focal + comparative DNA FASTAs:** launch `app.py` and select
   `Кандидатные гены + сравнительные FASTA`. In the focal file, use headers
   `>GENE_ID|FOCAL_SPECIES`, one target per ID. In the second file, use
   `>GENE_ID|COMPARISON_SPECIES|OPTIONAL_LOCUS`, with any number of
   comparisons per gene. The UI analyzes all candidate IDs automatically.
2. **Full assemblies:** upload focal FASTA/GFF3 and relative genomes in
   `app.py`, or use NCBI auto download via `auto_app.py`.
3. **Peng & Zhao comparative extract:** upload the three-section
   `supplementary_file_unannotated_alignments.fasta`; select any set of
   `FBgn...` IDs and comparative species. Its species data are already
   PRESELECTED by the paper, not an unbiased whole-genome search.

**New per-gene report (six case criteria):** input provenance, independent
sequence comparisons, explicit `insufficient_data` or preliminary evidence,
concrete observed evidence, stage trace, SHA-256 input checksums,
reproducibility, and separation of observations from evolutionary mechanisms.
Batch conclusions export as `JSON` and human-readable `Markdown` without
claiming unperformed ancestral reconstruction or phylogenetic supports.

CLI for one or many paper-derived genes:

```bash
python -m evogene.case_batch_cli \
  --input supplementary_file_unannotated_alignments.fasta \
  --gene FBgn0040343 --gene FBgn0262810 --gene FBgn0025616 \
  --method raw --output results/case2_multi_gene.json
```

For published-file multi-gene batches, the original 54MB source is read just
three times, then selected gene subsets are analyzed locally. This avoids
repeatedly scanning the entire file for each gene.

A matching coordinate or a local DNA alignment is not sufficient to determine
whether a gene arose by duplication, de novo emergence, or ancient divergence.
In such cases the correct reported outcome is **insufficient data**, including
which specific investigations are missing.


## NEW: one-input multi-gene comparative-genomics workflow

Use `streamlit run auto_app.py` (the canonical, redesigned interface):

1. Upload **FASTA of one or more candidate genes**, e.g.
   `>FBgn0040343` followed by valid A/C/G/T/N DNA (30–100000 nt).
2. Choose the **focal species** and any relevant **comparison species**
   from the 21-species catalogue. The published assembly accession is
   pinned when available; any NCBI fallback is recorded transparently.
   FlyBase r6.36 for D. melanogaster is NOT an NCBI RefSeq accession;
   the app warns if the actual assembly differs from the study's version.
3. Optional: add a **TE consensus library FASTA** to test candidate
   transposable element sequence similarities. This does not replace
   RepeatMasker/Dfam or whole-genome repeat annotation.
4. Run the project. The workflow calls the existing NCBI Datasets +
   minimap2 + GFF3 screening steps for each gene. It records additional
   mapped loci, preliminary protein matches, approximate flank-based
   synteny, candidate uninterrupted ORFs inside aligned syntenic
   fragments, and independent mechanism tests.
5. The pipeline extracts up to 21 observed nucleotide loci (one candidate
   hit per comparison species). It builds **an observed, unrooted tree**
   by MAFFT/IQ-TREE ML when installed. With MAFFT but no IQ-TREE, it
   computes NJ with 100 site-bootstrap replicates. With neither, it
   produces a clearly flagged exploratory k-mer NJ tree without
   branch support. For fewer than 3 suitable sequences **no tree is
   fabricated**. Even ML topology based on genomic loci is **not**
   verified gene-family orthology.
6. Inspect **per-gene hypotheses and argumentation**, candidate TE
   hits if searched, DNA alignment/ORF/synteny charts, UTC timestamped
   events, per-species assembly accessions, software versions,
   errors and limitations.
7. Download the entire project ZIP containing **JSON report, JSONL audit,
   genuine Newick and MSA where calculated**. Persistent outputs are saved
   in `results/evolutionary_pipeline/<run_id>/`.

External dependencies for the automatic workflow:
**NCBI datasets v2, minimap2** (required),
**MAFFT and IQ-TREE2** (recommended, needed for robust tree reconstruction).
Optional repeat library should be supplied for TE screening.
Python requirements: `pip install -r requirements.txt`.

### Scientific scope: paper analogue, not exact replication

Peng & Zhao (Nature Communications 2024) used **Progressive Cactus,
halLiftover, Genewise/Spaln calibrated coding predictions, all-vs-all
BLASTp and jackhmmer/TBLASTN** to evaluate ancestral noncoding syntenic
loci, distant homologs and candidate de novo age. The automatic EvoGene
screen presently uses minimap2/flank-based synteny instead of HAL and
a preliminary protein search, so it cannot independently establish
definitive *de novo* origin, historical loss, phylogenetic duplication
or the exact branch of origination. Its hypotheses are **diagnostic,
not probabilities**. Run the published author pipeline and/or additional
validated tools before making biological claims.
