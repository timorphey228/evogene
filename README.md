# EvoGene — evolutionary gene-origin screening

Upload a single-gene DNA FASTA, select *Drosophila* reference species and relatives, and run NCBI Datasets + minimap2. EvoGene downloads available genomes/GFF3 and produces **independent mechanistic evidence screens** for de novo, duplication, retroposition, relocation, and gene fusion, with conservation shown separately.

**Scientific status:** exploratory prototype, not a validated causal classifier. Missing hits and missing annotation are *unknown*, not evidence of absence. The report contains no invented origin probabilities. A preliminary candidate requires independent validation with protein homology, synteny, ancestral ORF and gene trees.

## Install and launch on macOS

```bash
cd ~/Desktop/evogene
git switch main
git pull --ff-only
python3 -m venv .venv  # skip if .venv already exists
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run auto_app.py
```

Install **NCBI Datasets v2** and **minimap2** separately. The UI checks their versions on startup and locates binaries in PATH, `~/bin`, `/opt/homebrew/bin`, `/usr/local/bin`. It disables analysis until both are available.

Use a *single gene* FASTA, not an unannotated whole-genome FASTA; choose its source species. Set 4–5 relatives for a richer comparison. A successful NCBI download does not ensure useful homology evidence. Protein-derived tests also require available annotations and detectable ORFs.

## Output

- **Mechanism evidence matrix:** status (`candidate`, `partial_signal`, `no_signal`, `not_assessed`), observed signals, caveats; **no forced percentages**.
- **Conservation:** separate list of relatives with aligned DNA, annotated CDS overlap, and syntenic matches.
- **Completeness:** which tests could not be carried out and why.
- **JSON:** automatically saved to `results/auto_ncbi/<gene>.json` and available via browser download.

```bash
python -m unittest discover -s tests -v
```

Biological benchmark validation using verified independently curated de novo genes, retrocopies, duplicated genes, and negative controls **is still pending**. Synthetic unit tests are not proof of empirical accuracy.

For details see `MECHANISM_V3.md`. Older modules and examples are retained for compatibility; `auto_app.py` and `evogene/ncbi_auto.py` are the current automatic pipeline.
