# EvoGene Auto-NCBI (third input mode)

**Only one uploaded file:** target genome FASTA, or a single-gene FASTA. 
The application downloads a **curated panel** of relative genomes from NCBI Datasets,
verifies each assembly's species and accession, caches genomes and available GFF3,
maps the DNA query with minimap2, and produces an explainable JSON report.

For whole-genome FASTA, specify the chosen **gene DNA** as pasted text or
enter contig, 1-based start/end and strand. FASTA does not contain gene IDs;
the ID input only labels the report. With a single-gene FASTA (one record
under 100 kb), no annotation or coordinates are required.

The automatic relative-species panel currently supports D. melanogaster,
D. simulans and D. yakuba; it is **not** a reconstructed species tree.
There is no scientific claim of orthology, duplication or de novo origin
from a simple DNA hit.

## Mac

    brew install ncbi-datasets-cli minimap2
    pip install -r requirements.txt
    streamlit run auto_app.py

Open http://localhost:8501 and upload FASTA.

CLI example (gene FASTA only):

    python -m evogene.auto_cli --target-fasta gene.fasta --gene-id FBgn0040343

CLI example (whole assembly; replace illustrative coordinates):

    python -m evogene.auto_cli --target-fasta dmel_genome.fna.gz --gene-id MY_GENE --chrom 2L --start 100000 --end 100600

Downloads are cached in genomes/ncbi_cache/ and *not* committed.
Reports are saved to results/auto_ncbi from CLI. Actual downloads must run
on a machine with internet access, NCBI Datasets CLI and minimap2.

Biological limitations: DNA similarity and CDS overlap are insufficient
to prove orthology, ancestral coding state or a de novo origin.
