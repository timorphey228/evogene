"""Reproducible DNA-locus phylogeny from *observed* comparative mappings.

Preferred: MAFFT multiple alignment -> IQ-TREE2 ML + UFBoot/SH-aLRT.
If IQ-TREE is absent but MAFFT works: NJ from MSA with site bootstrap.
If both absent: alignment-free 7-mer NJ *exploratory*, never a validated
gene-family phylogeny; no branch support is fabricated.

Unrooted topology; sequences are mapped genomic locus segments, not verified
orthologous CDS. This distinction is explicit in every result.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import math
import random
import shutil
import subprocess

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .genomes import fasta
from .tools import find_tool

MAX_SEQUENCES = 21
MAX_DNA_LENGTH = 12000
MIN_QUERY_COVERAGE = 0.65
MIN_IDENTITY = 0.60


def _clock():
    return datetime.now(timezone.utc).isoformat()


def safe_label(value):
    return "".join(c if c.isalnum() or c in "_-." else "_" for c in str(value))[:85]


def _reverse_complement(sequence):
    return str(Seq(sequence).reverse_complement())


def collect_locus_sequences(query, report, *, min_coverage=MIN_QUERY_COVERAGE,
                            min_identity=MIN_IDENTITY):
    """Extract one *candidate* locus per comparison species from genome FASTAs."""
    if not query or len(query) > MAX_DNA_LENGTH:
        return {}, [{"species": "focal", "state": "skipped",
                     "reason": "DNA query missing or exceeds 12 kb for tree inference"}]
    label = safe_label(report.get("target_species", "focal"))
    sequences = {label: query.upper()}
    details = [{"species": report.get("target_species", "focal"),
                "state": "focal", "length": len(query)}]
    assemblies = {entry["species"]: entry for entry in report.get("ncbi_assemblies", [])}
    for species in report.get("requested_relative_species", list(assemblies)):
        hits = report.get("relatives", {}).get(species, [])
        good = [hit for hit in hits if
                hit.get("query_coverage", 0) >= min_coverage and
                hit.get("identity_approx", 0) >= min_identity and
                hit.get("end", 0) > hit.get("start", 0)]
        good.sort(key=lambda h: (h.get("query_coverage", 0) *
                                 h.get("identity_approx", 0), h.get("mapq", 0)),
                  reverse=True)
        if not good:
            details.append({"species": species, "state": "excluded",
                            "reason": "No candidate locus passed coverage/identity thresholds"})
            continue
        best = good[0]
        if len(good) > 1 and (best.get("mapq") or 0) == 0:
            details.append({"species": species, "state": "excluded",
                            "reason": "Multiple candidate loci with ambiguous primary mapping"})
            continue
        if best["end"] - best["start"] > MAX_DNA_LENGTH:
            details.append({"species": species, "state": "excluded",
                            "reason": "Matched genomic interval exceeds 12 kb"})
            continue
        assembly = assemblies.get(species)
        if not assembly:
            details.append({"species": species, "state": "excluded",
                            "reason": "No verified assembly FASTA"})
            continue
        seq = None
        try:
            for contig, dna in fasta(assembly["fasta"]):
                if contig == best["chrom"]:
                    a, b = best["start"], best["end"]
                    if a >= 0 and b <= len(dna):
                        seq = dna[a:b]
                    break
        except (OSError, UnicodeError, ValueError) as exc:
            details.append({"species": species, "state": "error", "reason": str(exc)})
            continue
        if seq is None or len(seq) < 30 or not set(seq) <= set("ACGTN"):
            details.append({"species": species, "state": "excluded",
                            "reason": "Matched reference segment missing/invalid"})
            continue
        if best["strand"] == "-":
            seq = _reverse_complement(seq)
        if seq.count("N") / len(seq) > 0.10:
            details.append({"species": species, "state": "excluded",
                            "reason": "More than 10% undetermined bases"})
            continue
        clean = safe_label(species)
        if clean in sequences:
            clean += "_" + str(len(sequences))
        sequences[clean] = seq
        context = (report.get("synteny_evidence") or {}).get(species) or {}
        details.append({
            "species": species, "state": "selected", "label": clean,
            "length": len(seq), "identity_approx": best.get("identity_approx"),
            "coverage": best.get("query_coverage"), "mapq": best.get("mapq"),
            "chrom": best["chrom"], "start_0based": best["start"],
            "end_0based": best["end"], "strand": best["strand"],
            "cds_annotation": best.get("coding_overlap", "unknown"),
            "flank_synteny_status": context.get("status", "not_assessed"),
            "genome_accession": assembly.get("accession", "unknown"),
        })
        if len(sequences) >= MAX_SEQUENCES:
            break
    return sequences, details


def _to_fasta(sequences):
    return "".join(f">{name}\n" +
                   "\n".join(seq[i:i+80] for i in range(0, len(seq), 80)) +
                   "\n" for name, seq in sequences.items())


def _kmers(seq, size=7):
    seq = seq.upper()
    return {seq[i:i+size] for i in range(len(seq)-size+1)
            if "N" not in seq[i:i+size]}


def _kmer_distance(a, b, k=7):
    """7mer-Jaccard dissimilarity: screening distance, NOT substitutions/site."""
    x, y = _kmers(a, k), _kmers(b, k)
    if not x or not y:
        return 1.0
    return 1.0 - len(x & y) / len(x | y)


def _distance_matrix(sequences, mode="kmer", alignment=None):
    names = list(sequences)
    if mode == "msa":
        rows = {rec.id: str(rec.seq).upper() for rec in alignment}
        n = len(next(iter(rows.values())))
        distances = lambda a, b: (
            sum(rows[a][i] != rows[b][i] for i in range(n)) / n
            if n else 1.0
        )
    else:
        distances = lambda a, b: _kmer_distance(sequences[a], sequences[b])
    data = []
    for i, a in enumerate(names):
        data.append([max(0.0, min(1.0, distances(a, b))) if a != b else 0.0
                     for b in names[:i+1]])
    return DistanceMatrix(names, data)


def _splits(tree):
    tips = {t.name for t in tree.get_terminals()}
    result = set()
    for clade in tree.get_nonterminals():
        leaves = frozenset(t.name for t in clade.get_terminals())
        complement = frozenset(tips - leaves)
        if len(leaves) < 2 or len(complement) < 2:
            continue
        # Equivalent unrooted bipartitions share a canonical representation.
        result.add(min((leaves, complement), key=lambda x: tuple(sorted(x))))
    return result


def _bootstrap_support(tree, alignment, *, repeats=100, seed=2024):
    """Site bootstrap on real MSA; support is not posterior probability."""
    names = [rec.id for rec in alignment]
    rows = [str(rec.seq).upper() for rec in alignment]
    n = len(rows[0])
    if n < 30:
        return 0
    rand = random.Random(seed)
    counts = Counter()
    for _ in range(repeats):
        indices = [rand.randrange(n) for __ in range(n)]
        resampled = [
            SeqRecord(Seq("".join(row[i] for i in indices)), id=name)
            for name, row in zip(names, rows)
        ]
        dm = _distance_matrix(dict(zip(names, rows)), "msa",
                              MultipleSeqAlignment(resampled))
        rep_tree = DistanceTreeConstructor().nj(dm)
        counts.update(_splits(rep_tree))
    all_tips = {t.name for t in tree.get_terminals()}
    for clade in tree.get_nonterminals():
        leaves = frozenset(t.name for t in clade.get_terminals())
        other = frozenset(all_tips - leaves)
        if len(leaves) >= 2 and len(other) >= 2:
            split = min((leaves, other), key=lambda x: tuple(sorted(x)))
            clade.confidence = 100 * counts[split] / repeats
    return repeats


def infer_tree(sequences, *, try_external=True, bootstrap_replicates=100):
    """Return a real inferred tree or a diagnostic reason for impossibility."""
    if len(sequences) < 3:
        return {
            "status": "not_inferred", "newick": None,
            "reason": "Need at least three actual sequences (focal + two usable comparisons)",
            "method": None, "support_type": None,
        }
    if len(sequences) > MAX_SEQUENCES or any(
            len(x) < 30 or len(x) > MAX_DNA_LENGTH for x in sequences.values()):
        return {"status": "not_inferred", "newick": None,
                "reason": "Sequence count or length exceeds validated safety bounds",
                "method": None, "support_type": None}
    warnings = [
        "Unrooted tree; branch direction and ancestor are not established.",
        "Input is genomic DNA; locus similarities are not validated orthology.",
        "Unannotated or duplicated loci may distort topology; no evolutionary mechanism is proven.",
    ]
    mafft = find_tool("mafft") if try_external else None
    iqtree = (find_tool("iqtree2") or find_tool("iqtree")) if try_external else None
    aligned = None
    aligned_text = None
    diagnostics = []
    if mafft:
        try:
            with TemporaryDirectory(prefix="evogene_mafft_") as tmp:
                source = Path(tmp) / "input.fna"
                source.write_text(_to_fasta(sequences), encoding="utf8")
                command = [mafft, "--auto", "--thread", "2", str(source)]
                cp = subprocess.run(command, capture_output=True, text=True,
                                    timeout=240, check=False)
                diagnostics.append({"tool": "MAFFT", "command": command,
                                    "return_code": cp.returncode,
                                    "stderr_tail": cp.stderr[-350:]})
                if cp.returncode == 0:
                    candidate = AlignIO.read(StringIO(cp.stdout), "fasta")
                    if len(candidate) == len(sequences) and len(candidate[0]) <= 50000:
                        aligned, aligned_text = candidate, cp.stdout
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            diagnostics.append({"tool": "MAFFT", "error": str(exc)})
    if aligned is not None and iqtree:
        try:
            with TemporaryDirectory(prefix="evogene_iqtree_") as tmp:
                source = Path(tmp) / "alignment.fna"
                source.write_text(aligned_text, encoding="utf8")
                prefix = str(Path(tmp) / "gene_tree")
                command = [iqtree, "-s", str(source), "-m", "MFP",
                           "-B", "1000", "--alrt", "1000", "-T", "2",
                           "--seed", "2024", "--prefix", prefix, "--redo"]
                cp = subprocess.run(command, capture_output=True, text=True,
                                    timeout=1800, check=False)
                treefile = Path(prefix + ".treefile")
                diagnostics.append({"tool": "IQ-TREE", "command": command,
                                    "return_code": cp.returncode,
                                    "stderr_tail": cp.stderr[-700:]})
                if cp.returncode == 0 and treefile.is_file():
                    newick = treefile.read_text(encoding="utf8").strip()
                    tree = Phylo.read(StringIO(newick), "newick")
                    return {
                        "status": "inferred", "method": "maximum_likelihood_IQ-TREE2",
                        "newick": newick, "taxa": len(sequences),
                        "support_type": "UFBoot 1000 + SH-aLRT 1000 (IQ-TREE)",
                        "alignment": "MAFFT nucleotide MSA", "alignment_fasta": aligned_text,
                        "rooted": False, "warnings": warnings,
                        "diagnostics": diagnostics,
                        "description": "ML topology based on mapped nucleotide loci, not validated orthology.",
                    }
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            diagnostics.append({"tool": "IQ-TREE", "error": str(exc)})
    if aligned is not None:
        matrix = _distance_matrix(sequences, "msa", aligned)
        tree = DistanceTreeConstructor().nj(matrix)
        tree.rooted = False
        repeats = _bootstrap_support(tree, aligned, repeats=bootstrap_replicates)
        method = "NJ_from_MAFFT_nucleotide_MSA"
        unit = "fraction_variable_alignment_columns"
        support = "site_bootstrap_NJ" if repeats else None
        if not repeats:
            warnings.append("Not enough aligned columns to compute site-bootstrap support.")
    else:
        matrix = _distance_matrix(sequences)
        tree = DistanceTreeConstructor().nj(matrix)
        tree.rooted = False
        repeats = 0
        method = "EXPLORATORY_NJ_from_7mer_Jaccard"
        unit = "7mer_Jaccard_dissimilarity_not_substitutions_per_site"
        support = None
        warnings += [
            "No MAFFT multiple alignment was available.",
            "This is a preliminary alignment-free clustering, NOT a validated gene phylogeny.",
            "Node bootstrap supports are not available.",
        ]
    out = StringIO()
    Phylo.write(tree, out, "newick")
    return {
        "status": "inferred_exploratory" if aligned is None else "inferred",
        "method": method, "newick": out.getvalue().strip(),
        "taxa": len(sequences), "rooted": False, "branch_unit": unit,
        "support_type": support, "bootstrap_replicates": repeats,
        "alignment": "MAFFT nucleotide MSA" if aligned is not None else None,
        "alignment_fasta": aligned_text, "warnings": warnings,
        "diagnostics": diagnostics,
        "description": ("Measured DNA similarity tree; ancestry and evolutionary "
                        "mechanisms require independent synteny/orthology checks."),
    }


def infer_report_tree(query, report, *, try_external=True):
    sequences, records = collect_locus_sequences(query, report)
    result = infer_tree(sequences, try_external=try_external)
    result["selection"] = records
    result["query_length_nt"] = len(query)
    result["generated_at_utc"] = _clock()
    result["taxa_requested"] = len(report.get("requested_relative_species", [])) + 1
    return result
