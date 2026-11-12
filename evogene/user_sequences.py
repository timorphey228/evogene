"""Direct user scenario: focal-gene FASTA + comparative FASTA from other species.

Input:
  focal:       >gene_ID|optional_focal_species
               ACGT...
  comparisons: >gene_ID|species_code|optional_record_ID
               ACGT...

These are DNA sequences, *not* genome assemblies. Pairwise similarity cannot
establish orthology, synteny, ancestral coding states or gene origin.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
from io import StringIO
import re

from Bio import SeqIO

from .case_data import compare
from .scientific_report import SCHEMA, make_case_gene_result


def _read_fasta(contents, kind):
    if isinstance(contents, bytes):
        contents = contents.decode("utf-8-sig")
    if not isinstance(contents, str):
        raise ValueError("FASTA must be text or UTF-8 bytes")
    records = list(SeqIO.parse(StringIO(contents), "fasta"))
    if not records:
        raise ValueError(f"{kind}: FASTA is empty")
    out = []
    for record in records:
        ident = record.id.strip()
        parts = ident.split("|")
        if not parts or not re.fullmatch(r"[A-Za-z0-9_.-]{1,120}", parts[0]):
            raise ValueError(f"{kind}: invalid gene ID {ident!r}")
        if kind == "comparison" and (len(parts) < 2 or not parts[1]):
            raise ValueError(f"Comparison FASTA header must be >GENE|SPECIES: {ident!r}")
        if len(parts) >= 2 and not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", parts[1]):
            raise ValueError(f"{kind}: invalid species code in {ident!r}")
        sequence = str(record.seq).upper().replace("U", "T")
        if len(sequence) < 3 or any(c not in "ACGTN" for c in sequence):
            raise ValueError(f"{kind}: only unaligned A/C/G/T/N DNA is accepted ({ident!r})")
        out.append({"gene": parts[0],
                    "species": parts[1] if len(parts) >= 2 else "not_specified",
                    "header": ident, "sequence": sequence})
    return out


def analyze_user_fastas(focal_data, relative_data,
                        focal_filename="focal.fasta",
                        relative_filename="comparisons.fasta"):
    """Independent per-gene evidence with explicit 'insufficient_data' conclusion."""
    targets = _read_fasta(focal_data, "focal")
    relatives = _read_fasta(relative_data, "comparison")
    ids = [x["gene"] for x in targets]
    if len(set(ids)) != len(ids):
        raise ValueError("Only one focal sequence per gene ID is allowed")
    by_gene = defaultdict(list)
    for item in relatives:
        by_gene[item["gene"]].append(item)
    focal_bytes = focal_data.encode() if isinstance(focal_data, str) else focal_data
    rel_bytes = relative_data.encode() if isinstance(relative_data, str) else relative_data
    provenance = {
        "dataset": "User-provided candidate DNA and comparative DNA FASTA",
        "filename": f"{focal_filename} + {relative_filename}",
        "focal_sha256": sha256(focal_bytes).hexdigest(),
        "comparison_sha256": sha256(rel_bytes).hexdigest(),
        "sha256": sha256(focal_bytes + b"\0" + rel_bytes).hexdigest(),
        "method": "local_pairwise_dna",
        "software": {"analysis": SCHEMA, "biopython": version("biopython")},
        "format": "focal >GENE|SPECIES ; comparison >GENE|SPECIES|HIT",
    }
    raw_reports, results = {}, []
    errors = []
    for target in targets:
        gene_id = target["gene"]
        dna = target["sequence"]
        related = by_gene.get(gene_id, [])
        measurements = []
        for item in related:
            if target["species"] != "not_specified" and item["species"] == target["species"]:
                errors.append({"gene_id": gene_id, "status": "warning",
                               "reason": "Comparison record has focal species label: " + item["header"]})
            metrics = compare(dna, item["sequence"])
            measurements.append({
                "species_abbr": item["species"],
                "header": item["header"],
                "header_kind": "user_supplied_comparison",
                "length_nt": len(item["sequence"]),
                **metrics,
            })
        described_species = sorted({h["species_abbr"] for h in measurements})
        passing = sorted({
            h["species_abbr"] for h in measurements
            if h.get("comparison_status") == "aligned"
            and (h.get("identity_aligned") or 0) >= .65
            and (h.get("query_coverage") or 0) >= .6
        })
        raw = {
            "gene_id": gene_id,
            "input_type": "user_supplied_fasta_sequences",
            "focal_species": target["species"],
            "method": "local_pairwise_dna",
            "target_length_nt": len(dna),
            "hypothesis": "unresolved",
            "explanation": ("Найденные сходные последовательности — только "
                            "кандидаты в гомологи. Эволюционный механизм не доказан."),
            "hits": measurements,
            "summary": {
                "record_count": len(measurements),
                "all_species_in_gene": described_species,
                "annotated_hit_species": [],
                "genomic_coordinate_hit_species": [],
                "observed_dna_similarity_species": passing,
                "comparison_skipped": sum(
                    h.get("comparison_status") != "aligned" for h in measurements),
            },
            "method_comparison": [],
            "requested_species": described_species,
            "phylogeny": {
                "status": "not_assessed",
                "reason": "No verified orthologous MSA and ML tree; user supplied DNA only.",
            },
            "limitations": [
                "Unannotated DNA sequences do not demonstrate protein homology.",
                "No chromosome coordinates, genomic assembly accession or GFF3 were supplied.",
                "No syntenic ancestry, ancient coding status, duplication or loss was tested.",
                "Negative local alignment is not evidence of gene absence.",
            ],
        }
        raw_reports[gene_id] = raw
        results.append(make_case_gene_result(raw, provenance))
    unmatched = [item["header"] for item in relatives if item["gene"] not in set(ids)]
    summary = {
        "requested": len(targets), "completed": len(results), "failed": 0,
        "unmatched_comparison_headers": unmatched,
        "origin_categories": dict(Counter(x["conclusion"]["category"] for x in results)),
        "interpretation": "Similarity is measured; evolutionary mechanisms remain unconfirmed.",
    }
    return {
        "schema": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "requested_gene_ids": ids, "source": provenance,
        "gene_count": len(results), "results": results, "raw_reports": raw_reports,
        "errors": errors, "summary": summary,
    }
