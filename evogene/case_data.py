"""Parser for Peng & Zhao supplied three-section nonstandard FASTA data.

Sections are raw Cactus hits, Spaln spliced predictions, and Genewise
spliced predictions. These are PRESELECTED candidate regions; none is
automatically a verified ortholog, ancestral noncoding locus, or MSA.
"""
from __future__ import annotations

from collections import Counter
from statistics import median
from pathlib import Path
import re
from Bio.Align import PairwiseAligner

METHODS = {
    "raw": ("raw sequence for ", "Raw Cactus hits"),
    "spaln": ("spaln for ", "Spaln predictions"),
    "genewise": ("genewise for ", "Genewise predictions"),
}


def index_file(path):
    """Stream 50+MB once; no sequence retention. Discover every gene and species."""
    genes = set()
    species = set()
    record_counts = Counter()
    gene_counts = Counter()
    method = None
    with open(path, encoding="utf-8-sig") as stream:
        for line in stream:
            line = line.strip()
            if line.startswith("Section "):
                method = None
                continue
            matched = next(
                ((key, line[len(prefix):]) for key, (prefix, _) in METHODS.items()
                 if line.startswith(prefix)), None
            )
            if matched:
                method, gene = matched
                genes.add(gene.split("_", 1)[-1])
                gene_counts[method] += 1
            elif line.startswith(">") and method:
                record_counts[method] += 1
                species.add(line[1:].split("|", 1)[0])
    return {
        "genes": sorted(genes), "species": sorted(species),
        "gene_counts": dict(gene_counts),
        "record_counts": dict(record_counts),
        "methods": list(METHODS),
        "format": "three_section_preacted_cactus_spaln_genewise_not_msa",
    }


def sections(path, method="raw", gene_id=None):
    """Yield (dmel_FBgn..., [(header, sequence), ...]) for one valid section.

    Each next method or gene marker flushes prior records. This prevents the
    original bug: 'spaln for' or 'genewise for' being appended as DNA text.
    """
    if method not in METHODS:
        raise ValueError(f"Unknown method: {method}")
    active_method = None
    gene = None
    records = []
    header = None
    sequence = []
    collecting = False

    def add_record():
        nonlocal header, sequence
        if collecting and header is not None:
            records.append((header, "".join(sequence).upper()))
        header = None
        sequence = []

    def is_wanted(value):
        return gene_id is None or value == gene_id or value.endswith("_" + gene_id)

    with open(path, encoding="utf-8-sig") as stream:
        for line in stream:
            item = line.strip()
            if item.startswith("Section "):
                add_record()
                if collecting and gene is not None:
                    yield gene, records
                active_method = None
                gene = None
                records = []
                collecting = False
                continue
            marker = next(
                ((key, item[len(prefix):]) for key, (prefix, _) in METHODS.items()
                 if item.startswith(prefix)), None
            )
            if marker:
                add_record()
                if collecting and gene is not None:
                    yield gene, records
                active_method, gene = marker
                records = []
                collecting = active_method == method and is_wanted(gene)
                continue
            if not collecting:
                continue
            if item.startswith(">"):
                add_record()
                header = item[1:]
            elif item and header is not None:
                sequence.append(item)
        add_record()
    if collecting and gene is not None:
        yield gene, records


def header_info(header):
    fields = header.split("|")
    species = fields[0] if fields else "unknown"
    locus = fields[1] if len(fields) > 1 else ""
    kind = ("genomic_coordinates" if re.search(r"[+-]:\d+:\d+$", locus)
            else "gene_identifier")
    score = next((v.removeprefix("score_") for v in fields[2:]
                  if v.startswith("score_")), None)
    shift = next((v.removeprefix("shift_") for v in fields[2:]
                  if v.startswith("shift_")), None)
    return {"species_abbr": species, "locus": locus, "header_kind": kind,
            "prediction_score_raw": score, "shift_flag": shift}


def compare(a, b, max_cells=2_000_000):
    """Local sequence identity only; bounded to avoid huge quadratic alignment."""
    if not a or not b:
        return {"comparison_status": "empty_sequence"}
    if len(a) * len(b) > max_cells:
        return {"comparison_status": "skipped_quadratic_cost",
                "max_cells": max_cells}
    aligner = PairwiseAligner(mode="local", match_score=2, mismatch_score=-2,
                              open_gap_score=-5, extend_gap_score=-0.5)
    try:
        al = aligner.align(a, b)[0]
    except (IndexError, ValueError):
        return {"comparison_status": "alignment_failed"}
    equal = total = 0
    for (a0, a1), (b0, b1) in zip(*al.aligned):
        x, y = a[a0:a1], b[b0:b1]
        total += min(len(x), len(y))
        equal += sum(p == q and p in "ACGT" for p, q in zip(x, y))
    return {
        "comparison_status": "aligned" if total else "no_local_alignment",
        "identity_aligned": round(equal / total, 4) if total else None,
        "query_coverage": round(sum(end - start for start, end in al.aligned[0]) / len(a), 4),
        "aligned_columns_ungapped": total,
    }


def method_summary(path, gene_id):
    """Compare annotations predicted by all three methods, no artificial labels."""
    result = []
    for key, (_, name) in METHODS.items():
        records = next((r for _, r in sections(path, key, gene_id)), [])
        seq_lengths = [len(seq) for _, seq in records if seq]
        genomic = [header_info(h) for h, _ in records]
        result.append({
            "method": key, "label": name, "records": len(records),
            "species": len(set(item["species_abbr"] for item in genomic)),
            "coordinate_hits": sum(item["header_kind"] == "genomic_coordinates" for item in genomic),
            "gene_identifier_hits": sum(item["header_kind"] == "gene_identifier" for item in genomic),
            "empty_sequences": len(records) - len(seq_lengths),
            "shift_yes": sum((item["shift_flag"] or "").lower() == "yes" for item in genomic),
            "shift_unknown": sum(item["shift_flag"] is None for item in genomic),
            "median_nt": int(median(seq_lengths)) if seq_lengths else None,
        })
    return result


def analyze(path, gene_id, method="raw", selected_species=None):
    """Investigate one gene and one method, with cross-method completeness."""
    if method not in METHODS:
        raise ValueError("Available methods: " + ", ".join(METHODS))
    entries = next(((name, rows) for name, rows in sections(path, method, gene_id)), None)
    if entries is None:
        raise ValueError("Gene not found in " + method + ": " + gene_id)
    section, records = entries
    targets = [(h, s) for h, s in records if h.startswith("dmel|" + gene_id + "|") and s]
    if not targets:
        raise ValueError("Reference gene dmel|" + gene_id + " missing")
    header, query = max(targets, key=lambda x: len(x[1]))
    hits = []
    allow = set(selected_species) if selected_species is not None else None
    for h, s in records:
        if h == header or not s:
            continue
        metadata = header_info(h)
        if allow is not None and metadata["species_abbr"] not in allow:
            continue
        hits.append({
            **metadata, "header": h, "length_nt": len(s),
            **compare(query, s),
        })
    annotated = sorted({
        h["species_abbr"] for h in hits
        if h["header_kind"] == "gene_identifier" and
        h.get("comparison_status") == "aligned" and
        h.get("query_coverage", 0) >= .6 and
        h.get("identity_aligned", 0) >= .65
    })
    genomic = sorted({
        h["species_abbr"] for h in hits
        if h["header_kind"] == "genomic_coordinates" and
        h.get("comparison_status") == "aligned" and
        h.get("query_coverage", 0) >= .6 and
        h.get("identity_aligned", 0) >= .65
    })
    category = ("preselected_annotated_homology_candidates" if annotated else
                "preselected_genomic_region_matches" if genomic else
                "unresolved_preselected_hits")
    return {
        "gene_id": gene_id,
        "method": method,
        "requested_species": sorted(allow) if allow is not None else None,
        "source_section": section,
        "input_type": "published_preselected_hits_not_whole_genome_or_msa",
        "target_length_nt": len(query),
        "hypothesis": category,
        "explanation": (
            "Локальное ДНК-сходство в отобранных авторами Cactus/Spaln/Genewise "
            "фрагментах. Гомология белка, ортология, некодирующий предок и "
            "эволюционный механизм не подтверждены этим сравнением."
        ),
        "summary": {
            "annotated_hit_species": annotated,
            "genomic_coordinate_hit_species": genomic,
            "record_count": len(hits), "all_species_in_gene": sorted(
                {header_info(h)["species_abbr"] for h, _ in records}
            ),
            "comparison_skipped": sum(h.get("comparison_status") != "aligned" for h in hits),
        },
        "method_comparison": method_summary(path, gene_id),
        "hits": hits,
        "limitations": [
            "The source contains three independent prediction sections; they cannot be concatenated.",
            "Preselected hits do not replace unbiased whole-genome homology search.",
            "Coordinate-labelled sequences are not automatically verified noncoding ancestors.",
            "The method's score/shift label is metadata, not recalculated statistical significance.",
            "Local DNA identity does not establish ancestral coding state or paralogy.",
            "No multiple alignment, gene tree or phylogenetic support is inferred here.",
        ],
        "phylogeny": {"status": "not_assessed",
                      "reason": "No validated homologous protein MSA and no ML tree were computed."},
    }


def export_gene_fasta(path, gene_id, method="raw", species=None):
    """Export selected *unaligned* sequences in standard FASTA syntax."""
    selected = set(species) if species is not None else None
    entries = next((r for _, r in sections(path, method, gene_id)), None)
    if entries is None:
        raise ValueError("Gene not found")
    lines = []
    for header, sequence in entries:
        if selected is not None and header.split("|", 1)[0] not in selected:
            continue
        if not sequence:
            continue
        lines.append(">" + header)
        lines.extend(sequence[i:i + 80] for i in range(0, len(sequence), 80))
    if not lines:
        raise ValueError("No sequences selected")
    return "\n".join(lines) + "\n"
