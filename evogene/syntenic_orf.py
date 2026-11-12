"""Survey potential ORFs inside *observed* syntenic DNA matches.

This does not reconstruct ancestral ORFs: syntenic genomic sequences can
contain introns, and evolutionary polarization requires a species tree and
properly aligned orthologous sites.
"""
from __future__ import annotations

from .genomes import fasta
from .orf import find_orfs
from .phylogeny import _reverse_complement


def survey_syntenic_orfs(report, min_aa=25):
    assemblies = {r["species"]: r for r in report.get("ncbi_assemblies", [])}
    output={}
    for species, context in (report.get("synteny_evidence") or {}).items():
        if context.get("status")!="flanks_collinear":
            output[species]={"status":"not_assessed",
                             "reason":"Both flanking loci were not aligned in genomic order"}
            continue
        hits=context.get("syntenic_hits") or []
        if not hits:
            output[species]={"status":"not_assessed",
                             "reason":"No mapped candidate DNA segment in observed syntenic window"}
            continue
        best=max(hits,key=lambda h:h.get("end",0)-h.get("start",0))
        length=best.get("end",0)-best.get("start",0)
        if length < 75 or length>12000:
            output[species]={"status":"not_assessed",
                             "reason":"Syntenic candidate too short or >12kb for bounded ORF screen",
                             "length_nt":length}
            continue
        assembly=assemblies.get(species)
        if not assembly:
            output[species]={"status":"not_assessed","reason":"No comparable genome FASTA"}
            continue
        seq=None
        try:
            for chrom,contig in fasta(assembly["fasta"]):
                if chrom==best["chrom"]:
                    a,b=best["start"],best["end"]
                    if 0<=a<b<=len(contig):
                        seq=contig[a:b]
                    break
        except (OSError,UnicodeError,ValueError) as exc:
            output[species]={"status":"failed","reason":str(exc)}
            continue
        if seq is None or set(seq)-set("ACGTN"):
            output[species]={"status":"not_assessed",
                             "reason":"Candidate coordinates or bases unavailable"}
            continue
        if best.get("strand")=="-":
            seq=_reverse_complement(seq)
        try:
            frames=find_orfs(seq,min_aa=min_aa)
        except ValueError as exc:
            output[species]={"status":"failed","reason":str(exc)}
            continue
        output[species]={
            "status":"screened_genomic_fragment",
            "chrom":best["chrom"],"start_0based":best["start"],"end_0based":best["end"],
            "length_nt":len(seq),"candidate_complete_orf_count":len(frames),
            "longest_complete_orf_aa":frames[0].aa_length if frames else None,
            "longest_orfs":[item.to_dict() for item in frames[:5]],
            "coding_annotation_status":context.get("locus_status","unknown"),
            "interpretation":"Only potential uninterrupted ORFs inside this genomic hit; "
                             "not spliced CDS, ancestral state or demonstrated translation.",
        }
    return output
