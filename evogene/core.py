"""Transparent hypothesis ranking from comparative-genomics observations.

This is NOT an experimentally validated gene-origin classifier. In particular,
protein properties and labels from Peng & Zhao do not prove evolutionary origin.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from typing import Any

CATEGORIES = ('de_novo', 'duplication', 'ancient_divergence', 'other_mechanism', 'undetermined')


def truth(value: Any) -> bool | None:
    if value is None or str(value).strip().lower() in ('', 'na', 'nan', 'none', 'unknown', '?'):
        return None
    if str(value).strip().lower() in ('yes', 'true', '1', 'y', 'да'):
        return True
    if str(value).strip().lower() in ('no', 'false', '0', 'n', 'нет'):
        return False
    raise ValueError(f'Expected yes/no/unknown, got {value!r}')


@dataclass
class Evidence:
    gene_id: str
    species: str = 'Drosophila melanogaster'
    # Observations are tri-state; missing evidence is not negative evidence.
    paralog_detected: bool | None = None
    distant_coding_homolog: bool | None = None
    syntenic_outgroup_locus: bool | None = None
    outgroup_locus_noncoding: bool | None = None
    homolog_search_completed: bool | None = None
    rearrangement_evidence: bool | None = None
    outgroup_count: int = 0
    alignment_identity: float | None = None
    notes: str = ''

    @classmethod
    def from_mapping(cls, row: dict[str, Any]) -> 'Evidence':
        fields = ('paralog_detected', 'distant_coding_homolog', 'syntenic_outgroup_locus',
                  'outgroup_locus_noncoding', 'homolog_search_completed', 'rearrangement_evidence')
        kwargs = {k: truth(row.get(k)) for k in fields}
        count = row.get('outgroup_count', 0)
        kwargs['outgroup_count'] = int(float(count)) if str(count).strip().lower() not in ('', 'none', 'nan', 'na') else 0
        identity = row.get('alignment_identity')
        kwargs['alignment_identity'] = float(identity) if identity is not None and str(identity).strip().lower() not in ('', 'none', 'nan', 'na') else None
        kwargs['gene_id'] = str(row.get('gene_id') or row.get('FBID') or '').strip()
        kwargs['species'] = str(row.get('species') or 'Drosophila melanogaster')
        kwargs['notes'] = str(row.get('notes') or '')
        if not kwargs['gene_id']:
            raise ValueError('gene_id is required')
        if kwargs['outgroup_count'] < 0:
            raise ValueError('outgroup_count must be nonnegative')
        if kwargs['alignment_identity'] is not None and not 0 <= kwargs['alignment_identity'] <= 1:
            raise ValueError('alignment_identity must be between 0 and 1')
        return cls(**kwargs)


@dataclass
class Result:
    gene_id: str
    hypothesis: str
    evidence_strength: str
    scores: dict[str, float]
    supporting: list[str] = field(default_factory=list)
    conflicting: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    alternatives: list[str] = field(default_factory=list)
    research_label: str | None = None
    research_age: str | None = None
    properties: dict[str, Any] = field(default_factory=dict)
    provenance: str = 'Rule-based prototype; not calibrated probabilities'

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def classify(e: Evidence, properties: dict[str, Any] | None = None) -> Result:
    """Prioritize positive evidence, report contradictions; never infer from missing rows.

    Values are qualitative scores for display, *not* posterior probabilities.
    Evidence fields should be assembled from independent alignment/search outputs.
    """
    scores = {name: 0.0 for name in CATEGORIES if name != 'undetermined'}
    support: list[str] = []
    conflict: list[str] = []
    missing: list[str] = []
    p = properties or {}

    if e.paralog_detected is True:
        scores['duplication'] += 4
        support.append('Detected a candidate paralog (verify reciprocal alignment and gene tree).')
    elif e.paralog_detected is None:
        missing.append('Paralog search (BLAST/DIAMOND, including within-species copies)')

    if e.distant_coding_homolog is True:
        scores['ancient_divergence'] += 4
        support.append('Coding homolog in a distant lineage supports older origin.')
    elif e.distant_coding_homolog is None:
        missing.append('Distant protein/DNA homology search')

    if e.rearrangement_evidence is True:
        scores['other_mechanism'] += 4
        support.append('Genomic rearrangement/fusion/transposition evidence reported.')
    elif e.rearrangement_evidence is None:
        missing.append('Structural variant/rearrangement analysis')

    denovo_conditions = (e.syntenic_outgroup_locus is True,
                         e.outgroup_locus_noncoding is True,
                         e.homolog_search_completed is True,
                         e.outgroup_count >= 2)
    if all(denovo_conditions):
        scores['de_novo'] += 4
        support.append('Syntenic noncoding ancestral locus in >=2 outgroups and a completed homolog search.')
    else:
        if e.syntenic_outgroup_locus is not True:
            missing.append('Confirmed orthologous genomic locus (synteny) in outgroup')
        if e.outgroup_locus_noncoding is not True:
            missing.append('Evidence that the aligned ancestral locus was noncoding (ORF/annotation)')
        if e.homolog_search_completed is not True:
            missing.append('Completed search for homologs outside the syntenic locus')
        if e.outgroup_count < 2:
            missing.append('At least two independently sampled outgroup species')

    # Strong contradicting evidence blocks a confident de novo call.
    if e.paralog_detected is True or e.distant_coding_homolog is True:
        if scores['de_novo']:
            conflict.append('De novo evidence conflicts with observed coding homolog/paralog.')
        scores['de_novo'] = 0
    if e.outgroup_locus_noncoding is False and e.syntenic_outgroup_locus is True:
        scores['ancient_divergence'] += 2
        support.append('Outgroup syntenic locus is annotated as coding.')

    winners = sorted(scores, key=lambda k: (-scores[k], k))
    best = winners[0]
    if scores[best] <= 0 or (len(winners) > 1 and scores[winners[0]] == scores[winners[1]]):
        best = 'undetermined'
    if best == 'undetermined':
        strength = 'insufficient_or_conflicting'
    elif conflict or (sum(v > 0 for v in scores.values()) > 1):
        strength = 'conflicting'
    elif best == 'de_novo' and e.outgroup_count >= 3:
        strength = 'moderate'
    elif best in ('duplication', 'ancient_divergence', 'other_mechanism'):
        strength = 'preliminary'
    else:
        strength = 'preliminary'

    label = str(p.get('De novo or not?', '')).strip() or None
    age = str(p.get('de novo gene age', '')).strip() or None
    return Result(
        gene_id=e.gene_id, hypothesis=best, evidence_strength=strength,
        scores=scores, supporting=support, conflicting=conflict, missing=missing,
        alternatives=[w for w in winners if w != best and scores[w] > 0],
        research_label=label, research_age=age,
        properties={k: v for k, v in p.items() if k not in ('De novo or not?', 'de novo gene age')},
    )
