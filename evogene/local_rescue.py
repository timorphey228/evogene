"""Sensitive *bounded* local DNA rescue within already-established synteny.

Never search whole genomes with quadratic alignment. Restrict to a short
collinear flank-bounded interval and short queries. Reports tentative homology.
"""
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from .genomes import fasta, annotate


def rescue_in_interval(query, assembly, chrom, low, high, gff_path=None):
    size=high-low
    if not 60 <= len(query) <= 1800 or not 60 <= size <= 10000:
        return {'status':'skipped_size','query_nt':len(query),'window_nt':size}
    if len({query[i:i+5] for i in range(0,len(query)-4)}) < 15:
        return {'status':'skipped_low_complexity'}
    segment=None
    for name,sequence in fasta(assembly):
        if name==chrom:
            segment=sequence[low:high]
            break
    if not segment or len(segment)<60:
        return {'status':'skipped_missing_interval'}
    aligner=PairwiseAligner()
    aligner.mode='local'
    aligner.match_score=2
    aligner.mismatch_score=-1
    aligner.open_gap_score=-8
    aligner.extend_gap_score=-2
    best=None
    for strand,dna in (('+',query),('-',str(Seq(query).reverse_complement()))):
        alignments=aligner.align(dna,segment)
        if not alignments:continue
        alignment=alignments[0]
        paired=alignment.aligned
        aligned=0;matches=0
        for (a,b),(x,y) in zip(paired[0],paired[1]):
            left=dna[a:b];right=segment[x:y]
            aligned+=len(left)
            matches+=sum(p==q for p,q in zip(left,right))
        if aligned<60:continue
        coverage=aligned/len(dna)
        identity=matches/aligned
        if coverage<.45 or identity<.68:continue
        query_blocks=paired[0]
        ref_blocks=paired[1]
        start=low+int(ref_blocks[0][0])
        end=low+int(ref_blocks[-1][1])
        candidate={'chrom':chrom,'start':start,'end':end,'strand':strand,
            'query_coverage':round(coverage,4),'identity_approx':round(identity,4),
            'aligned_nt':aligned,'matching_nt':matches}
        if best is None or (matches,coverage)>(best['matching_nt'],best['query_coverage']):
            best=candidate
    if best is None:
        return {'status':'no_local_homology','query_nt':len(query),'window_nt':size}
    annotate([best],gff_path)
    return {'status':'candidate_local_homology','hit':best,'method':'bounded_local_pairwise_alignment',
            'warning':'Tentative local DNA homology; repeats and chance similarity require validation.'}
