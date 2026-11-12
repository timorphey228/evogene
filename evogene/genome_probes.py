"""Genomic neighborhood screens: local sequence anchors, not validated orthology."""
from .genomes import fasta, map_genome

def overlapping(a,b):
    return a['chrom']==b['chrom'] and min(a['end'],b['end'])>max(a['start'],b['start'])

def discover_target_locus(query,genome,all_hits=None):
    # Reuse reference mapping for every mechanism: a single minimap2 call.
    candidates=all_hits if all_hits is not None else map_genome(query,genome,'target')
    matches=[x for x in candidates if
             x.get('query_coverage',0)>=.5 and x.get('identity_approx',0)>=.7]
    if not matches:return None,[]
    main=matches[0]; extras=[]
    for x in matches[1:]:
        if not overlapping(x,main) and not any(overlapping(x,y) for y in extras):
            extras.append(x)
    return main,extras[:20]

def focal_flanks(genome,locus,size=4000):
    for chrom,sequence in fasta(genome):
        if chrom==locus['chrom']:
            left=sequence[max(0,locus['start']-size):locus['start']]
            right=sequence[locus['end']:locus['end']+size]
            return (left,right) if len(left)>=750 and len(right)>=750 else None
    return None

def assess_synteny(flanks,hits,relative_assembly,species):
    if not flanks:return {'status':'unknown_flanks_unavailable','locus_status':'unknown'}
    left=[h for h in map_genome(flanks[0],relative_assembly,species)
          if h['query_coverage']>=.55 and h['identity_approx']>=.65]
    right=[h for h in map_genome(flanks[1],relative_assembly,species)
           if h['query_coverage']>=.55 and h['identity_approx']>=.65]
    pairs=[]
    for a in left[:8]:
        for b in right[:8]:
            if a['chrom']!=b['chrom'] or a['strand']!=b['strand']:continue
            if a['strand']=='+' and a['start']>=b['start']:continue
            if a['strand']=='-' and b['start']>=a['start']:continue
            low=min(a['end'],b['end']); high=max(a['start'],b['start'])
            if low>=high or high-low>250000:continue
            score=a['query_coverage']*a['identity_approx']+b['query_coverage']*b['identity_approx']
            pairs.append((score,a,b,low,high))
    if not pairs:
        return {'status':'not_established','locus_status':'unknown','flank_hits':[len(left),len(right)]}
    _,a,b,low,high=max(pairs,key=lambda x:x[0])
    inside=[h for h in hits if h['chrom']==a['chrom'] and
            low <= (h['start']+h['end'])/2 <= high and h['query_coverage']>=.5 and
            h['identity_approx']>=.65]
    elsewhere=[h for h in hits if h not in inside and
               h['query_coverage']>=.5 and h['identity_approx']>=.65]
    if inside:
        if any(h.get('coding_overlap')=='overlaps_cds' for h in inside):
            locus='homolog_overlaps_cds'
        elif all(h.get('coding_overlap')=='no_cds_overlap_not_proof_of_noncoding' for h in inside):
            locus='homolog_no_cds_overlap'
        else:locus='homolog_annotation_unknown'
    elif elsewhere:locus='homolog_elsewhere'
    else:locus='no_detectable_gene_hit'
    return {'status':'flanks_collinear','locus_status':locus,
            'chrom':a['chrom'],'interval_0based':[low,high],
            'homologs_at_syntenic_locus':len(inside),'homologs_elsewhere':len(elsewhere)}
