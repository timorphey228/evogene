"""Read real LiZhaoLab-style raw Cactus hits: they are NOT whole genomes."""
from Bio.Align import PairwiseAligner

def sections(path):
    gene=None;records=[];header=None;seq=[]
    with open(path,encoding='utf8') as stream:
        for line in stream:
            line=line.strip()
            if line.startswith('raw sequence for '):
                if gene is not None:
                    if header:records.append((header,''.join(seq)))
                    yield gene,records
                gene=line.split('raw sequence for ',1)[1];records=[];header=None;seq=[]
            elif line.startswith('>'):
                if header:records.append((header,''.join(seq)))
                header=line[1:];seq=[]
            elif header and line:seq.append(line.upper())
    if gene is not None:
        if header:records.append((header,''.join(seq)))
        yield gene,records

def compare(a,b):
    aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-2,
                             open_gap_score=-5,extend_gap_score=-.5)
    al=aligner.align(a,b)[0]
    equal=total=0
    for (a0,a1),(b0,b1) in zip(*al.aligned):
        for x,y in zip(a[a0:a1],b[b0:b1]):
            equal+=int(x==y and x in 'ACGT');total+=1
    return {'identity_aligned':round(equal/total,4) if total else 0,
            'query_coverage':round(sum(b-a for a,b in al.aligned[0])/len(a),4)}

def analyze(path,gene_id):
    entries=None;name=None
    for section,records in sections(path):
        if section.endswith('_'+gene_id):
            entries=records;name=section;break
    if entries is None:raise ValueError('Gene not found: '+gene_id)
    targets=[(h,s) for h,s in entries if h.startswith('dmel|') and 'FBgn' in h]
    if not targets:raise ValueError('Missing target DNA')
    header,query=targets[0]
    annot=[];loci=[]
    for h,s in entries:
        if h==header or not s:continue
        item={'species_abbr':h.split('|',1)[0],'header':h,'length_nt':len(s),
              'header_kind':'genomic_coordinates' if ':' in h.split('|',1)[1] else 'gene_identifier',
              **compare(query,s)}
        if item['header_kind']=='gene_identifier':annot.append(item)
        elif len(loci)<8:loci.append(item)
    hits=annot+loci
    coded=sorted({h['species_abbr'] for h in annot if h['query_coverage']>=.6 and h['identity_aligned']>=.65})
    coords=sorted({h['species_abbr'] for h in loci if h['query_coverage']>=.6 and h['identity_aligned']>=.65})
    hyp=('widespread_annotated_homolog_candidates' if len(coded)>=2 else
         'near_lineage_annotated_homolog_candidate' if coded else
         'unannotated_genomic_locus_homology' if coords else 'unresolved')
    return {'gene_id':gene_id,'source_section':name,
            'input_type':'preselected_cactus_hits_not_whole_genomes','target_length_nt':len(query),
            'hypothesis':hyp,'explanation':'Similar DNA sequences in preselected Cactus hits. The evolutionary mechanism cannot yet be established.',
            'summary':{'annotated_hit_species':coded,'genomic_coordinate_hit_species':coords,'record_count':len(hits)},
            'hits':hits,'limitations':['Preselected hits cannot replace a whole-genome search.',
                'Gene-like identifiers do not independently establish protein orthology.',
                'Unannotated genomic coordinates do not prove noncoding ancestry.',
                'Pairwise local DNA similarity does not resolve gene duplication vs conserved ancestry.']}
