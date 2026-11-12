"""Optional local protein homology screen for intron-disrupted parent candidates.

Not a substitute for BLASTP/DIAMOND and phylogenetic direction. Uses a four-mer
prefilter and local protein alignment, reports candidate loci from GFF protein_id.
"""
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from .genomes import fasta, gff

def query_orfs(dna, minimum=55):
    sequence=''.join(c for c in dna.upper() if c in 'ACGTN')
    rc=str(Seq(sequence).reverse_complement())
    found=[]
    for strand,s in (('+',sequence),('-',rc)):
        for frame in range(3):
            translated=str(Seq(s[frame:frame+(len(s)-frame)//3*3]).translate())
            for segment in translated.split('*'):
                if len(segment)>=minimum and segment.count('X')/len(segment)<.10:
                    found.append({'aa':segment,'strand':strand,'frame':frame,'length_aa':len(segment)})
    found.sort(key=lambda e:e['length_aa'],reverse=True)
    return found[:3]

def protein_gene_models(gff_path, accessions):
    if not gff_path:return {}
    keys=set(accessions)
    records={}
    try:
        for row in gff(gff_path):
            if row['type']!='cds':continue
            protein=row['attrs'].get('protein_id','')
            if protein not in keys:continue
            record=records.setdefault(protein,{'chrom':row['chrom'],'start':row['start'],
                'end':row['end'],'exons':set(),'parents':set()})
            record['start']=min(record['start'],row['start'])
            record['end']=max(record['end'],row['end'])
            record['exons'].add((row['start'],row['end']))
            record['parents'].add(row['attrs'].get('Parent',''))
    except (OSError,UnicodeError,ValueError):
        return {}
    return {k:{'chrom':v['chrom'],'start':v['start'],'end':v['end'],
               'cds_exons':len(v['exons']),'transcript_parents':sorted(v['parents'])}
            for k,v in records.items()}

def target_protein_ids(gff_path, target_locus):
    """Protein accessions of CDS transcripts actually overlapping target locus."""
    if not gff_path or not target_locus:
        return []
    ids=set()
    for row in gff(gff_path):
        if row['type']!='cds' or row['chrom']!=target_locus['chrom']:
            continue
        if min(row['end'],target_locus['end'])<=max(row['start'],target_locus['start']):
            continue
        ids.update(x.strip() for x in row['attrs'].get('protein_id','').split(',') if x.strip())
    return sorted(ids)

def local_similarity(a,b):
    aligner=PairwiseAligner()
    aligner.mode='local'
    aligner.match_score=2
    aligner.mismatch_score=-1
    aligner.open_gap_score=-8
    aligner.extend_gap_score=-1
    result=aligner.align(a,b)
    if not result:return None
    blocks=result[0].aligned
    aligned=0;match=0
    for (a0,a1),(b0,b1) in zip(blocks[0],blocks[1]):
        x=a[a0:a1];y=b[b0:b1]
        aligned+=len(x)
        match+=sum(p==q for p,q in zip(x,y))
    if aligned==0:return None
    return {'identity':round(match/aligned,4),
            'query_coverage':round(aligned/len(a),4),
            'identical_aa':match,'aligned_aa':aligned}

def protein_homology(query_dna, protein_fasta, gff_path, target_locus=None, max_candidates=35):
    """Search protein family using annotated target protein first; six-frame ORF fallback.

    A genomic gene FASTA includes introns, so translating its entire DNA often
    produces no useful peptide. Query accessions from target GFF3 avoid this.
    """
    orfs=query_orfs(query_dna)
    ids=target_protein_ids(gff_path,target_locus)
    annotated=[]
    if ids:
        wanted=set(ids)
        for accession,sequence in fasta(protein_fasta):
            if accession in wanted and 40<=len(sequence)<=10000:
                annotated.append({'aa':sequence, 'length_aa':len(sequence),'source':'annotated_cds'})
    queries=(sorted(annotated,key=lambda x:x['length_aa'],reverse=True)[:3]
             if annotated else orfs)
    if not queries:
        return {'status':'not_assessed_no_query_protein','candidates':[],
                'query_source':'unavailable','annotated_target_protein_ids':ids,
                'warning':'No annotated target protein or ORF >=55 aa; protein homology not assessed.'}
    source='annotated_cds' if annotated else 'six_frame_orf'
    keys={q['aa'][i:i+4] for q in queries for i in range(len(q['aa'])-3)}
    pool=[]
    for accession,sequence in fasta(protein_fasta):
        if len(sequence)<45 or len(sequence)>10000:continue
        votes=len(keys.intersection({sequence[i:i+4] for i in range(len(sequence)-3)}))
        if votes>=3:pool.append((votes,accession,sequence))
    pool.sort(key=lambda x:x[0],reverse=True)
    checked=pool[:max_candidates*3]
    matches=[]
    for votes,accession,protein in checked:
        best=None
        for q in queries:
            score=local_similarity(q['aa'],protein)
            if score and (best is None or score['identical_aa']>best['identical_aa']):
                best={**score,'query_orf_length_aa':q['length_aa']}
        if best and best['identity']>=.30 and best['query_coverage']>=.35 and best['identical_aa']>=30:
            matches.append({'protein_accession':accession,'shared_4mers':votes,**best})
    matches.sort(key=lambda x:(x['identical_aa'],x['query_coverage']),reverse=True)
    matches=matches[:max_candidates]
    models=protein_gene_models(gff_path,[x['protein_accession'] for x in matches])
    for x in matches:
        model=models.get(x['protein_accession'])
        x['gene_model']=model
        x['overlaps_target']=bool(
            x['protein_accession'] in ids or
            (model and target_locus and model['chrom']==target_locus['chrom'] and
             min(model['end'],target_locus['end'])>max(model['start'],target_locus['start']))
        )
    return {'status':'searched','query_source':source,
            'annotated_target_protein_ids':ids,
            'query_protein_lengths_aa':[x['length_aa'] for x in queries],
            'screened_proteins_with_4mer_seed':len(pool),
            'aligned_candidates':len(checked),'candidates':matches,
            'warning':'Fast local protein screening, not BLASTP. Absence of hits is inconclusive.'}
