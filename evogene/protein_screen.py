"""Optional local protein homology screen for intron-disrupted parent candidates.

Not a substitute for BLASTP/DIAMOND and phylogenetic direction. Uses a four-mer
prefilter and local protein alignment, reports candidate loci from GFF protein_id.
"""
from collections import Counter
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

def protein_homology(query_dna,protein_fasta,gff_path,target_locus=None,max_candidates=35):
    orfs=query_orfs(query_dna)
    if not orfs:
        return {'status':'no_translatable_orf','candidates':[],
                'warning':'No ORF >=55 aa found; introns may interrupt genomic ORF.'}
    keys=set()
    for o in orfs:
        pep=o['aa']
        keys.update(pep[i:i+4] for i in range(len(pep)-3))
    if not keys:return {'status':'no_seed','candidates':[]}
    # Stream proteins; keep top candidates by shared 4-mers. We cannot claim
    # a negative BLAST result if this sensitivity-limited screen finds none.
    pool=[]
    for accession,sequence in fasta(protein_fasta):
        if len(sequence)<45 or len(sequence)>10000:continue
        kmers={sequence[i:i+4] for i in range(len(sequence)-3)}
        votes=len(keys.intersection(kmers))
        if votes>=3:pool.append((votes,accession,sequence))
    pool.sort(key=lambda x:x[0],reverse=True)
    checked=pool[:max_candidates*3]
    matches=[]
    for votes,acc,protein in checked:
        best=None
        for orf in orfs:
            score=local_similarity(orf['aa'],protein)
            if score and (best is None or score['identical_aa']>best['identical_aa']):
                best={**score,'query_orf_length_aa':orf['length_aa']}
        if best and best['identity']>=.30 and best['query_coverage']>=.35 and best['identical_aa']>=30:
            matches.append({'protein_accession':acc,'shared_4mers':votes,**best})
    matches.sort(key=lambda x:(x['identical_aa'],x['query_coverage']),reverse=True)
    matches=matches[:max_candidates]
    models=protein_gene_models(gff_path,[x['protein_accession'] for x in matches])
    for x in matches:
        model=models.get(x['protein_accession'])
        x['gene_model']=model
        x['overlaps_target']=bool(model and target_locus and
            model['chrom']==target_locus['chrom'] and
            min(model['end'],target_locus['end'])>max(model['start'],target_locus['start']))
    return {'status':'searched','orf_lengths_aa':[x['length_aa'] for x in orfs],
            'screened_proteins_with_4mer_seed':len(pool),
            'aligned_candidates':len(checked),
            'candidates':matches,
            'warning':'Local ORF/protein screen, not BLAST: absence of hit is inconclusive.'}
