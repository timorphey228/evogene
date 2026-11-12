"""Genome FASTA/GFF3 -> DNA locus mapping with minimap2 (not validated orthology)."""
from dataclasses import dataclass,asdict
import gzip, subprocess, shutil, tempfile
from pathlib import Path
from .support import distribution

def lines(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path,encoding='utf8')

def fasta(path):
    name=None; seq=[]
    with lines(path) as f:
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if name is not None: yield name,''.join(seq).upper()
                name=line[1:].split()[0]; seq=[]
            elif line and not line.startswith(';'):
                if name is None: raise ValueError('Invalid FASTA')
                seq.append(line)
    if name is not None: yield name,''.join(seq).upper()

def gff(path):
    from urllib.parse import unquote
    with lines(path) as f:
        for line in f:
            if line.startswith('#'): continue
            c=line.rstrip().split('\t')
            if len(c)!=9: continue
            try: a,b=int(c[3])-1,int(c[4])
            except ValueError: continue
            if a<0 or b<=a: continue
            attrs={k:unquote(v) for item in c[8].split(';') if '=' in item for k,v in [item.split('=',1)]}
            yield {'chrom':c[0],'type':c[2].lower(),'start':a,'end':b,'strand':c[6],'attrs':attrs}

def gene_region(fasta_path,gff_path,gene_id):
    gene=None
    for f in gff(gff_path):
        if f['type']!='gene': continue
        attrs=f['attrs']
        names=[attrs.get(k,'') for k in ('ID','Name','gene_id','locus_tag')]
        names+=attrs.get('Dbxref','').split(',')
        if gene_id in names or any(x.endswith(':'+gene_id) for x in names):
            gene=f;break
    if gene is None: raise ValueError('Gene not found in target GFF3: '+gene_id)
    for chrom,seq in fasta(fasta_path):
        if chrom==gene['chrom']:
            if gene['end']>len(seq): raise ValueError('GFF gene exceeds FASTA contig')
            dna=seq[gene['start']:gene['end']]
            if gene['strand']=='-':dna=dna.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
            return gene,dna
    raise ValueError('GFF chromosome not present in target FASTA')

def map_genome(query,reference,species):
    if shutil.which('minimap2') is None:
        raise RuntimeError('Install minimap2 first (macOS: brew install minimap2)')
    with tempfile.TemporaryDirectory() as temp:
        p=Path(temp)/'query.fa'
        p.write_text('>query\n'+query+'\n')
        cp=subprocess.run(['minimap2','-x','asm20','-c','--secondary=yes','-N','10',
                           str(reference),str(p)],text=True,capture_output=True,timeout=600)
    if cp.returncode: raise RuntimeError('minimap2: '+cp.stderr[-500:])
    hits=[]
    for line in cp.stdout.splitlines():
        v=line.split('\t')
        if len(v)<12:continue
        qlen=int(v[1]);aligned=int(v[10])
        hits.append({'species':species,'chrom':v[5],'strand':v[4],
                     'start':int(v[7]),'end':int(v[8]),
                     'query_coverage':round((int(v[3])-int(v[2]))/qlen,4),
                     'identity_approx':round(int(v[9])/aligned,4) if aligned else 0,
                     'mapq':int(v[11])})
    return sorted(hits,key=lambda h:h['identity_approx']*h['query_coverage'],reverse=True)

def annotate(hits,gff_path):
    if not gff_path:
        for hit in hits:hit['coding_overlap']='unknown_no_gff'
        return
    cds=[f for f in gff(gff_path) if f['type']=='cds']
    for hit in hits:
        hit['coding_overlap']=('overlaps_cds' if any(c['chrom']==hit['chrom'] and
            min(c['end'],hit['end'])>max(c['start'],hit['start']) for c in cds)
            else 'no_cds_overlap_not_proof_of_noncoding' if cds else 'unknown_no_cds_features')

def analyze(target_fasta,target_gff,gene_id,relatives,coverage=.5,identity=.65):
    feature,query=gene_region(target_fasta,target_gff,gene_id)
    if len(query)<30:raise ValueError('Target gene interval too short for minimap2')
    if not relatives: raise ValueError('Provide at least one relative assembly')
    out={}
    for r in relatives:
        hits=map_genome(query,r['fasta'],r['species'])
        annotate(hits,r.get('gff'))
        out[r['species']]=hits
    self_hits=map_genome(query,target_fasta,'target')
    extras=[h for h in self_hits if not(h['chrom']==feature['chrom'] and
            min(h['end'],feature['end'])>max(h['start'],feature['start'])) and
            h['query_coverage']>=coverage and h['identity_approx']>=identity]
    passing={s:[h for h in hs if h['query_coverage']>=coverage and h['identity_approx']>=identity] for s,hs in out.items()}
    coding=[s for s,hs in passing.items() if any(h['coding_overlap']=='overlaps_cds' for h in hs)]
    if extras: hyp='possible_duplication'; explanation='An additional target-genome copy was detected; paralogy needs verification.'
    elif len(coding)>=2: hyp='conserved_coding_candidate'; explanation='Similar loci overlap coding annotations in at least two relatives, but orthology is not verified.'
    elif any(passing.values()): hyp='conserved_dna_locus_origin_unresolved'; explanation='DNA similarities exist, but synteny and ancestral coding status are not established.'
    else: hyp='unresolved_no_detectable_hit'; explanation='No hits pass thresholds; missing hits do NOT prove de novo origin.'
    result={'gene_id':gene_id,'input_type':'whole_genome_fasta_gff3',
            'target':{'chrom':feature['chrom'],'start_0based':feature['start'],
            'end_0based':feature['end'],'length_nt':len(query)},
            'hypothesis':hyp,'explanation':explanation,'relatives':out,'extra_target_copies':extras,
            'coding_annotation_species':coding,
            'limitations':['DNA matches do not establish orthology or synteny; inspect flanking genes.',
                           'No CDS overlap is not proof of noncoding ancestry.',
                           'No hit is not proof of absence or de novo origin.',
                           'Gene tree, paralog search and protein analyses remain necessary.']}
    result['origin_support'] = distribution(result)
    return result
