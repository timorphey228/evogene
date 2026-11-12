"""NCBI automatic relative-assembly search. Requires datasets CLI and minimap2.

FASTA assemblies do not encode gene identifiers. Provide a gene sequence,
1-based coordinates, or a single-record gene FASTA as the input.
"""
import json, re, shutil, subprocess, zipfile
from pathlib import Path
from .genomes import fasta, map_genome, annotate

def fasta_iter(path):
    return fasta(path)

def get_contig(path, name):
    for header, sequence in fasta(path):
        if header == name:
            return sequence
    raise ValueError(f'Contig {name!r} not found in {path}')

def reverse_complement(sequence):
    return sequence.upper().translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


PANEL={
'Drosophila melanogaster':['Drosophila simulans','Drosophila sechellia','Drosophila yakuba','Drosophila erecta','Drosophila ananassae'],
'Drosophila simulans':['Drosophila melanogaster','Drosophila sechellia','Drosophila yakuba','Drosophila erecta','Drosophila ananassae'],
'Drosophila yakuba':['Drosophila erecta','Drosophila melanogaster','Drosophila simulans','Drosophila sechellia','Drosophila ananassae']}
CASE={'Drosophila simulans':'GCF_000754195.2','Drosophila sechellia':'GCF_004382195.1','Drosophila yakuba':'GCF_000005975.2','Drosophila erecta':'GCF_003286155.1','Drosophila ananassae':'GCF_003285975.2'}

def command(args,timeout=180):
    if not shutil.which('datasets'):raise RuntimeError('Install NCBI CLI: brew install ncbi-datasets-cli')
    p=subprocess.run(args,check=False,text=True,capture_output=True,timeout=timeout)
    if p.returncode:raise RuntimeError('NCBI datasets: '+p.stderr[-600:])
    return p.stdout

def find_assembly(species):
    def valid(s):
        r=json.loads(s).get('reports',[])
        return [x for x in r if re.fullmatch(r'GC[AF]_\d{9}\.\d+',str(x.get('accession',''))) and
                str(x.get('organism',{}).get('organism_name','')).lower().startswith(species.lower())]
    accession=CASE.get(species)
    if accession:
        try:
            matches=valid(command(['datasets','summary','genome','accession',accession]))
            if matches:return matches[0],'case_accession_verified'
        except (RuntimeError,ValueError):pass
    matches=valid(command(['datasets','summary','genome','taxon',species,'--annotated','--limit','20']))
    if not matches:matches=valid(command(['datasets','summary','genome','taxon',species,'--limit','20']))
    if not matches:raise LookupError('No valid NCBI assembly for '+species)
    def rank(r):
        i=r.get('assembly_info',{})
        return (int(bool(r.get('annotation_info'))),int(i.get('assembly_level')=='Chromosome'),
                int(r['accession'].startswith('GCF_')))
    return sorted(matches,key=rank,reverse=True)[0],'ncbi_taxon_search'

def download(species,cache='genomes/ncbi_cache'):
    metadata,method=find_assembly(species)
    accession=metadata['accession']
    folder=Path(cache)/accession; folder.mkdir(parents=True,exist_ok=True)
    manifest=folder/'manifest.json'
    if manifest.exists():
        data=json.loads(manifest.read_text())
        f=folder/data['fasta'];g=folder/data['gff'] if data.get('gff') else None
        if f.exists() and (not g or g.exists()):
            return {'species':species,'accession':accession,'fasta':str(f),'gff':str(g) if g else None,'cached':True}
    archive=folder/'download.zip'
    try:
        try:
            command(['datasets','download','genome','accession',accession,'--include','genome,gff3','--filename',str(archive),'--no-progressbar'],1200)
        except RuntimeError:
            archive.unlink(missing_ok=True)
            command(['datasets','download','genome','accession',accession,'--include','genome','--filename',str(archive),'--no-progressbar'],1200)
        if not archive.exists() or archive.stat().st_size>700*1024*1024:raise ValueError('Genome ZIP missing or >700 MB')
        with zipfile.ZipFile(archive) as z:
            members=[x for x in z.infolist() if x.filename.startswith('ncbi_dataset/data/') and not x.is_dir()]
            if sum(x.file_size for x in members)>3_000_000_000:raise ValueError('Uncompressed dataset too large')
            fasta=[x for x in members if x.filename.endswith(('.fna','.fna.gz'))]
            gffs=[x for x in members if x.filename.endswith(('.gff','.gff3','.gff.gz','.gff3.gz'))]
            if len(fasta)!=1:raise ValueError('Expected precisely one assembly FASTA')
            f=folder/('genome.fna.gz' if fasta[0].filename.endswith('.gz') else 'genome.fna')
            with z.open(fasta[0]) as inp,f.open('wb') as out:shutil.copyfileobj(inp,out)
            g=None
            if gffs:
                g=folder/('annotation.gff3.gz' if gffs[0].filename.endswith('.gz') else 'annotation.gff3')
                with z.open(gffs[0]) as inp,g.open('wb') as out:shutil.copyfileobj(inp,out)
        (folder/'metadata.json').write_text(json.dumps(metadata,indent=2))
        manifest.write_text(json.dumps({'species':species,'accession':accession,'source':method,'fasta':f.name,'gff':g.name if g else None},indent=2))
        return {'species':species,'accession':accession,'fasta':str(f),'gff':str(g) if g else None,'cached':False}
    finally:archive.unlink(missing_ok=True)

def target_gene(path,sequence='',chrom='',start=None,end=None,strand='+'):
    if sequence:
        dna=''.join(x.strip() for x in sequence.upper().splitlines() if not x.lstrip().startswith('>')).replace(' ','')
        mode='pasted_gene_sequence'
    elif chrom and start is not None and end is not None:
        if int(start)<1 or int(end)<int(start):raise ValueError('Coordinates must be 1-based inclusive')
        dna=get_contig(path,chrom)[int(start)-1:int(end)]
        if strand=='-':dna=reverse_complement(dna)
        mode='gene_coordinates'
    else:
        it=fasta_iter(path);first=next(it,None)
        if not first or len(first[1])>100000 or next(it,None) is not None:
            raise ValueError('Whole-genome FASTA does not identify the gene: enter coordinates or paste gene DNA.')
        dna=first[1];mode='single_record_gene_fasta'
    if len(dna)<30 or any(x not in 'ACGTN' for x in dna):raise ValueError('Gene query must be >=30 DNA nucleotides A/C/G/T/N')
    return dna,mode

def analyze(target_fasta,target_species,gene_id,*,gene_sequence='',chrom='',start=None,end=None,strand='+',count=2,cache='genomes/ncbi_cache',progress=None):
    if target_species not in PANEL:raise ValueError('Automatic panel not defined for '+target_species)
    if not 1<=count<=5:raise ValueError('Relative count must be 1..5')
    query,mode=target_gene(target_fasta,gene_sequence,chrom,start,end,strand)
    selected=[];failures=[];hits={}
    for species in PANEL[target_species][:count]:
        if progress:progress('NCBI: '+species)
        try:
            assembly=download(species,cache)
            selected.append(assembly)
            if progress:progress('Aligning: '+species)
            m=map_genome(query,assembly['fasta'],species)
            annotate(m,assembly['gff'])
            hits[species]=m
        except (OSError,RuntimeError,ValueError,LookupError,subprocess.TimeoutExpired,zipfile.BadZipFile) as exc:
            failures.append({'species':species,'error':str(exc)})
    if not hits:raise RuntimeError('No species analyzed successfully: '+str(failures))
    valid={s:[h for h in hs if h['query_coverage']>=.5 and h['identity_approx']>=.65] for s,hs in hits.items()}
    matching=[s for s,hs in valid.items() if hs]
    coding=[s for s,hs in valid.items() if any(h['coding_overlap']=='overlaps_cds' for h in hs)]
    category=('conserved_coding_candidate' if len(coding)>=2 else
              'conserved_dna_locus_origin_unresolved' if matching else 'unresolved_no_detectable_hit')
    return {'gene_id':gene_id,'target_species':target_species,'query_mode':mode,
            'hypothesis':category,'matching_species':matching,'coding_species':coding,
            'ncbi_assemblies':selected,'download_or_analysis_errors':failures,'relatives':hits,
            'limitations':['DNA mapping alone does not prove orthology, synteny or de novo origin.',
              'No CDS overlap is not proof of noncoding ancestry.',
              'No mapping hit is not proof of absence.',
              'Relative panel is curated, not inferred phylogenetically.']}
