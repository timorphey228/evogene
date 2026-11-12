"""NCBI automatic relative-assembly search. Requires datasets CLI and minimap2.

FASTA assemblies do not encode gene identifiers. Provide a gene sequence,
1-based coordinates, or a single-record gene FASTA as the input.
"""
import json, re, shutil, subprocess, zipfile
from datetime import datetime, timezone
from pathlib import Path
from .genomes import fasta, map_genome, annotate
from .mechanism_v3 import classify
from .gene_models import compare_retrocopy
from .protein_screen import protein_homology
from .tools import find_tool
from .genome_probes import discover_target_locus, focal_flanks, assess_synteny
from .orf import find_orfs

def fasta_iter(path):
    return fasta(path)

def get_contig(path, name):
    for header, sequence in fasta(path):
        if header == name:
            return sequence
    raise ValueError(f'Contig {name!r} not found in {path}')

def reverse_complement(sequence):
    return sequence.upper().translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


# All 20 fly species represented in the article-associated extract. These are
# a CHOICE PANEL, not a hard cap on scientific comparisons: users can select 1–19
# relatives. Availability and genome/annotation quality must be verified per NCBI.
SPECIES_BY_CODE = {
    "dmel": "Drosophila melanogaster",
    "dsim": "Drosophila simulans",
    "dsec": "Drosophila sechellia",
    "dyak": "Drosophila yakuba",
    "dere": "Drosophila erecta",
    "dfic": "Drosophila ficusphila",
    "drho": "Drosophila rhopaloa",
    "dele": "Drosophila elegans",
    "deug": "Drosophila eugracilis",
    "dtak": "Drosophila takahashii",
    "dbia": "Drosophila biarmipes",
    "dkik": "Drosophila kikkawai",
    "dana": "Drosophila ananassae",
    "dper": "Drosophila persimilis",
    "dpse": "Drosophila pseudoobscura",
    "dmir": "Drosophila miranda",
    "dwil": "Drosophila willistoni",
    "sleb": "Scaptodrosophila lebanonensis",
    "blat": "Bactrocera latifrons",
    "bdor": "Bactrocera dorsalis",
    "mdom": "Musca domestica",
}
SPECIES_ORDER = tuple(SPECIES_BY_CODE.values())
PANEL = {
    species: [relative for relative in SPECIES_ORDER if relative != species]
    for species in SPECIES_ORDER
}

# Pinned assembly versions supplied with Peng & Zhao (2024); the focal
# D. melanogaster FlyBase r6.36 release is NOT an NCBI accession.
REFERENCE_ASSEMBLIES = {
    "Drosophila melanogaster": "FlyBase r6.36",
    "Drosophila simulans": "GCF_000754195.2",
    "Drosophila sechellia": "GCF_004382195.1",
    "Drosophila yakuba": "GCF_000005975.2",
    "Drosophila erecta": "GCF_003286155.1",
    "Drosophila ficusphila": "GCF_000220665.1",
    "Drosophila rhopaloa": "GCF_000236305.1",
    "Drosophila elegans": "GCF_000224195.1",
    "Drosophila eugracilis": "GCF_000236325.1",
    "Drosophila biarmipes": "GCF_000233415.1",
    "Drosophila takahashii": "GCF_000224235.1",
    "Drosophila kikkawai": "GCF_000224215.1",
    "Drosophila ananassae": "GCF_003285975.2",
    "Drosophila miranda": "GCF_003369915.1",
    "Drosophila pseudoobscura": "GCF_009870125.1",
    "Drosophila persimilis": "GCF_003286085.1",
    "Drosophila willistoni": "GCF_000005925.1",
    "Scaptodrosophila lebanonensis": "GCF_003285725.1",
    "Bactrocera dorsalis": "GCF_000789215.1",
    "Bactrocera latifrons": "GCF_001853355.1",
    "Musca domestica": "GCF_000371365.1",
}
CASE={name:acc for name,acc in REFERENCE_ASSEMBLIES.items() if acc.startswith("GCF_")}

def command(args,timeout=180):
    executable=find_tool('datasets')
    if not executable:
        raise RuntimeError('NCBI datasets не обнаружен. EvoGene проверил PATH, ~/bin, /opt/homebrew/bin и /usr/local/bin. Установите datasets и проверьте ~/bin/datasets version.')
    p=subprocess.run([executable,*args[1:]],check=False,text=True,capture_output=True,timeout=timeout)
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
        # The fallback is explicitly flagged; never silently assert article reproducibility.
        fallback_note="case_accession_unavailable_fallback"
    else:
        fallback_note="not_pinned_flybase_ncbi_fallback" if species=="Drosophila melanogaster" else "ncbi_taxon_search"
    matches=valid(command(['datasets','summary','genome','taxon',species,'--annotated','--limit','20']))
    if not matches:matches=valid(command(['datasets','summary','genome','taxon',species,'--limit','20']))
    if not matches:raise LookupError('No valid NCBI assembly for '+species)
    def rank(r):
        i=r.get('assembly_info',{})
        return (int(bool(r.get('annotation_info'))),int(i.get('assembly_level')=='Chromosome'),
                int(r['accession'].startswith('GCF_')))
    return sorted(matches,key=rank,reverse=True)[0],fallback_note

def download(species,cache='genomes/ncbi_cache'):
    metadata,method=find_assembly(species)
    accession=metadata['accession']
    folder=Path(cache)/accession; folder.mkdir(parents=True,exist_ok=True)
    manifest=folder/'manifest.json'
    if manifest.exists():
        data=json.loads(manifest.read_text())
        f=folder/data['fasta'];g=folder/data['gff'] if data.get('gff') else None
        if f.exists() and (not g or g.exists()):
            return {'species':species,'accession':accession,'fasta':str(f),'gff':str(g) if g else None,'cached':True,'source':data.get('source'),'reference_expected':REFERENCE_ASSEMBLIES.get(species),'reference_matches_study':accession==REFERENCE_ASSEMBLIES.get(species)}
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
        return {'species':species,'accession':accession,'fasta':str(f),'gff':str(g) if g else None,'cached':False,'source':method,'reference_expected':REFERENCE_ASSEMBLIES.get(species),'reference_matches_study':accession==REFERENCE_ASSEMBLIES.get(species)}
    finally:archive.unlink(missing_ok=True)

def ensure_proteins(assembly,cache='genomes/ncbi_cache'):
    """Download focal proteome separately, with persistent cache."""
    accession=assembly['accession']
    folder=Path(cache)/accession
    folder.mkdir(parents=True,exist_ok=True)
    for suffix in ('.faa','.faa.gz'):
        path=folder/('proteins'+suffix)
        if path.is_file() and path.stat().st_size:
            return str(path)
    archive=folder/'protein_only.zip'
    try:
        command(['datasets','download','genome','accession',accession,
                 '--include','protein','--filename',str(archive),
                 '--no-progressbar'],1200)
        if not archive.exists() or archive.stat().st_size>300*1024*1024:
            raise ValueError('Protein NCBI archive missing or >300 MB')
        with zipfile.ZipFile(archive) as z:
            candidates=[v for v in z.infolist()
                        if v.filename.startswith('ncbi_dataset/data/') and
                        v.filename.endswith(('.faa','.faa.gz')) and not v.is_dir()]
            if len(candidates)!=1:
                raise ValueError('Expected exactly one protein FASTA (.faa)')
            entry=candidates[0]
            if entry.file_size>600*1024*1024:
                raise ValueError('Uncompressed protein FASTA too large')
            output=folder/('proteins.faa.gz' if entry.filename.endswith('.gz') else 'proteins.faa')
            with z.open(entry) as inp,output.open('wb') as out:
                shutil.copyfileobj(inp,out)
        return str(output)
    finally:
        archive.unlink(missing_ok=True)

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

def analyze(target_fasta,target_species,gene_id,*,gene_sequence='',chrom='',start=None,end=None,strand='+',count=2,relative_species=None,cache='genomes/ncbi_cache',progress=None):
    """Return evidence plus auditable, real-time pipeline events.

    Every step records completed/skipped/failed, never invented progress.
    """
    events=[]
    def emit(stage,state,detail='',species=None):
        entry={'stage':stage,'state':state,'detail':str(detail),
               'timestamp_utc':datetime.now(timezone.utc).isoformat()}
        if species:entry['species']=species
        events.append(entry)
        if progress:progress(entry)

    if target_species not in PANEL:
        raise ValueError('Automatic panel not defined for '+target_species)
    if relative_species is None:
        if not 1 <= count <= len(PANEL[target_species]):
            raise ValueError(f'Relative count must be 1..{len(PANEL[target_species])}')
        selected_species = PANEL[target_species][:count]
    else:
        selected_species = list(relative_species)
        if not selected_species:
            raise ValueError('Choose at least one comparison species')
        if len(set(selected_species)) != len(selected_species):
            raise ValueError('Comparison species must be unique')
        if any(species not in PANEL[target_species] for species in selected_species):
            raise ValueError('Comparison species must differ from target and belong to the supported panel')
    query,mode=target_gene(target_fasta,gene_sequence,chrom,start,end,strand)
    emit('input','completed',f'{len(query)} nt; input={mode}; target={target_species}')
    selected=[]; failures=[]; hits={}; synteny_evidence={}
    for species in selected_species:
        try:
            emit('download','running','Поиск генома и GFF3',species)
            assembly=download(species,cache)
            selected.append(assembly)
            emit('download','completed',
                 f"accession={assembly['accession']}; GFF3={'да' if assembly['gff'] else 'нет'}; cached={assembly.get('cached',False)}",species)
            emit('dna_alignment','running','minimap2: ген → геном',species)
            matches=map_genome(query,assembly['fasta'],species)
            emit('dna_alignment','completed',f'Выравниваний: {len(matches)}',species)
            emit('cds_annotation','running','Сопоставление выравниваний с CDS',species)
            annotate(matches,assembly['gff'])
            hits[species]=matches
            emit('cds_annotation','completed',
                 f"CDS-перекрытий: {sum(x['coding_overlap']=='overlaps_cds' for x in matches)}",species)
        except (OSError,RuntimeError,ValueError,LookupError,subprocess.TimeoutExpired,zipfile.BadZipFile) as exc:
            failures.append({'species':species,'error':str(exc)})
            emit('relative_analysis','failed',str(exc),species)
    if not hits:
        raise RuntimeError('No species analyzed successfully: '+str(failures))

    extra_target_copies=[]; origin_screen_errors=[]; focal_reference=None; main=None
    protein_evidence={'status':'not_run','candidates':[]}
    retrocopy_evidence={}
    try:
        emit('target_reference','running','Поиск исходного референсного генома')
        focal=download(target_species,cache)
        focal_reference={'species':target_species,'accession':focal['accession'],
                         'annotation_available':bool(focal['gff'])}
        emit('target_reference','completed',
             f"accession={focal['accession']}; GFF3={'да' if focal['gff'] else 'нет'}")
        emit('target_mapping','running','Картирование гена на собственный геном')
        self_hits=map_genome(query,focal['fasta'],'target')
        main,extra_target_copies=discover_target_locus(
            query,focal['fasta'],all_hits=self_hits)
        emit('target_mapping','completed',
             f"Совпадений: {len(self_hits)}; локус={'найден' if main else 'не найден'}; дополнительные копии={len(extra_target_copies)}")
        if main is None:
            origin_screen_errors.append('Target locus is not mapped to source reference; synteny and CDS structure unavailable')
        parent_candidates=([hit for hit in self_hits
            if hit['query_coverage']>=.25 and hit['identity_approx']>=.60 and
            not (main and hit['chrom']==main['chrom'] and
                 min(hit['end'],main['end'])>max(hit['start'],main['start']))][:30]
            if main else [])
        emit('gene_structure','running','Анализ CDS-экзонов исходного и похожих локусов')
        retrocopy_evidence=compare_retrocopy(main,parent_candidates,focal['gff'])
        target_cds=retrocopy_evidence.get('target_model',{})
        emit('gene_structure','completed',
             f"CDS={target_cds.get('status','unknown')}; сегментов={target_cds.get('exon_count')}; кандидаты={len(parent_candidates)}")
        try:
            emit('protein_homology','running','Белковый поиск, сначала аннотированный белок гена')
            proteins=ensure_proteins(focal,cache)
            protein_evidence=protein_homology(query,proteins,focal['gff'],main)
            emit('protein_homology','completed',
                 f"Источник={protein_evidence.get('query_source','unknown')}; найдено белковых совпадений={len(protein_evidence.get('candidates',[]))}")
        except (OSError,RuntimeError,ValueError,LookupError,subprocess.TimeoutExpired,zipfile.BadZipFile) as exc:
            origin_screen_errors.append('Protein homology: '+str(exc))
            protein_evidence={'status':'failed','candidates':[],'error':str(exc)}
            emit('protein_homology','failed',str(exc))
        if main:
            flanks=focal_flanks(focal['fasta'],main)
            if not flanks:
                emit('synteny','skipped','Не удалось выделить два фланкирующих участка')
            for assembly in selected:
                species=assembly['species']
                if species not in hits:continue
                try:
                    emit('synteny','running','Сравнение левого и правого флангов',species)
                    synteny_evidence[species]=assess_synteny(
                        flanks,hits[species],assembly['fasta'],species,
                        query=query,gff_path=assembly.get('gff'))
                    detail=(str(synteny_evidence[species].get('status'))+': '+
                            str(synteny_evidence[species].get('locus_status')))
                    rescued=synteny_evidence[species].get('local_rescue')
                    if rescued:detail+='; локальный поиск: '+rescued['status']
                    emit('synteny','completed',detail,species)
                except (OSError,RuntimeError,ValueError,subprocess.TimeoutExpired) as exc:
                    origin_screen_errors.append('Synteny '+species+': '+str(exc))
                    emit('synteny','failed',str(exc),species)
        else:
            emit('synteny','skipped','Не найден локус в исходном геноме')
    except (OSError,RuntimeError,ValueError,LookupError,subprocess.TimeoutExpired,zipfile.BadZipFile) as exc:
        origin_screen_errors.append('Target reference: '+str(exc))
        emit('target_reference','failed',str(exc))

    valid={species:[h for h in hs if h['query_coverage']>=.5 and
                   h['identity_approx']>=.65] for species,hs in hits.items()}
    matching=[s for s,hs in valid.items() if hs]
    coding=[s for s,hs in valid.items()
            if any(h['coding_overlap']=='overlaps_cds' for h in hs)]
    category=('conserved_coding_candidate' if len(coding)>=2 else
              'conserved_dna_locus_origin_unresolved' if matching else
              'unresolved_no_detectable_hit')
    # Genomic uninterrupted ORFs are diagnostic only; exon-spliced CDS
    # and ancestral ORF status require independent annotations/alignments.
    orf_skipped = len(query) > 100_000
    query_orfs = ([o.to_dict() for o in find_orfs(query, min_aa=25)[:80]]
                  if not orf_skipped else [])
    report={
        'query_orfs':query_orfs,
        'orf_note':('ORF scan skipped for DNA >100 kb. Supply spliced CDS for coding analysis.'
                   if orf_skipped else
                   'Candidate ATG-to-stop ORFs in the submitted contiguous DNA (both strands; >=25 aa). '
                   'NOT spliced CDS, ancestral coding state or verified translation.'),
        'phylogeny':{'status':'not_assessed',
                     'reason':'No verified homologous protein MSA, model selection or ML tree is computed.'},
        'gene_id':gene_id,'target_species':target_species,'query_mode':mode,
        'requested_relative_species':selected_species,
        'hypothesis':category,'query_length_nt':len(query),
        'matching_species':matching,'coding_species':coding,'ncbi_assemblies':selected,
        'download_or_analysis_errors':failures,'relatives':hits,
        'extra_target_copies':extra_target_copies,'synteny_evidence':synteny_evidence,
        'origin_screen_errors':origin_screen_errors,'focal_reference':focal_reference,
        'target_locus':main,'retrocopy_evidence':retrocopy_evidence,
        'protein_evidence':protein_evidence,
        'limitations':['No mapping hit is not proof of gene absence.',
                       'Noncoding GFF3 annotations do not establish ancestral noncoding status.',
                       'Protein hits and synteny are preliminary, not a proven origin mechanism.',
                       'Evidence strengths are not calibrated probabilities.'],
        'species_status':([
            {'species':species,'status':('matched' if species in matching else 'no_alignment'),
             'alignment_count':len(hits[species]),'passing_count':len(valid[species])}
            for species in hits
        ]+[
            {'species':f['species'],'status':'failed','error':f['error']}
            for f in failures
        ]),
    }
    emit('mechanism_tests','running','Независимые тесты: de novo, дупликация, ретропозиция, транспозиция, слияние')
    report['mechanism_report']=classify(report)
    for test in report['mechanism_report']['tests']:
        emit('test_'+test['mechanism'],'completed',
             f"{test['status']}; сигналов={test['signal_count']}")
    emit('mechanism_tests','completed',
         f"Кандидатов: {len(report['mechanism_report']['candidate_mechanisms'])}; непроверенных: {len(report['mechanism_report']['completeness']['not_assessed'])}")
    report['pipeline_events']=events
    return report
