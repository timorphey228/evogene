"""Evidence screens for hypotheses. These are candidates, not proofs."""
def passing(hits):
    return [h for h in hits if h.get('query_coverage',0)>=.5 and h.get('identity_approx',0)>=.65]

def screen(report):
    hs={s:passing(h) for s,h in (report.get('relatives') or {}).items()}
    matched=[s for s,h in hs.items() if h]
    coding=[s for s,h in hs.items() if any(x.get('coding_overlap')=='overlaps_cds' for x in h)]
    noncds=[s for s,h in hs.items() if h and all(x.get('coding_overlap')=='no_cds_overlap_not_proof_of_noncoding' for x in h)]
    extras=passing(report.get('extra_target_copies',[]))
    synteny=report.get('synteny_evidence',{})
    noncoding_syn=[s for s,v in synteny.items() if v.get('status')=='flanks_collinear' and v.get('locus_status')=='homolog_no_cds_overlap']
    coding_syn=[s for s,v in synteny.items() if v.get('status')=='flanks_collinear' and v.get('locus_status')=='homolog_overlaps_cds']
    displaced=[s for s,v in synteny.items() if v.get('status')=='flanks_collinear' and v.get('locus_status')=='homolog_elsewhere']
    denovo= len(noncoding_syn)>=2 and not extras and not coding
    return {
      'de_novo':{
        'status':'preliminary_candidate' if denovo else 'insufficient_evidence',
        'support_points':3*min(2,len(noncoding_syn)) if denovo else 0,
        'observations':{'syntenic_noncds_species':noncoding_syn,'noncds_without_synteny_species':noncds},
        'missing':['Verify ancestral ORF and unannotated homologs','Confirm appropriate phylogenetic outgroups'],
        'interpretation':'Syntenic DNA without annotated CDS suggests a candidate only; noncoding ancestry is not proven.'},
      'duplication':{
        'status':'preliminary_candidate' if extras else 'insufficient_evidence',
        'support_points':min(7,len(extras)*3),
        'observations':{'extra_target_loci':len(extras)},
        'missing':['Confirm paralogy with protein homology or gene trees','Exclude repetitive sequence matches'],
        'interpretation':'Additional non-overlapping target loci may indicate DNA duplication, not proven gene paralogy.'},
      'transposition':{
        'status':'preliminary_candidate' if displaced else 'insufficient_evidence',
        'support_points':2*min(2,len(displaced)),
        'observations':{'putative_displaced_homolog_species':displaced},
        'missing':['Check insertion signatures and mobile elements','Determine orthology and ancestral genomic location'],
        'interpretation':'Homology elsewhere despite conserved flanks is a relocation signal, not proof of transposition.'},
      'ancestral_conservation':{
        'status':'preliminary_candidate' if coding or coding_syn else ('dna_similarity_only' if matched else 'insufficient_evidence'),
        'support_points':min(6,len(matched))+min(6,2*len(coding))+min(4,2*len(coding_syn)),
        'observations':{'dna_match_species':matched,'cds_overlap_species':coding,'syntenic_coding_species':coding_syn},
        'missing':['Verify orthology and ancestral coding state'],
        'interpretation':'Coding-annotated cross-species DNA matches support a conserved candidate.'},
      'other_mechanism':{
        'status':'not_assessed','support_points':0,'observations':{},
        'missing':['Fusion, horizontal transfer, transposon domestication and exon shuffling require dedicated data'],
        'interpretation':'No dedicated screen for other mechanisms; they remain possible.'}
    }

def distribute(report):
    tests=screen(report)
    errors=len(report.get('download_or_analysis_errors',[]))
    w={key:v['support_points'] for key,v in tests.items()}
    w['unresolved']=8+2*errors
    order=('de_novo','duplication','transposition','ancestral_conservation','other_mechanism','unresolved')
    labels={'de_novo':'De novo','duplication':'Дупликация','transposition':'Транспозиция / ретротранспозиция',
            'ancestral_conservation':'Древнее сохранённое происхождение','other_mechanism':'Другой механизм',
            'unresolved':'Недостаточно данных'}
    raw={key:1000*w[key]/sum(w.values()) for key in order}
    units={key:int(raw[key]) for key in order}
    ranks=sorted(order,key=lambda key:(-(raw[key]-units[key]),order.index(key)))
    for key in ranks[:1000-sum(units.values())]:units[key]+=1
    return {'kind':'normalized_screening_evidence_not_probability','is_calibrated_probability':False,
       'scoring_method':'origin_screens_v2','tests':tests,
       'options':[{'mechanism':key,'label':labels[key],'support_percent':units[key]/10,
                   'evidence_weight':w[key],
                   'status':tests[key]['status'] if key in tests else 'uncertain',
                   'explanation':tests[key]['interpretation'] if key in tests else 'Mechanism unresolved.'}
                  for key in order],
       'warning':'Проценты — нормированные оценки диагностических признаков, НЕ вероятности. 0% не исключает механизм.',
       'conclusion':('Есть предварительные кандидаты по диагностическим признакам; происхождение требует проверки.'
                     if any(t['status']=='preliminary_candidate' for t in tests.values())
                     else 'Недостаточно данных, чтобы установить механизм происхождения.')}
