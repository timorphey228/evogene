"""Evidence-weighted origin MECHANISM screens, separate from cross-species conservation.

No claim of posterior probabilities: scores are rule weights and missing tests mean
an undetermined classification, not a zero biological probability.
"""
from .origin_tests import screen as legacy_screen
LABELS={'de_novo':'De novo','dna_duplication':'ДНК-дупликация',
        'retroposition':'Ретропозиция','transposition':'ДНК-транспозиция / релокация',
        'fusion':'Слияние / химерное происхождение','unresolved':'Не установлен'}
ORDER=tuple(LABELS)

def classify(report):
    s=legacy_screen(report)
    retro=report.get('retrocopy_evidence') or {}
    candidates=retro.get('parent_candidates') or []
    ref_target=report.get('target_locus') or report.get('target')
    proteins=report.get('protein_evidence') or {}
    protein_candidates=proteins.get('candidates') or []
    focal_target=retro.get('target_model') or {}
    eligible_protein_parents=[
        x for x in protein_candidates if x.get('gene_model') and
        not x.get('overlaps_target') and
        x['gene_model'].get('cds_exons',0)>=2 and
        focal_target.get('status')=='annotated' and
        focal_target.get('exon_count')==1 and
        ref_target and x['gene_model']['chrom']!=ref_target.get('chrom') and
        x.get('identity',0)>=.32 and x.get('query_coverage',0)>=.40
    ]
    protein_paralogs=[
        x for x in protein_candidates if x.get('gene_model') and ref_target and
        not x.get('overlaps_target') and
        x.get('identity',0)>=.35 and x.get('query_coverage',0)>=.50
    ]
    positive_retro=[x for x in candidates if x.get('parent_multi_exon_and_target_single_exon')]
    dup=report.get('extra_target_copies') or []
    synteny=report.get('synteny_evidence') or {}
    errors=report.get('download_or_analysis_errors') or []
    errors2=report.get('origin_screen_errors') or []
    outgroup_noncds=s['de_novo']['observations']['syntenic_noncds_species']
    cdss=s['ancestral_conservation']['observations']['cds_overlap_species']
    relocated=s['transposition']['observations']['putative_displaced_homolog_species']
    mosaic=s['other_mechanism']['observations']['split_homology_species']
    usable_target=retro.get('target_model',{}).get('status')=='annotated'
    target_reference_done=bool(report.get('focal_reference'))
    # Prevent DNA/CDS conservation from competing with mechanism hypotheses.
    w={'de_novo':0,'dna_duplication':0,'retroposition':0,
       'transposition':0,'fusion':0}
    # De novo: outgroup syntenic non-CDS evidence, not merely absent hit.
    if len(outgroup_noncds)>=2 and not dup and not cdss:
        w['de_novo']=7
    if dup or protein_paralogs:
        # Additional DNA or protein-coding loci indicate duplication broadly, but are also
        # compatible with retroposition; reduce nonspecific score when retro signal.
        w['dna_duplication']=min(6,2+2*len(dup)+2*len(protein_paralogs))
    if positive_retro or eligible_protein_parents:
        w['retroposition']=10
        w['dna_duplication']=min(w['dna_duplication'],2)
    if relocated:
        w['transposition']=min(6,2+2*len(relocated))
    if mosaic:
        w['fusion']=min(6,2+2*len(mosaic))
    # Meaningful unresolved mass; missing target annotation and tool errors
    # increase uncertainty, not competing mechanism scores.
    unresolved=8+2*len(errors)+2*len(errors2)+(3 if not usable_target else 0)
    if not target_reference_done:unresolved+=4
    w['unresolved']=unresolved
    summed=sum(w.values())
    raw={k:1000*w[k]/summed for k in ORDER}
    units={k:int(raw[k]) for k in ORDER}
    ranked=sorted(ORDER,key=lambda k:(-(raw[k]-units[k]),ORDER.index(k)))
    for k in ranked[:1000-sum(units.values())]:units[k]+=1
    statuses={}
    reasons={
        'de_novo':f"Синтенные участки без CDS у {len(outgroup_noncds)} видов; нужны предковый ORF и белковые гомологи.",
        'dna_duplication':f"ДНК-локусов: {len(dup)}; белковых гомологов-кандидатов: {len(protein_paralogs)}. Нужна филогения паралогов.",
        'retroposition':f"Кандидатов по структуре CDS и ДНК: {len(positive_retro)}; по белковому сходству и структуре CDS: {len(eligible_protein_parents)}. Нужно подтвердить направление происхождения.",
        'transposition':f"Признаки смещения гомолога при сохранённых фланках: {len(relocated)} видов; нужны признаки вставки.",
        'fusion':f"Два независимых выравнивания разных фрагментов гена обнаружены у {len(mosaic)} видов; требуется валидация.",
        'unresolved':'Часть необходимых тестов не выполнена или свидетельства неоднозначны.'
    }
    observations={
        'de_novo':{'syntenic_noncoding_candidate_species':outgroup_noncds,
                   'coding_hits':cdss},
        'dna_duplication':{'additional_loci':dup[:12], 'protein_paralogs':protein_paralogs[:10]},
        'retroposition':{'dna_architecture':retro,'protein_parent_candidates':eligible_protein_parents[:10],
                         'protein_homology_screen':proteins.get('status','not_run')},
        'transposition':{'displaced_loci_species':relocated},
        'fusion':{'mosaic_species':mosaic},
        'unresolved':{'errors':errors,'origin_screen_errors':errors2}
    }
    tests=[]
    for k in ORDER:
        state='candidate' if w[k] and k!='unresolved' else (
              'unresolved' if k=='unresolved' else 'not_supported_by_current_tests')
        tests.append({'mechanism':k,'label':LABELS[k],'status':state,
                      'support_percent':units[k]/10,
                      'evidence_weight':w[k],'interpretation':reasons[k],
                      'observations':observations[k]})
    possibilities=[k for k in ORDER if k!='unresolved' and w[k]]
    most=max(possibilities,key=lambda k:w[k]) if possibilities else None
    return {'kind':'normalized_mechanism_evidence_not_probabilities',
            'is_calibrated_probability':False,
            'tests':tests,'options':tests,
            'leading_candidate':most,
            'conclusion':('Предварительный кандидат: '+LABELS[most]+
                          '. Для подтверждения нужны независимые проверки.'
                          if most else 'Механизм возникновения гена не установлен.'),
            'conservation':{
                'dna_homology_species':s['ancestral_conservation']['observations']['dna_match_species'],
                'cds_overlap_species':cdss,
                'syntenic_coding_species':s['ancestral_conservation']['observations']['syntenic_coding_species'],
                'explanation':'Это отдельная характеристика сохранности, НЕ механизм возникновения.'},
            'warning':'Проценты отражают эвристические веса свидетельств, а не вероятности. Непроверенный механизм не исключён.'}
