"""Targeted CDS exon annotation for suspected parent/child genomic loci."""

def _overlap(a,b):
    return a['chrom']==b['chrom'] and min(a['end'],b['end'])>max(a['start'],b['start'])

def models_for_loci(gff_path,loci):
    """Two GFF passes regardless of the number of loci."""
    from .genomes import gff
    if not gff_path:
        return [{'status':'unavailable','exon_count':None,'reason':'No GFF3 annotation'} for _ in loci]
    parents=[set() for _ in loci]
    try:
        for row in gff(gff_path):
            if row['type']!='cds':continue
            for i,locus in enumerate(loci):
                if locus and _overlap(row,locus):
                    ids=row['attrs'].get('Parent') or row['attrs'].get('protein_id')
                    if ids:parents[i].update(ids.split(','))
        relevant=set().union(*parents) if parents else set()
        counts={p:set() for p in relevant}
        if relevant:
            for row in gff(gff_path):
                if row['type']!='cds':continue
                ids=row['attrs'].get('Parent') or row['attrs'].get('protein_id') or ''
                for p in ids.split(','):
                    if p in counts:
                        counts[p].add((row['chrom'],row['start'],row['end']))
    except (OSError,UnicodeError,ValueError) as exc:
        return [{'status':'unavailable','exon_count':None,'reason':str(exc)} for _ in loci]
    out=[]
    for matched in parents:
        if not matched:
            out.append({'status':'no_annotated_cds','exon_count':None,
                        'reason':'No CDS overlapping hit, not proof of noncoding DNA'})
            continue
        exon_counts=[len(counts[p]) for p in matched if p in counts]
        out.append({'status':'annotated','exon_count':min(exon_counts),
                    'max_exon_count':max(exon_counts),
                    'transcript_count':len(exon_counts),
                    'reason':'CDS segments across transcript; alternative isoforms may differ'})
    return out

def compare_retrocopy(target,candidates,target_gff):
    loci=([target] if target else [None])+list(candidates)
    models=models_for_loci(target_gff,loci)
    t=models[0]
    records=[]
    for h,model in zip(candidates,models[1:]):
        is_other=bool(target and h['chrom']!=target['chrom'])
        signal=(t['status']=='annotated' and t['exon_count']==1 and
                model['status']=='annotated' and model['exon_count']>=2 and
                is_other and h.get('query_coverage',0)>=.25 and
                h.get('identity_approx',0)>=.60)
        records.append({'chrom':h['chrom'],'start':h['start'],'end':h['end'],
            'query_coverage':h.get('query_coverage',0),
            'identity_approx':h.get('identity_approx',0),
            'other_chromosome':is_other,'candidate_parent_model':model,
            'parent_multi_exon_and_target_single_exon':signal})
    records.sort(key=lambda x:(x['parent_multi_exon_and_target_single_exon'],
                                x['query_coverage']*x['identity_approx']),reverse=True)
    return {'target_model':t,'parent_candidates':records[:20],
            'retrocopy_architecture_signal':any(x['parent_multi_exon_and_target_single_exon'] for x in records),
            'limitations':['Intron loss can occur without retroposition.',
                           'GFF absence is unknown, not intronless.',
                           'Parent/child direction needs gene tree and outgroup evidence.',
                           'DNA self-mapping is not sensitive enough for distant protein paralogs.']}
