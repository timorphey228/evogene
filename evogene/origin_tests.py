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
    # Screen for mosaic/fusion-like structure: disjoint query segments align to
    # distant genomic regions in an outgroup; can also reflect rearrangement or
    # fragmented assemblies and therefore is never a confirmed gene fusion.
    mosaic=[]
    for species, all_hits in (report.get('relatives') or {}).items():
        parts=[h for h in all_hits if h.get('query_length',0)>0 and
               .20<=h.get('query_coverage',0)<.80 and
               h.get('identity_approx',0)>=.65 and
               'query_start' in h and 'query_end' in h]
        for i,a in enumerate(parts):
            for b in parts[i+1:]:
                n=a['query_length']
                overlap=max(0,min(a['query_end'],b['query_end'])-
                            max(a['query_start'],b['query_start']))
                independent=(a.get('chrom')!=b.get('chrom') or
                             abs(a.get('start',0)-b.get('start',0))>100000)
                total=(a['query_end']-a['query_start']+
                       b['query_end']-b['query_start']-overlap)/n
                if independent and overlap/n<=.05 and total>=.55:
                    mosaic.append(species)
                    break
            if species in mosaic:break
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
        'status':'preliminary_candidate' if mosaic else 'insufficient_evidence',
        'support_points':min(4,2*len(mosaic)),
        'observations':{'split_homology_species':mosaic},
        'missing':['Validate gene fusion or mosaic ancestry by independent protein/domain alignments',
                   'Other mechanisms (HGT, TE domestication, exon shuffling) need dedicated tests'],
        'interpretation':'Disjoint query segments mapping to distant loci may suggest a chimeric/fusion origin, not proof.'}
    }
