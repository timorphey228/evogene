"""Run: python -m evogene.pipeline genomes ... OR python -m evogene.pipeline case ..."""
import argparse,csv,json
from pathlib import Path
from .genomes import analyze as analyze_genomes
from .case_data import analyze as analyze_case
def save(r,d):
    path=Path(d);path.mkdir(parents=True,exist_ok=True)
    name=r['gene_id']
    (path/(name+'.json')).write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    hits=r.get('hits') or [x for rows in r.get('relatives',{}).values() for x in rows]
    if hits:
        keys=list(dict.fromkeys(k for h in hits for k in h))
        with (path/(name+'_hits.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,keys);w.writeheader();w.writerows(hits)
    (path/(name+'.md')).write_text('# '+name+'\n\n**Hypothesis:** '+r['hypothesis']+
          '\n\n'+r['explanation']+'\n\n## Limitations\n'+
          '\n'.join('- '+x for x in r['limitations'])+'\n')
    return path/(name+'.json')
def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='mode',required=True)
    g=sub.add_parser('genomes')
    for arg in ('target-fasta','target-gff','gene'):g.add_argument('--'+arg,required=True)
    g.add_argument('--relative',action='append',required=True,help='SPECIES=FASTA or SPECIES=FASTA:GFF3')
    g.add_argument('--results',default='results')
    c=sub.add_parser('case');c.add_argument('--alignment-file',required=True)
    c.add_argument('--gene',action='append',required=True)
    c.add_argument('--results',default='results')
    a=p.parse_args()
    if a.mode=='case':
        for gene in a.gene: print(save(analyze_case(a.alignment_file,gene),a.results))
    else:
        relatives=[]
        for val in a.relative:
            species,paths=val.split('=',1)
            x=paths.rsplit(':',1)
            row={'species':species,'fasta':x[0]}
            if len(x)==2:row['gff']=x[1]
            relatives.append(row)
        print(save(analyze_genomes(a.target_fasta,a.target_gff,a.gene,relatives),a.results))
if __name__=='__main__':main()
