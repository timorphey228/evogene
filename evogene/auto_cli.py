import argparse,json
from pathlib import Path
from .ncbi_auto import analyze
p=argparse.ArgumentParser()
p.add_argument('--target-fasta',required=True)
p.add_argument('--species',default='Drosophila melanogaster')
p.add_argument('--gene-id',default='target_gene')
p.add_argument('--gene-sequence',default='')
p.add_argument('--chrom',default='')
p.add_argument('--start',type=int)
p.add_argument('--end',type=int)
p.add_argument('--strand',default='+',choices=['+','-'])
p.add_argument('--count',type=int,default=2)
p.add_argument('--results',default='results/auto_ncbi')
a=p.parse_args()
r=analyze(a.target_fasta,a.species,a.gene_id,gene_sequence=a.gene_sequence,
          chrom=a.chrom,start=a.start,end=a.end,strand=a.strand,count=a.count,
          progress=print)
d=Path(a.results);d.mkdir(parents=True,exist_ok=True)
p=d/(a.gene_id+'.json');p.write_text(json.dumps(r,indent=2,ensure_ascii=False))
print('\nМЕХАНИЗМ ВОЗНИКНОВЕНИЯ (независимые признаки, не вероятности):')
m=r.get('mechanism_report')
if m:
    for item in m['options']:
        print(f"{item['label']}: {item['status']} — признаков {item['signal_count']}")
        for signal in item['signals']:
            print('  +', signal)
    print('Сохранность ДНК у видов:', ', '.join(m['conservation']['dna_homology_species']))
    print(m['warning'])
else:
    print('Новый классификатор отсутствует в отчёте.')
print(r['hypothesis'],p)
