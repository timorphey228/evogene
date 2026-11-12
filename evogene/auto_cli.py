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
print('\nСВОДНАЯ ОЦЕНКА (не статистические вероятности):')
for item in r['origin_support']['options']:
    print(f"{item['label']}: {item['support_percent']:.1f}%")
print(r['origin_support']['warning'])
print(r['hypothesis'],p)
