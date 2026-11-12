"""Run: python -m evogene.cli --properties data/article_features.csv --evidence data/demo_evidence.csv"""
import argparse
from .core import classify
from .io import export_results, read_evidence, read_properties


def main():
    p = argparse.ArgumentParser(description='Evidence-aware Drosophila gene origin prototype')
    p.add_argument('--properties', help='Original Figure 2 supplementary .xlsx or exported CSV')
    p.add_argument('--evidence', required=True, help='Per-gene observations CSV; demo is synthetic')
    p.add_argument('--output', default='results.json')
    args = p.parse_args()
    props = read_properties(args.properties) if args.properties else {}
    results = [classify(e, props.get(e.gene_id)) for e in read_evidence(args.evidence)]
    export_results(results, args.output)
    for r in results:
        print(f'{r.gene_id}\t{r.hypothesis}\t{r.evidence_strength}\tarticle_label={r.research_label or "n/a"}')
    print(f'Saved {len(results)} results to {args.output}')

if __name__ == '__main__':
    main()
