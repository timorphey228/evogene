"""Command line entrypoint for the EvoGene comparative-genomics pipeline.

Example:
  python -m evogene.workflow_cli --fasta genes.fasta \
     --focal "Drosophila melanogaster" \
     --relative "Drosophila simulans" \
     --relative "Drosophila yakuba" \
     --relative "Drosophila ananassae" \
     --output results/evolutionary_pipeline
"""
from __future__ import annotations

import argparse
from pathlib import Path

from .evolutionary_workflow import run_project
from .ncbi_auto import PANEL


def main():
    p=argparse.ArgumentParser(
        description="EvoGene: reproducible multi-gene hypothesis + locus-tree workflow")
    p.add_argument("--fasta",required=True,help="FASTA of focal gene candidates")
    p.add_argument("--focal",required=True,choices=sorted(PANEL),help="Focal species")
    p.add_argument("--relative",action="append",required=True,dest="relatives",
                   help="Comparison species; repeat for each selected species")
    p.add_argument("--te-library",default=None,
                   help="Optional transposable-element library FASTA")
    p.add_argument("--no-tree",action="store_true",
                   help="Do not reconstruct a locus tree")
    p.add_argument("--cache",default="genomes/ncbi_cache")
    p.add_argument("--output",default="results/evolutionary_pipeline")
    args=p.parse_args()
    def progress(event):
        if event["state"]=="running":
            return
        print(f"[{event['timestamp_utc']}] {event['gene_id']} "
              f"{event['stage']} / {event['state']}: {event.get('detail','')}")
    report=run_project(
        Path(args.fasta).read_bytes(),args.focal,args.relatives,
        te_library=args.te_library,infer_phylogeny=not args.no_tree,
        cache=args.cache,output_folder=args.output,progress=progress)
    print("ID:",report["run_id"])
    print("Report:",Path(report["output_folder"])/"run_report.json")
    print("Audit:",Path(report["output_folder"])/"audit.jsonl")
    print("Genes:",report["summary"]["completed"],"completed /",
          report["summary"]["failed"],"failed")
    print("Trees:",report["summary"]["tree_computed"])
    if report["summary"]["failed"]:
        p.exit(1,"Some gene analyses failed. Inspect audit.jsonl.\n")


if __name__=="__main__":
    main()
