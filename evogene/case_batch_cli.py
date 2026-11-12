"""Repeatable multi-gene evidence report from article's three-section FASTA.

Example:
python -m evogene.case_batch_cli --input supplementary_file_unannotated_alignments.fasta \
 --gene FBgn0040343 --gene FBgn0262810 --gene FBgn0025616 \
 --method raw --output results/case2.json
"""
import argparse
import json
from pathlib import Path
from .scientific_report import analyze_case_batch, batch_as_markdown


def main():
    p = argparse.ArgumentParser(description="Auditable multi-gene EvoGene case 2 report")
    p.add_argument("--input", required=True, help="Three-section supplementary FASTA")
    p.add_argument("--gene", action="append", dest="genes", required=True,
                   help="FBgn ID; repeat for multiple candidates")
    p.add_argument("--method", choices=("raw", "spaln", "genewise"), default="raw")
    p.add_argument("--species", action="append",
                   help="Species code; repeat. Default is all comparison species.")
    p.add_argument("--output", default="results/case2_multi_gene.json")
    p.add_argument("--markdown-output", default=None)
    args = p.parse_args()
    report = analyze_case_batch(
        args.input, args.genes, method=args.method,
        selected_species=args.species, source_name=Path(args.input).name,
    )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n",
                   encoding="utf-8")
    md = Path(args.markdown_output) if args.markdown_output else out.with_suffix(".md")
    md.parent.mkdir(parents=True, exist_ok=True)
    md.write_text(batch_as_markdown(report), encoding="utf-8")
    print(f"Генов: {report['summary']['completed']} / "
          f"{report['summary']['requested']}; ошибок: {report['summary']['failed']}")
    print("Отчёт:", out)
    print("Аннотация:", md)
    if report["errors"]:
        p.exit(1, "Not all genes were analyzed; inspect errors in the JSON\n")


if __name__ == "__main__":
    main()
