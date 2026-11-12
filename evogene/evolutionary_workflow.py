"""One workflow: candidate-gene FASTA → specified genomes → synteny, homology,
locus phylogeny → alternative hypotheses → UTC audit log and downloadable files.

This implements a *screening analogue* of Peng & Zhao (2024), NOT their exact
Progressive Cactus/HAL/Genewise/Spaln/UniProt calibrated workflow.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import re
import shutil
import subprocess
import sys
import uuid

from Bio import SeqIO

from . import ncbi_auto
from .genomes import map_genome
from .phylogeny import infer_report_tree


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def candidate_records(data, *, max_genes=50):
    """Accept FASTA text/bytes; one unambiguous DNA record per candidate ID."""
    if isinstance(data, bytes):
        data = data.decode("utf-8-sig")
    if not isinstance(data, str):
        raise ValueError("Expected FASTA bytes or string")
    records = []
    for record in SeqIO.parse(StringIO(data), "fasta"):
        fields = record.id.split("|")
        gene_id = fields[0]
        if not re.fullmatch(r"[A-Za-z0-9_.:-]{1,120}", gene_id):
            raise ValueError("Invalid gene identifier: " + gene_id)
        dna = str(record.seq).upper().replace("U", "T")
        if not 30 <= len(dna) <= 100_000 or set(dna) - set("ACGTN"):
            raise ValueError("Each candidate must contain 30–100000 DNA bases A/C/G/T/N: " + gene_id)
        records.append({"gene_id": gene_id, "dna": dna,
                        "length_nt": len(dna), "description": record.description})
    if not records:
        raise ValueError("Candidate FASTA is empty")
    if len(records) > max_genes:
        raise ValueError(f"Maximum {max_genes} genes per submitted project")
    if len({x["gene_id"] for x in records}) != len(records):
        raise ValueError("Gene identifiers in a project must be unique")
    return records


def _versions():
    out={"python":sys.version.split()[0]}
    try:
        out["biopython"]=version("biopython")
    except PackageNotFoundError:
        out["biopython"]=None
    for name, command in (
        ("minimap2", ["minimap2", "--version"]),
        ("datasets", ["datasets", "version"]),
        ("mafft", ["mafft", "--version"]),
        ("iqtree", ["iqtree2", "--version"]),
    ):
        exe=shutil.which(command[0])
        out[name]={"path":exe,"version":None}
        if exe:
            try:
                cp=subprocess.run([exe]+command[1:], text=True,
                                  capture_output=True,timeout=8,check=False)
                out[name]["version"]=(cp.stdout+" "+cp.stderr).strip()[:250]
            except (OSError,subprocess.TimeoutExpired) as exc:
                out[name]["version_error"]=str(exc)
    return out


def _hypotheses(report, te, tree):
    """Differential evidence, no pseudo-posteriors or deterministic de novo."""
    tests=(report.get("mechanism_report") or {}).get("tests", [])
    t={item["mechanism"]:item for item in tests}
    evidence=[]
    duplication=t.get("dna_duplication",{})
    de_novo=t.get("de_novo",{})
    retro=t.get("retroposition",{})
    relocation=t.get("transposition",{})
    fusion=t.get("fusion",{})
    def emit(key,label,test,missing):
        signals=test.get("signals") or []
        status=test.get("status","not_assessed")
        evidence.append({
            "mechanism":key,"label":label,"status":status,
            "observations":test.get("observations",{}),
            "arguments":signals if isinstance(signals,list) else [str(signals)],
            "limitations":test.get("limitations",[]),
            "next_verifications":missing,
            "is_probability":False
        })
    emit("de_novo","Возникновение de novo",de_novo,[
        "Проверить, что синтенный предковый участок некодирующий (Genewise/Spaln, ORF).",
        "Исключить древних гомологов чувствительным BLASTp/jackhmmer/TBLASTN.",
    ])
    emit("duplication","Дупликация/паралогия",duplication,[
        "Подтвердить родство копий белковым выравниванием и деревом генов.",
        "Отличить сегментную/тандемную дупликацию от повторов и ошибок сборки.",
    ])
    emit("retroposition","Ретрокопирование",retro,[
        "Подтвердить структуру экзонов, интронов и происхождение родительской копии.",
    ])
    emit("relocation","Перемещение/транспозиция",relocation,[
        "Проверить структурные признаки мобильных элементов и геномные перестройки.",
    ])
    emit("fusion","Слияние генов",fusion,[
        "Подтвердить независимые гомологии доменов и их предковые локусы.",
    ])
    evidence.append({
        "mechanism":"transposable_element_origin",
        "label":"Происхождение от транспозонной последовательности",
        "status":te.get("status","not_assessed"),
        "observations":{"hits":te.get("hits",[])[:10]},
        "arguments":te.get("arguments",[]),
        "limitations":["Сходство с TE-библиотекой не доказывает механизм интеграции или экспрессию.",
                       "Отрицательный результат неполной библиотеки не исключает TE-происхождение."],
        "next_verifications":["Проверить RepeatMasker/Dfam, фланги, границы TE и родственные копии."],
        "is_probability":False,
    })
    annotated = report.get("coding_species",[])
    evidence.append({
        "mechanism":"ancient_origin",
        "label":"Более древнее происхождение",
        "status":"partial_signal" if len(annotated)>=2 else "not_assessed",
        "observations":{"CDS_annotation_overlaps":annotated},
        "arguments":[f"Кодирующие аннотации перекрываются у {len(annotated)} видов"]
                     if len(annotated)>=2 else [],
        "limitations":["Нельзя автоматически назвать локусы ортологами или оценить возраст.",
                       "Возраст гена устанавливается по древнейшему надёжному ортологу, дереву видов и внешним группам."],
        "next_verifications":["Проверить ортологию и более далёкие таксоны."],
        "is_probability":False,
    })
    candidates=[x for x in evidence if x["status"] in ("candidate","partial_signal")]
    tree_evidence=(
        "Имеется экспериментально вычисленная топология ML/дистанционная; "
        "она не подтверждает автоматически механизм происхождения."
        if tree.get("newick") else
        "Филогенетическая реконструкция невозможна при текущих данных."
    )
    return {
        "status":"preliminary_hypotheses" if candidates else "insufficient_data",
        "most_likely":"not_established",
        "argument":(
            "Сигналы не имеют откалиброванных вероятностей; биологическое происхождение "
            "должно быть проверено независимыми методами."
        ),
        "candidates":[x["mechanism"] for x in candidates],
        "hypotheses":evidence,
        "phylogeny_interpretation":tree_evidence,
        "not_assessed":[x["mechanism"] for x in evidence if x["status"]=="not_assessed"],
    }


def _te_screen(query, te_library):
    if te_library is None:
        return {"status":"not_assessed","arguments":[],
                "reason":"TE library not supplied; no transposon search performed"}
    try:
        hits=map_genome(query,te_library,"TE_consensus")
    except (OSError,RuntimeError,ValueError,subprocess.TimeoutExpired) as exc:
        return {"status":"failed","arguments":[],
                "reason":str(exc),"hits":[]}
    good=[h for h in hits if h.get("query_coverage",0)>=.35 and
          h.get("identity_approx",0)>=.75]
    return {
        "status":"partial_signal" if good else "no_signal",
        "arguments":[f"{len(good)} preliminary DNA matches to supplied TE library"]
                     if good else [],
        "hits":good[:15], "all_mapping_hits":len(hits),
        "reason":"BLAST/RepeatMasker validation and element boundaries remain unassessed",
    }


def _sha_file(path):
    dig=sha256()
    with open(path,"rb") as handle:
        for piece in iter(lambda:handle.read(1024*1024),b""):
            dig.update(piece)
    return dig.hexdigest()


def run_project(fasta_input, focal_species, comparison_species, *,
                cache="genomes/ncbi_cache", te_library=None,
                infer_phylogeny=True, progress=None, output_folder=None):
    """Process multiple target genes; incomplete species do not fabricate absence.

    progress callback gets stage events with gene_id, UTC timestamp, state.
    If output_folder is supplied, writes JSON, Newick and MSA if available.
    """
    if focal_species not in ncbi_auto.PANEL:
        raise ValueError("Unsupported focal species: "+focal_species)
    species=list(comparison_species)
    if not species or len(species)!=len(set(species)) or (
       set(species)-set(ncbi_auto.PANEL[focal_species])):
        raise ValueError("Choose unique supported comparison species different from focal")
    input_bytes=(fasta_input.encode("utf8") if isinstance(fasta_input,str) else fasta_input)
    genes=candidate_records(fasta_input)
    run_id=uuid.uuid4().hex[:12]
    started=utcnow()
    root=Path(output_folder) if output_folder is not None else None
    if root is not None:
        root.mkdir(parents=True,exist_ok=True)
    audit=[]
    versions=_versions()
    summary={
        "schema":"evogene-comparative-pipeline-v1",
        "run_id":run_id, "started_at_utc":started,
        "input_sha256":sha256(input_bytes).hexdigest(),
        "focal_species":focal_species,
        "requested_comparison_species":species,
        "expected_study_assemblies":{
            s:ncbi_auto.REFERENCE_ASSEMBLIES.get(s) for s in [focal_species]+species},
        "tools":versions,"settings":{
            "gene_tree_requested":infer_phylogeny,
            "tree_min_coverage":.65,"tree_min_identity":.60,
            "TE_library":str(te_library) if te_library is not None else None,
            "TE_library_sha256":_sha_file(te_library) if te_library is not None else None,
        },
        "genes":{},"gene_order":[x["gene_id"] for x in genes],"audit_log":audit,
    }
    def log(gene_id,stage,state,detail="",species_name=None,**kwargs):
        event={"timestamp_utc":utcnow(),"gene_id":gene_id,
               "stage":stage,"state":state,"detail":str(detail)}
        if species_name: event["species"]=species_name
        event.update(kwargs)
        audit.append(event)
        if progress:progress(event)
    for gene in genes:
        gene_id=gene["gene_id"]
        log(gene_id,"input","completed",
            f"FASTA validated; length={gene['length_nt']}nt; checksum="
            f"{sha256(gene['dna'].encode()).hexdigest()}")
        events=[]
        def callback(item):
            stamp=item.get("timestamp_utc",utcnow())
            entry={"timestamp_utc":stamp,"gene_id":gene_id,
                   "stage":item.get("stage","unknown"),"state":item.get("state","unknown"),
                   "detail":item.get("detail","")}
            if item.get("species"):entry["species"]=item["species"]
            audit.append(entry)
            events.append(entry)
            if progress:progress(entry)
        try:
            core=ncbi_auto.analyze(
                "",focal_species,gene_id,gene_sequence=gene["dna"],
                relative_species=species,cache=cache,progress=callback)
            log(gene_id,"phylogeny","running","Extract genomic locus sequences for available species")
            if infer_phylogeny:
                tree=infer_report_tree(gene["dna"],core)
            else:
                tree={"status":"not_assessed","newick":None,
                      "reason":"Phylogeny option disabled by user"}
            core["phylogeny"]=tree
            log(gene_id,"phylogeny","completed" if tree.get("newick") else "skipped",
                f"{tree.get('method')}; taxa={tree.get('taxa',0)}; status={tree['status']}",
                support_type=tree.get("support_type"))
            log(gene_id,"transposable_elements","running",
                "Check supplied TE consensus library" if te_library else "TE library missing")
            te=_te_screen(gene["dna"],te_library)
            log(gene_id,"transposable_elements",
                "skipped" if te["status"]=="not_assessed" else
                "failed" if te["status"]=="failed" else "completed",
                te.get("reason",str(te.get("arguments",[]))))
            inference=_hypotheses(core,te,tree)
            log(gene_id,"hypothesis_review","completed",
                f"{inference['status']}; signals={inference['candidates']}")
            record={
                "gene_id":gene_id,"query_length_nt":gene["length_nt"],
                "candidate_sha256":sha256(gene["dna"].encode()).hexdigest(),
                "outcome":"preliminary_hypotheses" if inference["candidates"] else "insufficient_data",
                "hypotheses":inference,
                "phylogeny":tree,
                "transposable_elements":te,
                "evidence":core,
                "limitations":[
                    "This is not a literal reproduction of Cactus/HAL/Spaln/Genewise.",
                    "A tree of DNA matches is not by itself proof of orthology or gene duplication.",
                    "Absence of GFF3/CDS or mapping hits is not proof of noncoding ancestry.",
                    "No global UniProt/jackhmmer/TBLASTN search was performed unless separately run.",
                ],
            }
            summary["genes"][gene_id]=record
            if root is not None:
                p=root/safe_file(gene_id)
                p.mkdir(exist_ok=True)
                (p/"evidence.json").write_text(
                    json.dumps(record,indent=2,ensure_ascii=False,default=str),encoding="utf8")
                if tree.get("newick"):
                    (p/"gene_tree.nwk").write_text(tree["newick"]+"\n",encoding="utf8")
                if tree.get("alignment_fasta"):
                    (p/"locus_alignment.fna").write_text(tree["alignment_fasta"],encoding="utf8")
        except Exception as exc:
            log(gene_id,"gene_analysis","failed",str(exc))
            summary["genes"][gene_id]={
                "gene_id":gene_id,"outcome":"failed",
                "error":str(exc),"phylogeny":{"status":"not_inferred","newick":None},
            }
    summary["finished_at_utc"]=utcnow()
    summary["summary"]={
        "requested":len(genes),
        "completed":sum(r.get("outcome")!="failed" for r in summary["genes"].values()),
        "failed":sum(r.get("outcome")=="failed" for r in summary["genes"].values()),
        "tree_computed":sum(bool(r.get("phylogeny",{}).get("newick"))
                            for r in summary["genes"].values()),
        "total_events":len(audit),
    }
    if root is not None:
        (root/"run_report.json").write_text(
            json.dumps(summary,indent=2,ensure_ascii=False,default=str),encoding="utf8")
        with (root/"audit.jsonl").open("w",encoding="utf8") as handle:
            for line in audit:
                handle.write(json.dumps(line,ensure_ascii=False)+"\n")
    return summary


def safe_file(value):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))[:120]
