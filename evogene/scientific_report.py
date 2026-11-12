"""Case 2: per-gene conclusions with auditable evidence, not invented mechanisms.

The published supplementary file contains *selected* Cactus/Spaln/Genewise
candidates and cannot alone establish de novo origin, duplication or losses.
This module handles one or MANY gene IDs and maps results to all six case
requirements without promoting missing data to a negative observation.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError
import json
import re
import sys

from .case_data import analyze_many as analyze_case_many

SCHEMA = "evogene-case2-auditable-v1"

# Five broad outcomes of the hackathon case; mechanisms are not inferred here.
OUTCOMES = ("de_novo", "duplication", "ancient_origin", "other_mechanism",
            "insufficient_data")


def _software_versions():
    from subprocess import run, TimeoutExpired
    try:
        biopython = version("biopython")
    except PackageNotFoundError:
        biopython = None
    try:
        repo = Path(__file__).resolve().parents[1]
        proc = run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                   text=True, capture_output=True, timeout=3, check=False)
        commit = proc.stdout.strip() if proc.returncode == 0 else None
    except (OSError, TimeoutError, TimeoutExpired):
        commit = None
    return {"python": sys.version.split()[0], "biopython": biopython,
            "git_commit": commit, "analysis": SCHEMA}


def hash_file(path, block_size=1024 * 1024):
    digest = sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def evidence_for_case(report):
    """Every item contains observed data + explicit interpretation boundary."""
    summary = report.get("summary") or {}
    hits = report.get("hits") or []
    annotated = summary.get("annotated_hit_species") or []
    coordinate = summary.get("genomic_coordinate_hit_species") or []
    complete = [h for h in hits if h.get("comparison_status") == "aligned"]
    skipped = [h for h in hits if h.get("comparison_status") != "aligned"]
    method_rows = report.get("method_comparison") or []
    evidence = [
        {
            "id": "observed_comparative_records",
            "type": "observation",
            "description": ("Пользовательские сравнительные последовательности"
                            if report.get("input_type") == "user_supplied_fasta_sequences"
                            else "Записи других видов для данного гена в выбранной секции"),
            "value": {
                "records": len(hits),
                "species": sorted({h.get("species_abbr") for h in hits
                                   if h.get("species_abbr")}),
                "annotated_identifier_species": annotated,
                "genomic_coordinate_species": coordinate,
            },
            "source": ("Заголовки пользовательского сравнительного FASTA"
                       if report.get("input_type") == "user_supplied_fasta_sequences"
                       else "FASTA headers in " + report.get("method", "unknown")),
            "interpretation": "Отобранные сравнительные последовательности доступны для проверки.",
            "does_not_prove": "Функцию белка, ортологию, паралогию или de novo происхождение.",
        },
        {
            "id": "observed_local_dna_similarity",
            "type": "observation",
            "description": "Локальные сравнения с последовательностью dmel",
            "value": {
                "aligned_records": len(complete),
                "not_aligned_or_skipped": len(skipped),
                "high_coverage_identity_records": sum(
                    (h.get("query_coverage") or 0) >= .6 and
                    (h.get("identity_aligned") or 0) >= .65 for h in complete),
                "parameters": {"query_coverage_min": .6, "identity_min": .65,
                               "comparison": "Biopython local PairwiseAligner",
                               "scoring": {"match": 2, "mismatch": -2,
                                           "gap_open": -5, "gap_extend": -0.5},
                               "alignment_cost_guard_cells": 2_000_000},
            },
            "source": "Попарное локальное ДНК-выравнивание, вычисленное EvoGene",
            "interpretation": "Есть измеримые совпадения между выбранными фрагментами.",
            "does_not_prove": "Независимый полно-геномный поиск гомологов и эволюционное событие.",
        },
        {
            "id": "three_independent_sections",
            "type": "observation",
            "description": "Сопоставление исходного Cactus и двух предсказаний",
            "value": [
                {
                    "method": row.get("method"),
                    "records": row.get("records"),
                    "species": row.get("species"),
                    "median_nt": row.get("median_nt"),
                    "shift_yes_header_count": row.get("shift_yes"),
                }
                for row in method_rows
            ],
            "source": "Raw Cactus / Spaln / Genewise sections, provided file",
            "interpretation": "Различия длины и флагов помогают выбрать локусы для проверки.",
            "does_not_prove": "Появление/потерю экзона или подтверждённый сдвиг рамки.",
        },
    ]
    if not method_rows:
        evidence = [e for e in evidence if e["id"] != "three_independent_sections"]
    return evidence


def make_case_gene_result(report, provenance):
    """Transparent conclusion: no definitive mechanism from these inputs alone."""
    gene_id = report["gene_id"]
    is_user = report.get("input_type") == "user_supplied_fasta_sequences"
    evidence = evidence_for_case(report)
    missing = [
        "Независимый поиск белковых и геномных гомологов за пределами заранее отобранных Cactus-хитов.",
        "Проверенная синтения, версии геномных сборок и статус кодирования соответствующих локусов.",
        "Реконструкция предковой ORF и оценка изменения рамки по ортологичным позициям.",
        "Проверка паралогичных копий и датировки событий на дереве видов.",
        "Выровненные подтверждённые ортологи/паралоги, ML-дерево и поддержка его ключевых ветвей.",
    ]
    observed = {
        "gene_id": gene_id,
        "focal_species": report.get("focal_species", "dmel"),
        "reference_sequence_length_nt": report.get("target_length_nt"),
        "selected_method": report.get("method"),
        "comparison_species_with_records": sorted({
            h.get("species_abbr") for h in report.get("hits", [])
            if h.get("species_abbr") and h.get("species_abbr") != "dmel"
        }),
        "selected_relative_species": report.get("requested_species"),
    }
    stages = [
        {"stage": "input", "state": "completed",
         "observation": (f"Выбран ген {gene_id} из FASTA пользователя" if is_user
                         else f"Выбран ген {gene_id} из набора Cactus/Spaln/Genewise")},
        {"stage": "comparative_sequence_retrieval", "state": "completed",
         "observation": f"Отобрано {len(report.get('hits', []))} сравнительных записей"},
        {"stage": "local_dna_alignment", "state": "partial" if
         (report.get("summary") or {}).get("comparison_skipped") else "completed",
         "observation": "Измерены локальная идентичность и покрытие; отдельные слишком большие сравнения могли быть пропущены"},
        {"stage": "cross_method_comparison", "state": ("not_assessed" if is_user else "completed"),
         "observation": ("В пользовательском FASTA нет трёх опубликованных методов"
                         if is_user else
                         "Сопоставлены Raw Cactus, Spaln и Genewise без объединения последовательностей")},
        {"stage": "protein_homology_validation", "state": "not_assessed",
         "observation": "Независимый BLASTp/TBLASTN/HMMER-поиск не выполнялся"},
        {"stage": "ancestral_state_reconstruction", "state": "not_assessed",
         "observation": "Предковые ORF, дупликации, потери и рамочные инделы не реконструированы"},
        {"stage": "gene_tree_inference", "state": "not_assessed",
         "observation": "ML-дерево и поддержка ветвей не рассчитаны"},
        {"stage": "origin_conclusion", "state": "completed",
         "observation": "По доступным свидетельствам механизм остаётся неопределённым"},
    ]
    conclusion = {
        "category": "insufficient_data",
        "label": "Недостаточно данных для установления механизма происхождения",
        "confidence": "not_assessed",
        "not_a_probability": True,
        "why": (
            ("Пользователь предоставил кандидатные гены и сравнительные ДНК-последовательности, "
             "но для определения происхождения не хватает проверенных ортологов, "
             "геномного окружения, предкового кодирующего состояния и дерева генов."
             if is_user else
             "В файле присутствуют заранее отобранные сравнительные участки, "
             "но отсутствует независимая проверка геномной и белковой гомологии, "
             "предковые кодирующие состояния и надёжное дерево генов. "
             "Данные позволяют исследовать кандидатов, но не доказать конкретный механизм.")
        ),
        "hypotheses_not_resolved": ["de_novo", "duplication", "ancient_origin",
                                     "other_mechanism"],
        "unsupported_claims": [
            "Отсутствие записи у вида означает эволюционную потерю.",
            "Метаданное shift_Yes доказывает исторический frameshift.",
            "Предсказание Spaln или Genewise доказывает трансляцию.",
            "Наличие геномного совпадения доказывает ортологию.",
        ],
    }
    criteria = [
        {"number": 1, "title": "Приём генов и сравнительных данных", "status": "met",
         "evidence": ("Приняты отдельные кандидатные гены и сравнительные FASTA"
                      if is_user else
                      "FASTA разделён на целевую последовательность dmel и сравнения других видов")},
        {"number": 2, "title": "Автоматизированный анализ нескольких генов", "status": "met_in_batch",
         "evidence": "Несколько генов обрабатываются отдельно с объединённым отчётом"},
        {"number": 3, "title": "Механизм происхождения либо недостаточность данных", "status": "met_with_uncertainty",
         "evidence": "insufficient_data без выдуманной вероятности происхождения"},
        {"number": 4, "title": "Конкретные свидетельства", "status": "met",
         "evidence": "Количество записей, виды, измеренные совпадения и метаданные трёх методов"},
        {"number": 5, "title": "Прослеживаемость и воспроизводимость", "status": "met_for_screen",
         "evidence": "Хеш входного файла, параметры методов и последовательность выполненных стадий"},
        {"number": 6, "title": "Разделение наблюдений и эволюционных гипотез",
         "status": "met",
         "evidence": "Наблюдения приведены отдельно от альтернатив и недостаточных данных"},
    ]
    return {
        "schema": SCHEMA, "gene_id": gene_id,
        "source": provenance,
        "input": observed,
        "conclusion": conclusion,
        "evidence": evidence,
        "analysis_steps": stages,
        "case_requirements": criteria,
        "missing_for_mechanistic_inference": missing,
        "limitations": report.get("limitations", []),
        "phylogeny": report.get("phylogeny", {"status": "not_assessed"}),
    }


def analyze_case_batch(path, genes, method="raw", selected_species=None, source_name=None):
    """Per-gene reports in one auditable, non-falsified batch.

    Three complete scans of the research file; each selected gene then uses a
    small extracted three-section input instead of rescanning the 54 MB source.
    """
    from .case_data import METHODS

    if method not in METHODS:
        raise ValueError("Unknown case method")
    genes = [str(g).strip() for g in genes]
    if not genes or len(set(genes)) != len(genes):
        raise ValueError("Select at least one unique candidate gene")
    if any(not re.fullmatch(r"FBgn[0-9]+", g) for g in genes):
        raise ValueError("Every candidate must have a valid FBgn gene identifier")
    relatives = None if selected_species is None else list(selected_species)
    if relatives is not None and (not relatives or len(set(relatives)) != len(relatives)):
        raise ValueError("Select at least one unique comparison species")
    path = Path(path)
    provenance = {
        "filename": source_name or path.name,
        "sha256": hash_file(path),
        "size_bytes": path.stat().st_size,
        "dataset": "Peng & Zhao supplementary preselected Cactus/Spaln/Genewise extracts",
        "method": method, "requested_species": relatives,
        "software": _software_versions(),
    }
    result = []
    raw_reports, per_gene_errors = analyze_case_many(
        path, genes, method=method, selected_species=relatives
    )
    errors = [
        {"gene_id": gene_id, "status": "failed", "reason": message}
        for gene_id, message in per_gene_errors.items()
    ]
    for gene_id in genes:
        if gene_id in raw_reports:
            result.append(make_case_gene_result(raw_reports[gene_id], provenance))
    category_counts = dict(Counter(g["conclusion"]["category"] for g in result))
    return {
        "schema": SCHEMA, "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "requested_gene_ids": genes,
        "source": provenance, "gene_count": len(result),
        "results": result, "errors": errors,
        "raw_reports": raw_reports,
        "summary": {
            "requested": len(genes), "completed": len(result),
            "failed": len(errors),
            "origin_categories": category_counts,
            "interpretation": (
                "Результаты описывают сравнительные последовательности и "
                "отсутствие подтверждения конкретного происхождения; "
                "это не биологически валидированная классификация."
            ),
        },
    }


def batch_as_markdown(batch):
    """User-readable 6-requirement report, suitable for presentation/export."""
    lines = [
        "# EvoGene — отчёт по кейсу № 2",
        "", "## Исходные данные и воспроизводимость",
        f"- Файл: {batch['source']['filename']}",
        f"- SHA-256: `{batch['source']['sha256']}`",
        f"- Метод анализа: {batch['source']['method']}",
        f"- Запрошено генов: {batch['summary']['requested']}",
        f"- Обработано: {batch['summary']['completed']}; ошибок: {batch['summary']['failed']}",
        "",
    ]
    for entry in batch["results"]:
        lines.extend([
            f"## {entry['gene_id']}",
            "", "**Вывод:** " + entry["conclusion"]["label"],
            "", "**Основание:** " + entry["conclusion"]["why"],
            "",
            "### Наблюдаемые данные",
        ])
        for e in entry["evidence"]:
            lines.append(f"- **{e['description']}**: " +
                         json.dumps(e["value"], ensure_ascii=False) +
                         f". Не доказывает: {e['does_not_prove']}")
        lines.append("")
        lines.append("### Этапы анализа")
        for step in entry["analysis_steps"]:
            lines.append(f"- {step['stage']} [{step['state']}]: {step['observation']}")
        lines.extend(["", "### Вывод, а не наблюдение", "",
                      "Допустимые альтернативы: " +
                      ", ".join(entry["conclusion"]["hypotheses_not_resolved"]),
                      "Статистические вероятности механизмов не рассчитывались.", "",
                      "### Соответствие шести требованиям кейса", ""])
        for criterion in entry["case_requirements"]:
            lines.append(
                f"- **{criterion['number']}. {criterion['title']}** "
                f"[{criterion['status']}]: {criterion['evidence']}"
            )
        lines.extend(["", "### Что необходимо дополнительно", ""])
        for missing in entry["missing_for_mechanistic_inference"]:
            lines.append("- " + missing)
        lines.append("")
    for error in batch["errors"]:
        lines.append(f"- **Ошибка {error['gene_id']}**: {error['reason']}")
    return "\n".join(lines) + "\n"
