"""Transparent, mechanism-agnostic evidence screens.

No forced percentages, no trained labels, no gene-name lookups.
A missing measurement is not negative evidence. Conservation is separate.
"""
from .origin_tests import screen as legacy_screen

LABELS = {
    "de_novo": "De novo",
    "dna_duplication": "ДНК-дупликация",
    "retroposition": "Ретропозиция",
    "transposition": "Транспозиция / релокация",
    "fusion": "Слияние / химерный ген",
}
KEYS = tuple(LABELS)


def _test(key, signals, limitations, observations, assessed, partial=False):
    if signals:
        status = "candidate" if not partial else "partial_signal"
    else:
        status = "no_signal" if assessed else "not_assessed"
    return {
        "mechanism": key,
        "label": LABELS[key],
        "status": status,
        "signals": list(signals),
        "limitations": list(limitations),
        "observations": observations,
        "assessed": bool(assessed),
        "signal_count": len(signals),
        # Deliberately no probability / support_percent / invented confidence.
        "interpretation": "; ".join(signals) if signals else
        ("Специфический признак не выявлен проведёнными проверками" if assessed
         else "Не хватает данных для оценки"),
    }


def classify(report):
    ancestry = legacy_screen(report)
    synteny = report.get("synteny_evidence") or {}
    protein = report.get("protein_evidence") or {}
    parent = report.get("retrocopy_evidence") or {}
    target = report.get("target_locus") or report.get("target") or {}
    model = parent.get("target_model") or {}
    focal_ok = bool(report.get("focal_reference"))
    synteny_assessed = any(
        item.get("status") == "flanks_collinear" for item in synteny.values()
    )
    annotated_target = model.get("status") == "annotated"
    matched = ancestry["ancestral_conservation"]["observations"]["dna_match_species"]
    coding = ancestry["ancestral_conservation"]["observations"]["cds_overlap_species"]
    noncoding = ancestry["de_novo"]["observations"]["syntenic_noncds_species"]
    displaced = ancestry["transposition"]["observations"]["putative_displaced_homolog_species"]
    mosaic = ancestry["other_mechanism"]["observations"]["split_homology_species"]
    extras = report.get("extra_target_copies") or []
    prot_matches = protein.get("candidates") or []
    prot_paralogs = [
        c for c in prot_matches
        if c.get("gene_model") and not c.get("overlaps_target")
        and c.get("identity", 0) >= .35 and
        c.get("query_coverage", 0) >= .50
    ]
    cds_parent_candidates = [
        c for c in prot_paralogs
        if c["gene_model"].get("cds_exons", 0) >= 2
    ]
    dna_parent_candidates = [
        c for c in parent.get("parent_candidates", [])
        if c.get("parent_multi_exon_and_target_single_exon")
    ]
    # Intron-containing targets may be impossible to translate in a single ORF.
    # Failure to identify the parent is a sensitivity limitation, not evidence
    # that parent proteins do not exist.
    retro_signals = []
    if annotated_target and model.get("exon_count") == 1:
        if dna_parent_candidates:
            retro_signals.append("Одноэкзонный CDS и многоэкзонный гомолог в референсе (ДНК)")
        if cds_parent_candidates:
            retro_signals.append("Одноэкзонный CDS и белковый гомолог с несколькими CDS-экзонами")
    elif cds_parent_candidates or dna_parent_candidates:
        retro_signals.append("Кандидат родительского гена найден, но структура CDS потомка неизвестна")

    duplicated = []
    if extras:
        duplicated.append(f"Дополнительных сходных ДНК-локусов: {len(extras)}")
    if prot_paralogs:
        duplicated.append(f"Белковых гомологов в другом локусе: {len(prot_paralogs)}")

    denovo_signals = []
    if noncoding:
        denovo_signals.append(
            f"Гомологичная ДНК без аннотированного CDS в синтенном окружении у {len(noncoding)} видов"
        )
    # One or two non-CDS homologs are suggestive but cannot prove that
    # an ancestral ORF was absent. Do not hard-require two species or silence
    # a weak signal because coding homologs exist elsewhere.
    relocation_signals = ([f"Похожие последовательности вне синтенного интервала у {len(displaced)} видов"]
                          if displaced else [])
    fusion_signals = ([f"Независимые совпадения разных участков запроса у {len(mosaic)} видов"]
                      if mosaic else [])

    items = [
        _test("de_novo", denovo_signals,
              ["Отсутствие аннотированного CDS не доказывает некодирующее прошлое",
               "Нужны реконструкция предковой ORF и поиск удалённых белковых гомологов"],
              {"syntenic_no_cds_species": noncoding, "coding_overlap_species": coding},
              synteny_assessed, partial=True),
        _test("dna_duplication", duplicated,
              ["ДНК-повторы и белковая гомология не определяют направление дупликации",
               "Нужны проверка геномного контекста и филогения паралогов"],
              {"extra_genomic_loci": extras[:10], "protein_paralogs": prot_paralogs[:10]},
              focal_ok and bool(target), partial=bool(prot_paralogs and not extras)),
        _test("retroposition", retro_signals,
              ["Необходимо независимо подтвердить родство и направление копирования",
               "Индивидуальный FASTA гена с интронами может не переводиться в полноценный белок",
               "Структура CDS и отсутствие интронов требуют аннотации"],
              {"target_cds": model, "dna_parent_candidates": dna_parent_candidates[:10],
               "protein_parent_candidates": cds_parent_candidates[:10],
               "protein_search_status": protein.get("status", "not_run")},
              annotated_target and protein.get("status") == "searched",
              partial=not (annotated_target and model.get("exon_count") == 1)),
        _test("transposition", relocation_signals,
              ["Смещение гомолога не доказывает транспозицию: возможны геномные перестройки"],
              {"displaced_species": displaced},
              synteny_assessed),
        _test("fusion", fusion_signals,
              ["Сборочные артефакты, повторы и перестройки могут имитировать слияние"],
              {"independent_query_segment_species": mosaic},
              bool(report.get("relatives"))),
    ]
    candidates = [x["mechanism"] for x in items if x["status"] in ("candidate", "partial_signal")]
    not_assessed = [x["mechanism"] for x in items if x["status"] == "not_assessed"]
    errors = (report.get("download_or_analysis_errors") or []) + (
        report.get("origin_screen_errors") or [])
    return {
        "kind": "independent_evidence_matrix_v4",
        "is_calibrated_probability": False,
        "options": items, "tests": items,
        "candidate_mechanisms": candidates,
        "leading_candidate": candidates[0] if len(candidates) == 1 else None,
        "completeness": {"assessed": len(items) - len(not_assessed),
                         "total": len(items), "not_assessed": not_assessed,
                         "processing_errors": errors},
        "conservation": {
            "dna_homology_species": matched,
            "cds_overlap_species": coding,
            "syntenic_coding_species": ancestry["ancestral_conservation"]["observations"]["syntenic_coding_species"],
            "explanation": "Межвидовая сохранность НЕ является механизмом возникновения.",
        },
        "conclusion": (
            "Обнаружены признаки нескольких механизмов; нужен разбор конкурирующих объяснений."
            if len(candidates) > 1 else
            "Есть предварительные признаки: " + LABELS[candidates[0]] + ". Нужна независимая проверка."
            if candidates else "Специфические признаки не установлены; это не подтверждение отсутствия механизмов."
        ),
        "warning": "Проценты удалены: нет валидированной вероятностной модели. Каждый механизм оценивается отдельно.",
    }
