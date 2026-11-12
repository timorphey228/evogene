"""Conservative normalized evidence support (NOT calibrated probabilities)."""
LABELS={
"de_novo":"De novo",
"duplication":"Дупликация",
"transposition":"Транспозиция / ретротранспозиция",
"ancestral_conservation":"Древнее сохранённое происхождение",
"other_mechanism":"Другой механизм",
"unresolved":"Недостаточно данных"}
ORDER=tuple(LABELS)

def distribution(report):
    """Normalize only observed evidence; missing evidence is not negative evidence."""
    relatives=report.get("relatives") or {}
    passed={species:[h for h in hits if h.get("query_coverage",0)>=.5 and h.get("identity_approx",0)>=.65]
            for species,hits in relatives.items()}
    matched=[s for s,hits in passed.items() if hits]
    coding=[s for s,hits in passed.items() if any(h.get("coding_overlap")=="overlaps_cds" for h in hits)]
    errors=len(report.get("download_or_analysis_errors",[]))
    extra=report.get("extra_target_copies",[])
    w={k:0 for k in ORDER}
    w["unresolved"]=8+2*errors
    w["ancestral_conservation"]=min(7,len(matched))+min(6,2*len(coding))
    if extra:w["duplication"]=min(8,2*len(extra))
    # No syntenic ancestry or transposition evidence is measured by current pipeline.
    total=sum(w.values())
    raw={k:1000*w[k]/total for k in ORDER}
    units={k:int(raw[k]) for k in ORDER}
    ranked=sorted(ORDER,key=lambda k:(-(raw[k]-units[k]),ORDER.index(k)))
    for k in ranked[:1000-sum(units.values())]:units[k]+=1
    explanations={
        "de_novo":"Не проверен синтенный предковый некодирующий участок.",
        "duplication":("Обнаружены дополнительные локусы, но паралогию нужно подтвердить."
                       if extra else "Систематический поиск паралогов не выполнен."),
        "transposition":"Нет данных о перемещении локуса или ретрокопии.",
        "ancestral_conservation":f"Сходные ДНК-участки у {len(matched)} видов, перекрытие CDS у {len(coding)}.",
        "other_mechanism":"Нет специальных свидетельств иных механизмов.",
        "unresolved":f"Не выполнены синтения, белковые деревья и проверка предковой ORF; ошибок сравнений: {errors}."}
    return {
        "kind":"normalized_evidence_support_not_probability",
        "is_calibrated_probability":False,
        "scoring_method":"conservative_rule_weights_v1",
        "comparison_species_analyzed":len(relatives),
        "comparison_species_failed":errors,
        "options":[{"mechanism":k,"label":LABELS[k],"support_percent":units[k]/10,
                    "weight":w[k],"explanation":explanations[k]}
                   for k in ORDER],
        "conclusion":("Обнаружено межвидовое сходство ДНК, но конкретный механизм происхождения не установлен."
                      if matched else "Недостаточно данных для установления механизма происхождения."),
        "warning":("Проценты — нормированные веса свидетельств, НЕ реальные вероятности. "
                   "0% означает, что диагностические данные не проверены, а не что механизм исключён.")
    }
