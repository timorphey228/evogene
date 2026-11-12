"""Evidence-first Streamlit visualizations for the ITMO case 2.

Plots only measured coordinates or user-supplied sequence alignments / trees.
A mapped DNA segment is NOT an evolutionary event or an inferred gene tree.
Rendering is optional: the core analysis remains usable without Streamlit.
"""
from __future__ import annotations

from io import StringIO
from collections import Counter
import json


def alignment_blocks(report, limit=400):
    """Query-coordinate PAF blocks. Not a base-by-base multiple alignment."""
    rows = []
    for species, hits in (report.get("relatives") or {}).items():
        for hit in hits:
            start = hit.get("query_start")
            end = hit.get("query_end")
            if start is None or end is None or end <= start:
                continue
            rows.append({
                "species": species, "start": int(start), "end": int(end),
                "chrom": str(hit.get("chrom", "?")),
                "genome_start": int(hit.get("start", 0)),
                "genome_end": int(hit.get("end", 0)),
                "strand": hit.get("strand", "?"),
                "identity": round(100 * float(hit.get("identity_approx", 0)), 2),
                "coverage": round(100 * float(hit.get("query_coverage", 0)), 2),
                "annotation": hit.get("coding_overlap", "annotation_unknown"),
                "mapq": hit.get("mapq"),
                "alignment_insertion_nt": hit.get("alignment_insertion_nt"),
                "alignment_deletion_nt": hit.get("alignment_deletion_nt"),
            })
    return sorted(rows, key=lambda x: (x["species"], x["start"], -x["coverage"]))[:limit]


def synteny_blocks(report):
    """Flank locations are shown relative to their *observed* genomic window."""
    rows = []
    for species, result in (report.get("synteny_evidence") or {}).items():
        if result.get("status") != "flanks_collinear":
            continue
        left = result.get("left_anchor")
        right = result.get("right_anchor")
        if not left or not right:
            continue
        origin = min(left["start"], right["start"])
        for key, hit in (("left flank", left), ("right flank", right)):
            rows.append({
                "species": species, "feature": key,
                "start": hit["start"] - origin, "end": hit["end"] - origin,
                "chrom": hit.get("chrom", result.get("chrom", "?")),
                "genomic_start": hit["start"], "genomic_end": hit["end"],
                "status": result.get("locus_status", "unknown"),
            })
        for hit in result.get("syntenic_hits", [])[:20]:
            rows.append({
                "species": species, "feature": "candidate DNA homolog",
                "start": hit["start"] - origin, "end": hit["end"] - origin,
                "chrom": hit.get("chrom", result.get("chrom", "?")),
                "genomic_start": hit["start"], "genomic_end": hit["end"],
                "status": result.get("locus_status", "unknown"),
            })
    return rows


def orf_blocks(report):
    """Reading frames of the submitted contiguous DNA, not spliced CDS."""
    rows = []
    for item in report.get("query_orfs", []):
        rows.append({
            "frame": f"{item['strand']}{item['frame'] + 1}",
            "start": int(item["start"]), "end": int(item["end"]),
            "aa_length": int(item["aa_length"]), "strand": item["strand"],
        })
    return rows


def alignment_gap_rows(report):
    """CIGAR I/D base counts are mapper operations, not inferred mutations."""
    rows = []
    for species, hits in (report.get("relatives") or {}).items():
        for hit in hits[:20]:
            for feature, key in (
                ("I (query relative to reference)", "alignment_insertion_nt"),
                ("D (query relative to reference)", "alignment_deletion_nt"),
            ):
                count = hit.get(key)
                if count is not None:
                    rows.append({
                        "species": species, "operation": feature, "bases": int(count),
                        "chrom": hit.get("chrom", "?"), "mapq": hit.get("mapq"),
                    })
    return rows


def mechanism_rows(report):
    m = report.get("mechanism_report") or {}
    return [{
        "mechanism": t.get("label", t.get("mechanism", "")),
        "key": t.get("mechanism", ""),
        "status": t.get("status", "not_assessed"),
        "signals": t.get("signal_count", 0),
        "assessed": bool(t.get("assessed", False)),
        "interpretation": t.get("interpretation", ""),
    } for t in m.get("tests", [])]


def parse_aligned_fasta(text, *, max_sequences=80, max_columns=20000):
    """Validate a real aligned FASTA; do not align unaligned sequences silently."""
    from Bio import SeqIO
    sequences = [(record.id, str(record.seq).upper()) for record in SeqIO.parse(StringIO(text), "fasta")]
    if not 2 <= len(sequences) <= max_sequences:
        raise ValueError(f"Expected 2–{max_sequences} aligned sequences")
    length = len(sequences[0][1])
    if not length or length > max_columns or any(len(seq) != length for _, seq in sequences):
        raise ValueError("FASTA must contain equally long aligned sequences (with gaps)")
    if any(ch not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ-.X*" for _, seq in sequences for ch in seq):
        raise ValueError("Unexpected residue/character in alignment")
    if len({name for name, _ in sequences}) != len(sequences):
        raise ValueError("Sequence identifiers must be unique")
    return sequences


def alignment_cells(sequences, start=0, width=60, max_sequences=30):
    if not sequences:
        return []
    end = min(len(sequences[0][1]), start + width)
    return [
        {"sequence": name, "column": i + 1, "base": seq[i],
         "kind": "gap" if seq[i] in "-." else seq[i]}
        for name, seq in sequences[:max_sequences] for i in range(start, end)
    ]


def newick_to_dot(newick, *, max_tips=200):
    """Draw the provided Newick topology; support is shown only if supplied."""
    from Bio import Phylo
    tree = Phylo.read(StringIO(newick), "newick")
    tips = tree.get_terminals()
    if not tips or all(t.name is None for t in tips) or len(tips) > max_tips:
        raise ValueError(f"Newick needs 1–{max_tips} tips")
    lines = ["digraph Tree {", 'graph [rankdir=LR, bgcolor="transparent"];',
             'node [fontname="Arial",fontsize=10];', 'edge [color="#94A3B8"];']
    counter = [0]

    def visit(node):
        node_id = f"n{counter[0]}"
        counter[0] += 1
        if node.is_terminal():
            label = node.name or "unnamed tip"
            shape = "box"
        else:
            # Confidence belongs to the node in Newick; its statistical
            # interpretation depends on the program that built the tree.
            label = (str(node.name) if node.name else "") if node.confidence is None else f"{node.confidence:g}"
            shape = "circle"
        lines.append(f'{node_id} [label={json.dumps(str(label))},shape={shape}];')
        for child in node.clades:
            child_id = visit(child)
            length = "" if child.branch_length is None else f"{child.branch_length:.4g}"
            lines.append(
                f'{node_id} -> {child_id} [label={json.dumps(length)}];'
            )
        return node_id

    visit(tree.root)
    lines.append("}")
    return "\n".join(lines)


def render_report_visualizations(report):
    """Interactive evidence panels. No inferred phylogeny is invented."""
    import streamlit as st
    import altair as alt

    st.header("Визуальный анализ доказательств")
    st.caption("Изображены только результаты измерений или отдельно загруженные данные. "
               "DNA-сходство не доказывает ортологию или эволюционное происхождение.")
    tab_align, tab_synteny, tab_orf, tab_gap, tab_mechanisms, tab_tree = st.tabs([
        "Выравнивания", "Синтения", "ORF", "Вставки/делеции", "Гипотезы и неопределённость", "Филогения"
    ])
    with tab_align:
        hits = alignment_blocks(report)
        if hits:
            st.caption("Покрытие координат запроса из PAF/minimap2; это не посимвольное MSA. "
                       "Наведи курсор на блок для идентичности, MAPQ и координат.")
            chart = alt.Chart(alt.Data(values=hits)).mark_bar(size=12).encode(
                x=alt.X("start:Q", title="Положение на ДНК-запросе, нт"),
                x2="end:Q",
                y=alt.Y("species:N", title="Вид"),
                color=alt.Color("annotation:N", title="CDS-аннотация"),
                tooltip=["species:N", "chrom:N", "genome_start:Q", "genome_end:Q",
                         "strand:N", "identity:Q", "coverage:Q", "mapq:Q", "annotation:N"],
            ).properties(height=max(150, 34 * len({h["species"] for h in hits}))).interactive()
            st.altair_chart(chart, use_container_width=True)
            st.dataframe(hits, use_container_width=True, hide_index=True)
        elif report.get("hits"):
            st.caption("Предварительно отобранные Cactus-фрагменты: локальная "
                       "идентичность и покрытие. Это не независимый поиск "
                       "гомологов в целых геномах.")
            rows = [{
                "species": hit.get("species_abbr", "?"),
                "identity": round(100 * float(hit.get("identity_aligned", 0)), 2),
                "coverage": round(100 * float(hit.get("query_coverage", 0)), 2),
                "kind": hit.get("header_kind", "unknown"),
                "header": hit.get("header", ""),
            } for hit in report["hits"]]
            chart = alt.Chart(alt.Data(values=rows)).mark_circle(size=95).encode(
                x=alt.X("coverage:Q", title="Покрытие запроса, %", scale=alt.Scale(domain=[0, 100])),
                y=alt.Y("identity:Q", title="Идентичность выровненных позиций, %",
                        scale=alt.Scale(domain=[0, 100])),
                color=alt.Color("kind:N", title="Тип Cactus-записи"),
                tooltip=["species:N", "header:N", "coverage:Q", "identity:Q"],
            ).properties(height=300).interactive()
            st.altair_chart(chart, use_container_width=True)
            st.dataframe(rows, use_container_width=True, hide_index=True)
        else:
            st.info("Блоков выравнивания нет или формат не содержит координат запроса.")
        st.markdown("**Посимвольное множественное выравнивание (отдельный файл)**")
        msa_upload = st.file_uploader("Загрузить готовое выравнивание FASTA с пропусками",
                                      type=["fa", "fasta", "fas", "aln"], key="evo_msa")
        if msa_upload is not None:
            try:
                sequences = parse_aligned_fasta(msa_upload.getvalue().decode("utf-8-sig"))
                length = len(sequences[0][1])
                first = st.slider("Начальная колонка", 1, length, 1, key="msa_offset")
                width = min(70, length - first + 1)
                cells = alignment_cells(sequences, first - 1, width)
                grid = alt.Chart(alt.Data(values=cells)).mark_rect(stroke="#E2E8F0").encode(
                    x=alt.X("column:O", title="Колонка выравнивания"),
                    y=alt.Y("sequence:N", title="Последовательность",
                            sort=[s[0] for s in sequences[:30]]),
                    color=alt.Color("kind:N", title="Нуклеотид/пропуск"),
                    tooltip=["sequence:N", "column:O", "base:N"],
                ).properties(height=max(110, len(sequences[:30]) * 21))
                st.altair_chart(grid, use_container_width=True)
                st.caption(f"Показаны {width} колонок и не более 30 последовательностей; "
                           "файл должен быть действительно выровнен вне приложения.")
            except (UnicodeDecodeError, ValueError) as exc:
                st.error(f"Не удалось прочитать MSA: {exc}")
    with tab_synteny:
        blocks = synteny_blocks(report)
        if blocks:
            st.caption("Расположение фланков и возможных гомологов по координатам "
                       "генома каждого вида; нулевая позиция отдельная для каждого локуса.")
            chart = alt.Chart(alt.Data(values=blocks)).mark_bar(size=14).encode(
                x=alt.X("start:Q", title="Смещение от начала локального окна, нт"),
                x2="end:Q",
                y=alt.Y("species:N", title="Вид"),
                color=alt.Color("feature:N", title="Геномный элемент"),
                tooltip=["species:N", "feature:N", "chrom:N", "genomic_start:Q",
                         "genomic_end:Q", "status:N"],
            ).properties(height=max(150, 32 * len({b["species"] for b in blocks}))).interactive()
            st.altair_chart(chart, use_container_width=True)
        else:
            st.info("Для схемы нужны координаты обоих сохранённых флангов. "
                    "Синтения может оставаться неопределённой.")
        statuses = [
            {"species": name, "status": value.get("status", "unknown"),
             "locus": value.get("locus_status", "unknown"),
             "interval": value.get("interval_0based")}
            for name, value in (report.get("synteny_evidence") or {}).items()
        ]
        if statuses:
            st.dataframe(statuses, use_container_width=True, hide_index=True)
    with tab_orf:
        rows = orf_blocks(report)
        if rows:
            chart = alt.Chart(alt.Data(values=rows)).mark_bar(size=12).encode(
                x=alt.X("start:Q", title="Координата исходной ДНК, 0-based"),
                x2="end:Q", y=alt.Y("frame:N", title="Цепь и рамка"),
                color=alt.Color("strand:N", title="Цепь"),
                tooltip=["frame:N", "start:Q", "end:Q", "aa_length:Q"],
            ).properties(height=230).interactive()
            st.altair_chart(chart, use_container_width=True)
            st.dataframe(rows, use_container_width=True, hide_index=True)
        else:
            st.info("ORF не были найдены при заданном пороге или анализ не выполнялся.")
        st.caption(report.get("orf_note",
            "ORF показывает потенциальное кодирование непрерывной ДНК. "
            "Не является подтверждением транскрипции, трансляции или предкового ORF."))
    with tab_gap:
        rows = alignment_gap_rows(report)
        st.caption("Это операции I и D в CIGAR выравнивания minimap2. "
                   "Они не являются сами по себе эволюционными инделами, "
                   "доказательствами потери гена или сдвига рамки считывания. "
                   "Для этого нужна реконструкция ортологичной CDS и предкового состояния.")
        if rows:
            chart = alt.Chart(alt.Data(values=rows)).mark_bar().encode(
                x=alt.X("species:N", title="Вид"),
                y=alt.Y("sum(bases):Q", title="Сумма нуклеотидов в операциях"),
                color=alt.Color("operation:N", title="CIGAR"),
                tooltip=["species:N", "operation:N", "sum(bases):Q"],
            ).properties(height=280)
            st.altair_chart(chart, use_container_width=True)
            st.dataframe(rows, use_container_width=True, hide_index=True)
        else:
            st.info("CIGAR операции не предоставлены. Геномные инделы не реконструированы.")
    with tab_mechanisms:
        rows = mechanism_rows(report)
        if rows:
            st.caption("Цвет и статус отражают выполненные диагностические проверки, "
                       "но НЕ вычисленные вероятности происхождения.")
            state_order = ["candidate", "partial_signal", "no_signal", "not_assessed"]
            color = alt.Scale(domain=state_order,
                range=["#16815D", "#D29D3C", "#8394A8", "#CBD5E1"])
            chart = alt.Chart(alt.Data(values=rows)).mark_rect(
                stroke="white", strokeWidth=1).encode(
                x=alt.X("status:N", sort=state_order, title="Статус проверки"),
                y=alt.Y("mechanism:N", title="Механизм"),
                color=alt.Color("status:N", scale=color, legend=None),
                tooltip=["mechanism:N", "status:N", "signals:Q", "interpretation:N"],
            ).properties(height=max(160, 38 * len(rows)))
            st.altair_chart(chart, use_container_width=True)
            st.dataframe(rows, use_container_width=True, hide_index=True)
            completeness = (report.get("mechanism_report") or {}).get("completeness") or {}
            assessed = completeness.get("assessed", 0)
            total = completeness.get("total", 0)
            st.metric("Оценено диагностических механизмов", f"{assessed}/{total}")
            if completeness.get("not_assessed"):
                st.warning("Не оценены: " + ", ".join(completeness["not_assessed"]))
            if completeness.get("processing_errors"):
                st.error("Ошибки обработки: " + str(completeness["processing_errors"]))
        else:
            st.info("Матрица диагностических проверок отсутствует для этого режима.")
        events = report.get("pipeline_events") or []
        if events:
            final = {}
            for event in events:
                final[(event.get("stage"), event.get("species"))] = event
            st.markdown("**Фактически выполненные этапы**")
            st.dataframe(list(final.values()), use_container_width=True, hide_index=True)
    with tab_tree:
        st.caption("Дерево строится не из таблицы совпадений: загрузите Newick "
                   "после настоящей реконструкции (например, IQ-TREE). "
                   "Числа в узлах показываются без переименования в bootstrap, "
                   "поскольку их смысл задаёт программа-источник.")
        tree_file = st.file_uploader("Newick дерево (.nwk, .treefile, .newick)",
                                     type=["nwk", "newick", "treefile", "tree", "txt"],
                                     key="evo_tree")
        if tree_file is None:
            st.info("Филогенетическое дерево не рассчитано текущим анализом. "
                    "Нет основания рисовать вымышленную топологию или поддержку.")
        else:
            try:
                dot = newick_to_dot(tree_file.getvalue().decode("utf-8-sig"))
                st.graphviz_chart(dot, use_container_width=True)
                st.download_button("Экспортировать Newick", tree_file.getvalue(),
                    file_name="uploaded_gene_tree.nwk", mime="text/plain")
            except (UnicodeDecodeError, ValueError, TypeError, AttributeError, IndexError) as exc:
                st.error(f"Некорректное дерево: {exc}")
