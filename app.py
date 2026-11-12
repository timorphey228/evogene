"""Streamlit GUI for local comparative genomes and published case extract."""
import json
import tempfile
from pathlib import Path

import streamlit as st
from evogene.genomes import analyze as analyze_genomes
from evogene.case_data import (
    analyze as analyze_case, index_file, export_gene_fasta, METHODS,
)
from evogene.visualization import render_report_visualizations
from evogene.analysis_notes import render_analysis_annotations

st.set_page_config(page_title="EvoGene · ITMO case 2", layout="wide")
st.title("EvoGene · эволюционное происхождение белок-кодирующих генов")
st.warning(
    "ДНК-сходство или предсказание CDS не доказывают ортологию, синтению, "
    "дупликацию либо происхождение de novo."
)
mode = st.radio("Источник данных", [
    "Локальные геномы FASTA + GFF3",
    "Данные статьи Peng & Zhao: Cactus / Spaln / Genewise",
])


def save_upload(upload, folder):
    if upload is None:
        return None
    path = Path(folder) / (str(len(list(Path(folder).iterdir()))) +
                           "".join(Path(upload.name).suffixes))
    path.write_bytes(upload.getvalue())
    return str(path)


@st.cache_data(show_spinner="Индексирую три раздела файла без загрузки всех последовательностей в память")
def index_upload(data: bytes):
    with tempfile.NamedTemporaryFile(suffix=".fasta") as temp:
        temp.write(data)
        temp.flush()
        return index_file(temp.name)


if mode.startswith("Локальные"):
    st.caption("Здесь выполняется проверка двух вручную заданных сравнительных геномов. "
               "Для более широкого набора до 19 видов используйте auto_app.py; "
               "для полного Cactus-набора — второй режим.")
    gene = st.text_input("Идентификатор гена из GFF3", value="FBgn0040343")
    tf = st.file_uploader("Целевой геном FASTA", key="target_fasta")
    tg = st.file_uploader("GFF3 целевого генома", key="target_gff")
    s1 = st.text_input("Первый родственный вид", value="Drosophila simulans")
    r1 = st.file_uploader("FASTA первого вида", key="relative1")
    g1 = st.file_uploader("GFF3 первого вида (необязательно)", key="relative1_gff")
    s2 = st.text_input("Второй родственный вид", value="Drosophila yakuba")
    r2 = st.file_uploader("FASTA второго вида", key="relative2")
    g2 = st.file_uploader("GFF3 второго вида (необязательно)", key="relative2_gff")
    if st.button("Анализировать локальные сборки"):
        if not all((tf, tg, r1)):
            st.error("Нужны целевые FASTA/GFF3 и как минимум один сравнительный FASTA.")
        else:
            try:
                with tempfile.TemporaryDirectory() as folder:
                    relatives = [{
                        "species": s1,
                        "fasta": save_upload(r1, folder),
                        "gff": save_upload(g1, folder),
                    }]
                    if r2:
                        relatives.append({
                            "species": s2,
                            "fasta": save_upload(r2, folder),
                            "gff": save_upload(g2, folder),
                        })
                    st.session_state["evo_result"] = analyze_genomes(
                        save_upload(tf, folder),
                        save_upload(tg, folder), gene, relatives,
                    )
                    st.session_state.pop("case_export", None)
            except Exception as exc:
                st.error(f"Ошибка анализа: {exc}")
else:
    st.markdown(
        "Файл `supplementary_file_unannotated_alignments.fasta` содержит "
        "**три самостоятельных раздела**, а не одно готовое множественное "
        "выравнивание. Исходные Cactus-хиты, Spaln и Genewise сравниваются "
        "отдельно. 20 видов доступны без ограничения в пять."
    )
    file = st.file_uploader(
        "Загрузить supplementary_file_unannotated_alignments.fasta",
        type=["fasta", "fa", "fas", "txt"], key="published_extract",
    )
    if file is not None:
        raw = file.getvalue()
        try:
            index = index_upload(raw)
            st.caption(
                f"Обнаружено {len(index['genes'])} генов, "
                f"{len(index['species'])} обозначений видов. "
                f"Записей по методам: {index['record_counts']}. "
                "Не следует считать эти данные независимой полной проверкой генома."
            )
            default = index["genes"].index("FBgn0040343") if "FBgn0040343" in index["genes"] else 0
            gene = st.selectbox("Выбрать любой ген из файла", index["genes"], index=default)
            method = st.selectbox(
                "Тип последовательностей", options=list(METHODS),
                format_func=lambda key: METHODS[key][1],
            )
            chosen_species = st.multiselect(
                "Виды для сравнения (все 20 доступны)",
                options=index["species"],
                default=index["species"],
                help="Фильтр распространяется на таблицу совпадений и FASTA-экспорт. "
                     "Полнота по всем видам отображается отдельно.",
            )
            st.caption("Raw = геномная ДНК; Spaln/Genewise = предсказания "
                       "сплайсированных кодирующих фрагментов. Score/shift — "
                       "авторская метаинформация; без опубликованной калибровки "
                       "её нельзя переводить в p-value.")
            if st.button("Исследовать выбранный ген", type="primary",
                         disabled=not chosen_species):
                with tempfile.NamedTemporaryFile(suffix=".fasta") as temp:
                    temp.write(raw)
                    temp.flush()
                    try:
                        result = analyze_case(
                            temp.name, gene, method=method,
                            selected_species=chosen_species,
                        )
                        export = export_gene_fasta(
                            temp.name, gene, method,
                            species=set(chosen_species) | {"dmel"},
                        )
                        st.session_state["evo_result"] = result
                        st.session_state["case_export"] = export
                    except (OSError, ValueError, UnicodeError) as exc:
                        st.error(f"Не удалось обработать файл: {exc}")
        except (OSError, ValueError, UnicodeError) as exc:
            st.error(f"Не удалось индексировать файл: {exc}")
    else:
        st.info("Загрузите файл из исследования; встроенная демо-таблица "
                "не заменяет реальных последовательностей.")

result = st.session_state.get("evo_result")
if result is not None:
    st.divider()
    st.subheader(result["gene_id"] + " — " + result["hypothesis"])
    st.write(result["explanation"])
    if result.get("method_comparison"):
        st.subheader("Сравнение трёх разделов исследования")
        st.dataframe(result["method_comparison"],
                     use_container_width=True, hide_index=True)
        st.caption("Показатели описывают входные записи и не являются "
                   "вероятностями кодирования или происхождения.")
    render_report_visualizations(result)
    render_analysis_annotations(result)
    with st.expander("Первичные данные и научные ограничения"):
        if "hits" in result:
            st.dataframe(result["hits"], use_container_width=True, hide_index=True)
        else:
            st.json(result.get("relatives", {}))
        for item in result.get("limitations", []):
            st.write("• " + str(item))
    payload = json.dumps(result, ensure_ascii=False, indent=2)
    st.download_button(
        "Скачать JSON с доказательствами", payload,
        file_name=result["gene_id"] + "_evidence.json",
        mime="application/json",
    )
    if result.get("input_type") == "published_preselected_hits_not_whole_genome_or_msa":
        export = st.session_state.get("case_export")
        if export:
            st.download_button(
                "Скачать выбранные последовательности FASTA (НЕ выровнены)",
                export, file_name=result["gene_id"] + "_" + result["method"] + ".fasta",
                mime="text/plain",
            )
