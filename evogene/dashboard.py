"""Streamlit research dashboard for EvoGene candidate-gene projects.

Launch: streamlit run auto_app.py
The current pipeline is an evidence-first screening analogue of the 2024
paper. It never claims that the fallback NJ tree proves gene ancestry.
"""
from __future__ import annotations

from datetime import datetime
from io import BytesIO
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile
import json

import streamlit as st

from .evolutionary_workflow import candidate_records, run_project
from .ncbi_auto import PANEL, REFERENCE_ASSEMBLIES
from .tools import check_tools
from .visualization import newick_to_dot, render_report_visualizations
from .analysis_notes import render_analysis_annotations

STAGES = {
    "input": "Подготовка данных",
    "download": "Загрузка геномов",
    "dna_alignment": "Поиск геномных совпадений",
    "cds_annotation": "Проверка CDS",
    "target_reference": "Исходная сборка",
    "target_mapping": "Дополнительные копии",
    "gene_structure": "Экзоны и интроны",
    "protein_homology": "Белковое сходство",
    "synteny": "Синтения",
    "mechanism_tests": "Эволюционные гипотезы",
    "phylogeny": "Филогенетика",
    "transposable_elements": "Транспозонные элементы",
    "hypothesis_review": "Заключение",
}
STATUS = {"candidate": "Обнаружен диагностический признак",
          "partial_signal": "Предварительный признак",
          "no_signal": "При выполненной проверке признак не найден",
          "not_assessed": "Не проверено",
          "failed": "Проверка не выполнена"}


def _style():
    st.markdown("""
<style>
  .block-container {max-width: 1390px; padding-top: 1.8rem; padding-bottom: 5rem;}
  h1,h2,h3 {letter-spacing:-0.035em;}
  h1 {font-weight:750 !important;}
  [data-testid="stMetric"] {border:1px solid rgba(128,145,164,.24);
      background:rgba(124,138,164,.055);padding:1rem;border-radius:16px;}
  .evo-intro{padding:26px 30px;border:1px solid rgba(115,141,184,.28);
      border-radius:22px;background:linear-gradient(125deg,rgba(53,113,172,.16),rgba(68,157,124,.10));
      margin-bottom:20px;}
  .evo-overline{font-size:.75rem;letter-spacing:.12em;color:#3F9E84;
      font-weight:700;text-transform:uppercase;}
  .evo-intro h1{font-size:2.5rem;margin:6px 0 5px;}
  .evo-intro p{opacity:.75;margin:0;}
  .evo-phase{font-weight:650;letter-spacing:.02em;font-size:.9rem;
      padding:.65rem 0;border-bottom:1px solid rgba(120,140,155,.2);}
  .evo-stamp{border-left:3px solid #2D9B82;padding-left:1rem;margin:1rem 0;}
  div[data-testid="stTabs"] button {border-radius:8px;}
</style>
""", unsafe_allow_html=True)


def _bundle(report):
    stream = BytesIO()
    with ZipFile(stream,"w",compression=ZIP_DEFLATED) as archive:
        archive.writestr("run_report.json",json.dumps(report,ensure_ascii=False,indent=2,default=str))
        archive.writestr("audit.jsonl","\n".join(
            json.dumps(x,ensure_ascii=False) for x in report.get("audit_log",[]))+"\n")
        for gene_id, result in report.get("genes",{}).items():
            safe="".join(c if c.isalnum() or c in "_.-" else "_" for c in gene_id)
            if (result.get("phylogeny") or {}).get("newick"):
                archive.writestr(f"{safe}/gene_tree.nwk",result["phylogeny"]["newick"]+"\n")
            if (result.get("phylogeny") or {}).get("alignment_fasta"):
                archive.writestr(f"{safe}/msa.fasta",result["phylogeny"]["alignment_fasta"])
    return stream.getvalue()


def _input_panel():
    st.subheader("01 · Входные данные")
    st.caption("Предоставьте ДНК последовательности одного или нескольких кандидатов. "
               "Система скачает сравнительные геномы по указанным референсам и сохранит версии сборок.")
    left, right = st.columns([1.1,1],gap="large")
    with left:
        file = st.file_uploader("FASTA кандидатных генов",type=["fa","fna","fasta"],
                                help="Каждая FASTA-запись — отдельный кандидат. Заголовок >FBgn... или >gene_id.")
        st.caption("Пример: `>FBgn0040343` и далее его нуклеотидная последовательность. "
                   "Для генов с интронами допускается геномная последовательность, но поиск ORF будет предварительным.")
        if file is not None:
            try:
                entries=candidate_records(file.getvalue())
                st.success(f"Генов: {len(entries)}; суммарная длина: "
                           f"{sum(e['length_nt'] for e in entries):,} нт")
                st.dataframe([{"Ген":e["gene_id"],"Длина, нт":e["length_nt"]}
                              for e in entries],hide_index=True,use_container_width=True)
            except (ValueError,UnicodeError) as exc:
                st.error("FASTA: "+str(exc))
    with right:
        focal=st.selectbox("Исходный вид",list(PANEL),index=list(PANEL).index("Drosophila melanogaster"))
        default=[s for s in ("Drosophila simulans","Drosophila yakuba",
                            "Drosophila ananassae") if s in PANEL[focal]]
        chosen=st.multiselect("Геномы для сравнения",PANEL[focal],default=default,
                              help="Добавьте близкие виды и внешние группы. "
                                   "Наличие гомолога вне синтенного локуса проверяется отдельно.")
        with st.expander("Конкретные версии геномных сборок",expanded=False):
            st.dataframe([{"Вид":s,"Версия из статьи":REFERENCE_ASSEMBLIES.get(s,"не указана")}
                          for s in [focal]+chosen],hide_index=True,use_container_width=True)
            st.caption("NCBI Datasets может заменить недоступную сборку; это отдельно фиксируется в журнале. "
                       "Для D. melanogaster FlyBase r6.36 не является NCBI accession.")
    with st.expander("Дополнительные настройки",expanded=False):
        te=st.file_uploader("Библиотека транспозонов FASTA (необязательно)",
                           type=["fa","fna","fasta"],key="te_library",
                           help="Используется как скрининг сходства; RepeatMasker/Dfam этим не заменяются.")
        phylogeny=st.checkbox("Автоматически построить дерево из доступных последовательностей",
                             value=True,help="При наличии MAFFT/IQ-TREE: ML или NJ с поддержкой. "
                             "Без них предварительная NJ-кластеризация 7-мерами, без поддержки.")
        cache=st.text_input("Каталог кэша геномов",value="genomes/ncbi_cache",
                            help="Сборки скачиваются один раз и повторно используются.")
    return file,focal,chosen,te,phylogeny,cache


def _run(file,focal,chosen,te,phylogeny,cache):
    ready=check_tools()
    checks=st.columns(2)
    for col,name in zip(checks,("datasets","minimap2")):
        col.caption(f"{name}: {'готов' if ready[name]['ok'] else 'не найден'}")
    missing=[k for k in ("datasets","minimap2") if not ready[k]["ok"]]
    if missing:
        st.warning("Для автоматической загрузки и выравнивания установите: "+", ".join(missing))
    disabled=(not file or not chosen or bool(missing))
    if st.button("Запустить исследование",type="primary",use_container_width=True,disabled=disabled):
        st.session_state.pop("evo_project",None)
        with st.status("Выполняется анализ выбранных генов",expanded=True) as status:
            lines=st.empty()
            track=[]
            def progress(event):
                if event["state"] != "running":
                    track.append(event)
                    if len(track)>9:
                        track.pop(0)
                    lines.markdown("\n".join(
                        "- **"+(x.get("gene_id") or "")+"** · "+
                        STAGES.get(x["stage"],x["stage"])+" · "+
                        x["state"]+" — "+x.get("detail","")[:130]
                        for x in track
                    ))
            try:
                from tempfile import TemporaryDirectory
                with TemporaryDirectory(prefix="evogene_te_") as work:
                    te_path=None
                    if te:
                        te_path=Path(work)/"transposable_elements.fna"
                        te_path.write_bytes(te.getvalue())
                    report=run_project(
                        file.getvalue(),focal,chosen,cache=cache,
                        te_library=te_path,infer_phylogeny=phylogeny,
                        output_folder="results/evolutionary_pipeline",
                        progress=progress)
                st.session_state["evo_project"]=report
                status.update(label="Анализ завершён. Откройте результаты ниже.",
                              state="complete",expanded=False)
            except Exception as exc:
                status.update(label="Не удалось запустить исследование",state="error")
                st.error(str(exc))
    return


def _hypothesis_panel(entry):
    result=entry.get("hypotheses",{})
    st.subheader("Эволюционные гипотезы")
    st.caption("Диагностические наблюдения не являются статистически откалиброванной "
               "вероятностью происхождения. Категория не присваивается искусственно.")
    st.info("Итог: "+("предварительные кандидаты" if result.get("candidates") else
                       "недостаточно подтверждений для выбора механизма"))
    for h in result.get("hypotheses",[]):
        with st.expander(h["label"]+" · "+STATUS.get(h["status"],h["status"]),
                         expanded=h["status"] in ("candidate","partial_signal")):
            for argument in h.get("arguments",[]):
                st.write("**Аргумент:**",argument)
            if not h.get("arguments"):
                st.caption("Положительного диагностического признака не установлено.")
            st.json(h.get("observations",{}))
            st.markdown("**Что ещё проверить**")
            for item in h.get("next_verifications",[]):
                st.write("• "+item)
            st.caption("Ограничения: "+"; ".join(h.get("limitations",[])))


def _tree_panel(entry):
    phy=entry.get("phylogeny") or {}
    st.subheader("Филогенетическое дерево")
    if not phy.get("newick"):
        st.warning("Дерево не построено: "+phy.get("reason","недостаточно последовательностей"))
        if phy.get("selection"):
            st.dataframe(phy["selection"],use_container_width=True,hide_index=True)
        return
    if phy.get("status")=="inferred_exploratory":
        st.warning("Это предварительная NJ-кластеризация сходства 7-меров, "
                   "не подтверждённая филогения генов. Поддержка ветвей отсутствует.")
    else:
        st.success("Дерево рассчитано: "+str(phy.get("method")))
    st.caption("Дерево не укоренено. Горизонтальное расположение ветвей "
               "не указывает направление эволюции. Геномные совпадения не гарантируют ортологию.")
    try:
        st.graphviz_chart(newick_to_dot(phy["newick"]),use_container_width=True)
    except (ValueError,TypeError) as exc:
        st.warning("Newick доступен для скачивания, но графическое отображение не удалось: "+str(exc))
    c1,c2,c3=st.columns(3)
    c1.metric("Последовательностей",phy.get("taxa",0))
    c2.metric("Метод",str(phy.get("method","—")).replace("_"," ")[:35])
    c3.metric("Поддержка ветвей",phy.get("support_type") or "Не вычислена")
    if phy.get("warnings"):
        with st.expander("Научные ограничения"):
            for warning in phy["warnings"]:
                st.write("• "+warning)
    if phy.get("selection"):
        with st.expander("Отобранные локусы и качество картирования"):
            st.dataframe(phy["selection"],hide_index=True,use_container_width=True)
    st.download_button("Скачать дерево Newick",phy["newick"],
                       file_name=entry["gene_id"]+".nwk",mime="text/plain")
    if phy.get("alignment_fasta"):
        st.download_button("Скачать множественное выравнивание",
                           phy["alignment_fasta"],file_name=entry["gene_id"]+"_msa.fasta",
                           mime="text/plain")


def _results(report):
    overall=report.get("summary",{})
    a,b,c,d=st.columns(4)
    a.metric("Кандидатных генов",overall.get("requested",0))
    b.metric("Завершено",overall.get("completed",0))
    c.metric("Деревьев построено",overall.get("tree_computed",0))
    d.metric("Событий аудита",overall.get("total_events",0))
    st.caption("Время запуска (UTC): "+report.get("started_at_utc","—")+
               " · завершение: "+report.get("finished_at_utc","—"))
    overview=[{"Ген":k,"Итог":v.get("outcome"),
               "Дерево":(v.get("phylogeny") or {}).get("method") or "не построено",
               "Сигналов":len((v.get("hypotheses") or {}).get("candidates",[]))}
              for k,v in report["genes"].items()]
    st.dataframe(overview,hide_index=True,use_container_width=True)
    gene_id=st.selectbox("Выберите ген для подробного анализа",report["gene_order"],
                         key="evo_result_gene")
    entry=report["genes"][gene_id]
    if entry.get("outcome")=="failed":
        st.error(entry.get("error","Исследование гена не завершено"))
        return
    tabs=st.tabs(["Дерево","Гипотезы","Геномные данные","Технические проверки"])
    with tabs[0]:_tree_panel(entry)
    with tabs[1]:_hypothesis_panel(entry)
    with tabs[2]:
        st.subheader("Наблюдения по выбранным геномам")
        core=entry["evidence"]
        st.dataframe(core.get("species_status",[]),use_container_width=True,hide_index=True)
        st.caption("Отсутствие найденного совпадения не доказывает утрату гена.")
        orf_rows=entry.get("syntenic_orf_evidence") or {}
        with st.expander("Потенциальные ORF в синтенных локусах"):
            st.caption("Это полные ATG–stop рамки в отдельных геномных фрагментах. "
                       "Без сплайсированного сравнения и реконструкции предковых состояний "
                       "они не доказывают ни кодирование, ни происхождение de novo.")
            st.dataframe([
                {"Вид":species,"Статус":item.get("status"),
                 "Длина фрагмента, нт":item.get("length_nt"),
                 "Кандидатных ORF":item.get("candidate_complete_orf_count"),
                 "Максимальная ORF, а.о.":item.get("longest_complete_orf_aa"),
                 "Ограничение":item.get("reason",item.get("interpretation"))}
                for species,item in orf_rows.items()
            ],hide_index=True,use_container_width=True)
        render_report_visualizations(core)
    with tabs[3]:
        render_analysis_annotations(entry["evidence"])
        st.subheader("Методы, качество и независимые проверки")
        st.json(entry.get("limitations",[]))
        st.json((entry["evidence"].get("mechanism_report") or {}).get("completeness",{}))


def _logs(report):
    st.subheader("Воспроизводимость и журнал проверок")
    st.caption("Каждая запись снабжена временем UTC, геном, этапом и результатом. "
               "Ошибки не скрываются. Внешние программы и версии включены в метаданные.")
    events=report.get("audit_log",[])
    genes=st.multiselect("Фильтр генов",report["gene_order"],default=report["gene_order"])
    filtered=[e for e in events if e.get("gene_id") in genes]
    st.dataframe([{
        "Время UTC":e.get("timestamp_utc",""),
        "Ген":e.get("gene_id",""),
        "Вид":e.get("species","—"),
        "Этап":STAGES.get(e.get("stage"),e.get("stage")),
        "Статус":e.get("state"),
        "Результат":e.get("detail",""),
    } for e in filtered],hide_index=True,use_container_width=True)
    with st.expander("Версии инструментов и референсные сборки"):
        st.json(report.get("tools",{}))
        st.json(report.get("expected_study_assemblies",{}))
        st.code("SHA-256 входных генов: "+report["input_sha256"])
    st.download_button("Скачать журнал JSONL",
                       "\n".join(json.dumps(x,ensure_ascii=False) for x in events)+"\n",
                       file_name="evogene_audit.jsonl",mime="application/x-ndjson")
    st.download_button("Скачать весь проект ZIP",_bundle(report),
                       file_name="evogene_project_"+report["run_id"]+".zip",
                       mime="application/zip")


def main():
    st.set_page_config(page_title="EvoGene · происхождение генов",page_icon="🧬",
                       layout="wide",initial_sidebar_state="collapsed")
    _style()
    st.markdown('<div class="evo-intro"><div class="evo-overline">ITMO · СИНТЕЗ · Evolutionary genomics</div>'
                '<h1>EvoGene</h1><p>От генетической последовательности — к проверяемым '
                'гипотезам происхождения, дереву и воспроизводимому отчёту.</p></div>',
                unsafe_allow_html=True)
    st.info("Научная модель: Peng & Zhao (2024). Текущий конвейер — аналог для "
            "скрининга, не дословное воспроизведение Cactus/HAL/Spaln/Genewise.")
    input_tab, result_tab, log_tab, method_tab=st.tabs([
        "01 · Загрузить и исследовать", "02 · Результаты", "03 · Журнал",
        "04 · Методика и ограничения"
    ])
    with input_tab:
        file,focal,chosen,te,phylogeny,cache=_input_panel()
        from hashlib import sha256
        fingerprint=(
            sha256(file.getvalue()).hexdigest() if file else None,
            focal,tuple(chosen),sha256(te.getvalue()).hexdigest() if te else None,
            phylogeny,cache,
        )
        if st.session_state.get("evo_run_input_fingerprint") != fingerprint:
            st.session_state.pop("evo_project",None)
            st.session_state["evo_run_input_fingerprint"]=fingerprint
        _run(file,focal,chosen,te,phylogeny,cache)
    report=st.session_state.get("evo_project")
    with result_tab:
        if report:_results(report)
        else:st.info("Загрузите FASTA и запустите исследование в первой вкладке.")
    with log_tab:
        if report:_logs(report)
        else:st.info("Здесь появятся временные отметки каждого выполненного теста.")
    with method_tab:
        st.markdown("""
### Почему два разных поиска?

**Синтенные локусы:** проверяются родственные геномные участки; отсутствие CDS
не является доказательством некодирующего предка.

**Весь геном:** отдельно ищутся иные копии, возможные паралоги и перемещённые
последовательности; один локус недостаточен для гипотезы de novo.

**Транспозоны:** при предоставлении TE-библиотеки выполняется отдельный
скрининг её сходства с геном. Это не заменяет RepeatMasker и Dfam.

**Филогенетика:** предпочтительно MAFFT + IQ-TREE (ML с UFBoot и SH-aLRT).
Если IQ-TREE недоступен, но есть MAFFT, строится NJ с site-bootstrap.
Если ни один инструмент не установлен, создаётся только предварительная
7-mer NJ-кластеризация — без статистической поддержки и без права
называть её доказанной историей гена.

**Не выполнено автоматически:** полная Progressive Cactus/HAL-реконструкция,
калибровка Genewise/Spaln, глобальный UniProt/jackhmmer/TBLASTN,
реконструкция предкового ORF и экспериментальная проверка трансляции.
Система отмечает эти проверки как не выполненные.
""")
        st.markdown("[Peng & Zhao (2024)](https://www.nature.com/articles/s41467-024-45028-1)")
        st.markdown("[Оригинальные скрипты](https://github.com/LiZhaoLab/DrosophilaDenovoGene)")


if __name__=="__main__":
    main()
