"""Run with: streamlit run auto_app.py"""
import json,tempfile,re
from pathlib import Path
import streamlit as st
from evogene.ncbi_auto import PANEL,analyze
from evogene.tools import check_tools
st.set_page_config(page_title='EvoGene Auto NCBI',layout='wide')
st.title('EvoGene — автоматическая загрузка геномов NCBI')
st.caption('Вы загружаете только свой FASTA. Сравнительные геномы и доступные GFF3 скачиваются автоматически.')
st.warning('FASTA полного генома не содержит названия и координат генов. Нужно ввести последовательность гена либо координаты. Это не доказательство de novo.')
st.subheader('Проверка инструментов на этом Mac')
tools_status=check_tools()
for name, status in tools_status.items():
    if status['ok']:
        st.success(f"{name}: {status['path']} — {status['detail'][:100]}")
    else:
        st.error(f"{name}: {status['detail']}")
if not all(status['ok'] for status in tools_status.values()):
    st.info('Анализ недоступен, пока оба инструмента не запускаются. PATH не требуется настраивать вручную, если программы находятся в стандартных папках.')
f=st.file_uploader('Целевой геном FASTA или FASTA одного гена',type=['fa','fna','fasta','gz'])
species=st.selectbox('Целевой вид',list(PANEL))
gene=st.text_input('ID выбранного гена (подпись отчёта)',value='FBgn0040343')
count=st.slider('Количество родственных видов',1,5,2)
method=st.radio('Как выбрать ген?', ['Вставить последовательность ДНК','Указать координаты','FASTA содержит только один ген'])
dna='';chrom='';start=end=None;strand='+'
if method=='Вставить последовательность ДНК':
    dna=st.text_area('ДНК гена (A/C/G/T/N), можно FASTA-текстом')
elif method=='Указать координаты':
    chrom=st.text_input('Название хромосомы в FASTA')
    c1,c2,c3=st.columns(3)
    start=c1.number_input('Начало (с 1)',min_value=1,value=1)
    end=c2.number_input('Конец (включительно)',min_value=1,value=1000)
    strand=c3.selectbox('Цепь',['+','-'])
if st.button('Скачать геномы NCBI и запустить анализ',type='primary',
             disabled=not all(status['ok'] for status in tools_status.values())):
    if not f:st.error('Загрузите FASTA')
    else:
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/('target'+''.join(Path(f.name).suffixes))
            path.write_bytes(f.getvalue())
            status=st.empty()
            with st.spinner('NCBI download and minimap2...'):
                try:
                    st.session_state['auto_report']=analyze(path,species,gene,gene_sequence=dna,
                        chrom=chrom,start=start,end=end,strand=strand,count=count,
                        progress=lambda message:status.info(message))
                except Exception as exc:st.error(str(exc))
            status.empty()
if 'auto_report' in st.session_state:
    r=st.session_state['auto_report']
    st.subheader('Проверка всех выбранных видов')
    if r.get('species_status'):
        labels={'matched':'Найдено сходство','no_alignment':'Выравниваний нет','failed':'Ошибка обработки'}
        st.dataframe([{'Вид':entry['species'],
                       'Статус':labels.get(entry['status'],entry['status']),
                       'Всего попаданий':entry.get('alignment_count',0),
                       'Прошли порог':entry.get('passing_count',0),
                       'Ошибка':entry.get('error','')}
                      for entry in r['species_status']],
                     use_container_width=True,hide_index=True)
    st.caption('Пустая таблица по виду означает отсутствие выравниваний при текущих настройках; '
               'это не доказательство отсутствия гена.')
    st.subheader('Итог: варианты эволюционного происхождения')
    support = r.get('origin_support')
    if support:
        st.warning(support['warning'])
        st.write(support['conclusion'])
        st.dataframe([{
            'Механизм': item['label'],
            'Поддержка (%)': item['support_percent'],
            'Пояснение': item['explanation']
        } for item in support['options']], use_container_width=True, hide_index=True)
        st.bar_chart([
            {'Гипотеза': item['label'], 'Поддержка (%)': item['support_percent']}
            for item in support['options']
        ], x='Гипотеза', y='Поддержка (%)', horizontal=True)
        st.caption('Проценты НЕ являются статистически калиброванными вероятностями.')
        st.subheader('Проверки каждой гипотезы')
        for mechanism, outcome in support['tests'].items():
            with st.expander(mechanism.replace('_',' ').title()+' — '+outcome['status']):
                st.write(outcome['interpretation'])
                st.write('Свидетельства:',outcome['observations'])
                st.write('Для проверки ещё нужно:',outcome['missing'])
        if r.get('origin_screen_errors'):
            st.warning('Не удалось выполнить часть проверок: '+str(r['origin_screen_errors']))
    st.subheader('Технические результаты: '+r['hypothesis'])
    st.write('NCBI assembly accessions:',[(x['species'],x['accession']) for x in r['ncbi_assemblies']])
    st.write('Matching species:',r['matching_species'])
    st.write('Download errors:',r['download_or_analysis_errors'])
    st.write('Limitations:',r['limitations'])
    for s,hits in r['relatives'].items():
        with st.expander(s):
            if hits:st.dataframe(hits)
            else:st.info('Сравнение выполнено, но minimap2 не нашёл выравниваний '
                         'при текущих настройках. Это не доказывает отсутствие гена.')
    st.subheader('Полный JSON-отчёт')
    payload=json.dumps(r,indent=2,ensure_ascii=False).encode('utf-8')
    safe_name=re.sub(r'[^A-Za-z0-9_.-]+','_',str(r.get('gene_id','gene')))[:100]
    st.download_button('Скачать JSON-отчёт',data=payload,
                       file_name=safe_name+'.json',mime='application/json',
                       key='auto_json_download',use_container_width=True)
    st.caption('Если скачивание в браузере не происходит, локальная копия сохранена по пути ниже.')
    local_path=Path('results/auto_ncbi')/(safe_name+'.json')
    try:
        local_path.parent.mkdir(parents=True,exist_ok=True)
        local_path.write_bytes(payload)
        st.code(str(local_path.resolve()))
    except OSError as exc:
        st.warning('Не удалось сохранить JSON локально: '+str(exc))
