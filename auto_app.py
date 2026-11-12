"""Run with: streamlit run auto_app.py"""
import json,tempfile
from pathlib import Path
import streamlit as st
from evogene.ncbi_auto import PANEL,analyze
st.set_page_config(page_title='EvoGene Auto NCBI',layout='wide')
st.title('EvoGene — автоматическая загрузка геномов NCBI')
st.caption('Вы загружаете только свой FASTA. Сравнительные геномы и доступные GFF3 скачиваются автоматически.')
st.warning('FASTA полного генома не содержит названия и координат генов. Нужно ввести последовательность гена либо координаты. Это не доказательство de novo.')
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
if st.button('Скачать геномы NCBI и запустить анализ',type='primary'):
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
    st.subheader(r['hypothesis'])
    st.write('NCBI assembly accessions:',[(x['species'],x['accession']) for x in r['ncbi_assemblies']])
    st.write('Matching species:',r['matching_species'])
    st.write('Download errors:',r['download_or_analysis_errors'])
    st.write('Limitations:',r['limitations'])
    for s,hits in r['relatives'].items():
        with st.expander(s):st.dataframe(hits)
    st.download_button('Download JSON report',json.dumps(r,indent=2,ensure_ascii=False),gene+'.json')
