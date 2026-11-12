"""Streamlit frontend for genome FASTA/GFF3 and uploaded case alignment extracts."""
import json,tempfile
from pathlib import Path
import streamlit as st
from evogene.genomes import analyze as analyze_genomes
from evogene.case_data import analyze as analyze_case

st.set_page_config(page_title='EvoGene v2',layout='wide')
st.title('EvoGene · происхождение белок-кодирующих генов')
st.warning('Сходство ДНК не доказывает ортологию, синтению, дупликацию или происхождение de novo.')
mode=st.radio('Режим', ['Геном FASTA + GFF3','Предварительно найденные Cactus-фрагменты'])
def saved(upload,folder):
    if upload is None:return None
    path=Path(folder)/(str(len(list(Path(folder).iterdir())))+''.join(Path(upload.name).suffixes))
    path.write_bytes(upload.getvalue())
    return str(path)
if mode.startswith('Геном'):
    st.caption('Нужен minimap2 (macOS: brew install minimap2). Для больших геномов используйте CLI.')
    gene=st.text_input('Ген из GFF3',value='FBgn0040343')
    tf=st.file_uploader('Целевой геном FASTA',key='target')
    tg=st.file_uploader('GFF3 целевого генома',key='targetg')
    s1=st.text_input('Первый родственный вид',value='dsim')
    r1=st.file_uploader('FASTA первого родственного вида',key='r1')
    g1=st.file_uploader('GFF3 первого вида (необязательно)',key='g1')
    s2=st.text_input('Второй родственный вид',value='dyak')
    r2=st.file_uploader('FASTA второго родственного вида',key='r2')
    g2=st.file_uploader('GFF3 второго вида (необязательно)',key='g2')
    if st.button('Анализировать FASTA'):
        if not all((tf,tg,r1)):
            st.error('Нужны целевой FASTA, целевой GFF3 и хотя бы один родственный FASTA.')
        else:
            try:
                with tempfile.TemporaryDirectory() as folder:
                    relatives=[{'species':s1,'fasta':saved(r1,folder),'gff':saved(g1,folder)}]
                    if r2: relatives.append({'species':s2,'fasta':saved(r2,folder),'gff':saved(g2,folder)})
                    st.session_state['evo_result']=analyze_genomes(saved(tf,folder),saved(tg,folder),gene,relatives)
            except Exception as ex:st.error(str(ex))
else:
    st.caption('TXT содержит уже выбранные фрагменты, а не собранные геномы.')
    file=st.file_uploader('supplementary_file_unannotated_alignments.txt',type='txt')
    gene=st.selectbox('Ген',['FBgn0262810','FBgn0025616','FBgn0029535','FBgn0040343'])
    if st.button('Анализировать фрагменты'):
        if not file:st.error('Загрузите TXT')
        else:
            try:
                with tempfile.NamedTemporaryFile(suffix='.txt') as temp:
                    temp.write(file.getvalue());temp.flush()
                    st.session_state['evo_result']=analyze_case(temp.name,gene)
            except Exception as ex:st.error(str(ex))
if 'evo_result' in st.session_state:
    result=st.session_state['evo_result']
    st.subheader(result['gene_id']+' — '+result['hypothesis'])
    st.write(result['explanation'])
    if 'hits' in result:st.dataframe(result['hits'],use_container_width=True)
    else:st.json(result['relatives'])
    st.write('Научные ограничения:',result['limitations'])
    st.download_button('Скачать JSON',json.dumps(result,ensure_ascii=False,indent=2),
                       file_name=result['gene_id']+'.json',mime='application/json')
