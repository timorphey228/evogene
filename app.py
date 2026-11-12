"""Streamlit UI. Upload original supplementary XLSX and manually curated evidence CSV."""
import tempfile
from pathlib import Path
import streamlit as st
from evogene.core import classify
from evogene.io import read_evidence, read_properties
from evogene.orf import find_orfs

st.set_page_config(page_title='EvoGene', layout='wide')
st.title('EvoGene · Evolutionary origin of genes')
st.caption('Prototype grounded in LiZhaoLab workflow and Peng & Zhao (2024). Hypotheses are not verified ancestry claims.')
st.info('Important: article labels and protein characteristics are reference metadata, not independently verified evidence of gene origin.')

left, right = st.columns(2)
with left:
    evidence_upload = st.file_uploader('Evidence CSV (homology, synteny, outgroups)', type=['csv'])
with right:
    features_upload = st.file_uploader('Article supplementary XLSX or features CSV', type=['xlsx', 'csv'])

def save_temp(upload):
    with tempfile.NamedTemporaryFile(delete=False, suffix=Path(upload.name).suffix) as out:
        out.write(upload.getvalue())
        return Path(out.name)

if evidence_upload:
    paths = []
    try:
        ep = save_temp(evidence_upload)
        paths.append(ep)
        properties = {}
        if features_upload:
            pp = save_temp(features_upload)
            paths.append(pp)
            properties = read_properties(pp)
        results = [classify(e, properties.get(e.gene_id)) for e in read_evidence(ep)]
        st.metric('Genes analyzed', len(results))
        st.dataframe([{'Gene': r.gene_id, 'Hypothesis': r.hypothesis,
                       'Evidence strength': r.evidence_strength, 'Article label': r.research_label,
                       'Age branch': r.research_age} for r in results], use_container_width=True)
        import json
        st.download_button('Download JSON', json.dumps([r.to_dict() for r in results], indent=2), 'evogene_results.json', 'application/json')
        for r in results:
            with st.expander(f'{r.gene_id} — {r.hypothesis}'):
                st.write('**Supporting evidence:**', r.supporting or 'None')
                st.write('**Contradictions:**', r.conflicting or 'None')
                st.write('**Missing data:**', r.missing or 'None')
                st.write('**Alternatives:**', r.alternatives or 'None')
                st.write('**Qualitative scores (not probabilities):**', r.scores)
                if r.properties:
                    st.json(r.properties)
    except Exception as exc:
        st.error(f'Invalid file: {exc}')
    finally:
        for path in paths:
            path.unlink(missing_ok=True)
else:
    st.write('Upload a manually checked evidence CSV to start. Example: ' + 'data/demo_evidence.csv' + '.')

st.divider()
st.subheader('ORF survey (optional)')
seq = st.text_area('Paste DNA from syntenic ancestral locus (A/C/G/T/N)', height=100)
min_aa = st.number_input('Minimum amino acid length', min_value=1, value=30)
if st.button('Find ORFs') and seq.strip():
    try:
        st.dataframe([o.to_dict() for o in find_orfs(seq, min_aa)], use_container_width=True)
        st.caption('ORFs indicate coding potential only; they do not establish ancestral transcription or translation.')
    except ValueError as exc:
        st.error(str(exc))
