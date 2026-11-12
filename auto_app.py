"""Canonical Streamlit entry point for EvoGene comparative-genomics workflow.

Start with:
    streamlit run auto_app.py

The legacy single-gene UI has been replaced by an auditable, multi-gene
research dashboard. Other modes remain available in app.py.
"""
from evogene.dashboard import main

if __name__ == "__main__":
    main()
