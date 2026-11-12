"""Import the article's supplementary CSV/XLSX and independently produced evidence CSV.

XLSX is read via Python standard-library ZIP/XML, so no Excel dependency is needed.
"""
from __future__ import annotations
import csv
import io
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from .core import Evidence

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


def _column_index(cell: str) -> int:
    letters = re.match(r'[A-Z]+', cell)
    if not letters:
        raise ValueError(f'Invalid cell reference: {cell}')
    out = 0
    for c in letters.group():
        out = out * 26 + ord(c) - 64
    return out - 1


def read_xlsx_sheet(path: str | Path, sheet_name: str = 'Figure 2 a,c,d,e,f,g,h') -> list[list[str]]:
    with zipfile.ZipFile(path) as archive:
        strings: list[str] = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            root = ET.fromstring(archive.read('xl/sharedStrings.xml'))
            strings = [''.join(t.text or '' for t in si.iter(NS + 't')) for si in root.findall(NS + 'si')]
        workbook = ET.fromstring(archive.read('xl/workbook.xml'))
        rels = ET.fromstring(archive.read('xl/_rels/workbook.xml.rels'))
        rmap = {r.attrib['Id']: r.attrib['Target'] for r in rels}
        sheet = next((s for s in workbook.find(NS + 'sheets') if s.attrib['name'] == sheet_name), None)
        if sheet is None:
            raise ValueError(f'Sheet {sheet_name!r} not found')
        rid = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
        target = rmap[rid].lstrip('/')
        entry = target if target.startswith('xl/') else 'xl/' + target
        root = ET.fromstring(archive.read(entry))
        rows: list[list[str]] = []
        for row in root.iter(NS + 'row'):
            values: dict[int, str] = {}
            for cell in row.findall(NS + 'c'):
                v = cell.find(NS + 'v')
                inline = cell.find(NS + 'is')
                if v is not None:
                    raw = v.text or ''
                    val = strings[int(raw)] if cell.attrib.get('t') == 's' else raw
                elif inline is not None:
                    val = ''.join(t.text or '' for t in inline.iter(NS + 't'))
                else:
                    val = ''
                values[_column_index(cell.attrib['r'])] = val
            if values:
                rows.append([values.get(i, '') for i in range(max(values) + 1)])
        return rows


def read_properties(path: str | Path) -> dict[str, dict[str, str]]:
    path = Path(path)
    if path.suffix.lower() == '.xlsx':
        rows = read_xlsx_sheet(path)
    else:
        with path.open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.reader(stream))
    header_index = next((i for i, r in enumerate(rows) if r and r[0].strip() == 'FBID'), None)
    if header_index is None:
        raise ValueError('Expected FBID header (Figure 2 a,c,d,e,f,g,h)')
    headers = rows[header_index]
    result = {}
    for row in rows[header_index + 1:]:
        if not row or not re.fullmatch(r'FBgn\d+', row[0].strip()):
            continue
        result[row[0].strip()] = {key: row[i].strip() if i < len(row) else '' for i, key in enumerate(headers) if key}
    return result


def read_evidence(path: str | Path) -> list[Evidence]:
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return [Evidence.from_mapping(row) for row in csv.DictReader(stream)]


def export_results(results: list, output: str | Path) -> None:
    import json
    Path(output).write_text(json.dumps([r.to_dict() for r in results], indent=2, ensure_ascii=False), encoding='utf-8')
