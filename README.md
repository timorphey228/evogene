# EvoGene v2 — анализ геномных последовательностей

Два режима анализа:

1. **genomes**: FASTA целевого генома + GFF3 (с координатами выбранного гена) + FASTA родственных видов и, желательно, их GFF3. Реальное сопоставление ДНК с помощью minimap2, оценка покрытия и примерной идентичности, поиск дополнительных локусов и перекрытия аннотированных CDS.
2. **case**: анализ загруженного в кейсе файла `supplementary_file_unannotated_alignments.txt`, содержащего заранее выделенные Cactus-фрагменты для генов. Он **не является набором полных геномов**.

## Установка на Mac

```bash
brew install minimap2
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Запуск по целым геномам

```bash
python -m evogene.pipeline genomes \
  --target-fasta /data/dmel_genome.fna.gz \
  --target-gff /data/dmel_annotations.gff3.gz \
  --gene FBgn0040343 \
  --relative dsim=/data/dsim_genome.fna.gz:/data/dsim_annotations.gff3.gz \
  --relative dyak=/data/dyak_genome.fna.gz:/data/dyak_annotations.gff3.gz \
  --results results/genome_analysis
```

GFF3 должен содержать запись `gene` с соответствующим ID. Сопоставление выполнено для ДНК всего генного интервала, а не для выделенного белка или отдельных экзонов. Сходный участок не означает ортологию, а отсутствие аннотированного CDS не доказывает некодирующее происхождение.

## Примеры по данным кейса

```bash
python -m evogene.pipeline case \
  --alignment-file /path/to/supplementary_file_unannotated_alignments.txt \
  --gene FBgn0262810 --gene FBgn0025616 \
  --gene FBgn0029535 --gene FBgn0040343 \
  --results results/case_examples
```

Результаты четырёх реальных примеров сохранены в `results/case_examples/`. Они основаны на *предварительно отобранных* фрагментах исходного файла, а не на автоматическом поиске по полным сборкам. Для полного комплекта JSON/CSV и тестов также доступен архив EvoGene_genome_analysis_v2.zip.

Сборки из приложенной таблицы `supplementary_file_genome_assembly.csv` перечислены в `data/assembly_accessions.csv`. Это идентификаторы сборок, **не FASTA**.

## Научные ограничения

- Поиск сходства ДНК не определяет ортологию или механизм возникновения сам по себе.
- Необходимо отдельно подтвердить синтению по фланкирующим генам, гомологию белков, дерево генов и предковое состояние ORF.
- Отсутствие попадания может объясняться потерей гена, дивергенцией или качеством сборки.
- Алгоритм не выдаёт калиброванную вероятность de novo или точный возраст.
- Результаты для примеров — **наблюдения о сходстве последовательностей, а не доказанные механизмы эволюции**.

Биоинформатическая основа: [LiZhaoLab/DrosophilaDenovoGene](https://github.com/LiZhaoLab/DrosophilaDenovoGene), стадии Cactus, ортологи, проверка скрытых гомологов. Код прототипа реализован отдельно.
