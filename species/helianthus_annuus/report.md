# Отчет по построению генетической карты подсолнечника

## Цель

Для подсолнечника (*Helianthus annuus*) нужно было подготовить генетическую карту на актуальный референсный геном:

`GCF_002127325.2_HanXRQr2.0-SUNRISE`

Финальный ожидаемый формат карты:

~~~text
chr	pos	cM
~~~

## Исходные данные

В качестве источника использовался supplementary Dataset S6 из статьи PLOS ONE по SFP-карте подсолнечника.

Исходный файл:

~~~text
species/helianthus_annuus/data/raw/pone.0051360.s006.xlsx
~~~

Из таблицы были извлечены 25-bp SFP features. Для каждого feature были сохранены:

- `feature_id`
- `chip_x`
- `chip_y`
- `unigene`
- `linkage_group`
- `cM`
- `plants_scored`
- `mismatches_to_template`
- `sequence`

После парсинга Dataset S6 было получено:

~~~text
Total feature rows: 67846
Unique feature IDs: 67846
Unique sequences: 66606
Unique unigenes: 22487
Bad sequences: 0
~~~

Все последовательности имели длину 25 bp.

Созданные файлы:

~~~text
species/helianthus_annuus/data/metadata/sunflower_features_metadata.tsv
species/helianthus_annuus/data/markers/sunflower_features_25bp.fasta
species/helianthus_annuus/results/qc/sunflower_features_dataset_summary.txt
~~~

## Особенность данных

В статье использовались SFP-маркеры — single-feature polymorphisms.
Физическую координату маркера на актуальном референсе нужно было получить через выравнивание самой 25-bp последовательности на геном.

Так как последовательности очень короткие, был выбран строгий подход, учитывать только идеальные и уникальные попадания.

## Референсный геном

Для выравнивания использовался серверный референс:

~~~text
/mnt/reference/genomes/heliantus_annuus/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.unmasked.fna
~~~

Обычный файл `*_genomic.fna`, указанный в таблице путей, на сервере фактически отсутствовал. В папке был доступен только `*_genomic.unmasked.fna`, поэтому для работы был использован unmasked reference genome.

На сервере рядом с FASTA уже были подготовлены BWA-индексы, поэтому дополнительная BLAST-база не создавалась.

## Выравнивание

Выравнивание 25-bp SFP features выполнялось на сервере `novaplant1`.

Использованные инструменты:

~~~text
bwa 0.7.17
samtools
python3
~~~

Команда выравнивания:

~~~bash
bwa aln \
  -n 0 \
  -o 0 \
  -l 25 \
  -t 8 \
  /mnt/reference/genomes/heliantus_annuus/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.unmasked.fna \
  species/helianthus_annuus/data/markers/sunflower_features_25bp.fasta \
  > species/helianthus_annuus/results/intermediate/sunflower_features_25bp.bwa_exact.sai

bwa samse \
  /mnt/reference/genomes/heliantus_annuus/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.unmasked.fna \
  species/helianthus_annuus/results/intermediate/sunflower_features_25bp.bwa_exact.sai \
  species/helianthus_annuus/data/markers/sunflower_features_25bp.fasta \
  > species/helianthus_annuus/results/intermediate/sunflower_features_25bp.bwa_exact.sam
~~~

Параметры:

- `-n 0` — не допускаются mismatch;
- `-o 0` — не допускаются gap openings;
- `-l 25` — seed length равен полной длине feature;
- `-t 8` — использовано 8 потоков.

## Результат выравнивания

BWA обработал все 67 846 feature-последовательностей.

Лог выравнивания:

~~~text
67846 sequences have been processed.
SAM alignment rows: 67846
~~~

Основные файлы выравнивания:

~~~text
species/helianthus_annuus/results/intermediate/sunflower_features_25bp.bwa_exact.sai
species/helianthus_annuus/results/intermediate/sunflower_features_25bp.bwa_exact.sam
species/helianthus_annuus/logs/02_align_sunflower_features_bwa_exact.log
~~~

## Фильтрация SAM и построение карты

Для построения карты использовался скрипт:

~~~text
species/helianthus_annuus/scripts/03_build_sunflower_map_from_bwa_exact.py
~~~

Критерии включения feature в финальную карту:

1. SFP выравнивались на референсные хромосомные последовательности;
2. референсная последовательность соответствует одной из хромосом `NC_035433.2`–`NC_035449.2`;
3. выравнивание полноразмерное: `25M`;
4. совпадение идеальное: `NM:i:0`;
5. попадание уникальное: `X0:i:1`;
6. опубликованная linkage group совпадает с физической хромосомой;
7. удалены точные дубли `chr-pos-cM`;
8. удалены позиции, где один и тот же `chr-pos` имел конфликтующие значения `cM`.

Соответствие accession-хромосом и номеров хромосом:

~~~text
NC_035433.2 -> chr1
NC_035434.2 -> chr2
NC_035435.2 -> chr3
NC_035436.2 -> chr4
NC_035437.2 -> chr5
NC_035438.2 -> chr6
NC_035439.2 -> chr7
NC_035440.2 -> chr8
NC_035441.2 -> chr9
NC_035442.2 -> chr10
NC_035443.2 -> chr11
NC_035444.2 -> chr12
NC_035445.2 -> chr13
NC_035446.2 -> chr14
NC_035447.2 -> chr15
NC_035448.2 -> chr16
NC_035449.2 -> chr17
~~~

Физическая координата `pos` рассчитывалась как середина 25-bp выравнивания.

## QC фильтрации

Итоговая статистика:

~~~text
Input features in metadata: 67846
SAM alignment rows: 67846
Unmapped features: 31930
Exact unique chromosome-level hits before LG filtering: 26411
LG/chr conflict rows excluded: 4248
Rows after LG/chr filtering: 22163
Rows after exact chr-pos-cM deduplication: 21708
chr-pos coordinates with conflicting cM values: 122
Rows excluded due to chr-pos cM conflicts: 266
Final unique chr-pos-cM rows: 21442
~~~

Распределение BWA `X0` tag показало, что часть features имела несколько возможных мест выравнивания. Такие мультимапперы были исключены из финальной карты.

~~~text
X0=1: 26442
X0=2: 5811
X0=3: 1747
X0=4: 893
X0=5: 357
~~~

Это ожидаемо для коротких 25-bp последовательностей и большого повторного генома подсолнечника.

## Итоговая карта

Финальный файл:

~~~text
species/helianthus_annuus/results/final/sunflower_genetic_map.bwa_exact_unique.tsv
~~~

Размер:

~~~text
21442 маркера
~~~


## Распределение финальных маркеров по хромосомам

~~~text
chr1    1053
chr2    603
chr3    1449
chr4    1306
chr5    1600
chr6    693
chr7    648
chr8    1327
chr9    1676
chr10   2235
chr11   1008
chr12   1367
chr13   1473
chr14   1206
chr15   1065
chr16   1171
chr17   1562
~~~

Все 17 хромосом подсолнечника представлены в финальной карте.



## Визуализация финальной карты

График:

~~~text
species/helianthus_annuus/results/figures/sunflower_genetic_map_coverage.svg
~~~

![Покрытие генетической карты по cM](results/figures/sunflower_genetic_map_coverage.svg)

Финальная карта содержит 21,442 маркера, распределенных по всем 17 хромосомам подсолнечника.


### Покрытие физической карты

График:

~~~text
species/helianthus_annuus/results/figures/sunflower_physical_map_coverage.svg
~~~
![Покрытие физической карты](results/figures/sunflower_physical_map_coverage.svg)

На этом графике те же финальные маркеры показаны уже не в генетических координатах cM, а в физических координатах референсного генома, то есть в bp/Mb.



