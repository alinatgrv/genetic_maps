# Zea mays — построение генетической карты с координатами на NAM-5.0

## Цель

Для кукурузы требовалось получить таблицу генетической карты в формате:

~~~text
chr    pos    cM
~~~

где:

- `chr` — номер хромосомы в актуальной сборке кукурузы;
- `pos` — физическая позиция на референсном геноме;
- `cM` — генетическая координата в сантиморганах.

В качестве источника генетической карты использовались данные MaizeGDB для карты **Genetic 1**, начиная с record id `1203637`.

Целевой референсный геном:

~~~text
GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0
~~~

Именно на эту сборку должны были быть перенесены физические координаты маркеров.

## Исходные данные

Данные были получены с MaizeGDB:

~~~text
https://www.maizegdb.org/data_center/map?id=1203637
~~~

Карта **Genetic 1** представлена отдельными linkage groups / chromosomes. Для загрузки использовались map ids:

~~~text
1203637 — LG1 / chr1
1203638 — LG2 / chr2
1203639 — LG3 / chr3
1203640 — LG4 / chr4
1203641 — LG5 / chr5
1203642 — LG6 / chr6
1203643 — LG7 / chr7
1203644 — LG8 / chr8
1203645 — LG9 / chr9
1203646 — LG10 / chr10
~~~

Попытка скачать `/map_text?id=...` напрямую через `curl` не сработала, потому что сайт вернул Cloudflare HTML-заглушку `Just a moment...`. Поэтому данные были скачаны через браузерную JavaScript-консоль, где страница уже была загружена и имела доступ к MaizeGDB.

Скачанные файлы сохранены в:

~~~text
species/zea_mays/data/raw/maizegdb_map_1203637/map_text_browser/
~~~

Файлы:

~~~text
maizegdb_genetic1_1203637_map_text_BROWSER.txt
maizegdb_genetic1_1203638_map_text_BROWSER.txt
maizegdb_genetic1_1203639_map_text_BROWSER.txt
maizegdb_genetic1_1203640_map_text_BROWSER.txt
maizegdb_genetic1_1203641_map_text_BROWSER.txt
maizegdb_genetic1_1203642_map_text_BROWSER.txt
maizegdb_genetic1_1203643_map_text_BROWSER.txt
maizegdb_genetic1_1203644_map_text_BROWSER.txt
maizegdb_genetic1_1203645_map_text_BROWSER.txt
maizegdb_genetic1_1203646_map_text_BROWSER.txt
~~~

Файлы `map_text_BROWSER.txt` содержат генетические координаты маркеров, а также физические координаты на нескольких сборках генома кукурузы. Для выполнения задания использовались только колонки целевой сборки:

~~~text
Zm-B73-REFERENCE-NAM-5.0_gene_model
Zm-B73-REFERENCE-NAM-5.0_chr
Zm-B73-REFERENCE-NAM-5.0_start
Zm-B73-REFERENCE-NAM-5.0_end
~~~

Колонки, относящиеся к более старым или другим сборкам, например `B73 RefGen_v3` или `Zm-B73-REFERENCE-GRAMENE-4.0`, не использовались для финальной карты, потому что требуемая система координат — именно `GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0`.

## Метод

Для обработки данных был написан скрипт:

~~~text
species/zea_mays/scripts/02_build_maize_map_from_maizegdb_map_text.py
~~~

Скрипт выполняет следующие шаги:

1. Читает все `map_text_BROWSER.txt` файлы для map ids `1203637–1203646`.
2. Объединяет строки MaizeGDB в одну таблицу.
3. Извлекает:
   - `Locus`;
   - `Coordinate`, который соответствует генетической координате `cM`;
   - `Zm-B73-REFERENCE-NAM-5.0_chr`;
   - `Zm-B73-REFERENCE-NAM-5.0_start`;
   - `Zm-B73-REFERENCE-NAM-5.0_end`.
4. Оставляет только строки, где есть валидная генетическая координата и валидные физические координаты на NAM-5.0.
5. Приводит chromosome к числовому виду: `chr1 → 1`, `chr2 → 2` и так далее.
6. Рассчитывает физическую позицию маркера как середину интервала:

~~~text
pos = round((start + end) / 2)
~~~

7. Формирует финальную таблицу только с тремя колонками:

~~~text
chr    pos    cM
~~~

8. Удаляет полностью дублирующиеся строки `chr-pos-cM`.

## Результаты построения карты

Общая статистика обработки:

~~~text
Input files: 10
All parsed MaizeGDB rows: 23683
Rows with numeric cM: 23683
Rows with NAM-5.0 chr/start/end/cM: 21355
Final unique chr-pos-cM rows: 21288
~~~

Финальная карта сохранена в файл:

~~~text
species/zea_mays/results/final/zea_mays_genetic_map.tsv
~~~

Первые строки финальной карты:

~~~text
chr    pos     cM
1      37410   0.0
1      43988   0.01
1      111468  0.04
1      189070  0.09
1      194512  0.09
1      201828  0.1
1      208171  0.12
1      246832  0.13
1      315532  0.14
~~~

## Распределение маркеров по хромосомам

После фильтрации и удаления дубликатов финальная карта содержит **21,288** уникальных строк `chr-pos-cM`.

Распределение по хромосомам:

~~~text
chr1     3381
chr2     2575
chr3     2266
chr4     2173
chr5     2373
chr6     1781
chr7     1681
chr8     1904
chr9     1611
chr10    1543
~~~

Полная таблица с дополнительной статистикой сохранена в:

~~~text
species/zea_mays/results/qc/maize_marker_counts_by_chr.tsv
~~~

## Контроль качества

Для проверки финальной карты был написан QC-скрипт:

~~~text
species/zea_mays/scripts/03_qc_maize_nam5_map.py
~~~

Он проверяет:

1. количество строк с валидными NAM-5.0 координатами;
2. количество строк, не имеющих координат на NAM-5.0;
3. дубликаты физических позиций;
4. соответствие linkage group и chromosome;
5. корреляцию между физической позицией и генетической координатой.

QC summary сохранён в:

~~~text
species/zea_mays/results/qc/maize_nam5_map_qc_summary.txt
~~~

Основные QC-результаты:

~~~text
Rows in intermediate table: 23683
Rows with valid NAM-5.0 position: 21355
Rows without valid NAM-5.0 position: 2328
Rows in final chr-pos-cM table: 21288
Duplicate chr-pos rows before final deduplication: 230 rows
~~~

### Соответствие linkage group и chromosome

Проверка `chromosome vs linkage group` показала идеальное диагональное соответствие:

~~~text
LG1  → chr1
LG2  → chr2
LG3  → chr3
LG4  → chr4
LG5  → chr5
LG6  → chr6
LG7  → chr7
LG8  → chr8
LG9  → chr9
LG10 → chr10
~~~

Это подтверждает, что координаты были корректно извлечены для соответствующих хромосом, без заметных межхромосомных конфликтов.

### Корреляция физической и генетической координаты

Для каждой хромосомы была рассчитана корреляция между физической позицией `pos` и генетической координатой `cM`.

Spearman correlation:

~~~text
chr1     0.999997
chr2     0.998497
chr3     0.999984
chr4     0.999790
chr5     1.000000
chr6     0.999999
chr7     0.999999
chr8     0.999997
chr9     0.999999
chr10    0.996949
~~~

Значения Spearman correlation близки к 1 для всех хромосом. Это означает, что порядок маркеров по физической координате практически полностью согласуется с порядком по генетической координате.

Небольшое число локальных уменьшений `cM` при сортировке по физической позиции наблюдалось на всех хромосомах, но оно невелико относительно общего числа маркеров. Такие случаи могут быть связаны с локальными особенностями курирования карты, одинаковыми или близкими генетическими координатами, а также с тем, что часть позиций получена по gene model intervals.

QC-файлы:

~~~text
species/zea_mays/results/qc/maize_nam5_duplicate_chr_pos.tsv
species/zea_mays/results/qc/maize_nam5_chr_vs_linkage_group.tsv
species/zea_mays/results/qc/maize_nam5_rows_without_valid_position.tsv
species/zea_mays/results/qc/maize_nam5_physical_vs_genetic_correlation.tsv
~~~

## Визуализация

Для визуализации финальной карты был написан скрипт:

~~~text
species/zea_mays/scripts/04_visualize_maize_final_map.py
~~~

Графики построены в том же стиле, что и для капусты: в виде SVG-схем с распределением маркеров по хромосомам.

Получены три основных графика:

~~~text
species/zea_mays/results/figures/maize_genetic_map_coverage.svg
species/zea_mays/results/figures/maize_physical_map_coverage.svg
species/zea_mays/results/figures/maize_marker_density_10cM.svg
~~~

### `maize_genetic_map_coverage.svg`

Показывает покрытие каждой хромосомы маркерами по генетической координате `cM`.

### `maize_physical_map_coverage.svg`

Показывает распределение маркеров по физической координате на сборке NAM-5.0.

### `maize_marker_density_10cM.svg`

Показывает плотность маркеров в интервалах по 10 cM.

Также была сохранена таблица плотности маркеров:

~~~text
species/zea_mays/results/qc/maize_marker_density_10cM_bins.tsv
~~~

## Итог

Для кукурузы была построена генетическая карта в требуемом формате `chr-pos-cM` на актуальной сборке:

~~~text
GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0
~~~

Исходная генетическая информация была взята из MaizeGDB Genetic 1, а физические координаты были извлечены из колонок, относящихся именно к сборке `Zm-B73-REFERENCE-NAM-5.0`.

Финальная карта содержит:

~~~text
21288 unique chr-pos-cM rows
~~~

QC показал:

- полное соответствие linkage groups и chromosomes;
- высокую согласованность физического и генетического порядка маркеров;
- отсутствие межхромосомных конфликтов;
- плотное покрытие всех 10 хромосом кукурузы.

Финальный файл:

~~~text
species/zea_mays/results/final/zea_mays_genetic_map.tsv
~~~

Эта карта может быть использована как финальный результат для пункта задания по кукурузе.
