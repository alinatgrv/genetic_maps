# Genetic maps

## Описание проекта

Цель проекта — получить для каждой культуры таблицу генетической карты в формате:

~~~text
chr    pos    cM
~~~

где:

- `chr` — хромосома в актуальной референсной сборке;
- `pos` — физическая позиция маркера;
- `cM` — генетическая координата маркера.

Для следующих культур: 
- кукуруза — *Zea mays*;
- подсолнечник — *Helianthus annuus*;
- капуста белокочанная — *Brassica oleracea* var. *capitata*.
- рис - *Oryza sativa*
- картофель - *Solanum tuberosum*
- рапс - *Brassica napus*

## Структура проекта

~~~text
genetic_maps_test_task/
├── README.md
└── species/
    ├── zea_mays/
    ├── helianthus_annuus/
    └── brassica_oleracea_capitata/
~~~

Для каждой культуры структура однотипная:

~~~text
species/<species>/
├── data/             # исходные данные, маркеры, метаданные
├── scripts/          # скрипты анализа
├── results/
│   ├── intermediate/ # промежуточные таблицы
│   ├── final/        # финальные карты
│   ├── qc/           # QC-таблицы
│   └── figures/      # графики
└── logs/             # логи
~~~

## Краткая сводка по культурам

### 1. Кукуруза — *Zea mays*

Для кукурузы использовались данные MaizeGDB Genetic 1 для linkage groups LG1–LG10.  
Целевой референс: `GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0`.

Основная задача состояла в том, чтобы привести данные MaizeGDB к актуальной сборке и получить финальную таблицу с колонками `chr`, `pos`, `cM`.

Финальные результаты лежат в:

~~~text
species/zea_mays/results/final/
~~~

График покрытия генетической карты:

![Maize genetic map coverage](species/zea_mays/results/figures/maize_genetic_map_coverage.svg)

---

### 2. Подсолнечник — *Helianthus annuus*

Для подсолнечника использовались 25-bp SFP-маркеры из опубликованных данных.  
Целевой референс: `GCF_002127325.2_HanXRQr2.0-SUNRISE`.

Физические координаты маркеров были получены через выравнивание 25-bp последовательностей на референсный геном. В финальную карту вошли только однозначные хиты, соответствующие ожидаемой хромосоме.

Финальный результат: `21,442` markers.

Финальная карта:

~~~text
species/helianthus_annuus/results/final/sunflower_genetic_map.bwa_exact_unique.tsv
~~~

График покрытия генетической карты:

![Sunflower genetic map coverage](species/helianthus_annuus/results/figures/sunflower_genetic_map_coverage.svg)

---

### 3. Капуста белокочанная — *Brassica oleracea* var. *capitata*

Для капусты использовалась опубликованная генетическая карта из SSR- и SNP/SNAP-маркеров.  
Целевой референс: `GCA_018177695.1_Cabbage_OX-heart_923_BVRC`.

Физические координаты определялись через выравнивание праймеров на референс и реконструкцию ампликонов. По итогам были собраны две версии карты: более строгая карта с однозначными хитами и карта с более мягким порогом.

Финальные результаты:

- strict / high-confidence map: `196` markers;
- relaxed mm1 collinear map: `222` markers.

Финальные карты лежат в:

~~~text
species/brassica_oleracea_capitata/results/final/
~~~

График покрытия генетической карты:

![Cabbage genetic map coverage](species/brassica_oleracea_capitata/results/figures/cabbage_genetic_map_coverage.svg)
