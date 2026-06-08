# Cabbage primer mismatch spectrum QC

## Primer-level minimum mismatch distribution

```tsv
min_mismatch_bin	n_primers	percent_primers
0	1435	58.47595762021189
1	147	5.990220048899755
2	50	2.037489812550937
3	10	0.4074979625101874
4	1	0.0407497962510187
5	0	0.0
>5	0	0.0
no_full_length_no_gap	811	33.0480847595762
```

## Marker rescue potential by threshold

```tsv
threshold_per_primer	n_markers_with_both_primers_min_mm_le_threshold	new_markers_vs_previous_threshold	n_markers_with_any_multihit_at_min_mismatch	note
0	383	383	94	Оба праймера имеют full-length no-gap hit с min mismatch <= 0
1	467	84	128	Оба праймера имеют full-length no-gap hit с min mismatch <= 1
2	496	29	143	Оба праймера имеют full-length no-gap hit с min mismatch <= 2
3	504	8	146	Оба праймера имеют full-length no-gap hit с min mismatch <= 3
4	505	1	146	Оба праймера имеют full-length no-gap hit с min mismatch <= 4
5	505	0	146	Оба праймера имеют full-length no-gap hit с min mismatch <= 5
```

## Marker rescue potential by threshold and marker type

```tsv
threshold_per_primer	marker_type	n_markers
0	SSR	381
0	SNP	2
1	SSR	461
1	SNP	6
2	SSR	486
2	SNP	10
3	SSR	489
3	SNP	15
4	SSR	489
4	SNP	16
5	SSR	489
5	SNP	16
```

## Output files

- `cabbage_primer_mismatch_spectrum.tsv`
- `cabbage_primer_min_mismatch_distribution.tsv`
- `cabbage_marker_min_mismatch_spectrum.tsv`
- `cabbage_marker_rescue_potential_by_threshold.tsv`
- `cabbage_marker_rescue_potential_by_threshold_and_type.tsv`
- `cabbage_marker_rescue_candidates_mm2_mm3.tsv`
